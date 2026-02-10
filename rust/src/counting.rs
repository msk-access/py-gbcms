use pyo3::prelude::*;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{self, Read, Record};
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;

use crate::stats::fisher_strand_bias;
use crate::types::{BaseCounts, Variant};

use rayon::prelude::*;

use anyhow::{Context, Result};
use log::debug;

/// Evidence accumulated for a single fragment (read pair) at a variant site.
/// Tracks the best base quality seen for each allele across both reads,
/// enabling quality-weighted consensus to resolve R1-vs-R2 conflicts.
#[derive(Debug, Clone)]
struct FragmentEvidence {
    /// Best base quality seen supporting REF across reads in this fragment
    best_ref_qual: u8,
    /// Best base quality seen supporting ALT across reads in this fragment
    best_alt_qual: u8,
    /// Orientation of read 1 (true = forward, false = reverse), if seen
    read1_orientation: Option<bool>,
    /// Orientation of read 2 (true = forward, false = reverse), if seen
    read2_orientation: Option<bool>,
}

impl FragmentEvidence {
    fn new() -> Self {
        FragmentEvidence {
            best_ref_qual: 0,
            best_alt_qual: 0,
            read1_orientation: None,
            read2_orientation: None,
        }
    }

    /// Record a read's allele call and orientation into this fragment's evidence.
    fn observe(&mut self, is_ref: bool, is_alt: bool, base_qual: u8, is_read1: bool, is_forward: bool) {
        if is_ref && base_qual > self.best_ref_qual {
            self.best_ref_qual = base_qual;
        }
        if is_alt && base_qual > self.best_alt_qual {
            self.best_alt_qual = base_qual;
        }
        if is_read1 {
            self.read1_orientation = Some(is_forward);
        } else {
            self.read2_orientation = Some(is_forward);
        }
    }

    /// Resolve this fragment's allele call using quality-weighted consensus.
    /// Returns (is_ref, is_alt) — exactly one will be true.
    fn resolve(&self, qual_diff_threshold: u8) -> (bool, bool) {
        let has_ref = self.best_ref_qual > 0;
        let has_alt = self.best_alt_qual > 0;

        match (has_ref, has_alt) {
            (true, false) => (true, false),   // Only REF evidence
            (false, true) => (false, true),   // Only ALT evidence
            (true, true) => {
                // Conflict: both alleles seen across reads in this fragment.
                // Use quality-weighted consensus: higher quality wins.
                // If quality difference is within threshold, call REF (conservative).
                if self.best_ref_qual > self.best_alt_qual + qual_diff_threshold {
                    (true, false)  // REF wins by quality margin
                } else if self.best_alt_qual > self.best_ref_qual + qual_diff_threshold {
                    (false, true)  // ALT wins by quality margin
                } else {
                    // Within threshold — conservative: call REF
                    (true, false)
                }
            }
            (false, false) => (false, false),  // Should not happen (filtered earlier)
        }
    }

    /// Get the fragment orientation. Prefers read 1 orientation when both reads present.
    fn orientation(&self) -> Option<bool> {
        // Prefer read 1 orientation (consistent with original behavior)
        self.read1_orientation.or(self.read2_orientation)
    }
}

/// Hash a QNAME to u64 for memory-efficient fragment tracking.
/// Using DefaultHasher for speed — collision probability is negligible
/// for typical variant-level read counts (~1000 fragments).
#[inline]
fn hash_qname(qname: &[u8]) -> u64 {
    let mut hasher = DefaultHasher::new();
    qname.hash(&mut hasher);
    hasher.finish()
}

/// Count bases for a list of variants in a BAM file.
#[pyfunction]
pub fn count_bam(
    py: Python<'_>,
    bam_path: String,
    variants: Vec<Variant>,
    min_mapq: u8,
    min_baseq: u8,
    filter_duplicates: bool,
    filter_secondary: bool,
    filter_supplementary: bool,
    filter_qc_failed: bool,
    filter_improper_pair: bool,
    filter_indel: bool,
    threads: usize,
) -> PyResult<Vec<BaseCounts>> {
    // We cannot share a single IndexedReader across threads because it's not Sync.
    // Instead, we use rayon's map_init to initialize a reader for each thread.
    // This is efficient because map_init reuses the thread-local state (the reader)
    // for multiple items processed by that thread.

    // Configure thread pool
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Failed to build thread pool: {}", e)))?;

    // Release GIL for parallel execution
    #[allow(deprecated)]
    let results: Result<Vec<BaseCounts>, anyhow::Error> = py.allow_threads(move || {
        pool.install(|| {
            variants
                .par_iter()
                .map_init(
                    || {
                        // Initialize thread-local BAM reader
                        bam::IndexedReader::from_path(&bam_path).map_err(|e| {
                            anyhow::anyhow!("Failed to open BAM: {}", e)
                        })
                    },
                    |bam_result, variant| {
                        // Get the reader or return error if initialization failed
                        let bam = match bam_result {
                            Ok(b) => b,
                            Err(e) => return Err(anyhow::anyhow!("BAM init failed: {}", e)),
                        };

                        count_single_variant(
                            bam,
                            variant,
                            min_mapq,
                            min_baseq,
                            filter_duplicates,
                            filter_secondary,
                            filter_supplementary,
                            filter_qc_failed,
                            filter_improper_pair,
                            filter_indel,
                        )
                    },
                )
                .collect()
        })
    });

    // Map anyhow::Error back to PyErr
    match results {
        Ok(r) => Ok(r),
        Err(e) => Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("{}", e))),
    }
}

fn count_single_variant(
    bam: &mut bam::IndexedReader,
    variant: &Variant,
    min_mapq: u8,
    min_baseq: u8,
    filter_duplicates: bool,
    filter_secondary: bool,
    filter_supplementary: bool,
    filter_qc_failed: bool,
    filter_improper_pair: bool,
    filter_indel: bool,
) -> Result<BaseCounts> {
    let tid = bam.header().tid(variant.chrom.as_bytes()).ok_or_else(|| {
        anyhow::anyhow!("Chromosome not found in BAM: {}", variant.chrom)
    })?;

    // Fetch region. Add buffer for indels.
    // Variant pos is 0-based.
    let start = variant.pos;
    let end = variant.pos + 1; // Fetch at least one base (anchor)

    bam.fetch((tid, start, end)).context("Failed to fetch region")?;

    let mut counts = BaseCounts::default();

    // Fragment tracking: QNAME hash -> FragmentEvidence
    // Using u64 hash keys instead of String for memory efficiency.
    let mut fragments: HashMap<u64, FragmentEvidence> = HashMap::new();

    // Quality threshold for fragment consensus tiebreaking.
    // When R1 and R2 disagree, the allele with higher quality wins
    // only if the quality difference exceeds this threshold.
    let qual_diff_threshold: u8 = 10;

    for result in bam.records() {
        let record = result.context("Error reading BAM record")?;

        // Apply filters
        if filter_duplicates && record.is_duplicate() {
            continue;
        }
        if filter_secondary && record.is_secondary() {
            continue;
        }
        if filter_supplementary && record.is_supplementary() {
            continue;
        }
        if filter_qc_failed && record.is_quality_check_failed() {
            continue;
        }
        if filter_improper_pair && !record.is_proper_pair() {
            continue;
        }

        // Indel filter: check CIGAR for Ins or Del
        if filter_indel {
            let has_indel = record.cigar().iter().any(|op| matches!(op, Cigar::Ins(_) | Cigar::Del(_)));
            if has_indel {
                continue;
            }
        }

        if record.mapq() < min_mapq {
            continue;
        }

        // Determine allele status and base quality at variant position
        let (is_ref, is_alt, base_qual) = check_allele_with_qual(&record, variant, min_baseq);

        if !is_ref && !is_alt {
            continue;
        }

        // Update per-read counts (unchanged — these are read-level, not fragment-level)
        counts.dp += 1;
        let is_reverse = record.is_reverse();

        if is_reverse {
            counts.dp_rev += 1;
        } else {
            counts.dp_fwd += 1;
        }

        if is_ref {
            counts.rd += 1;
            if is_reverse {
                counts.rd_rev += 1;
            } else {
                counts.rd_fwd += 1;
            }
        } else if is_alt {
            counts.ad += 1;
            if is_reverse {
                counts.ad_rev += 1;
            } else {
                counts.ad_fwd += 1;
            }
        }

        // Accumulate fragment evidence using QNAME hash
        let qname_hash = hash_qname(record.qname());
        let is_read1 = record.is_first_in_template();
        let is_forward = !is_reverse;

        let evidence = fragments.entry(qname_hash).or_insert_with(FragmentEvidence::new);
        evidence.observe(is_ref, is_alt, base_qual, is_read1, is_forward);
    }

    // Resolve fragment-level counts using quality-weighted consensus.
    // Each fragment contributes exactly ONE allele call (REF xor ALT),
    // preventing the double-counting bug where R1=REF + R2=ALT
    // inflated both rdf and adf.
    for (_qname_hash, evidence) in &fragments {
        let (frag_ref, frag_alt) = evidence.resolve(qual_diff_threshold);

        if !frag_ref && !frag_alt {
            continue; // Should not happen, but guard defensively
        }

        // Get fragment orientation (prefer read 1)
        let orientation = match evidence.orientation() {
            Some(o) => o,
            None => continue, // No orientation data (should not happen)
        };

        counts.dpf += 1;

        if frag_ref {
            counts.rdf += 1;
            if orientation {
                counts.rdf_fwd += 1;
            } else {
                counts.rdf_rev += 1;
            }
        } else if frag_alt {
            counts.adf += 1;
            if orientation {
                counts.adf_fwd += 1;
            } else {
                counts.adf_rev += 1;
            }
        }
    }

    // Calculate stats
    let (sb_pval, sb_or) =
        fisher_strand_bias(counts.rd_fwd, counts.rd_rev, counts.ad_fwd, counts.ad_rev);
    counts.sb_pval = sb_pval;
    counts.sb_or = sb_or;

    let (fsb_pval, fsb_or) = fisher_strand_bias(
        counts.rdf_fwd,
        counts.rdf_rev,
        counts.adf_fwd,
        counts.adf_rev,
    );
    counts.fsb_pval = fsb_pval;
    counts.fsb_or = fsb_or;

    Ok(counts)
}

/// Check if a read supports the reference or alternate allele.
/// Returns (is_ref, is_alt, base_quality) where base_quality is the
/// quality score at the variant position (used for fragment consensus).
fn check_allele_with_qual(record: &Record, variant: &Variant, min_baseq: u8) -> (bool, bool, u8) {
    let variant_type = &variant.variant_type;
    debug!("check_allele type={} pos={} ref={} alt={}", variant_type, variant.pos, variant.ref_allele, variant.alt_allele);

    let (is_ref, is_alt) = match variant_type.as_str() {
        "SNP" => check_snp(record, variant, min_baseq),
        "INSERTION" => check_insertion(record, variant),
        "DELETION" => check_deletion(record, variant),
        "MNP" => {
            debug!("Dispatching to check_mnp");
            check_mnp(record, variant, min_baseq)
        },
        "COMPLEX" => check_complex(record, variant, min_baseq),
        _ => {
            // Auto-detect MNP if type is not explicit but looks like one
            if variant.ref_allele.len() == variant.alt_allele.len()
                && variant.ref_allele.len() > 1
            {
                debug!("Auto-detected MNP");
                check_mnp(record, variant, min_baseq)
            } else {
                // Fallback to complex check for anything else (e.g. DelIns)
                check_complex(record, variant, min_baseq)
            }
        }
    };

    // Extract base quality at variant position for fragment consensus.
    // For SNPs/MNPs: quality at the variant position.
    // For indels: quality of the anchor base (position before indel).
    // For complex: minimum quality across the reconstructed haplotype.
    let base_qual = if is_ref || is_alt {
        get_variant_base_qual(record, variant)
    } else {
        0
    };

    (is_ref, is_alt, base_qual)
}

/// Get the base quality at the variant position for fragment consensus scoring.
/// Returns the quality of the most relevant base for this variant type.
fn get_variant_base_qual(record: &Record, variant: &Variant) -> u8 {
    match variant.variant_type.as_str() {
        "SNP" => {
            // Quality at the variant position
            match find_read_pos(record, variant.pos) {
                Some(p) => record.qual()[p],
                None => 0,
            }
        }
        "MNP" => {
            // Minimum quality across MNP positions
            let len = variant.ref_allele.len();
            let start = match find_read_pos(record, variant.pos) {
                Some(p) => p,
                None => return 0,
            };
            let quals = record.qual();
            (0..len)
                .filter_map(|i| {
                    let p = start + i;
                    if p < quals.len() { Some(quals[p]) } else { None }
                })
                .min()
                .unwrap_or(0)
        }
        "INSERTION" | "DELETION" => {
            // Quality of the anchor base (position before the indel)
            match find_read_pos(record, variant.pos) {
                Some(p) => record.qual()[p],
                None => 0,
            }
        }
        _ => {
            // Complex / unknown: quality at variant start position
            match find_read_pos(record, variant.pos) {
                Some(p) => record.qual()[p],
                None => 0,
            }
        }
    }
}

fn check_complex(record: &Record, variant: &Variant, min_baseq: u8) -> (bool, bool) {
    // Haplotype Reconstruction Strategy
    // 1. Identify genomic region of REF allele: [start, end)
    // 2. Traverse CIGAR to reconstruct the sequence the read has for this region.
    //    - Match/Mismatch: Append base.
    //    - Insertion: Append inserted bases (if occurring within/at start of region).
    //    - Deletion: Skip (do not append).
    // 3. Compare reconstructed sequence to ALT allele.

    let start_pos = variant.pos;
    let end_pos = variant.pos + variant.ref_allele.len() as i64; // Exclusive end

    let cigar = record.cigar();
    let mut ref_pos = record.pos();
    let mut read_pos = 0;
    let seq = record.seq();
    let quals = record.qual();

    let mut reconstructed_seq = Vec::new();
    let mut min_qual_observed = 255u8;

    debug!("check_complex start: pos={} ref={} alt={}", start_pos, variant.ref_allele, variant.alt_allele);

    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len_i64 = *len as i64;
                let len_usize = *len as usize;

                // Intersection of [ref_pos, ref_pos + len) and [start_pos, end_pos)
                let overlap_start = std::cmp::max(ref_pos, start_pos);
                let overlap_end = std::cmp::min(ref_pos + len_i64, end_pos);

                debug!("Match op len={} ref_pos={} overlap=[{}, {})", len, ref_pos, overlap_start, overlap_end);

                if overlap_start < overlap_end {
                    // We have overlap
                    let offset_in_op = (overlap_start - ref_pos) as usize;
                    let overlap_len = (overlap_end - overlap_start) as usize;

                    let current_read_pos = read_pos + offset_in_op;
                    
                    for i in 0..overlap_len {
                        let p = current_read_pos + i;
                        if p >= seq.len() { break; } // Should not happen if BAM is valid
                        
                        let base = seq[p];
                        let qual = quals[p];
                        reconstructed_seq.push(base);
                        if qual < min_qual_observed {
                            min_qual_observed = qual;
                        }
                    }
                }
                ref_pos += len_i64;
                read_pos += len_usize;
            }
            Cigar::Ins(len) => {
                let len_usize = *len as usize;
                debug!("Ins op len={} ref_pos={}", len, ref_pos);
                
                if ref_pos >= start_pos && ref_pos <= end_pos {
                     for i in 0..len_usize {
                        let p = read_pos + i;
                        if p >= seq.len() { break; }
                        let base = seq[p];
                        let qual = quals[p];
                         reconstructed_seq.push(base);
                         if qual < min_qual_observed {
                            min_qual_observed = qual;
                        }
                     }
                }
                
                read_pos += len_usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                debug!("Del op len={} ref_pos={}", len, ref_pos);
                ref_pos += *len as i64;
            }
            Cigar::SoftClip(len) => {
                read_pos += *len as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }

    let reconstructed_str = String::from_utf8_lossy(&reconstructed_seq);
    debug!("Reconstructed: '{}' (len={})", reconstructed_str, reconstructed_seq.len());

    // Compare reconstructed sequence to ALT
    let alt_bytes = variant.alt_allele.as_bytes();
    
    // Convert reconstructed to bytes for comparison (it's already Vec<u8>)
    // But we need to handle case sensitivity
    if reconstructed_seq.len() != alt_bytes.len() {
        debug!("Len mismatch ALT: {} != {}", reconstructed_seq.len(), alt_bytes.len());
    } else {
        let mut matches_alt = true;
        for (i, &b) in reconstructed_seq.iter().enumerate() {
            if b.to_ascii_uppercase() != alt_bytes[i].to_ascii_uppercase() {
                matches_alt = false;
                break;
            }
        }
        
        if matches_alt {
            if min_qual_observed < min_baseq {
                debug!("Matches ALT but low qual: {} < {}", min_qual_observed, min_baseq);
                return (false, false);
            }
            debug!("Matches ALT");
            return (false, true); // Matches ALT
        }
    }
    
    // Check REF?
    let ref_bytes = variant.ref_allele.as_bytes();
    debug!("Checking REF. ReconLen={} RefLen={}", reconstructed_seq.len(), ref_bytes.len());
    
    if reconstructed_seq.len() == ref_bytes.len() {
        let mut matches_ref = true;
        for (i, &b) in reconstructed_seq.iter().enumerate() {
            let rb = ref_bytes[i].to_ascii_uppercase();
            let bb = b.to_ascii_uppercase();
            if bb != rb {
                debug!("Mismatch at {}: {} != {}", i, bb as char, rb as char);
                matches_ref = false;
                break;
            }
        }
        if matches_ref {
             if min_qual_observed < min_baseq {
                debug!("Matches REF but low qual: {} < {}", min_qual_observed, min_baseq);
                return (false, false);
            }
            debug!("Matches REF");
            return (true, false);
        }
    }

    debug!("No match");
    (false, false)
}

fn check_snp(record: &Record, variant: &Variant, min_baseq: u8) -> (bool, bool) {
    let read_pos = match find_read_pos(record, variant.pos) {
        Some(p) => p,
        None => return (false, false),
    };

    let qual = record.qual()[read_pos];
    if qual < min_baseq {
        return (false, false);
    }

    let base = record.seq()[read_pos] as char;
    let ref_char = variant
        .ref_allele
        .chars()
        .next()
        .unwrap()
        .to_ascii_uppercase();
    let alt_char = variant
        .alt_allele
        .chars()
        .next()
        .unwrap()
        .to_ascii_uppercase();
    let base_upper = base.to_ascii_uppercase();

    (base_upper == ref_char, base_upper == alt_char)
}

fn check_insertion(record: &Record, variant: &Variant) -> (bool, bool) {
    // For insertion, variant.pos is the anchor base.
    // We check if there is an insertion AFTER this base.
    // VCF: REF=A, ALT=AT. Insertion of T.

    // Re-implement with index access for lookahead
    let cigar_view = record.cigar();
    let mut ref_pos = record.pos();
    let mut read_pos = 0;

    for (i, op) in cigar_view.iter().enumerate() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len_i64 = *len as i64;

                // Check if this block contains the anchor
                if variant.pos >= ref_pos && variant.pos < ref_pos + len_i64 {
                    // Found anchor block.
                    // Is anchor the last base?
                    if variant.pos == ref_pos + len_i64 - 1 {
                        // Anchor is at end of match. Check next op.
                        if let Some(Cigar::Ins(ins_len)) = cigar_view.get(i + 1) {
                            // Found insertion!
                            // Check length and sequence
                            let ins_len_usize = *ins_len as usize;
                            let expected_ins_len = variant.alt_allele.len() - 1; // VCF ALT includes anchor

                            if ins_len_usize == expected_ins_len {
                                // Check sequence
                                // Read pos at start of Ins is current read_pos + len of match
                                let ins_start = read_pos + *len as usize;
                                let ins_seq = &record.seq().as_bytes()
                                    [ins_start..ins_start + ins_len_usize];

                                // Expected: ALT without anchor
                                let expected_seq = &variant.alt_allele.as_bytes()[1..];

                                if ins_seq == expected_seq {
                                    return (false, true); // ALT
                                }
                            }
                        }
                        // If no insertion follows, it's REF
                        return (true, false);
                    } else {
                        // Anchor is in middle of match -> REF
                        return (true, false);
                    }
                }
                ref_pos += len_i64;
                read_pos += *len as usize;
            }
            Cigar::Ins(len) => {
                read_pos += *len as usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_pos += *len as i64;
            }
            Cigar::SoftClip(len) => {
                read_pos += *len as usize;
            }
            _ => {}
        }
    }

    (false, false)
}

fn check_deletion(record: &Record, variant: &Variant) -> (bool, bool) {
    // Deletion: VCF REF=AT, ALT=A. Deletion of T.
    // Anchor is A.

    let cigar_view = record.cigar();
    let mut ref_pos = record.pos();

    for (i, op) in cigar_view.iter().enumerate() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let len_i64 = *len as i64;

                if variant.pos >= ref_pos && variant.pos < ref_pos + len_i64 {
                    if variant.pos == ref_pos + len_i64 - 1 {
                        // Anchor at end of match. Check next op.
                        if let Some(Cigar::Del(del_len)) = cigar_view.get(i + 1) {
                            // Found deletion!
                            let expected_del_len = variant.ref_allele.len() - 1;
                            if *del_len as usize == expected_del_len {
                                return (false, true); // ALT
                            }
                        }
                        return (true, false); // REF
                    } else {
                        return (true, false); // REF
                    }
                }
                ref_pos += len_i64;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                ref_pos += *len as i64;
            }
            _ => {}
        }
    }

    (false, false)
}

fn check_mnp(record: &Record, variant: &Variant, min_baseq: u8) -> (bool, bool) {
    // MNP: REF=AT, ALT=CG. Lengths equal, > 1.
    let len = variant.ref_allele.len();
    
    // Find read position of the first base
    let start_read_pos = match find_read_pos(record, variant.pos) {
        Some(p) => p,
        None => return (false, false),
    };

    // Check if the read covers the entire MNP
    if start_read_pos + len > record.seq().len() {
        return (false, false);
    }

    // Check qualities and sequence
    let quals = record.qual();
    let seq = record.seq();
    let seq_bytes = seq.as_bytes();
    
    let ref_bytes = variant.ref_allele.as_bytes();
    let alt_bytes = variant.alt_allele.as_bytes();

    let mut matches_ref = true;
    let mut matches_alt = true;

    for i in 0..len {
        let pos = start_read_pos + i;
        
        // Check quality
        if quals[pos] < min_baseq {
            return (false, false);
        }

        let base = seq_bytes[pos].to_ascii_uppercase();
        let r = ref_bytes[i].to_ascii_uppercase();
        let a = alt_bytes[i].to_ascii_uppercase();

        if base != r {
            matches_ref = false;
        }
        if base != a {
            matches_alt = false;
        }
    }

    // Ensure no indels in the MNP region?
    // Simple check: Re-verify read pos for the last base of MNP
    let end_read_pos = match find_read_pos(record, variant.pos + len as i64 - 1) {
        Some(p) => p,
        None => return (false, false),
    };

    // If contiguous, end - start should be len - 1
    if end_read_pos - start_read_pos != len - 1 {
        // println!("MNP Debug: Indel detected");
        return (false, false); // Indel detected within MNP
    }

    (matches_ref, matches_alt)
}

/// Find the read index corresponding to a genomic position.
fn find_read_pos(record: &Record, target_pos: i64) -> Option<usize> {
    let cigar = record.cigar();
    let mut ref_pos = record.pos();
    let mut read_pos = 0;

    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                if target_pos >= ref_pos && target_pos < ref_pos + *len as i64 {
                    return Some(read_pos + (target_pos - ref_pos) as usize);
                }
                ref_pos += *len as i64;
                read_pos += *len as usize;
            }
            Cigar::Ins(len) => {
                read_pos += *len as usize;
            }
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                if target_pos >= ref_pos && target_pos < ref_pos + *len as i64 {
                    return None; // Position is deleted
                }
                ref_pos += *len as i64;
            }
            Cigar::SoftClip(len) => {
                read_pos += *len as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    None
}
