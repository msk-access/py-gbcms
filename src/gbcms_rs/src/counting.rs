use pyo3::prelude::*;
use rust_htslib::bam::record::{Cigar, CigarStringView};
use rust_htslib::bam::{self, Read, Record};
use std::collections::{HashMap, HashSet};

use crate::stats::fisher_strand_bias;
use crate::types::{BaseCounts, Variant};

/// Count bases for a list of variants in a BAM file.
#[pyfunction]
pub fn count_bam(
    bam_path: String,
    variants: Vec<Variant>,
    min_mapq: u8,
    min_baseq: u8,
    filter_duplicates: bool,
    filter_secondary: bool,
    filter_supplementary: bool,
) -> PyResult<Vec<BaseCounts>> {
    let mut bam = bam::IndexedReader::from_path(&bam_path).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to open BAM: {}", e))
    })?;

    let mut results = Vec::with_capacity(variants.len());

    for variant in variants {
        let counts = count_single_variant(
            &mut bam,
            &variant,
            min_mapq,
            min_baseq,
            filter_duplicates,
            filter_secondary,
            filter_supplementary,
        )?;
        results.push(counts);
    }

    Ok(results)
}

fn count_single_variant(
    bam: &mut bam::IndexedReader,
    variant: &Variant,
    min_mapq: u8,
    min_baseq: u8,
    filter_duplicates: bool,
    filter_secondary: bool,
    filter_supplementary: bool,
) -> PyResult<BaseCounts> {
    let tid = bam.header().tid(variant.chrom.as_bytes()).ok_or_else(|| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "Chromosome not found in BAM: {}",
            variant.chrom
        ))
    })?;

    // Fetch region. Add buffer for indels.
    // Variant pos is 0-based.
    let start = variant.pos;
    let end = variant.pos + 1; // Fetch at least one base (anchor)

    bam.fetch((tid, start, end)).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Failed to fetch region: {}", e))
    })?;

    let mut counts = BaseCounts::default();

    // Fragment tracking
    // fragment_name -> (read1_orientation, read2_orientation)
    // orientation: true = forward, false = reverse
    let mut fragment_orientations: HashMap<String, HashMap<u8, bool>> = HashMap::new();

    // fragment_name -> (has_ref, has_alt)
    let mut fragment_alleles: HashMap<String, (bool, bool)> = HashMap::new();

    for result in bam.records() {
        let record = result.map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
                "Error reading BAM record: {}",
                e
            ))
        })?;

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
        if record.mapq() < min_mapq {
            continue;
        }

        // Determine allele status
        let (is_ref, is_alt) = check_allele(&record, variant, min_baseq);

        if !is_ref && !is_alt {
            continue;
        }

        // Update basic counts
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

        // Track fragment info
        let qname = String::from_utf8_lossy(record.qname()).to_string();
        let is_read1 = record.is_first_in_template();
        let read_idx = if is_read1 { 1 } else { 2 };

        fragment_orientations
            .entry(qname.clone())
            .or_default()
            .insert(read_idx, !is_reverse); // Store orientation (forward = true)

        let entry = fragment_alleles.entry(qname).or_insert((false, false));
        if is_ref {
            entry.0 = true;
        }
        if is_alt {
            entry.1 = true;
        }
    }

    // Process fragments (Majority Rule)
    for (qname, orientations) in fragment_orientations {
        if !fragment_alleles.contains_key(&qname) {
            continue;
        }
        let (has_ref, has_alt) = fragment_alleles[&qname];

        // Determine orientation
        let orientation = if orientations.len() == 2 {
            // Both reads present
            let o1 = orientations.get(&1);
            let o2 = orientations.get(&2);
            if o1 == o2 {
                *o1.unwrap() // Consistent orientation
            } else {
                *o1.unwrap() // Fallback to read 1 (or handle discordance)
            }
        } else {
            // Single read
            *orientations.values().next().unwrap()
        };

        counts.dpf += 1;

        if has_ref {
            counts.rdf += 1;
            if orientation {
                counts.rdf_fwd += 1;
            } else {
                counts.rdf_rev += 1;
            }
        }

        if has_alt {
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
fn check_allele(record: &Record, variant: &Variant, min_baseq: u8) -> (bool, bool) {
    let variant_type = &variant.variant_type;

    match variant_type.as_str() {
        "SNP" => check_snp(record, variant, min_baseq),
        "INSERTION" => check_insertion(record, variant),
        "DELETION" => check_deletion(record, variant),
        _ => (false, false), // Complex not fully supported yet
    }
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

    let cigar = record.cigar();
    let mut ref_pos = record.pos();
    let mut read_pos = 0;

    for op in cigar.iter() {
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                // Check if we are at the anchor position
                let len_i64 = *len as i64;
                if variant.pos >= ref_pos && variant.pos < ref_pos + len_i64 {
                    // We are inside a match block.
                    // If the insertion is immediately after this block?
                    // Or if we are at the end of this block and the next op is Ins?

                    // Calculate offset in this block
                    let offset = variant.pos - ref_pos;

                    // If anchor is the LAST base of this match block
                    if offset == len_i64 - 1 {
                        // Check next operation
                        // This requires peeking, which is hard with iterator.
                        // Simpler: Check if we found the anchor.
                        // If we are here, we matched the anchor.
                        // Now we need to see if the NEXT op is an insertion matching our ALT.
                        // This logic is tricky inside a single pass loop.
                        // Let's use a specialized finder.
                    }
                }
                ref_pos += len_i64;
                read_pos += *len as usize;
            }
            Cigar::Ins(len) => {
                // If we just passed the anchor (ref_pos is now anchor + 1),
                // then this insertion is AT the anchor.
                // Wait, ref_pos is the position on the reference.
                // If we had a Match(10) starting at 0. ref_pos goes 0->10.
                // If anchor is at 9.
                // Next op is Ins(1).
                // The insertion is between 9 and 10.
                read_pos += *len as usize;
            }
            _ => {
                // Handle other ops
                match op {
                    Cigar::Del(len) | Cigar::RefSkip(len) => ref_pos += *len as i64,
                    Cigar::SoftClip(len) => read_pos += *len as usize,
                    _ => {}
                }
            }
        }
    }

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
                        if i + 1 < cigar_view.len() {
                            if let Cigar::Ins(ins_len) = cigar_view[i + 1] {
                                // Found insertion!
                                // Check length and sequence
                                let ins_len_usize = ins_len as usize;
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
                        if i + 1 < cigar_view.len() {
                            if let Cigar::Del(del_len) = cigar_view[i + 1] {
                                // Found deletion!
                                let expected_del_len = variant.ref_allele.len() - 1;
                                if del_len as usize == expected_del_len {
                                    return (false, true); // ALT
                                }
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
