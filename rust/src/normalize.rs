//! Variant normalization: left-alignment, MAF anchor resolution, REF validation.
//!
//! Consolidates all FASTA-dependent variant preparation into a single pass:
//! 1. MAF→VCF anchor base fetch (if `is_maf`)
//! 2. REF allele validation against reference
//! 3. bcftools-style left-alignment (`realign_left`)
//! 4. `ref_context` fetch for Smith-Waterman haplotype alignment
//!
//! Uses rayon `par_iter().map_init()` with thread-local FASTA readers,
//! matching the pattern in `counting::count_bam()`.

use std::fs::File;

use bio::io::fasta;
use log::{debug, info, warn};
use pyo3::prelude::*;
use rayon::prelude::*;

use crate::types::Variant;

// ---------------------------------------------------------------------------
// PreparedVariant — return type from prepare_variants()
// ---------------------------------------------------------------------------

/// Result of variant preparation: normalized coords, ref_context, and validation info.
///
/// Every input variant produces exactly one `PreparedVariant`, even if validation
/// fails — this ensures the output always has the same row count as input.
#[pyclass]
#[derive(Debug, Clone)]
pub struct PreparedVariant {
    /// Ready-to-count variant (normalized coords + ref_context populated).
    /// For invalid variants, contains best-effort coords with no ref_context.
    #[pyo3(get)]
    pub variant: Variant,

    /// Validation result: "PASS", "PASS_WARN_REF_CORRECTED",
    /// "PASS_WARN_HOMOPOLYMER_DECOMP", "REF_MISMATCH", or "FETCH_FAILED".
    #[pyo3(get, set)]
    pub validation_status: String,

    /// True if left-alignment changed the variant's coordinates.
    #[pyo3(get)]
    pub was_normalized: bool,

    /// Original 0-based position before any transformation.
    #[pyo3(get)]
    pub original_pos: i64,

    /// Original REF allele before any transformation.
    #[pyo3(get)]
    pub original_ref: String,

    /// Original ALT allele before any transformation.
    #[pyo3(get)]
    pub original_alt: String,

    /// Corrected variant for homopolymer decomposition dual-counting.
    /// When a complex variant spans a homopolymer and appears to be a
    /// miscollapsed D(n)+SNV event, this holds the corrected allele
    /// (e.g., CCCCCC→CCCCT instead of CCCCCC→T).
    /// `None` for normal variants where no decomposition is detected.
    #[pyo3(get)]
    pub decomposed_variant: Option<Variant>,
}

// ---------------------------------------------------------------------------
// check_homopolymer_decomp — detect miscollapsed homopolymer events
// ---------------------------------------------------------------------------

/// Check if a complex variant looks like a miscollapsed homopolymer event.
///
/// Returns the corrected ALT allele if the pattern is detected, `None` otherwise.
///
/// **Pattern**: The REF allele is a homopolymer run (all same base, ≥3bp),
/// the variant is a net deletion (ref_len > alt_len), and the ALT's last base
/// matches the reference base immediately AFTER the REF span — suggesting
/// D(n) + SNV were merged into a larger deletion than the reads actually support.
///
/// # Example
/// ```text
/// ref = CCCCCC, alt = T, next_ref_base = T
/// Homopolymer base = C, alt ends with T == next_ref_base, T ≠ C
/// → corrected alt = CCCCT (delete 1 C, change last C→T)
/// ```
fn check_homopolymer_decomp(
    ref_al: &str,
    alt_al: &str,
    next_ref_base: u8,
) -> Option<String> {
    let ref_bytes = ref_al.as_bytes();
    let alt_bytes = alt_al.as_bytes();

    // Must be a net deletion with a remaining allele
    if ref_bytes.len() <= alt_bytes.len() || alt_bytes.is_empty() {
        return None;
    }

    // REF must be at least 3bp (meaningful homopolymer)
    if ref_bytes.len() < 3 {
        return None;
    }

    // REF must be a homopolymer (all same base)
    let homo_base = ref_bytes[0].to_ascii_uppercase();
    if !ref_bytes
        .iter()
        .all(|&b| b.to_ascii_uppercase() == homo_base)
    {
        return None;
    }

    // ALT's last base must match the next reference base after the run
    let alt_last = alt_bytes.last()?.to_ascii_uppercase();
    if alt_last != next_ref_base.to_ascii_uppercase() {
        return None;
    }

    // ALT's last base must differ from the homopolymer base (it's the "SNV")
    if alt_last == homo_base {
        return None;
    }

    // Construct the corrected ALT allele:
    // Keep (ref_len - net_del) bases from the homopolymer, then replace
    // the last one with the SNV base.
    //
    // For SOX2:  ref=CCCCCC(6), alt=T(1), net_del = 6-1 = 5
    //   keep = 6 - 5 = 1 homopolymer base → but we need alt_len bases total
    //   Actually: corrected = ref[net_del..] with last base → alt_last
    //   ref[5..] = "C", change last to T → "T"... that gives alt_len=1
    //   That's wrong. The correct decomposition:
    //   D(1) removes 1 C → CCCCC remains → change last C→T → CCCCT
    //   So corrected_alt = homopolymer[..ref_len-1] + alt_last_char
    //   = "CCCCC" + "T" = "CCCCT" (5bp) — BUT: what's the corrected ref?
    //   ref stays CCCCCC(6bp), alt becomes CCCCT(5bp) → 6bp→5bp = D(1)+SNV
    //
    // General: corrected alt = ref[0..ref_len-1] as string + alt's last char
    // This assumes D(1) + SNV at the boundary. This is the minimal decomposition.
    let mut corrected = String::with_capacity(ref_bytes.len());
    // All but the last base of the homopolymer
    for &b in &ref_bytes[..ref_bytes.len() - 1] {
        corrected.push(b as char);
    }
    // Replace the last position with the SNV base
    corrected.push(alt_last as char);

    debug!(
        "Homopolymer decomp detected: ref={} alt={} → corrected alt={}",
        ref_al, alt_al, corrected
    );

    Some(corrected)
}

// ---------------------------------------------------------------------------
// left_align_variant — pure bcftools realign_left() algorithm
// ---------------------------------------------------------------------------

/// Left-align a variant using the bcftools `realign_left()` algorithm.
///
/// Three-phase algorithm:
///
/// **Phase 1 — Extend left:**
/// Prepend up to `aln_win` upstream reference bases to both alleles,
/// giving Phase 2 the full repeat context for proper right-trimming.
///
/// **Phase 2 — Right-trim:**
/// While the rightmost base is identical in ref and alt, and both alleles
/// have length > 1 (preserving VCF anchor), remove the rightmost base.
///
/// **Phase 3 — Left-trim:**
/// While the leftmost base is identical and both alleles have length > 1,
/// remove the leftmost base and advance `pos`.
///
/// # Arguments
/// * `pos` — 0-based genomic position of the variant
/// * `ref_allele` — Reference allele bytes (e.g., b"ACTCT")
/// * `alt_allele` — Alternate allele bytes (e.g., b"A")
/// * `reference` — Reference sequence covering the region
/// * `ref_offset` — 0-based genomic start position of `reference`
/// * `aln_win` — Maximum number of bases to extend left (bcftools default: 100)
///
/// # Returns
/// `(new_pos, new_ref, new_alt, was_modified)` — normalized coordinates
pub fn left_align_variant(
    pos: i64,
    ref_allele: &[u8],
    alt_allele: &[u8],
    reference: &[u8],
    ref_offset: i64,
    aln_win: usize,
) -> (i64, Vec<u8>, Vec<u8>, bool) {
    let original_pos = pos;
    let original_ref = ref_allele.to_vec();
    let original_alt = alt_allele.to_vec();

    // SNPs and same-length substitutions: cannot be left-aligned
    if ref_allele.len() == 1 && alt_allele.len() == 1 {
        return (pos, original_ref, original_alt, false);
    }

    // Phase 1: Extend both alleles to the left by up to aln_win bases.
    // This gives Phase 2 the full repeat context to right-trim through.
    let pos_in_ref = (pos - ref_offset) as usize;
    let actual_extend = aln_win.min(pos_in_ref);

    let mut r = Vec::with_capacity(actual_extend + ref_allele.len());
    let mut a = Vec::with_capacity(actual_extend + alt_allele.len());

    // Prepend upstream reference bases to both alleles
    let prefix_start = pos_in_ref - actual_extend;
    for &base in reference.iter().take(pos_in_ref).skip(prefix_start) {
        r.push(base);
        a.push(base);
    }
    let mut cur_pos = pos - actual_extend as i64;

    // Append original alleles
    r.extend_from_slice(ref_allele);
    a.extend_from_slice(alt_allele);

    // Phase 2: Right-trim — remove matching trailing bases
    // Keep at least 1 base in each allele (VCF anchor requirement)
    while r.len() > 1 && a.len() > 1
        && r.last().unwrap().eq_ignore_ascii_case(a.last().unwrap())
    {
        r.pop();
        a.pop();
    }

    // Phase 3: Left-trim — remove matching leading bases
    // Keep at least 1 base in each allele
    while r.len() > 1 && a.len() > 1
        && r[0].eq_ignore_ascii_case(&a[0])
    {
        r.remove(0);
        a.remove(0);
        cur_pos += 1;
    }

    let was_modified = cur_pos != original_pos || r != original_ref || a != original_alt;
    (cur_pos, r, a, was_modified)
}


// ---------------------------------------------------------------------------
// resolve_maf_anchor — MAF-style indel → VCF-style with anchor base
// ---------------------------------------------------------------------------

/// Convert a MAF-style indel to VCF-style by fetching the anchor base.
///
/// Handles three cases based on the `is_maf` flag:
/// - `ref = "-"` (insertion): Anchor is at `start_pos` (1-based), prepend to ALT
/// - `alt = "-"` (deletion): Anchor is at `start_pos - 1` (1-based), prepend to REF
/// - Complex (different-length, non-dash): Anchor at `start_pos - 1`, prepend to both
///
/// # Arguments
/// * `reader` — Thread-local FASTA reader
/// * `chrom` — Chromosome name (tries both as-is and with "chr" prefix)
/// * `start_pos` — MAF Start_Position (1-based)
/// * `ref_allele` — MAF Reference_Allele (may be "-")
/// * `alt_allele` — MAF Tumor_Seq_Allele2 (may be "-")
///
/// # Returns
/// `(pos_0based, vcf_ref, vcf_alt, variant_type)` or error if FASTA fetch fails.
fn resolve_maf_anchor(
    reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    start_pos: i64,
    ref_allele: &str,
    alt_allele: &str,
) -> anyhow::Result<(i64, String, String, String)> {
    let is_insertion = ref_allele == "-";
    let is_deletion = alt_allele == "-";
    let is_snp = ref_allele.len() == 1
        && alt_allele.len() == 1
        && !is_insertion
        && !is_deletion;

    if is_snp {
        // SNPs: MAF Start_Position == VCF POS (1-based), convert to 0-based
        return Ok((
            start_pos - 1,
            ref_allele.to_string(),
            alt_allele.to_string(),
            "SNP".to_string(),
        ));
    }

    // Determine anchor position (0-based) and variant type
    let (anchor_pos_0based, vtype) = if is_insertion {
        // MAF insertion: Start_Position is the base BEFORE the insertion
        // Anchor is at start_pos (1-based) → start_pos - 1 (0-based)
        (start_pos - 1, "INSERTION")
    } else if is_deletion {
        // MAF deletion: Start_Position is the FIRST DELETED base
        // Anchor is one base before → start_pos - 2 (0-based)
        (start_pos - 2, "DELETION")
    } else {
        // Complex: different-length, non-dash alleles
        // Anchor is one base before start → start_pos - 2 (0-based)
        (start_pos - 2, "COMPLEX")
    };

    // Fetch anchor base, trying both chrom names
    let anchor_base = fetch_single_base(reader, chrom, anchor_pos_0based)?;
    let anchor_upper = (anchor_base as char).to_uppercase().to_string();

    // Build VCF-style alleles
    let (vcf_ref, vcf_alt) = if is_insertion {
        (anchor_upper.clone(), format!("{}{}", anchor_upper, alt_allele))
    } else if is_deletion {
        (format!("{}{}", anchor_upper, ref_allele), anchor_upper.clone())
    } else {
        // Complex
        (
            format!("{}{}", anchor_upper, ref_allele),
            format!("{}{}", anchor_upper, alt_allele),
        )
    };

    // VCF POS (0-based) = anchor position
    Ok((anchor_pos_0based, vcf_ref, vcf_alt, vtype.to_string()))
}

/// Fetch a single base from the reference, delegating to `fetch_region`.
fn fetch_single_base(
    reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    pos_0based: i64,
) -> anyhow::Result<u8> {
    let pos = pos_0based as u64;
    let buf = fetch_region(reader, chrom, pos, pos + 1)?;
    buf.first()
        .map(|b| b.to_ascii_uppercase())
        .ok_or_else(|| anyhow::anyhow!(
            "FASTA fetch returned empty for {}:{}", chrom, pos_0based
        ))
}

// ---------------------------------------------------------------------------
// validate_ref — check REF allele against reference genome
// ---------------------------------------------------------------------------

/// Validate that the REF allele matches the reference genome.
///
/// Case-insensitive comparison with tolerance for partial mismatches.
/// Tries both the given chromosome name and with/without "chr" prefix.
///
/// # Returns
/// Tuple of `(status, Option<fasta_ref>)`:
/// - `("PASS", None)` — exact match
/// - `("PASS_WARN_REF_CORRECTED", Some(fasta_ref))` — ≥90% match, corrected
/// - `("REF_MISMATCH", None)` — <90% match, rejected
/// - `("FETCH_FAILED", None)` — region could not be fetched
fn validate_ref(
    reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    pos_0based: i64,
    ref_allele: &str,
) -> (String, Option<String>) {
    let ref_len = ref_allele.len() as u64;
    let pos = pos_0based as u64;

    let fetch_result = fetch_region(reader, chrom, pos, pos + ref_len);
    match fetch_result {
        Ok(ref_seq) => {
            // Fast path: exact match
            if ref_seq.eq_ignore_ascii_case(ref_allele.as_bytes()) {
                ("PASS".to_string(), None)
            } else {
                // Compute similarity: count matching bases (case-insensitive)
                let max_len = ref_seq.len().max(ref_allele.len());
                if max_len == 0 {
                    return ("REF_MISMATCH".to_string(), None);
                }
                let matches = ref_seq
                    .iter()
                    .zip(ref_allele.as_bytes())
                    .filter(|(a, b)| a.eq_ignore_ascii_case(b))
                    .count();
                let similarity = matches as f64 / max_len as f64;

                if similarity >= 0.90 {
                    let fasta_ref = String::from_utf8_lossy(&ref_seq).to_uppercase();
                    warn!(
                        "REF partially mismatched at {}:{} — {}/{} bases match ({:.1}%), \
                         correcting to FASTA REF",
                        chrom,
                        pos + 1,
                        matches,
                        max_len,
                        similarity * 100.0,
                    );
                    ("PASS_WARN_REF_CORRECTED".to_string(), Some(fasta_ref))
                } else {
                    ("REF_MISMATCH".to_string(), None)
                }
            }
        }
        Err(_) => ("FETCH_FAILED".to_string(), None),
    }
}

/// Fetch a region from FASTA, trying chrom variants (with/without chr prefix).
fn fetch_region(
    reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    start: u64,
    end: u64,
) -> anyhow::Result<Vec<u8>> {
    let mut buf = Vec::new();

    // Try as-is
    if reader.fetch(chrom, start, end).is_ok() {
        buf.clear();
        reader.read(&mut buf)?;
        if !buf.is_empty() {
            return Ok(buf);
        }
    }

    // Try with chr prefix
    let chr_name = format!("chr{}", chrom);
    if reader.fetch(&chr_name, start, end).is_ok() {
        buf.clear();
        reader.read(&mut buf)?;
        if !buf.is_empty() {
            return Ok(buf);
        }
    }

    // Try stripping chr prefix
    if let Some(stripped) = chrom.strip_prefix("chr") {
        if reader.fetch(stripped, start, end).is_ok() {
            buf.clear();
            reader.read(&mut buf)?;
            if !buf.is_empty() {
                return Ok(buf);
            }
        }
    }

    anyhow::bail!("FASTA fetch failed for {}:{}-{}", chrom, start, end)
}

// ---------------------------------------------------------------------------
// Tandem repeat detection — for adaptive context padding
// ---------------------------------------------------------------------------

/// Detect the longest tandem repeat (motif 1–6bp) touching a given position.
///
/// Scans the provided sequence `seq` for tandem repeats whose footprint
/// includes `pos_in_seq`.  Returns `(motif_length, total_span)` for the
/// longest repeat found, or `(1, 1)` if no repeat is detected.
///
/// # Algorithm
///
/// For each candidate motif length *k* (1..=6):
/// 1. Extract the motif starting at `pos_in_seq`
/// 2. Extend left while preceding *k* bases match the motif
/// 3. Extend right while following *k* bases match the motif
/// 4. If the span covers ≥ 2 full copies of the motif, record it
///
/// The longest span (breaking ties in favour of larger motif) wins.
fn find_tandem_repeat(seq: &[u8], pos_in_seq: usize) -> (usize, usize) {
    let len = seq.len();
    let mut best_unit = 1;
    let mut best_span = 1;

    for k in 1..=6usize {
        if pos_in_seq + k > len {
            break;
        }

        let motif = &seq[pos_in_seq..pos_in_seq + k];

        // Extend left: how many full motif copies precede pos_in_seq?
        let mut left_copies = 0;
        let mut cursor = pos_in_seq;
        while cursor >= k {
            if &seq[cursor - k..cursor] == motif {
                left_copies += 1;
                cursor -= k;
            } else {
                break;
            }
        }

        // Extend right: how many full motif copies follow pos_in_seq?
        let mut right_copies = 0;
        let mut cursor = pos_in_seq + k;
        while cursor + k <= len {
            if &seq[cursor..cursor + k] == motif {
                right_copies += 1;
                cursor += k;
            } else {
                break;
            }
        }

        let total_copies = left_copies + 1 + right_copies; // 1 = the copy at pos_in_seq
        let span = total_copies * k;

        // Must have ≥2 copies to count as a repeat
        if total_copies >= 2 && span > best_span {
            best_unit = k;
            best_span = span;
        }
    }

    (best_unit, best_span)
}

/// Compute adaptive context padding based on nearby tandem repeats.
///
/// Strategy:
/// 1. Fetch a wide scan window (±`scan_radius` bp) around the variant
/// 2. Run `find_tandem_repeat()` at the variant position within the window
/// 3. `effective = max(default_pad, repeat_span / 2 + 3)`, capped at `max_pad`
/// 4. Falls back to `default_pad` on FASTA fetch failure
///
/// # Arguments
/// * `reader` — FASTA reader (mutable)
/// * `chrom` — Chromosome name
/// * `pos` — 0-based variant position
/// * `ref_len` — Length of the REF allele
/// * `default_pad` — Minimum padding (from `--context-padding`)
/// * `max_pad` — Hard cap on adaptive padding
fn compute_adaptive_padding(
    reader: &mut fasta::IndexedReader<File>,
    chrom: &str,
    pos: i64,
    ref_len: usize,
    default_pad: i64,
    max_pad: i64,
) -> i64 {
    let scan_radius: i64 = 30;
    let scan_start = (pos - scan_radius).max(0);
    let scan_end = pos + ref_len as i64 + scan_radius;

    let scan_seq = match fetch_region(reader, chrom, scan_start as u64, scan_end as u64) {
        Ok(s) => s,
        Err(_) => {
            debug!(
                "Adaptive scan fetch failed for {}:{}-{}, using default padding {}",
                chrom, scan_start, scan_end, default_pad
            );
            return default_pad;
        }
    };

    // Position within the scan window
    let pos_in_scan = (pos - scan_start) as usize;

    let (motif_len, repeat_span) = find_tandem_repeat(&scan_seq, pos_in_scan);
    let adaptive = (repeat_span as i64) / 2 + 3;
    let effective = default_pad.max(adaptive).min(max_pad);

    if effective > default_pad {
        debug!(
            "Adaptive padding for {}:{}: repeat motif={}bp span={}bp → padding {} (default {})",
            chrom, pos + 1, motif_len, repeat_span, effective, default_pad
        );
    }

    effective
}

// ---------------------------------------------------------------------------
// prepare_variants — PyO3 entry point, consolidates all FASTA access
// ---------------------------------------------------------------------------

/// Prepare variants for counting in a single pass over the reference FASTA.
///
/// For each variant, this function performs (in order):
/// 1. **MAF anchor fetch** — if `is_maf`, resolve `-` alleles to VCF-style
/// 2. **REF validation** — check REF against the reference genome
/// 3. **Left-alignment** — bcftools `realign_left()` for indels
/// 4. **ref_context fetch** — flanking sequence for Smith-Waterman alignment
///
/// Uses rayon `par_iter().map_init()` with thread-local FASTA readers,
/// matching the established pattern in `count_bam()`.
///
/// # Arguments
/// * `variants` — Input variants (raw MAF or VCF coords, 0-based)
/// * `fasta_path` — Path to indexed reference FASTA
/// * `context_padding` — Minimum flanking bases for ref_context (e.g. 5)
/// * `is_maf` — If true, perform MAF→VCF anchor resolution for indels
/// * `threads` — Number of rayon worker threads
/// * `adaptive_context` — If true, dynamically increase padding in repeat regions
///
/// # Returns
/// One `PreparedVariant` per input variant, in the same order.
#[pyfunction]
#[pyo3(signature = (variants, fasta_path, context_padding, is_maf, threads=1, adaptive_context=true))]
pub fn prepare_variants(
    py: Python<'_>,
    variants: Vec<Variant>,
    fasta_path: String,
    context_padding: i64,
    is_maf: bool,
    threads: usize,
    adaptive_context: bool,
) -> PyResult<Vec<PreparedVariant>> {
    info!(
        "prepare_variants: {} variants, is_maf={}, context_padding={}, adaptive={}, threads={}",
        variants.len(),
        is_maf,
        context_padding,
        adaptive_context,
        threads,
    );

    // Build rayon thread pool (same pattern as count_bam)
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| {
            pyo3::exceptions::PyRuntimeError::new_err(format!(
                "Failed to build thread pool: {}",
                e
            ))
        })?;

    let fasta_path_clone = fasta_path.clone();

    // Release GIL for parallel execution
    #[allow(deprecated)]
    let results: Result<Vec<PreparedVariant>, anyhow::Error> =
        py.allow_threads(move || {
            pool.install(|| {
                variants
                    .par_iter()
                    .map_init(
                        || {
                            // Thread-local FASTA reader initialization
                            fasta::IndexedReader::from_file(&fasta_path_clone)
                                .map_err(|e| {
                                    anyhow::anyhow!(
                                        "Failed to open FASTA: {}",
                                        e
                                    )
                                })
                        },
                        |reader_result, variant| {
                            prepare_single_variant(
                                reader_result,
                                variant,
                                context_padding,
                                is_maf,
                                adaptive_context,
                            )
                        },
                    )
                    .collect()
            })
        });

    match results {
        Ok(r) => {
            // Log summary
            let valid = r.iter().filter(|p| p.validation_status == "PASS").count();
            let normalized = r.iter().filter(|p| p.was_normalized).count();
            info!(
                "prepare_variants complete: {}/{} valid, {} normalized",
                valid,
                r.len(),
                normalized,
            );
            Ok(r)
        }
        Err(e) => Err(pyo3::exceptions::PyRuntimeError::new_err(format!(
            "prepare_variants failed: {}",
            e
        ))),
    }
}

/// Process a single variant through the full preparation pipeline.
///
/// Steps: MAF anchor → validate REF → left-align → adaptive context → fetch ref_context.
fn prepare_single_variant(
    reader_result: &mut Result<fasta::IndexedReader<File>, anyhow::Error>,
    variant: &Variant,
    context_padding: i64,
    is_maf: bool,
    adaptive_context: bool,
) -> Result<PreparedVariant, anyhow::Error> {
    let reader = reader_result.as_mut().map_err(|e| {
        anyhow::anyhow!("FASTA reader not available: {}", e)
    })?;

    // Save originals before any transformation
    let original_pos = variant.pos;
    let original_ref = variant.ref_allele.clone();
    let original_alt = variant.alt_allele.clone();

    // Step 1: MAF anchor resolution (if needed)
    let (mut pos, mut ref_al, mut alt_al, mut vtype) = if is_maf
        && (variant.ref_allele == "-"
            || variant.alt_allele == "-"
            || variant.ref_allele.len() != variant.alt_allele.len())
    {
        // MAF indel/complex: resolve anchor base
        // variant.pos is 0-based (from maf_to_internal), start_pos is 1-based
        let start_pos_1based = variant.pos + 1;
        match resolve_maf_anchor(
            reader,
            &variant.chrom,
            start_pos_1based,
            &variant.ref_allele,
            &variant.alt_allele,
        ) {
            Ok(result) => result,
            Err(e) => {
                warn!(
                    "MAF anchor fetch failed for {}:{} {}>{}:  {}",
                    variant.chrom,
                    variant.pos + 1,
                    variant.ref_allele,
                    variant.alt_allele,
                    e,
                );
                return Ok(PreparedVariant {
                    variant: variant.clone(),
                    validation_status: "FETCH_FAILED".to_string(),
                    was_normalized: false,
                    original_pos,
                    original_ref,
                    original_alt,
                    decomposed_variant: None,
                });
            }
        }
    } else {
        // VCF-style or MAF SNP: use coords as-is
        (
            variant.pos,
            variant.ref_allele.clone(),
            variant.alt_allele.clone(),
            variant.variant_type.clone(),
        )
    };

    // Step 2: REF validation (with tolerance for partial mismatches)
    let (status, corrected_ref) = validate_ref(reader, &variant.chrom, pos, &ref_al);

    // If REF was partially mismatched but ≥90% similar, correct it to the FASTA REF
    if let Some(fasta_ref) = corrected_ref {
        ref_al = fasta_ref;
    }

    if !status.starts_with("PASS") {
        debug!(
            "REF validation {}: {}:{} {}>{} ({})",
            status,
            variant.chrom,
            pos + 1,
            ref_al,
            alt_al,
            status,
        );
        return Ok(PreparedVariant {
            variant: Variant {
                chrom: variant.chrom.clone(),
                pos,
                ref_allele: ref_al,
                alt_allele: alt_al,
                variant_type: vtype,
                ref_context: None,
                ref_context_start: 0,
            },
            validation_status: status,
            was_normalized: false,
            original_pos,
            original_ref,
            original_alt,
            decomposed_variant: None,
        });
    }

    // Step 3: Left-alignment (only for indels/complex)
    let mut was_normalized = false;
    let is_indel = ref_al.len() != alt_al.len()
        || (ref_al.len() > 1 && alt_al.len() > 1);

    if is_indel {
        let norm_window: i64 = 100; // bcftools default
        let pad = context_padding.max(norm_window);
        let wide_start = (pos - pad).max(0);
        let wide_end = pos + ref_al.len() as i64 + pad;

        if let Ok(wide_ref) = fetch_region(
            reader,
            &variant.chrom,
            wide_start as u64,
            wide_end as u64,
        ) {
            let (new_pos, new_ref, new_alt, modified) = left_align_variant(
                pos,
                ref_al.as_bytes(),
                alt_al.as_bytes(),
                &wide_ref,
                wide_start,
                norm_window as usize,
            );

            if modified {
                debug!(
                    "Left-aligned: {}:{} {}>{} → {}:{} {}>{}",
                    variant.chrom,
                    pos + 1,
                    ref_al,
                    alt_al,
                    variant.chrom,
                    new_pos + 1,
                    String::from_utf8_lossy(&new_ref),
                    String::from_utf8_lossy(&new_alt),
                );
                pos = new_pos;
                ref_al = String::from_utf8(new_ref)
                    .unwrap_or_else(|_| ref_al.clone());
                alt_al = String::from_utf8(new_alt)
                    .unwrap_or_else(|_| alt_al.clone());
                was_normalized = true;

                // Re-determine variant type after normalization
                vtype = if ref_al.len() == 1 && alt_al.len() == 1 {
                    "SNP".to_string()
                } else if ref_al.len() == 1 && alt_al.len() > 1 {
                    "INSERTION".to_string()
                } else if ref_al.len() > 1 && alt_al.len() == 1 {
                    "DELETION".to_string()
                } else {
                    "COMPLEX".to_string()
                };
            }
        } else {
            debug!(
                "Wide ref fetch failed for {}:{}-{}, skipping normalization",
                variant.chrom, wide_start, wide_end,
            );
        }
    }

    // Step 4: Fetch ref_context at (possibly normalized) position
    //         With adaptive_context, padding is increased in repeat regions.
    let (ref_context, ref_context_start) = if is_indel || vtype == "COMPLEX" {
        let effective_padding = if adaptive_context {
            compute_adaptive_padding(
                reader,
                &variant.chrom,
                pos,
                ref_al.len(),
                context_padding,
                50,  // max cap
            )
        } else {
            context_padding
        };

        let ctx_start = (pos - effective_padding).max(0);
        let ctx_end = pos + ref_al.len() as i64 + effective_padding;

        match fetch_region(
            reader,
            &variant.chrom,
            ctx_start as u64,
            ctx_end as u64,
        ) {
            Ok(ctx_bytes) => (
                Some(String::from_utf8_lossy(&ctx_bytes).to_string()),
                ctx_start,
            ),
            Err(_) => {
                warn!(
                    "ref_context fetch failed for {}:{}-{}, SW alignment will be skipped",
                    variant.chrom, ctx_start, ctx_end,
                );
                (None, 0)
            }
        }
    } else {
        // SNPs don't need ref_context
        (None, 0)
    };

    // Step 5: Homopolymer decomposition detection
    // Check if the variant looks like a miscollapsed D(n)+SNV in a homopolymer.
    // If detected, build a corrected Variant for dual-counting.
    let decomposed_variant = if ref_al.len() > alt_al.len() && ref_al.len() >= 3 {
        // Fetch the next reference base after the ref span
        let next_pos = (pos + ref_al.len() as i64) as u64;
        let next_base = fetch_region(reader, &variant.chrom, next_pos, next_pos + 1)
            .ok()
            .and_then(|b| b.first().copied());

        next_base.and_then(|nb| {
            check_homopolymer_decomp(&ref_al, &alt_al, nb).map(|corrected_alt| {
                let decomp_vtype = if ref_al.len() == corrected_alt.len() {
                    if ref_al.len() == 1 { "SNP" } else { "COMPLEX" }
                } else if ref_al.len() > corrected_alt.len() {
                    "DELETION"
                } else {
                    "INSERTION"
                }
                .to_string();

                Variant {
                    chrom: variant.chrom.clone(),
                    pos,
                    ref_allele: ref_al.clone(),
                    alt_allele: corrected_alt,
                    variant_type: decomp_vtype,
                    ref_context: ref_context.clone(),
                    ref_context_start,
                }
            })
        })
    } else {
        None
    };

    Ok(PreparedVariant {
        variant: Variant {
            chrom: variant.chrom.clone(),
            pos,
            ref_allele: ref_al,
            alt_allele: alt_al,
            variant_type: vtype,
            ref_context,
            ref_context_start,
        },
        validation_status: "PASS".to_string(),
        was_normalized,
        original_pos,
        original_ref,
        original_alt,
        decomposed_variant,
    })
}

// ===========================================================================
// Unit Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // -- left_align_variant tests --

    #[test]
    fn test_snp_passthrough() {
        // SNPs should never be modified
        let reference = b"ATCGATCG";
        let (pos, r, a, modified) =
            left_align_variant(3, b"G", b"T", reference, 0, 100);
        assert_eq!(pos, 3);
        assert_eq!(r, b"G");
        assert_eq!(a, b"T");
        assert!(!modified);
    }

    #[test]
    fn test_deletion_in_homopolymer() {
        // Reference: GAAAAAAG (poly-A at pos 1-6)
        //            01234567
        // Deletion at pos 5: REF=AA, ALT=A (delete one A)
        // Should normalize to pos 0: REF=GA, ALT=G
        let reference = b"GAAAAAAG";
        let (pos, r, a, modified) =
            left_align_variant(5, b"AA", b"A", reference, 0, 100);

        assert!(modified, "Deletion in homopolymer should be left-shifted");
        assert_eq!(pos, 0, "Should shift to leftmost position");
        assert_eq!(r, b"GA");
        assert_eq!(a, b"G");
    }

    #[test]
    fn test_insertion_in_homopolymer() {
        // Reference: GAAAAAAG (poly-A at pos 1-6)
        //            01234567
        // Insertion at pos 5: REF=A, ALT=AA (insert one A)
        // Should normalize to pos 0: REF=G, ALT=GA
        let reference = b"GAAAAAAG";
        let (pos, r, a, modified) =
            left_align_variant(5, b"A", b"AA", reference, 0, 100);

        assert!(modified, "Insertion in homopolymer should be left-shifted");
        assert_eq!(pos, 0, "Should shift to leftmost position");
        assert_eq!(r, b"G");
        assert_eq!(a, b"GA");
    }

    #[test]
    fn test_deletion_in_dinuc_repeat() {
        // Reference: CTCTCTCTCTCT (CT repeat)
        //            0 1 2 3 4 5 6 7 8 9 10 11
        // Deletion at pos 4: REF=CTCT, ALT=CT (delete one CT unit)
        // ref[4..8] = C,T,C,T ✓
        // After extend left 4 + right-trim + left-trim:
        //   Extended: r=CTCTCTCT, a=CTCTCT, cur_pos=0
        //   Right-trim: removes matching pairs until a.len()==1 → r=CTC, a=C
        //   Left-trim: a.len()==1 → no trim
        // Result: pos=0, REF=CTC, ALT=C (minimal representation)
        let reference = b"CTCTCTCTCTCT";
        let (pos, r, a, modified) =
            left_align_variant(4, b"CTCT", b"CT", reference, 0, 100);

        assert!(modified, "Deletion in dinuc repeat should be left-shifted");
        assert_eq!(pos, 0, "Should shift to leftmost position");
        assert_eq!(r, b"CTC");
        assert_eq!(a, b"C");
    }

    #[test]
    fn test_insertion_in_dinuc_repeat() {
        // Reference: CTCTCTCTCTCT (CT repeat)
        //            0 1 2 3 4 5 6 7 8 9 10 11
        // Insertion at pos 4: REF=C, ALT=CCT (insert one CT unit)
        // ref[4] = C ✓
        // After extend left 4: r=CTCTC, a=CTCTCCT
        // Right-trim: C vs T → no match (immediate stop since last bases differ)
        //
        // The algorithm correctly can't shift this because the trailing bases
        // are different. We need to ensure the alleles share a trailing context.
        // Instead, let's use: pos=4, REF=CT, ALT=CTCT (insert CT before next CT)
        // ref[4..6] = C,T ✓
        // After extend left 4: r=CTCTCT, a=CTCTCTCT
        // Right-trim: T==T, C==C, T==T, C==C → trim 4 → r=CT, a=CTCT
        //   Keep trimming: T==T → r=C, a=CTC. C==C → r="" stop (len=1 for r)
        //   Actually r.len()>1 fails. So r=C, a=CTC.
        // Wait, let me retrace. After extend 4: r=[C,T,C,T,C,T], a=[C,T,C,T,C,T,C,T]
        // Right: T==T → pop → r=[C,T,C,T,C], a=[C,T,C,T,C,T,C]
        // Right: C==C → pop → r=[C,T,C,T], a=[C,T,C,T,C,T]
        // Right: T==T → pop → r=[C,T,C], a=[C,T,C,T,C]
        // Right: C==C → pop → r=[C,T], a=[C,T,C,T]
        // Right: T==T → pop → r=[C], a=[C,T,C] — stop (r.len()==1)
        // Left: C==C → r.len()==1, stop
        // Result: pos=0+3=3? No, cur_pos = 4-4=0 (from extend), then left-trim doesn't happen.
        // Final: pos=0, r=[C], a=[C,T,C]  = REF=C, ALT=CTC
        //
        // That's correct! pos=0, REF=C, ALT=CTC is the leftmost representation.
        let reference = b"CTCTCTCTCTCT";
        let (pos, r, a, modified) =
            left_align_variant(4, b"CT", b"CTCT", reference, 0, 100);

        assert!(modified, "Insertion in dinuc repeat should be left-shifted");
        assert_eq!(pos, 0, "Should shift to leftmost position");
        assert_eq!(r, b"C");
        assert_eq!(a, b"CTC");
    }

    #[test]
    fn test_no_shift_needed() {
        // Reference: ATCGATCG — non-repetitive
        // Deletion at pos=1: REF=TC, ALT=T
        // ref[1..3] = T,C ✓
        // After extend 1: r=[A,T,C], a=[A,T]
        // Right-trim: C vs T → no match
        // Left-trim: A==A → r=[T,C], a=[T], pos=1. T==T → r.len()>1 but a.len()==1, stop
        // But now pos=1 and r=[T,C], a=[T] — same as original!
        // Actually this IS already leftmost. Let me verify:
        // Original: pos=1, REF=TC, ALT=T. After normalize: pos=1, REF=TC, ALT=T.
        let reference = b"ATCGATCG";
        let (pos, r, a, modified) =
            left_align_variant(1, b"TC", b"T", reference, 0, 100);
        assert_eq!(pos, 1);
        assert_eq!(r, b"TC");
        assert_eq!(a, b"T");
        assert!(!modified);
    }

    #[test]
    fn test_pos_at_start() {
        // Variant at position 0 — no room to shift left
        // Reference: AATCGATCG
        //            012345678
        // pos=0, REF=AA, ALT=A — deletion of one A
        let reference = b"AATCGATCG";
        let (pos, r, a, _modified) =
            left_align_variant(0, b"AA", b"A", reference, 0, 100);
        assert_eq!(pos, 0, "Cannot shift past position 0");
        assert_eq!(r, b"AA");
        assert_eq!(a, b"A");
    }

    #[test]
    fn test_with_ref_offset() {
        // Reference window starts at offset 100 (simulates fetching a window)
        // Poly-A: reference = GAAAAAAG at positions 100-107
        // Deletion at pos=105: REF=AA, ALT=A
        // ref_in_window[5..7] = A,A ✓
        let reference = b"GAAAAAAG";
        let (pos, r, a, modified) =
            left_align_variant(105, b"AA", b"A", reference, 100, 100);

        assert!(modified);
        assert_eq!(pos, 100, "Should shift to start of window");
        assert_eq!(r, b"GA");
        assert_eq!(a, b"G");
    }

    #[test]
    fn test_complex_no_normalize() {
        // Complex variant (MNP-like): REF=AG, ALT=TC — same length, no shared bases
        let reference = b"AAAGTTTT";
        let (pos, r, a, modified) =
            left_align_variant(2, b"AG", b"TC", reference, 0, 100);
        assert_eq!(pos, 2);
        assert_eq!(r, b"AG");
        assert_eq!(a, b"TC");
        assert!(!modified);
    }

    // -- find_tandem_repeat tests --

    #[test]
    fn test_find_tandem_repeat_homopolymer() {
        // GATCAAAAAAGCTT — 6 A's at pos 4-9
        let seq = b"GATCAAAAAAGCTT";
        let (unit, span) = find_tandem_repeat(seq, 4);
        assert_eq!(unit, 1, "Should detect homopolymer (unit=1)");
        assert_eq!(span, 6, "6 consecutive A's");
    }

    #[test]
    fn test_find_tandem_repeat_dinuc() {
        // GATCACACACCGTT — CA repeat at pos 4-9 (3 copies)
        let seq = b"GATCACACACGTT";
        let (unit, span) = find_tandem_repeat(seq, 4);
        assert_eq!(unit, 2, "Should detect dinucleotide repeat");
        assert_eq!(span, 6, "3 copies of CA = 6bp span");
    }

    #[test]
    fn test_find_tandem_repeat_trinuc() {
        // CAGCAGCAGTTTT — 3 copies of CAG at pos 0
        let seq = b"CAGCAGCAGTTTT";
        let (unit, span) = find_tandem_repeat(seq, 0);
        assert_eq!(unit, 3, "Should detect trinucleotide repeat");
        assert_eq!(span, 9, "3 copies of CAG = 9bp span");
    }

    #[test]
    fn test_find_tandem_repeat_no_repeat() {
        // ATCGATCG — no repeat at pos 2
        let seq = b"ATCGATCG";
        let (unit, span) = find_tandem_repeat(seq, 2);
        assert_eq!(span, 1, "No repeat detected, span should be 1");
        assert_eq!(unit, 1);
    }

    #[test]
    fn test_find_tandem_repeat_at_edge() {
        // AAAAAAA — homopolymer at pos 0
        let seq = b"AAAAAAA";
        let (unit, span) = find_tandem_repeat(seq, 0);
        assert_eq!(unit, 1);
        assert_eq!(span, 7, "Entire sequence is a homopolymer");
    }

    #[test]
    fn test_find_tandem_repeat_middle_of_run() {
        // TTAAAAATT — 5 A's, queried at pos 4 (middle of run)
        let seq = b"TTAAAAATT";
        let (unit, span) = find_tandem_repeat(seq, 4);
        assert_eq!(unit, 1);
        assert_eq!(span, 5, "5 A's even when querying from the middle");
    }

    // -- compute_adaptive_padding (formula) tests --

    #[test]
    fn test_adaptive_padding_default() {
        // Non-repeat: span=1 → adaptive = 1/2+3 = 3 → max(5,3) = 5
        let span = 1;
        let default_pad: i64 = 5;
        let max_pad: i64 = 50;
        let adaptive = (span as i64) / 2 + 3;
        let effective = default_pad.max(adaptive).min(max_pad);
        assert_eq!(effective, 5, "Non-repeat should use default padding");
    }

    #[test]
    fn test_adaptive_padding_homopoly() {
        // Poly-A span=14 → adaptive = 14/2+3 = 10 → max(5,10) = 10
        let span = 14;
        let default_pad: i64 = 5;
        let max_pad: i64 = 50;
        let adaptive = (span as i64) / 2 + 3;
        let effective = default_pad.max(adaptive).min(max_pad);
        assert_eq!(effective, 10, "Homopolymer should increase padding");
    }

    #[test]
    fn test_adaptive_padding_dinuc() {
        // CA-repeat span=20 → adaptive = 20/2+3 = 13 → max(5,13) = 13
        let span = 20;
        let default_pad: i64 = 5;
        let max_pad: i64 = 50;
        let adaptive = (span as i64) / 2 + 3;
        let effective = default_pad.max(adaptive).min(max_pad);
        assert_eq!(effective, 13, "Dinucleotide repeat should increase padding");
    }

    #[test]
    fn test_adaptive_padding_capped() {
        // Very long repeat span=120 → adaptive = 120/2+3 = 63 → min(63,50) = 50
        let span = 120;
        let default_pad: i64 = 5;
        let max_pad: i64 = 50;
        let adaptive = (span as i64) / 2 + 3;
        let effective = default_pad.max(adaptive).min(max_pad);
        assert_eq!(effective, 50, "Should be capped at max_pad");
    }
}
