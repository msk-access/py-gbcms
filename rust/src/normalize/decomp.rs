//! Homopolymer decomposition detection.
//!
//! Detects miscollapsed D(n)+SNV events in complex variants spanning
//! homopolymer runs, and produces a corrected ALT allele for dual-counting.

use log::debug;

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
pub(crate) fn check_homopolymer_decomp(
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
