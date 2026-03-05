//! Base counting engine for variant classification.
//!
//! This module orchestrates read-level allele classification for a list of
//! variants in a BAM file. It coordinates the variant-type-specific checkers,
//! fragment-level evidence tracking, and alignment backends.
//!
//! ## Submodules
//!
//! - [`engine`] — Core orchestration: `count_bam`, `count_single_variant`, `check_allele_with_qual`
//! - [`fragment`] — Fragment evidence tracking and QNAME hashing
//! - [`alignment`] — Smith-Waterman alignment backend (Phase 3)
//! - [`pairhmm`] — PairHMM alignment backend (probabilistic Phase 3 alternative)
//! - [`variant_checks`] — Per-variant-type classification (SNP, MNP, Ins, Del, Complex)
//! - [`utils`] — Shared utility functions (position lookup, quality, masked comparison)
//! - [`mfsd`] — Mutant Fragment Size Distribution statistics (KS test, LLR, mean)

mod engine;
mod fragment;
mod alignment;
pub mod pairhmm;
mod variant_checks;
mod utils;
pub(crate) mod mfsd;

// Re-export the PyO3 entry point so lib.rs can call counting::count_bam
pub use engine::count_bam;

// Re-export AlignmentBackend for sibling modules (used by pairhmm tests)
pub use engine::AlignmentBackend;
