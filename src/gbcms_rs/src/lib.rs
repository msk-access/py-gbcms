#![allow(unsafe_op_in_unsafe_fn)]
use pyo3::prelude::*;

mod counting;
mod stats;
mod types;

/// A Python module implemented in Rust.
#[pymodule]
fn gbcms_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    pyo3_log::init();
    m.add_function(wrap_pyfunction!(counting::count_bam, m)?)?;
    m.add_class::<types::Variant>()?;
    m.add_class::<types::BaseCounts>()?;
    Ok(())
}
