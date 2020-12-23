use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

pub mod bin;
pub mod grid;
pub mod lagrange_subgrid;
pub mod lumi;
pub mod subgrid;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

/// A Python module implemented in Rust.
/// NOTE: this name has to match the one in Cargo.toml 'lib.name'
#[pymodule]
fn pineappl(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_class::<bin::PyBinRemapper>()?;

    Ok(())
}
