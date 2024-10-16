//! Generate PyO3 interface.

#![allow(unsafe_op_in_unsafe_fn)]

use pyo3::prelude::*;

pub mod bin;
pub mod boc;
pub mod convolutions;
pub mod evolution;
pub mod fk_table;
pub mod grid;
pub mod interpolation;
pub mod pids;
pub mod subgrid;

/// PyO3 Python module that contains all exposed classes from Rust.
#[pymodule]
fn vector_length_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    bin::register(m)?;
    boc::register(m)?;
    evolution::register(m)?;
    convolutions::register(m)?;
    grid::register(m)?;
    interpolation::register(m)?;
    pids::register(m)?;
    fk_table::register(m)?;
    subgrid::register(m)?;
    m.add("version", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}
