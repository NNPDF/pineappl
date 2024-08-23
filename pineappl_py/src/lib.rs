//! Generate PyO3 interface.

#![allow(unsafe_op_in_unsafe_fn)]

use pyo3::prelude::*;

pub mod bin;
pub mod boc;
pub mod evolution;
pub mod fk_table;
pub mod grid;
pub mod import_only_subgrid;
pub mod pids;
pub mod subgrid;

/// PyO3 Python module that contains all exposed classes from Rust.
#[pymodule]
fn pineappl(m: &Bound<'_, PyModule>) -> PyResult<()> {
    bin::register(m)?;
    boc::register(m)?;
    grid::register(m)?;
    import_only_subgrid::register(m)?;
    evolution::register(m)?;
    fk_table::register(m)?;
    pids::register(m)?;
    subgrid::register(m)?;
    m.add("version", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}
