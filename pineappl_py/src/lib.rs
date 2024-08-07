// this is needed for PyO3 to work
#![allow(unsafe_op_in_unsafe_fn)]

use pyo3::prelude::*;
use pyo3::wrap_pymodule;

pub mod bin;
pub mod evolution;
pub mod fk_table;
pub mod grid;
pub mod import_only_subgrid;
pub mod lumi;
pub mod subgrid;

/// PyO3 Python module that contains all exposed classes from Rust.
///
/// NOTE: this name has to match the one in Cargo.toml 'lib.name'
#[pymodule]
fn pineappl(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(bin::bin))?;
    m.add_wrapped(wrap_pymodule!(grid::grid))?;
    m.add_wrapped(wrap_pymodule!(import_only_subgrid::import_only_subgrid))?;
    m.add_wrapped(wrap_pymodule!(evolution::evolution))?;
    m.add_wrapped(wrap_pymodule!(lumi::lumi))?;
    m.add_wrapped(wrap_pymodule!(fk_table::fk_table))?;
    m.add_wrapped(wrap_pymodule!(subgrid::subgrid))?;
    m.add("version", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}
