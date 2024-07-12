// this is needed for PyO3 to work
#![allow(unsafe_op_in_unsafe_fn)]

use pyo3::prelude::*;

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
    m.add_class::<bin::PyBinRemapper>()?;
    m.add_class::<evolution::PyEvolveInfo>()?;
    m.add_class::<grid::PyGrid>()?;
    m.add_class::<grid::PyOrder>()?;
    m.add_class::<grid::PyOperatorSliceInfo>()?;
    m.add_class::<grid::PyPidBasis>()?;
    m.add_class::<lumi::PyLumiEntry>()?;
    m.add_class::<import_only_subgrid::PyImportOnlySubgridV1>()?;
    m.add_class::<import_only_subgrid::PyImportOnlySubgridV2>()?;
    m.add_class::<fk_table::PyFkTable>()?;
    m.add_class::<fk_table::PyFkAssumptions>()?;
    m.add_class::<subgrid::PySubgridEnum>()?;
    m.add_class::<subgrid::PySubgridParams>()?;
    m.add_class::<subgrid::PyMu2>()?;
    m.add("version", env!("CARGO_PKG_VERSION"))?;

    Ok(())
}
