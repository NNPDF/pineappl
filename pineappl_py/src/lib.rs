use pyo3::prelude::*;

pub mod bin;
pub mod fk_table;
pub mod grid;
pub mod import_only_subgrid;
pub mod lumi;
pub mod subgrid;

/// A Python module implemented in Rust.
/// NOTE: this name has to match the one in Cargo.toml 'lib.name'
#[pymodule]
fn pineappl(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<bin::PyBinRemapper>()?;
    m.add_class::<grid::PyGrid>()?;
    m.add_class::<grid::PyOrder>()?;
    m.add_class::<lumi::PyLumiEntry>()?;
    m.add_class::<import_only_subgrid::PyImportOnlySubgridV1>()?;
    m.add_class::<fk_table::PyFkTable>()?;
    m.add_class::<subgrid::PySubgridEnum>()?;
    m.add_class::<subgrid::PySubgridParams>()?;

    Ok(())
}
