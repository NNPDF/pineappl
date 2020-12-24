use pyo3::prelude::*;

pub mod bin;
pub mod grid;
pub mod lagrange_subgrid;
pub mod lumi;
pub mod subgrid;

/// A Python module implemented in Rust.
/// NOTE: this name has to match the one in Cargo.toml 'lib.name'
#[pymodule]
fn pineappl(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<bin::PyBinRemapper>()?;
    m.add_class::<grid::PyGrid>()?;
    m.add_class::<grid::PyOrder>()?;
    m.add_class::<lagrange_subgrid::PyLagrangeSubgridV2>()?;
    m.add_class::<subgrid::PyExtraSubgridParams>()?;
    m.add_class::<subgrid::PySubgridParams>()?;

    Ok(())
}
