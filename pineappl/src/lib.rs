#![warn(clippy::all, clippy::cargo, clippy::nursery, clippy::pedantic)]
#![warn(missing_docs)]
#![allow(clippy::module_name_repetitions)]

//! `PineAPPL` is not an extension of `APPLgrid`.

use pyo3::prelude::*;

mod convert;

pub mod bin;
pub mod grid;
pub mod lagrange_subgrid;
pub mod lumi;
pub mod ntuple_subgrid;
pub mod sparse_array3;
pub mod subgrid;

/// A Python module implemented in Rust.
/// NOTE: this name has to match the one in Cargo.toml 'lib.name'
#[pymodule]
fn pineappl(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<bin::BinRemapper>()?;
    //m.add_class::<grid::Grid>()?;
    //m.add_class::<grid::Order>()?;
    //m.add_class::<lagrange_subgrid::LagrangeSubgridV2>()?;
    //m.add_class::<lumi::LumiEntry>()?;
    //m.add_class::<subgrid::ExtraSubgridParams>()?;
    //m.add_class::<subgrid::SubgridParams>()?;

    Ok(())
}
