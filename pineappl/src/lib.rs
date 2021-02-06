#![warn(clippy::all, clippy::cargo, clippy::nursery, clippy::pedantic)]
#![warn(missing_docs)]
#![allow(clippy::module_name_repetitions)]

//! `PineAPPL` is not an extension of `APPLgrid`.

mod convert;

pub mod bin;
pub mod grid;
pub mod lagrange_subgrid;
pub mod lumi;
pub mod ntuple_subgrid;
pub mod read_only_sparse_subgrid;
pub mod sparse_array3;
pub mod subgrid;
