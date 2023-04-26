//! `PineAPPL` is not an extension of `APPLgrid`.

#![warn(clippy::all, clippy::cargo, clippy::nursery, clippy::pedantic)]
#![warn(missing_docs)]
#![allow(clippy::module_name_repetitions)]
#![allow(clippy::similar_names)]

mod convert;

pub mod bin;
pub mod empty_subgrid;
pub mod evolution;
pub mod fk_table;
pub mod grid;
pub mod import_only_subgrid;
pub mod lagrange_subgrid;
pub mod lumi;
pub mod ntuple_subgrid;
pub mod pids;
pub mod sparse_array3;
pub mod subgrid;
