//! `PineAPPL` is not an extension of `APPLgrid`.
//!
//! # Overview
//!
//! The main type of this crate is [`Grid`], which represents the interpolation grids that
//! `PineAPPL` implements. Roughly speaking, a `Grid` is a three-dimensional array of [`Subgrid`]
//! objects together with metadata. The three dimensions are
//! 1. (perturbative) orders, represented by the type [`Order`] and accessible by
//!    [`Grid::orders()`],
//! 2. bins, whose limits can be accessed by [`Grid::bin_info()`], and
//! 3. channels, whose definition is returned by [`Grid::channels()`].
//!
//! `Subgrid` is a `trait` and objects that implement it are of the type [`SubgridEnum`]. The
//! latter is an `enum` of different types that are optimized to different scenarios: fast event
//! filling, small storage profile, etc.
//!
//! [`Grid`]: grid::Grid
//! [`Grid::bin_info()`]: grid::Grid::bin_info
//! [`Grid::channels()`]: grid::Grid::channels
//! [`Grid::orders()`]: grid::Grid::orders
//! [`Subgrid`]: subgrid::Subgrid
//! [`SubgridEnum`]: subgrid::SubgridEnum
//! [`Order`]: order::Order
//!
//! ## Metadata
//!
//! Metadata is a collection of key--value pairs, in which both keys and values are `String`
//! objects. In metadata anything a user whishes can be stored. However, there are [special keys],
//! which have meaning to `PineAPPL` and/or its CLI `pineappl`. This metadata enables the CLI to
//! automatically generate plots that are correctly labeled, for instance. For more applications
//! see also the [CLI tutorial].
//!
//! [special keys]: https://nnpdf.github.io/pineappl/docs/metadata.html
//! [CLI tutorial]: https://nnpdf.github.io/pineappl/docs/cli-tutorial.html

mod convert;

pub mod bin;
pub mod boc;
pub mod convolutions;
pub mod empty_subgrid;
pub mod evolution;
pub mod fk_table;
pub mod grid;
pub mod import_only_subgrid;
pub mod lagrange_subgrid;
pub mod ntuple_subgrid;
pub mod packed_array;
pub mod pids;
pub mod sparse_array3;
pub mod subgrid;
