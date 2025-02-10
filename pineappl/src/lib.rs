//! `PineAPPL` is not an extension of `APPLgrid`.
//!
//! # Overview
//!
//! The main type of this crate is [`Grid`], which represents the interpolation grids that
//! `PineAPPL` implements. Roughly speaking, a `Grid` is a three-dimensional array of [`Subgrid`]
//! objects together with metadata. The three dimensions are
//! 1. bins, whose limits can be accessed by [`Grid::bwfl()`], and
//! 2. (perturbative) orders, represented by the type [`Order`] and accessible by
//!    [`Grid::orders()`],
//! 3. channels, whose definition is returned by [`Grid::channels()`].
//!
//! `Subgrid` is a `trait` and objects that implement it are of the type [`SubgridEnum`]. The
//! latter is an `enum` of different types that are optimized to different scenarios: fast event
//! filling, small storage profile, etc.
//!
//! [`Grid`]: grid::Grid
//! [`Grid::bwfl()`]: grid::Grid::bwfl
//! [`Grid::channels()`]: grid::Grid::channels
//! [`Grid::orders()`]: grid::Grid::orders
//! [`Subgrid`]: subgrid::Subgrid
//! [`SubgridEnum`]: subgrid::SubgridEnum
//! [`Order`]: boc::Order
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
mod v0;

pub mod boc;
pub mod convolutions;
pub mod error;
pub mod evolution;
pub mod fk_table;
pub mod grid;
pub mod interpolation;
pub mod packed_array;
pub mod pids;
pub mod reference;
pub mod subgrid;
