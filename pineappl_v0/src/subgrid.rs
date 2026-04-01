//! Module containing the trait `Subgrid` and supporting structs.

use super::empty_subgrid::EmptySubgridV1;
use super::import_only_subgrid::{ImportOnlySubgridV1, ImportOnlySubgridV2};
use super::lagrange_subgrid::{LagrangeSparseSubgridV1, LagrangeSubgridV1, LagrangeSubgridV2};
use super::ntuple_subgrid::NtupleSubgridV1;
use enum_dispatch::enum_dispatch;
use serde::Deserialize;
use std::borrow::Cow;

/// Enum which lists all possible `Subgrid` variants possible.
#[enum_dispatch(Subgrid)]
#[derive(Deserialize)]
pub enum SubgridEnum {
    // WARNING: never change the order or content of this enum, only add to the end of it
    /// Lagrange-interpolation subgrid.
    LagrangeSubgridV1,
    /// N-tuple subgrid.
    NtupleSubgridV1,
    /// Lagrange-interpolation subgrid.
    LagrangeSparseSubgridV1,
    /// Lagrange-interpolation subgrid with possibly different x1 and x2 bins.
    LagrangeSubgridV2,
    /// Import-only sparse subgrid with possibly different x1 and x2 bins.
    ImportOnlySubgridV1,
    /// Empty subgrid.
    EmptySubgridV1,
    /// Same as [`ImportOnlySubgridV1`], but with support for different renormalization and
    /// factorization scales choices.
    ImportOnlySubgridV2,
}

/// Structure denoting renormalization and factorization scale values.
#[derive(Deserialize, Clone)]
pub struct Mu2 {
    /// The (squared) renormalization scale value.
    pub ren: f64,
    /// The (squared) factorization scale value.
    pub fac: f64,
}

/// Trait each subgrid must implement.
#[enum_dispatch]
pub trait Subgrid {
    /// Return a slice of [`Mu2`] values corresponding to the (squared) renormalization and
    /// factorization values of the grid. If the subgrid does not use a grid, this method should
    /// return an empty slice.
    fn mu2_grid(&self) -> Cow<'_, [Mu2]>;

    /// Return a slice of values of `x1`. If the subgrid does not use a grid, this method should
    /// return an empty slice.
    fn x1_grid(&self) -> Cow<'_, [f64]>;

    /// Return a slice of values of `x2`. If the subgrid does not use a grid, this method should
    /// return an empty slice.
    fn x2_grid(&self) -> Cow<'_, [f64]>;

    /// Returns true if `fill` was never called for this grid.
    fn is_empty(&self) -> bool;

    /// Return an iterator over all non-zero elements of the subgrid.
    fn indexed_iter(&self) -> SubgridIndexedIter<'_>;
}

/// Type to iterate over the non-zero contents of a subgrid. The tuple contains the indices of the
/// `mu2_grid`, the `x1_grid` and finally the `x2_grid`.
pub type SubgridIndexedIter<'a> = Box<dyn Iterator<Item = ((usize, usize, usize), f64)> + 'a>;

/// Subgrid creation parameters for subgrids that perform interpolation.
#[derive(Deserialize)]
pub struct SubgridParams {
    _q2_bins: usize,
    _q2_max: f64,
    _q2_min: f64,
    _q2_order: usize,
    _reweight: bool,
    _x_bins: usize,
    _x_max: f64,
    _x_min: f64,
    _x_order: usize,
}
