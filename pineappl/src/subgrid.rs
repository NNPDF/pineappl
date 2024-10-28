//! Module containing the trait `Subgrid` and supporting structs.

use super::empty_subgrid::EmptySubgridV1;
use super::import_subgrid::ImportSubgridV1;
use super::interp_subgrid::InterpSubgridV1;
use enum_dispatch::enum_dispatch;
// use float_cmp::approx_eq;
// use ndarray::Array3;
// use super::evolution::EVOLVE_INFO_TOL_ULPS;
use super::interpolation::Interp;
use serde::{Deserialize, Serialize};

/// Enum which lists all possible `Subgrid` variants possible.
#[enum_dispatch(Subgrid)]
#[derive(Clone, Deserialize, Serialize)]
pub enum SubgridEnum {
    // WARNING: never change the order or content of this enum, only add to the end of it
    /// Subgrid type that supports filling.
    InterpSubgridV1,
    /// Empty subgrid.
    EmptySubgridV1,
    /// TODO
    ImportSubgridV1,
}

/// Structure denoting renormalization and factorization scale values.
#[derive(Debug, Deserialize, Clone, PartialEq, PartialOrd, Serialize)]
pub struct Mu2 {
    /// The (squared) renormalization scale value.
    pub ren: f64,
    /// The (squared) factorization scale value.
    pub fac: f64,
    /// The (squared) fragmentation scale value.
    pub frg: f64,
}

/// Size-related statistics for a subgrid.
#[derive(Debug, Eq, PartialEq)]
pub struct Stats {
    /// Number of possible total entries for a subgrid. This number is the product of the lengths
    /// of the slices returned by [`Subgrid::mu2_grid`], [`Subgrid::x1_grid`] and
    /// [`Subgrid::x2_grid`].
    pub total: usize,
    /// Number of allocated entries for a subgrid. This number is always smaller or equal than
    /// [`Self::total`].
    pub allocated: usize,
    /// Number of allocated zero entries for a subgrid. This number is always smaller or equal than
    /// [`Self::allocated`] and contributes to [`Self::overhead`].
    pub zeros: usize,
    /// The overhead of a [`Subgrid`] is the size of internal data not used to store grid values.
    pub overhead: usize,
    /// This value multiplied with any other member of this struct gives an approximate size in
    /// bytes.
    pub bytes_per_value: usize,
}

/// Trait each subgrid must implement.
#[enum_dispatch]
pub trait Subgrid {
    /// TODO
    fn node_values(&self) -> Vec<Vec<f64>>;

    /// Fill the subgrid with `weight` that is being interpolated with `interps` using the
    /// kinematic information in `ntuple`. The parameter `ntuple` assumes the same ordering given
    /// by `kinematics` in [`Grid::new`] that was used to create the grid.
    fn fill(&mut self, interps: &[Interp], ntuple: &[f64], weight: f64);

    /// Returns true if `fill` was never called for this grid.
    fn is_empty(&self) -> bool;

    /// Merge `other` into this subgrid, possibly transposing the two dimensions given by
    /// `transpose`.
    fn merge(&mut self, other: &SubgridEnum, transpose: Option<(usize, usize)>);

    /// Scale the subgrid by `factor`.
    fn scale(&mut self, factor: f64);

    /// Assume that the convolution functions for indices `a` and `b` for this grid are the same
    /// and use this to optimize the size of the grid.
    fn symmetrize(&mut self, a: usize, b: usize);

    /// Return an iterator over all non-zero elements of the subgrid.
    fn indexed_iter(&self) -> SubgridIndexedIter;

    /// Return statistics for this subgrid.
    fn stats(&self) -> Stats;

    /// TODO
    fn optimize_static_nodes(&mut self);
}

/// Type to iterate over the non-zero contents of a subgrid. The tuple contains the indices of the
/// `mu2_grid`, the `x1_grid` and finally the `x2_grid`.
pub type SubgridIndexedIter<'a> = Box<dyn Iterator<Item = (Vec<usize>, f64)> + 'a>;
