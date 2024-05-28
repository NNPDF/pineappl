//! Module containing the trait `Subgrid` and supporting structs.

use super::empty_subgrid::EmptySubgridV1;
use super::grid::Ntuple;
use super::import_only_subgrid::{ImportOnlySubgridV1, ImportOnlySubgridV2};
use super::lagrange_subgrid::{LagrangeSparseSubgridV1, LagrangeSubgridV1, LagrangeSubgridV2};
use super::ntuple_subgrid::NtupleSubgridV1;
use enum_dispatch::enum_dispatch;
use ndarray::Array3;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;

/// Enum which lists all possible `Subgrid` variants possible.
#[enum_dispatch(Subgrid)]
#[derive(Clone, Deserialize, Serialize)]
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
#[derive(Debug, Deserialize, Clone, PartialEq, PartialOrd, Serialize)]
pub struct Mu2 {
    /// The (squared) renormalization scale value.
    pub ren: f64,
    /// The (squared) factorization scale value.
    pub fac: f64,
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
    /// Return a slice of [`Mu2`] values corresponding to the (squared) renormalization and
    /// factorization values of the grid. If the subgrid does not use a grid, this method should
    /// return an empty slice.
    fn mu2_grid(&self) -> Cow<[Mu2]>;

    /// Return a slice of values of `x1`. If the subgrid does not use a grid, this method should
    /// return an empty slice.
    fn x1_grid(&self) -> Cow<[f64]>;

    /// Return a slice of values of `x2`. If the subgrid does not use a grid, this method should
    /// return an empty slice.
    fn x2_grid(&self) -> Cow<[f64]>;

    /// Convolute the subgrid with a luminosity function, which takes indices as arguments that
    /// correspond to the entries given in the slices `x1`, `x2` and `mu2`.
    fn convolve(
        &self,
        x1: &[f64],
        x2: &[f64],
        mu2: &[Mu2],
        lumi: &mut dyn FnMut(usize, usize, usize) -> f64,
    ) -> f64;

    /// Fills the subgrid with `weight` for the parton momentum fractions `x1` and `x2`, and the
    /// scale `q2`. Filling is currently only support where both renormalization and factorization
    /// scale have the same value.
    fn fill(&mut self, ntuple: &Ntuple<f64>);

    /// Returns true if `fill` was never called for this grid.
    fn is_empty(&self) -> bool;

    /// Merges `other` into this subgrid.
    fn merge(&mut self, other: &mut SubgridEnum, transpose: bool);

    /// Scale the subgrid by `factor`.
    fn scale(&mut self, factor: f64);

    /// Assumes that the initial states for this grid are the same and uses this to optimize the
    /// grid by getting rid of almost half of the entries.
    fn symmetrize(&mut self);

    /// Returns an empty copy of the current subgrid.
    fn clone_empty(&self) -> SubgridEnum;

    /// Return an iterator over all non-zero elements of the subgrid.
    fn indexed_iter(&self) -> SubgridIndexedIter;

    /// Return statistics for this subgrid.
    fn stats(&self) -> Stats;

    /// Return the static (single) scale, if this subgrid has one.
    fn static_scale(&self) -> Option<Mu2>;
}

// this is needed in the Python interface
impl From<&SubgridEnum> for Array3<f64> {
    fn from(subgrid: &SubgridEnum) -> Self {
        let mut result = Self::zeros((
            subgrid.mu2_grid().len(),
            subgrid.x1_grid().len(),
            subgrid.x2_grid().len(),
        ));

        for ((imu2, ix1, ix2), value) in subgrid.indexed_iter() {
            result[[imu2, ix1, ix2]] = value;
        }

        result
    }
}

/// Type to iterate over the non-zero contents of a subgrid. The tuple contains the indices of the
/// `mu2_grid`, the `x1_grid` and finally the `x2_grid`.
pub type SubgridIndexedIter<'a> = Box<dyn Iterator<Item = ((usize, usize, usize), f64)> + 'a>;

/// Subgrid creation parameters for subgrids that perform interpolation.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct SubgridParams {
    q2_bins: usize,
    q2_max: f64,
    q2_min: f64,
    q2_order: usize,
    reweight: bool,
    x_bins: usize,
    x_max: f64,
    x_min: f64,
    x_order: usize,
}

impl Default for SubgridParams {
    fn default() -> Self {
        Self {
            q2_bins: 40,
            q2_max: 1e8,
            q2_min: 1e2,
            q2_order: 3,
            reweight: true,
            x_bins: 50,
            x_max: 1.0,
            x_min: 2e-7,
            x_order: 3,
        }
    }
}

impl SubgridParams {
    /// Returns the number of bins for the $Q^2$ axis.
    #[must_use]
    pub const fn q2_bins(&self) -> usize {
        self.q2_bins
    }

    /// Returns the upper limit of the $Q^2$ axis.
    #[must_use]
    pub const fn q2_max(&self) -> f64 {
        self.q2_max
    }

    /// Returns the lower limit of the $Q^2$ axis.
    #[must_use]
    pub const fn q2_min(&self) -> f64 {
        self.q2_min
    }

    /// Returns the interpolation order for the $Q^2$ axis.
    #[must_use]
    pub const fn q2_order(&self) -> usize {
        self.q2_order
    }

    /// Returns whether reweighting is enabled or not.
    #[must_use]
    pub const fn reweight(&self) -> bool {
        self.reweight
    }

    /// Sets the number of bins for the $Q^2$ axis.
    pub fn set_q2_bins(&mut self, q2_bins: usize) {
        self.q2_bins = q2_bins;
    }

    /// Sets the upper limit of the $Q^2$ axis.
    pub fn set_q2_max(&mut self, q2_max: f64) {
        self.q2_max = q2_max;
    }

    /// Sets the lower limit of the $Q^2$ axis.
    pub fn set_q2_min(&mut self, q2_min: f64) {
        self.q2_min = q2_min;
    }

    /// Sets the interpolation order for the $Q^2$ axis.
    pub fn set_q2_order(&mut self, q2_order: usize) {
        self.q2_order = q2_order;
    }

    /// Sets the reweighting parameter.
    pub fn set_reweight(&mut self, reweight: bool) {
        self.reweight = reweight;
    }

    /// Sets the number of bins for the $x$ axes.
    pub fn set_x_bins(&mut self, x_bins: usize) {
        self.x_bins = x_bins;
    }

    /// Sets the upper limit of the $x$ axes.
    pub fn set_x_max(&mut self, x_max: f64) {
        self.x_max = x_max;
    }

    /// Sets the lower limit of the $x$ axes.
    pub fn set_x_min(&mut self, x_min: f64) {
        self.x_min = x_min;
    }

    /// Sets the interpolation order for the $x$ axes.
    pub fn set_x_order(&mut self, x_order: usize) {
        self.x_order = x_order;
    }

    /// Returns the number of bins for the $x$ axes.
    #[must_use]
    pub const fn x_bins(&self) -> usize {
        self.x_bins
    }

    /// Returns the upper limit of the $x$ axes.
    #[must_use]
    pub const fn x_max(&self) -> f64 {
        self.x_max
    }

    /// Returns the lower limit of the $x$ axes.
    #[must_use]
    pub const fn x_min(&self) -> f64 {
        self.x_min
    }

    /// Returns the interpolation order for the $x$ axes.
    #[must_use]
    pub const fn x_order(&self) -> usize {
        self.x_order
    }
}

/// Extra grid creation parameters when the limits for `x1` and `x2` are different.
#[derive(Deserialize, Serialize)]
pub struct ExtraSubgridParams {
    reweight2: bool,
    x2_bins: usize,
    x2_max: f64,
    x2_min: f64,
    x2_order: usize,
}

impl Default for ExtraSubgridParams {
    fn default() -> Self {
        Self {
            reweight2: true,
            x2_bins: 50,
            x2_max: 1.0,
            x2_min: 2e-7,
            x2_order: 3,
        }
    }
}

impl From<&SubgridParams> for ExtraSubgridParams {
    fn from(subgrid_params: &SubgridParams) -> Self {
        Self {
            reweight2: subgrid_params.reweight(),
            x2_bins: subgrid_params.x_bins(),
            x2_max: subgrid_params.x_max(),
            x2_min: subgrid_params.x_min(),
            x2_order: subgrid_params.x_order(),
        }
    }
}

impl ExtraSubgridParams {
    /// Returns whether reweighting is enabled for the `x2` axis or not.
    #[must_use]
    pub const fn reweight2(&self) -> bool {
        self.reweight2
    }

    /// Sets the reweighting parameter for the `x2` axis.
    pub fn set_reweight2(&mut self, reweight2: bool) {
        self.reweight2 = reweight2;
    }

    /// Sets the number of bins for the `x2` axes.
    pub fn set_x2_bins(&mut self, x_bins: usize) {
        self.x2_bins = x_bins;
    }

    /// Sets the upper limit of the `x2` axes.
    pub fn set_x2_max(&mut self, x_max: f64) {
        self.x2_max = x_max;
    }

    /// Sets the lower limit of the `x2` axes.
    pub fn set_x2_min(&mut self, x_min: f64) {
        self.x2_min = x_min;
    }

    /// Sets the interpolation order for the `x2` axes.
    pub fn set_x2_order(&mut self, x_order: usize) {
        self.x2_order = x_order;
    }

    /// Returns the number of bins for the `x2` axes.
    #[must_use]
    pub const fn x2_bins(&self) -> usize {
        self.x2_bins
    }

    /// Returns the upper limit of the `x2` axes.
    #[must_use]
    pub const fn x2_max(&self) -> f64 {
        self.x2_max
    }

    /// Returns the lower limit of the `x2` axes.
    #[must_use]
    pub const fn x2_min(&self) -> f64 {
        self.x2_min
    }

    /// Returns the interpolation order for the `x2` axes.
    #[must_use]
    pub const fn x2_order(&self) -> usize {
        self.x2_order
    }
}
