//! Module containing the trait `Subgrid` and supporting structs.

use super::empty_subgrid::EmptySubgridV1;
use super::grid::Ntuple;
use super::import_only_subgrid::ImportOnlySubgridV1;
use super::lagrange_subgrid::{LagrangeSparseSubgridV1, LagrangeSubgridV1, LagrangeSubgridV2};
use super::ntuple_subgrid::NtupleSubgridV1;
use either::Either;
use enum_dispatch::enum_dispatch;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::ops::Range;

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
}

/// Trait each subgrid must implement.
#[enum_dispatch]
pub trait Subgrid {
    /// Return a `Vec` of values of `q2`. If the subgrid does not use a grid, this method should
    /// return an empty `Vec`.
    fn q2_grid(&self) -> Cow<[f64]>;

    /// Return a `Vec` of values of `x1`. If the subgrid does not use a grid, this method should
    /// return an empty `Vec`.
    fn x1_grid(&self) -> Cow<[f64]>;

    /// Return a `Vec` of values of `x2`. If the subgrid does not use a grid, this method should
    /// return an empty `Vec`.
    fn x2_grid(&self) -> Cow<[f64]>;

    /// Convolute the subgrid with a luminosity function, which either takes indices as arguments,
    /// in which case the `x1`, `x2` and `q2` values can be read from the given slices, or takes
    /// the usual values `x1`, `x2`, and `q2`. If the method `x1_grid` and `x2_grid` return a
    /// non-empty vector, this method must use the indexed luminosity function.
    fn convolute(
        &self,
        x1: &[f64],
        x2: &[f64],
        q2: &[f64],
        lumi: Either<&dyn Fn(usize, usize, usize) -> f64, &dyn Fn(f64, f64, f64) -> f64>,
    ) -> f64;

    /// Fills the subgrid with `weight` for the parton momentum fractions `x1` and `x2`, and the
    /// scale `q2`.
    fn fill(&mut self, ntuple: &Ntuple<f64>);

    /// Returns true if `fill` was never called for this grid.
    fn is_empty(&self) -> bool;

    /// Merges `other` into this subgrid.
    fn merge(&mut self, other: &mut SubgridEnum, transpose: bool);

    /// Scale the subgrid by `factor`.
    fn scale(&mut self, factor: f64);

    /// Returns the half-open interval of indices of filled q2 slices.
    fn q2_slice(&self) -> Range<usize>;

    // TODO: rename the function to export_q2_slice

    /// Fill the q2-slice with index `q2_slice` into `grid`.
    fn fill_q2_slice(&self, q2_slice: usize, grid: &mut [f64]);

    /// Assumes that the initial states for this grid are the same and uses this to optimize the
    /// grid by getting rid of almost half of the entries.
    fn symmetrize(&mut self);

    /// Returns an empty copy of the current subgrid.
    fn clone_empty(&self) -> SubgridEnum;
}

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
    /// Returns whether reweighting is enabled for the $x_2$ axis or not.
    #[must_use]
    pub const fn reweight2(&self) -> bool {
        self.reweight2
    }

    /// Sets the reweighting parameter for the $x_2$ axis.
    pub fn set_reweight2(&mut self, reweight2: bool) {
        self.reweight2 = reweight2;
    }

    /// Sets the number of bins for the $x_2$ axes.
    pub fn set_x2_bins(&mut self, x_bins: usize) {
        self.x2_bins = x_bins;
    }

    /// Sets the upper limit of the $x_2$ axes.
    pub fn set_x2_max(&mut self, x_max: f64) {
        self.x2_max = x_max;
    }

    /// Sets the lower limit of the $x_2$ axes.
    pub fn set_x2_min(&mut self, x_min: f64) {
        self.x2_min = x_min;
    }

    /// Sets the interpolation order for the $x_2$ axes.
    pub fn set_x2_order(&mut self, x_order: usize) {
        self.x2_order = x_order;
    }

    /// Returns the number of bins for the $x_2$ axes.
    #[must_use]
    pub const fn x2_bins(&self) -> usize {
        self.x2_bins
    }

    /// Returns the upper limit of the $x_2$ axes.
    #[must_use]
    pub const fn x2_max(&self) -> f64 {
        self.x2_max
    }

    /// Returns the lower limit of the $x_2$ axes.
    #[must_use]
    pub const fn x2_min(&self) -> f64 {
        self.x2_min
    }

    /// Returns the interpolation order for the $x_2$ axes.
    #[must_use]
    pub const fn x2_order(&self) -> usize {
        self.x2_order
    }
}
