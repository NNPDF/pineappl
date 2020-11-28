//! Module containing the trait `Subgrid` and supporting structs.

use super::grid::Ntuple;
use super::lagrange_subgrid::{LagrangeSparseSubgridV1, LagrangeSubgridV1};
use super::ntuple_subgrid::NtupleSubgridV1;
use either::Either;
use enum_dispatch::enum_dispatch;
use serde::{Deserialize, Serialize};

/// Enum which lists all possible `Subgrid` variants possible.
#[enum_dispatch(Subgrid)]
#[derive(Deserialize, Serialize)]
pub enum SubgridEnum {
    // WARNING: never change the order or content of this enum, only add to the end of it
    /// Lagrange-interpolation subgrid.
    LagrangeSubgridV1,
    /// N-tuple subgrid.
    NtupleSubgridV1,
    /// Lagrange-interpolation subgrid.
    LagrangeSparseSubgridV1,
}

/// Trait each subgrid must implement.
#[enum_dispatch]
pub trait Subgrid {
    /// Return a `Vec` of values of `q2`. If the subgrid does not use a grid, this method should
    /// return an empty `Vec`.
    fn grid_q2(&self) -> Vec<f64>;

    /// Return a `Vec` of values of `x`. If the subgrid does not use a grid, this method should
    /// return an empty `Vec`.
    fn grid_x(&self) -> Vec<f64>;

    /// Convolute the subgrid with a luminosity function, which is either takes indices as
    /// arguments, in which case the `x` and `q2` values can be read from the given slices, or
    /// takes the usual values `x1`, `x2`, and `q2`. If the method `grid_x` returns a non-empty
    /// vector, this method must use the indexed luminosity function.
    fn convolute(
        &self,
        x: &[f64],
        q2: &[f64],
        lumi: Either<&dyn Fn(usize, usize, usize) -> f64, &dyn Fn(f64, f64, f64) -> f64>,
    ) -> f64;

    /// Fills the subgrid with `weight` for the parton momentum fractions `x1` and `x2`, and the
    /// scale `q2`.
    fn fill(&mut self, ntuple: &Ntuple<f64>);

    /// Returns true if `fill` was never called for this grid.
    fn is_empty(&self) -> bool;

    /// Merges `other` into this subgrid.
    fn merge(&mut self, other: &mut SubgridEnum);

    /// Scale the subgrid by `factor`.
    fn scale(&mut self, factor: f64);

    // TODO: the following should be a Range

    /// Returns the half-open interval of indices of filled q2 slices.
    fn q2_slice(&self) -> (usize, usize);

    // TODO: rename the function to export_q2_slice

    /// Fill the q2-slice with index `q2_slice` into `grid`.
    fn fill_q2_slice(&self, q2_slice: usize, grid: &mut [f64]);

    // TODO: rename the function to import_applgrid_f2_q2_slice

    /// Writes into subgrid.
    fn write_q2_slice(&mut self, q2_slice: usize, grid: &[f64]);
}

/// Subgrid creation parameters for subgrids that perform interpolation.
#[derive(Deserialize, Serialize)]
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
