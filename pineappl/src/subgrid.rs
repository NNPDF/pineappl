//! Module containing the trait `Subgrid` and supporting structs.

use super::empty_subgrid::EmptySubgridV1;
use super::import_subgrid::ImportSubgridV1;
use super::interp_subgrid::InterpSubgridV1;
use enum_dispatch::enum_dispatch;
// use float_cmp::approx_eq;
// use ndarray::Array3;
// use super::evolution::EVOLVE_INFO_TOL_ULPS;
use super::interpolation::Interp;
use float_cmp::approx_eq;
use serde::{Deserialize, Serialize};

/// TODO
#[must_use]
pub fn node_value_eq(lhs: f64, rhs: f64) -> bool {
    approx_eq!(f64, lhs, rhs, ulps = 4096)
}

/// TODO
#[must_use]
pub fn node_value_eq_ref_mut(lhs: &mut f64, rhs: &mut f64) -> bool {
    node_value_eq(*lhs, *rhs)
}

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
    fn optimize_nodes(&mut self);
}

/// Type to iterate over the non-zero contents of a subgrid. The tuple contains the indices of the
/// `mu2_grid`, the `x1_grid` and finally the `x2_grid`.
pub type SubgridIndexedIter<'a> = Box<dyn Iterator<Item = (Vec<usize>, f64)> + 'a>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::interpolation::{Interp, InterpMeth, Map, ReweightMeth};

    #[test]
    fn check_old_x_values() {
        let new_50_x = Interp::new(
            2e-7,
            1.0,
            50,
            3,
            ReweightMeth::ApplGridX,
            Map::ApplGridF2,
            InterpMeth::Lagrange,
        )
        .node_values();
        let old_50_x = [
            1.0,
            9.309440808717544e-1,
            8.627839323906108e-1,
            7.956242522922756e-1,
            7.295868442414312e-1,
            6.648139482473823e-1,
            6.01472197967335e-1,
            5.397572337880445e-1,
            4.798989029610255e-1,
            4.221667753589648e-1,
            3.668753186482242e-1,
            3.1438740076927585e-1,
            2.651137041582823e-1,
            2.195041265003886e-1,
            1.7802566042569432e-1,
            1.4112080644440345e-1,
            1.0914375746330703e-1,
            8.228122126204893e-2,
            6.0480028754447364e-2,
            4.341491741702269e-2,
            3.0521584007828916e-2,
            2.108918668378717e-2,
            1.4375068581090129e-2,
            9.699159574043399e-3,
            6.496206194633799e-3,
            4.328500638820811e-3,
            2.8738675812817515e-3,
            1.9034634022867384e-3,
            1.2586797144272762e-3,
            8.314068836488144e-4,
            5.487795323670796e-4,
            3.6205449638139736e-4,
            2.3878782918561914e-4,
            1.5745605600841445e-4,
            1.0381172986576898e-4,
            6.843744918967897e-5,
            4.511438394964044e-5,
            2.97384953722449e-5,
            1.9602505002391748e-5,
            1.292101569074731e-5,
            8.516806677573355e-6,
            5.613757716930151e-6,
            3.7002272069854957e-6,
            2.438943292891682e-6,
            1.607585498470808e-6,
            1.0596094959101024e-6,
            6.984208530700364e-7,
            4.6035014748963906e-7,
            3.034304765867952e-7,
            1.9999999999999954e-7,
        ];

        // check that the old x-grid values are 'equal' to the new ones
        for (old, new) in old_50_x.into_iter().zip(new_50_x) {
            assert!(node_value_eq(old, new), "{old} {new}");
        }
    }
}
