//! Module containing the trait `Subgrid` and supporting structs.

use super::empty_subgrid::EmptySubgridV1;
use super::lagrange_subgrid::LagrangeSubgridV2;
use super::packed_subgrid::PackedQ1X2SubgridV1;
use enum_dispatch::enum_dispatch;
// use float_cmp::approx_eq;
// use ndarray::Array3;
// use super::evolution::EVOLVE_INFO_TOL_ULPS;
use super::interpolation::Interp;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;

/// TODO
#[derive(Clone)]
pub enum NodeValues {
    /// TODO
    UseFromGrid,
    /// TODO
    UseThese(Vec<f64>),
}

impl NodeValues {
    /// TODO
    pub fn extend(&mut self, other: &Self) {
        match (self, other) {
            (NodeValues::UseFromGrid, NodeValues::UseFromGrid) => (),
            (NodeValues::UseThese(a), NodeValues::UseThese(b)) => {
                a.extend_from_slice(b);
                a.sort_by(|lhs, rhs| lhs.partial_cmp(rhs).unwrap());
                // TODO: use some tolerance
                a.dedup();
            }
            _ => unimplemented!(),
        }
    }

    /// TODO
    pub fn len(&self) -> usize {
        match self {
            NodeValues::UseFromGrid => unimplemented!(),
            NodeValues::UseThese(a) => a.len(),
        }
    }

    /// TODO
    pub fn find(&self, value: f64) -> Option<usize> {
        match self {
            NodeValues::UseFromGrid => unimplemented!(),
            NodeValues::UseThese(a) => a.iter().position(|&x|
                // approx_eq!(f64, x, value, ulps = EVOLVE_INFO_TOL_ULPS)
                x == value),
        }
    }

    /// TODO
    pub fn get(&self, index: usize) -> f64 {
        match self {
            NodeValues::UseFromGrid => unimplemented!(),
            NodeValues::UseThese(a) => a[index],
        }
    }
}

impl PartialEq for NodeValues {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (NodeValues::UseFromGrid, NodeValues::UseFromGrid) => true,
            (NodeValues::UseThese(a), NodeValues::UseThese(b)) => a.iter().zip(b).all(
                // TODO: use some tolerance
                |(&a, &b)| a == b,
                //approx_eq!(f64, a, b, ulps = EVOLVE_INFO_TOL_ULPS)
            ),
            // TODO: the remaining cases could still be the same, but we don't know the values from `UseFromGrid`.
            _ => false,
        }
    }
}

/// Enum which lists all possible `Subgrid` variants possible.
#[enum_dispatch(Subgrid)]
#[derive(Clone, Deserialize, Serialize)]
pub enum SubgridEnum {
    // WARNING: never change the order or content of this enum, only add to the end of it
    /// Lagrange-interpolation subgrid with possibly different x1 and x2 bins.
    LagrangeSubgridV2,
    /// Empty subgrid.
    EmptySubgridV1,
    /// TODO
    PackedQ1X2SubgridV1,
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

    /// TODO
    fn node_values(&self) -> Vec<NodeValues>;

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

    /// Return the static (single) scale, if this subgrid has one.
    fn static_scale(&self) -> Option<Mu2>;
}

// // this is needed in the Python interface
// impl From<&SubgridEnum> for Array3<f64> {
//     fn from(subgrid: &SubgridEnum) -> Self {
//         let mut result = Self::zeros((
//             subgrid.mu2_grid().len(),
//             subgrid.x1_grid().len(),
//             subgrid.x2_grid().len(),
//         ));
//
//         for ((imu2, ix1, ix2), value) in subgrid.indexed_iter() {
//             result[[imu2, ix1, ix2]] = value;
//         }
//
//         result
//     }
// }

/// Type to iterate over the non-zero contents of a subgrid. The tuple contains the indices of the
/// `mu2_grid`, the `x1_grid` and finally the `x2_grid`.
pub type SubgridIndexedIter<'a> = Box<dyn Iterator<Item = (Vec<usize>, f64)> + 'a>;
