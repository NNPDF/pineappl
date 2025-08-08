//! Module containing the trait `Subgrid` and supporting structs.

use super::interpolation::{self, Interp};
use super::packed_array::PackedArray;
use enum_dispatch::enum_dispatch;
use float_cmp::approx_eq;
use itertools::izip;
use serde::{Deserialize, Serialize};
use std::{iter, mem};

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

/// A subgrid type that is always empty.
#[derive(Clone, Default, Deserialize, Serialize)]
pub struct EmptySubgridV1;

impl Subgrid for EmptySubgridV1 {
    fn fill(&mut self, _: &[Interp], _: &[f64], _: f64) {
        panic!("EmptySubgridV1 doesn't support the fill operation");
    }

    fn node_values(&self) -> Vec<Vec<f64>> {
        Vec::new()
    }

    fn shape(&self) -> &[usize] {
        panic!("EmptySubgridV1 doesn't have a shape");
    }

    fn is_empty(&self) -> bool {
        true
    }

    fn merge_impl(&mut self, subgrid: &SubgridEnum, _: Option<(usize, usize)>) {
        assert!(
            subgrid.is_empty(),
            "EmptySubgridV1 doesn't support the merge operation for non-empty subgrids"
        );
    }

    fn scale(&mut self, _: f64) {}

    fn symmetrize(&mut self, _: usize, _: usize) {}

    fn indexed_iter(&self) -> SubgridIndexedIter {
        Box::new(iter::empty())
    }

    fn stats(&self) -> Stats {
        Stats {
            total: 0,
            allocated: 0,
            zeros: 0,
            overhead: 0,
            bytes_per_value: 0,
        }
    }

    fn optimize_nodes(&mut self) {}

    fn repair(&mut self) -> bool { false }
}

/// TODO
#[derive(Clone, Deserialize, Serialize)]
pub struct ImportSubgridV1 {
    array: PackedArray<f64>,
    node_values: Vec<Vec<f64>>,
}

impl Subgrid for ImportSubgridV1 {
    fn fill(&mut self, _: &[Interp], _: &[f64], _: f64) {
        panic!("ImportSubgridV1 doesn't support the fill operation");
    }

    fn node_values(&self) -> Vec<Vec<f64>> {
        self.node_values.clone()
    }

    fn is_empty(&self) -> bool {
        self.array.is_empty()
    }

    fn repair(&mut self) -> bool { self.array.clear_if_empty() }

    fn merge_impl(&mut self, other: &SubgridEnum, transpose: Option<(usize, usize)>) {
        let lhs_node_values = self.node_values();
        let mut rhs_node_values = other.node_values();
        let mut new_node_values = lhs_node_values.clone();
        if let Some((a, b)) = transpose {
            rhs_node_values.swap(a, b);
        }

        if new_node_values != rhs_node_values {
            for (new, rhs) in new_node_values.iter_mut().zip(&rhs_node_values) {
                new.extend(rhs);
                new.sort_by(f64::total_cmp);
                new.dedup_by(node_value_eq_ref_mut);
            }

            let mut array = PackedArray::new(new_node_values.iter().map(Vec::len).collect());

            for (indices, value) in self.array.indexed_iter() {
                let target: Vec<_> = izip!(indices, &new_node_values, &lhs_node_values)
                    .map(|(index, new, lhs)| {
                        new.iter()
                            .position(|&value| node_value_eq(value, lhs[index]))
                            // UNWRAP: must succeed, `new_node_values` is the union of
                            // `lhs_node_values` and `rhs_node_values`
                            .unwrap()
                    })
                    .collect();

                array[target.as_slice()] = value;
            }

            self.array = array;
            self.node_values.clone_from(&new_node_values);
        }

        for (mut indices, value) in other.indexed_iter() {
            if let Some((a, b)) = transpose {
                indices.swap(a, b);
            }

            let target: Vec<_> = izip!(indices, &new_node_values, &rhs_node_values)
                .map(|(index, new, rhs)| {
                    new.iter()
                        .position(|&value| node_value_eq(value, rhs[index]))
                        // UNWRAP: must succeed, `new_node_values` is the union of
                        // `lhs_node_values` and `rhs_node_values`
                        .unwrap()
                })
                .collect();

            self.array[target.as_slice()] += value;
        }
    }

    fn scale(&mut self, factor: f64) {
        self.array *= factor;
    }

    fn symmetrize(&mut self, a: usize, b: usize) {
        let mut new_array = PackedArray::new(self.array.shape().to_vec());

        for (mut index, sigma) in self.array.indexed_iter() {
            // TODO: why not the other way around?
            if index[b] < index[a] {
                index.swap(a, b);
            }

            new_array[index.as_slice()] += sigma;
        }

        self.array = new_array;
    }

    fn indexed_iter(&self) -> SubgridIndexedIter {
        Box::new(self.array.indexed_iter())
    }

    fn shape(&self) -> &[usize] {
        self.array.shape()
    }

    fn stats(&self) -> Stats {
        Stats {
            total: self.array.shape().iter().product(),
            allocated: self.array.non_zeros() + self.array.explicit_zeros(),
            zeros: self.array.explicit_zeros(),
            overhead: self.array.overhead(),
            bytes_per_value: mem::size_of::<f64>(),
        }
    }

    fn optimize_nodes(&mut self) {}
}

impl ImportSubgridV1 {
    /// Constructor.
    #[must_use]
    pub const fn new(array: PackedArray<f64>, node_values: Vec<Vec<f64>>) -> Self {
        Self { array, node_values }
    }
}

impl From<&SubgridEnum> for ImportSubgridV1 {
    fn from(subgrid: &SubgridEnum) -> Self {
        // find smallest ranges
        let ranges: Vec<_> = subgrid.indexed_iter().fold(
            subgrid
                .node_values()
                .iter()
                .map(|values| values.len()..0)
                .collect(),
            |mut prev, (indices, _)| {
                for (i, index) in indices.iter().enumerate() {
                    prev[i].start = prev[i].start.min(*index);
                    prev[i].end = prev[i].end.max(*index + 1);
                }
                prev
            },
        );

        let node_values: Vec<_> = subgrid
            .node_values()
            .iter()
            .zip(&ranges)
            .map(|(values, range)| values[range.clone()].to_vec())
            .collect();

        let mut array = PackedArray::new(node_values.iter().map(Vec::len).collect());

        for (mut indices, value) in subgrid.indexed_iter() {
            for (index, range) in indices.iter_mut().zip(&ranges) {
                *index -= range.start;
            }

            array[indices.as_slice()] += value;
        }

        Self::new(array, node_values)
    }
}

/// Subgrid that uses interpolation.
#[derive(Clone, Deserialize, Serialize)]
pub struct InterpSubgridV1 {
    array: PackedArray<f64>,
    interps: Vec<Interp>,
    static_nodes: Vec<Option<f64>>,
}

impl InterpSubgridV1 {
    /// Constructor.
    #[must_use]
    pub fn new(interps: &[Interp]) -> Self {
        Self {
            array: PackedArray::new(interps.iter().map(Interp::nodes).collect()),
            interps: interps.to_vec(),
            static_nodes: vec![Some(-1.0); interps.len()],
        }
    }
}

impl Subgrid for InterpSubgridV1 {
    fn fill(&mut self, interps: &[Interp], ntuple: &[f64], weight: f64) {
        debug_assert_eq!(interps.len(), ntuple.len());

        if interpolation::interpolate(interps, ntuple, weight, &mut self.array) {
            for (value, previous_node) in ntuple.iter().zip(&mut self.static_nodes) {
                if let Some(previous_value) = previous_node {
                    if *previous_value < 0.0 {
                        *previous_value = *value;
                    } else if !node_value_eq(*previous_value, *value) {
                        *previous_node = None;
                    }
                }
            }
        }
    }

    fn node_values(&self) -> Vec<Vec<f64>> {
        self.interps.iter().map(Interp::node_values).collect()
    }

    fn is_empty(&self) -> bool {
        self.array.is_empty()
    }

    fn repair(&mut self) -> bool { self.array.clear_if_empty() }

    fn shape(&self) -> &[usize] {
        self.array.shape()
    }

    fn merge_impl(&mut self, other: &SubgridEnum, transpose: Option<(usize, usize)>) {
        // we cannot use `Self::indexed_iter` because it multiplies with `reweight`
        if let SubgridEnum::InterpSubgridV1(other) = other {
            // TODO: make sure `other` has the same interpolation as `self`
            for (mut index, value) in other.array.indexed_iter() {
                if let Some((a, b)) = transpose {
                    index.swap(a, b);
                }
                self.array[index.as_slice()] += value;
            }
        } else {
            unimplemented!();
        }
    }

    fn scale(&mut self, factor: f64) {
        self.array *= factor;
    }

    fn symmetrize(&mut self, a: usize, b: usize) {
        let mut new_array = PackedArray::new(self.array.shape().to_vec());

        for (mut index, sigma) in self.array.indexed_iter() {
            // TODO: why not the other way around?
            if index[b] < index[a] {
                index.swap(a, b);
            }

            new_array[index.as_slice()] += sigma;
        }

        self.array = new_array;
    }

    fn indexed_iter(&self) -> SubgridIndexedIter {
        let nodes: Vec<_> = self.interps.iter().map(Interp::node_values).collect();

        Box::new(self.array.indexed_iter().map(move |(indices, weight)| {
            let reweight = self
                .interps
                .iter()
                .enumerate()
                .map(|(i, interp)| interp.reweight(nodes[i][indices[i]]))
                .product::<f64>();
            (indices, weight * reweight)
        }))
    }

    fn stats(&self) -> Stats {
        Stats {
            total: self.array.shape().iter().product(),
            allocated: self.array.non_zeros() + self.array.explicit_zeros(),
            zeros: self.array.explicit_zeros(),
            overhead: self.array.overhead(),
            bytes_per_value: mem::size_of::<f64>(),
        }
    }

    fn optimize_nodes(&mut self) {
        // find the optimal ranges in which the nodes are used
        let ranges: Vec<_> = self.array.indexed_iter().fold(
            self.node_values()
                .iter()
                .map(|values| values.len()..0)
                .collect(),
            |mut prev, (indices, _)| {
                for (i, index) in indices.iter().enumerate() {
                    prev[i].start = prev[i].start.min(*index);
                    prev[i].end = prev[i].end.max(*index + 1);
                }
                prev
            },
        );

        let mut new_array = PackedArray::new(
            ranges
                .iter()
                .zip(&self.static_nodes)
                .map(|(range, static_node)| {
                    if static_node.is_some() {
                        1
                    } else {
                        range.clone().count()
                    }
                })
                .collect(),
        );

        for (mut index, value) in self.array.indexed_iter() {
            for (idx, range, static_node) in izip!(&mut index, &ranges, &self.static_nodes) {
                if static_node.is_some() {
                    *idx = 0;
                } else {
                    *idx -= range.start;
                }
            }
            new_array[index.as_slice()] += value;
        }

        self.array = new_array;

        for (interp, static_node, range) in izip!(&mut self.interps, &mut self.static_nodes, ranges)
        {
            *interp = if let &mut Some(value) = static_node {
                Interp::new(
                    value,
                    value,
                    1,
                    0,
                    interp.reweight_meth(),
                    interp.map(),
                    interp.interp_meth(),
                )
            } else {
                interp.sub_interp(range)
            };
        }
    }
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

impl SubgridEnum {
    /// TODO
    pub fn merge(&mut self, other: &Self, transpose: Option<(usize, usize)>) {
        if other.is_empty() {
            return;
        }
        if let Self::EmptySubgridV1(_) = self {
            if transpose.is_none() {
                *self = other.clone();
            } else {
                todo!();
            }
        } else {
            self.merge_impl(other, transpose);
        }
    }
}

/// Size-related statistics for a subgrid.
#[derive(Debug, Eq, PartialEq)]
pub struct Stats {
    /// Number of possible total entries for a subgrid. This number is the product of the lengths
    /// of the containers returned by [`Subgrid::node_values`].
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
    /// by `kinematics` in [`Grid::new`](super::grid::Grid::new) that was used to create the grid.
    fn fill(&mut self, interps: &[Interp], ntuple: &[f64], weight: f64);

    /// Returns true if `fill` was never called for this grid.
    fn is_empty(&self) -> bool;

    /// Merge `other` into this subgrid, possibly transposing the two dimensions given by
    /// `transpose`.
    fn merge_impl(&mut self, other: &SubgridEnum, transpose: Option<(usize, usize)>);

    /// Scale the subgrid by `factor`.
    fn scale(&mut self, factor: f64);

    /// Return the shape of the subgrid
    fn shape(&self) -> &[usize];

    /// Assume that the convolution functions for indices `a` and `b` for this grid are the same
    /// and use this to optimize the size of the grid.
    fn symmetrize(&mut self, a: usize, b: usize);

    /// Return an iterator over all non-zero elements of the subgrid.
    fn indexed_iter(&self) -> SubgridIndexedIter;

    /// Return statistics for this subgrid.
    fn stats(&self) -> Stats;

    /// TODO
    fn optimize_nodes(&mut self);

    /// Repair subgrid if necessary.
    fn repair (&mut self) -> bool;
}

/// Type to iterate over the non-zero contents of a subgrid. The tuple contains the indices of the
/// `mu2_grid`, the `x1_grid` and finally the `x2_grid`.
pub type SubgridIndexedIter<'a> = Box<dyn Iterator<Item = (Vec<usize>, f64)> + 'a>;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::interpolation::{InterpMeth, Map, ReweightMeth};
    use crate::v0;

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

    #[test]
    fn empty_subgrid_v1_new() {
        let mut subgrid = EmptySubgridV1;
        assert!(subgrid.is_empty());
        subgrid.merge_impl(&EmptySubgridV1.into(), None);
        subgrid.scale(2.0);
        subgrid.symmetrize(1, 2);
        subgrid.optimize_nodes();
        assert_eq!(
            subgrid.stats(),
            Stats {
                total: 0,
                allocated: 0,
                zeros: 0,
                overhead: 0,
                bytes_per_value: 0,
            }
        );
    }

    #[test]
    #[should_panic(expected = "EmptySubgridV1 doesn't support the fill operation")]
    fn empty_subgrid_v1_fill() {
        let mut subgrid = EmptySubgridV1;
        subgrid.fill(&v0::default_interps(false, 2), &[0.0; 3], 0.0);
    }

    #[test]
    #[should_panic(
        expected = "EmptySubgridV1 doesn't support the merge operation for non-empty subgrids"
    )]
    fn empty_subgrid_v1_merge_non_empty() {
        let mut subgrid_lhs = EmptySubgridV1;

        let mut array = PackedArray::new(vec![1, 1]);
        array[0] = 1.0;
        let node_values = vec![vec![1.0]; 2];
        let subgrid_rhs = ImportSubgridV1::new(array, node_values).into();

        subgrid_lhs.merge_impl(&subgrid_rhs, None);
    }

    #[test]
    fn empty_subgrid_v1_node_values() {
        assert!(EmptySubgridV1.node_values().is_empty());
    }

    #[test]
    #[should_panic(expected = "EmptySubgridV1 doesn't have a shape")]
    fn empty_subgrid_v1_shape() {
        let _ = EmptySubgridV1.shape();
    }

    #[test]
    #[should_panic(expected = "ImportSubgridV1 doesn't support the fill operation")]
    fn import_subgrid_v1_fill() {
        let mut subgrid =
            ImportSubgridV1::new(PackedArray::new(vec![0, 0, 0]), vec![Vec::new(); 3]);
        subgrid.fill(&v0::default_interps(false, 2), &[0.0; 3], 0.0);
    }

    #[test]
    fn import_subgrid_v1_test() {
        let x = vec![
            0.015625, 0.03125, 0.0625, 0.125, 0.1875, 0.25, 0.375, 0.5, 0.75, 1.0,
        ];
        let mut grid1: SubgridEnum = ImportSubgridV1::new(
            PackedArray::new(vec![1, 10, 10]),
            vec![vec![0.0], x.clone(), x.clone()],
        )
        .into();

        assert_eq!(grid1.node_values(), vec![vec![0.0], x.clone(), x.clone()]);

        assert!(grid1.is_empty());

        // only use exactly representable numbers here so that we can avoid using approx_eq
        if let SubgridEnum::ImportSubgridV1(ref mut x) = grid1 {
            x.array[[0, 1, 2]] = 1.0;
            x.array[[0, 1, 3]] = 2.0;
            x.array[[0, 4, 3]] = 4.0;
            x.array[[0, 7, 1]] = 8.0;
        }

        assert!(!grid1.is_empty());

        assert_eq!(grid1.indexed_iter().next(), Some((vec![0, 1, 2], 1.0)));
        assert_eq!(grid1.indexed_iter().nth(1), Some((vec![0, 1, 3], 2.0)));
        assert_eq!(grid1.indexed_iter().nth(2), Some((vec![0, 4, 3], 4.0)));
        assert_eq!(grid1.indexed_iter().nth(3), Some((vec![0, 7, 1], 8.0)));

        // create grid with transposed entries, but different q2
        let mut grid2: SubgridEnum = ImportSubgridV1::new(
            PackedArray::new(vec![1, 10, 10]),
            vec![vec![1.0], x.clone(), x],
        )
        .into();
        if let SubgridEnum::ImportSubgridV1(ref mut x) = grid2 {
            x.array[[0, 2, 1]] = 1.0;
            x.array[[0, 3, 1]] = 2.0;
            x.array[[0, 3, 4]] = 4.0;
            x.array[[0, 1, 7]] = 8.0;
        }

        assert_eq!(grid2.indexed_iter().next(), Some((vec![0, 1, 7], 8.0)));
        assert_eq!(grid2.indexed_iter().nth(1), Some((vec![0, 2, 1], 1.0)));
        assert_eq!(grid2.indexed_iter().nth(2), Some((vec![0, 3, 1], 2.0)));
        assert_eq!(grid2.indexed_iter().nth(3), Some((vec![0, 3, 4], 4.0)));

        grid1.merge(&grid2, None);

        // the luminosity function is symmetric, so after symmetrization the result must be
        // unchanged
        grid1.symmetrize(1, 2);

        grid1.scale(2.0);

        assert_eq!(
            grid1.stats(),
            Stats {
                total: 200,
                allocated: 8,
                zeros: 0,
                overhead: 12,
                bytes_per_value: 8,
            }
        );
    }

    #[test]
    fn interp_subgrid_v1_fill_zero() {
        let interps = v0::default_interps(false, 2);
        let mut subgrid = InterpSubgridV1::new(&interps);

        subgrid.fill(&interps, &[1000.0, 0.5, 0.5], 0.0);

        assert!(subgrid.is_empty());
        assert_eq!(subgrid.indexed_iter().count(), 0);
        assert_eq!(
            subgrid.stats(),
            Stats {
                total: 100000,
                allocated: 0,
                zeros: 0,
                overhead: 0,
                bytes_per_value: mem::size_of::<f64>()
            }
        );
    }

    #[test]
    fn interp_subgrid_v1_fill_outside_range() {
        let interps = v0::default_interps(false, 2);
        let mut subgrid = InterpSubgridV1::new(&interps);

        subgrid.fill(&interps, &[1000.0, 1e-10, 0.5], 0.0);

        assert!(subgrid.is_empty());
        assert_eq!(subgrid.indexed_iter().count(), 0);
        assert_eq!(
            subgrid.stats(),
            Stats {
                total: 100000,
                allocated: 0,
                zeros: 0,
                overhead: 0,
                bytes_per_value: mem::size_of::<f64>()
            }
        );
    }

    #[test]
    fn interp_subgrid_v1_fill() {
        let interps = v0::default_interps(false, 2);
        let mut subgrid = InterpSubgridV1::new(&interps);

        subgrid.fill(&interps, &[1000.0, 0.5, 0.5], 1.0);

        assert!(!subgrid.is_empty());
        assert_eq!(subgrid.indexed_iter().count(), 4 * 4 * 4);
        assert_eq!(
            subgrid.stats(),
            Stats {
                total: 100000,
                allocated: 64,
                zeros: 0,
                overhead: 32,
                bytes_per_value: mem::size_of::<f64>()
            }
        );

        subgrid.fill(&interps, &[1000000.0, 0.5, 0.5], 1.0);

        assert!(!subgrid.is_empty());
        assert_eq!(subgrid.indexed_iter().count(), 2 * 4 * 4 * 4);
        assert_eq!(
            subgrid.stats(),
            Stats {
                total: 100000,
                allocated: 128,
                zeros: 0,
                overhead: 64,
                bytes_per_value: mem::size_of::<f64>()
            }
        );

        subgrid.optimize_nodes();

        let node_values = subgrid.node_values();

        assert_eq!(node_values[0].len(), 23);
        assert_eq!(node_values[1].len(), 1);
        assert_eq!(node_values[2].len(), 1);

        assert_eq!(subgrid.shape(), [23, 1, 1]);
    }
}
