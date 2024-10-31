//! Module containing the Lagrange-interpolation subgrid.

use super::interpolation::{self, Interp};
use super::packed_array::PackedArray;
use super::subgrid::{self, Stats, Subgrid, SubgridEnum, SubgridIndexedIter};
use itertools::izip;
use serde::{Deserialize, Serialize};
use std::mem;

/// Subgrid which uses Lagrange-interpolation.
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
                    } else if !subgrid::node_value_eq(*previous_value, *value) {
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

    fn merge(&mut self, other: &SubgridEnum, transpose: Option<(usize, usize)>) {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::v0;

    #[test]
    fn fill_zero() {
        let interps = v0::default_interps(2);
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
    fn fill_outside_range() {
        let interps = v0::default_interps(2);
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
    fn fill() {
        let interps = v0::default_interps(2);
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
    }
}
