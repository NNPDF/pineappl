//! Module containing the Lagrange-interpolation subgrid.

use super::interpolation::{self, Interp};
use super::packed_array::PackedArray;
use super::subgrid::{Mu2, NodeValues, Stats, Subgrid, SubgridEnum, SubgridIndexedIter};
use float_cmp::approx_eq;
use serde::{Deserialize, Serialize};
use std::mem;

/// Subgrid which uses Lagrange-interpolation.
#[derive(Clone, Deserialize, Serialize)]
pub struct LagrangeSubgridV2 {
    array: PackedArray<f64>,
    interps: Vec<Interp>,
    pub(crate) static_q2: f64,
}

impl LagrangeSubgridV2 {
    /// Constructor.
    #[must_use]
    pub fn new(interps: &[Interp]) -> Self {
        Self {
            array: PackedArray::new(interps.iter().map(Interp::nodes).collect()),
            interps: interps.to_vec(),
            static_q2: 0.0,
        }
    }
}

impl Subgrid for LagrangeSubgridV2 {
    fn fill(&mut self, interps: &[Interp], ntuple: &[f64], weight: f64) {
        debug_assert_eq!(interps.len(), ntuple.len());

        if interpolation::interpolate(interps, ntuple, weight, &mut self.array) {
            // TODO: make this more general
            let q2 = ntuple[0];
            if self.static_q2 == 0.0 {
                self.static_q2 = q2;
            } else if (self.static_q2 != -1.0) && !approx_eq!(f64, self.static_q2, q2, ulps = 4) {
                self.static_q2 = -1.0;
            }
        }
    }

    fn node_values(&self) -> Vec<NodeValues> {
        self.interps
            .iter()
            .map(|interp| NodeValues::UseThese(interp.node_values()))
            .collect()
    }

    fn is_empty(&self) -> bool {
        self.array.is_empty()
    }

    fn merge(&mut self, other: &SubgridEnum, transpose: Option<(usize, usize)>) {
        // we cannot use `Self::indexed_iter` because it multiplies with `reweight`
        if let SubgridEnum::LagrangeSubgridV2(other) = other {
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

    fn static_scale(&self) -> Option<Mu2> {
        (self.static_q2 > 0.0).then_some(Mu2 {
            ren: self.static_q2,
            fac: self.static_q2,
            frg: -1.0,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::v0;

    #[test]
    fn fill_zero() {
        let interps = v0::default_interps(2);
        let mut subgrid = LagrangeSubgridV2::new(&interps);

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
        let mut subgrid = LagrangeSubgridV2::new(&interps);

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
        let mut subgrid = LagrangeSubgridV2::new(&interps);

        subgrid.fill(&interps, &[1000.0, 0.5, 0.5], 1.0);

        assert!(!subgrid.is_empty());
        assert_eq!(subgrid.indexed_iter().count(), 4 * 4 * 4);
        assert_eq!(
            subgrid.static_scale(),
            Some(Mu2 {
                ren: 1000.0,
                fac: 1000.0,
                frg: -1.0
            })
        );
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
        assert_eq!(subgrid.static_scale(), None);
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
    }
}
