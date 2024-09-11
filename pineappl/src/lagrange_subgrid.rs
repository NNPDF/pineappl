//! Module containing the Lagrange-interpolation subgrid.

use super::interpolation::{self, Interp};
use super::packed_array::PackedArray;
use super::subgrid::{Mu2, Stats, Subgrid, SubgridEnum, SubgridIndexedIter};
use float_cmp::approx_eq;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::mem;

/// Subgrid which uses Lagrange-interpolation.
#[derive(Clone, Deserialize, Serialize)]
pub struct LagrangeSubgridV2 {
    grid: PackedArray<f64, 3>,
    interps: Vec<Interp>,
    pub(crate) static_q2: f64,
}

impl LagrangeSubgridV2 {
    /// Constructor.
    #[must_use]
    pub fn new(interps: &[Interp]) -> Self {
        debug_assert_eq!(interps.len(), 3);
        Self {
            grid: PackedArray::new([interps[0].nodes(), interps[1].nodes(), interps[2].nodes()]),
            interps: interps.to_vec(),
            static_q2: 0.0,
        }
    }
}

impl Subgrid for LagrangeSubgridV2 {
    fn fill(&mut self, interps: &[Interp], ntuple: &[f64], weight: f64) {
        debug_assert_eq!(interps.len(), ntuple.len());

        // TODO: lift the restriction to exactly three kinematics variables
        debug_assert_eq!(interps.len(), 3);

        if interpolation::interpolate(interps, &ntuple, weight, &mut self.grid) {
            let q2 = ntuple[0];
            if self.static_q2 == 0.0 {
                self.static_q2 = q2;
            } else if (self.static_q2 != -1.0) && !approx_eq!(f64, self.static_q2, q2, ulps = 4) {
                self.static_q2 = -1.0;
            }
        }
    }

    fn mu2_grid(&self) -> Cow<[Mu2]> {
        self.interps[0]
            .node_values()
            .iter()
            .map(|&q2| Mu2 {
                ren: q2,
                fac: q2,
                frg: -1.0,
            })
            .collect()
    }

    fn x1_grid(&self) -> Cow<[f64]> {
        self.interps[1].node_values().into()
    }

    fn x2_grid(&self) -> Cow<[f64]> {
        self.interps[2].node_values().into()
    }

    fn is_empty(&self) -> bool {
        self.grid.is_empty()
    }

    fn merge(&mut self, other: &mut SubgridEnum, transpose: bool) {
        // we cannot use `Self::indexed_iter` because it multiplies with `reweight`
        if let SubgridEnum::LagrangeSubgridV2(other) = other {
            // TODO: make sure `other` has the same interpolation as `self`
            for (mut index, value) in other.grid.indexed_iter() {
                if transpose {
                    index.swap(1, 2);
                }
                self.grid[index] += value;
            }
        } else {
            unimplemented!();
        }
    }

    fn scale(&mut self, factor: f64) {
        self.grid *= factor;
    }

    fn symmetrize(&mut self) {
        let mut new_array = PackedArray::new([
            self.mu2_grid().len(),
            self.x1_grid().len(),
            self.x2_grid().len(),
        ]);

        for ([i, j, k], sigma) in self.grid.indexed_iter().filter(|([_, j, k], _)| k >= j) {
            new_array[[i, j, k]] = sigma;
        }
        // do not change the diagonal entries (k==j)
        for ([i, j, k], sigma) in self.grid.indexed_iter().filter(|([_, j, k], _)| k < j) {
            new_array[[i, k, j]] += sigma;
        }

        self.grid = new_array;
    }

    fn clone_empty(&self) -> SubgridEnum {
        self.clone().into()
    }

    fn indexed_iter(&self) -> SubgridIndexedIter {
        let nodes: Vec<_> = self.interps.iter().map(Interp::node_values).collect();

        Box::new(self.grid.indexed_iter().map(move |(index, v)| {
            (
                (index[0], index[1], index[2]),
                v * self
                    .interps
                    .iter()
                    .enumerate()
                    .map(|(i, interp)| interp.reweight(nodes[i][index[i]]))
                    .product::<f64>(),
            )
        }))
    }

    fn stats(&self) -> Stats {
        Stats {
            total: self.mu2_grid().len() * self.x1_grid().len() * self.x2_grid().len(),
            allocated: self.grid.non_zeros() + self.grid.explicit_zeros(),
            zeros: self.grid.explicit_zeros(),
            overhead: self.grid.overhead(),
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
        let interps = v0::default_interps();
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
        let interps = v0::default_interps();
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
        let interps = v0::default_interps();
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
