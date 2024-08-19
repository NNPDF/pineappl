//! TODO

use super::packed_array::PackedArray;
use super::subgrid::{Mu2, Stats, Subgrid, SubgridEnum, SubgridIndexedIter};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::mem;

/// TODO
#[derive(Clone, Default, Deserialize, Serialize)]
pub struct PackedQ1X2SubgridV1 {
    array: PackedArray<f64, 3>,
    mu2_grid: Vec<Mu2>,
    x1_grid: Vec<f64>,
    x2_grid: Vec<f64>,
}

impl PackedQ1X2SubgridV1 {
    /// Constructor.
    #[must_use]
    pub const fn new(
        array: PackedArray<f64, 3>,
        mu2_grid: Vec<Mu2>,
        x1_grid: Vec<f64>,
        x2_grid: Vec<f64>,
    ) -> Self {
        Self {
            array,
            mu2_grid,
            x1_grid,
            x2_grid,
        }
    }

    /// Return the array containing the numerical values of the grid.
    pub fn array_mut(&mut self) -> &mut PackedArray<f64, 3> {
        &mut self.array
    }
}

impl Subgrid for PackedQ1X2SubgridV1 {
    fn fill(&mut self, _: &[f64], _: f64) {
        panic!("PackedQ1X2SubgridV1 doesn't support the fill operation");
    }

    fn mu2_grid(&self) -> Cow<[Mu2]> {
        Cow::Borrowed(&self.mu2_grid)
    }

    fn x1_grid(&self) -> Cow<[f64]> {
        Cow::Borrowed(&self.x1_grid)
    }

    fn x2_grid(&self) -> Cow<[f64]> {
        Cow::Borrowed(&self.x2_grid)
    }

    fn is_empty(&self) -> bool {
        self.array.is_empty()
    }

    fn merge(&mut self, other: &mut SubgridEnum, transpose: bool) {
        if self.is_empty() && !transpose {
            if let SubgridEnum::PackedQ1X2SubgridV1(other) = other {
                *self = mem::take(other);
                return;
            }
        }

        let rhs_mu2 = other.mu2_grid().into_owned();
        let rhs_x1 = if transpose {
            other.x2_grid()
        } else {
            other.x1_grid()
        };
        let rhs_x2 = if transpose {
            other.x1_grid()
        } else {
            other.x2_grid()
        };

        if (self.mu2_grid != rhs_mu2) || (self.x1_grid() != rhs_x1) || (self.x2_grid() != rhs_x2) {
            let mut mu2_grid = self.mu2_grid.clone();
            let mut x1_grid = self.x1_grid.clone();
            let mut x2_grid = self.x2_grid.clone();

            mu2_grid.extend_from_slice(&rhs_mu2);
            mu2_grid.sort_by(|a, b| a.partial_cmp(b).unwrap());
            mu2_grid.dedup();
            x1_grid.extend_from_slice(&rhs_x1);
            x1_grid.sort_by(|a, b| a.partial_cmp(b).unwrap());
            x1_grid.dedup();
            x2_grid.extend_from_slice(&rhs_x2);
            x2_grid.sort_by(|a, b| a.partial_cmp(b).unwrap());
            x2_grid.dedup();

            let mut array = PackedArray::new([mu2_grid.len(), x1_grid.len(), x2_grid.len()]);

            for ([i, j, k], value) in self.array.indexed_iter() {
                let target_i = mu2_grid
                    .iter()
                    .position(|mu2| *mu2 == self.mu2_grid[i])
                    .unwrap_or_else(|| unreachable!());
                let target_j = x1_grid
                    .iter()
                    .position(|&x| x == self.x1_grid[j])
                    .unwrap_or_else(|| unreachable!());
                let target_k = x2_grid
                    .iter()
                    .position(|&x| x == self.x2_grid[k])
                    .unwrap_or_else(|| unreachable!());

                array[[target_i, target_j, target_k]] = value;
            }

            self.array = array;
            self.mu2_grid = mu2_grid;
            self.x1_grid = x1_grid;
            self.x2_grid = x2_grid;
        }

        for ((i, j, k), value) in other.indexed_iter() {
            let (j, k) = if transpose { (k, j) } else { (j, k) };
            let target_i = self
                .mu2_grid
                .iter()
                .position(|x| *x == rhs_mu2[i])
                .unwrap_or_else(|| unreachable!());
            let target_j = self
                .x1_grid
                .iter()
                .position(|&x| x == rhs_x1[j])
                .unwrap_or_else(|| unreachable!());
            let target_k = self
                .x2_grid
                .iter()
                .position(|&x| x == rhs_x2[k])
                .unwrap_or_else(|| unreachable!());

            self.array[[target_i, target_j, target_k]] += value;
        }
    }

    fn scale(&mut self, factor: f64) {
        self.array *= factor;
    }

    fn symmetrize(&mut self) {
        let mut new_array =
            PackedArray::new([self.mu2_grid.len(), self.x1_grid.len(), self.x2_grid.len()]);

        for ([i, j, k], sigma) in self.array.indexed_iter().filter(|([_, j, k], _)| k >= j) {
            new_array[[i, j, k]] = sigma;
        }
        // do not change the diagonal entries (k==j)
        for ([i, j, k], sigma) in self.array.indexed_iter().filter(|([_, j, k], _)| k < j) {
            new_array[[i, k, j]] += sigma;
        }

        mem::swap(&mut self.array, &mut new_array);
    }

    fn clone_empty(&self) -> SubgridEnum {
        Self {
            array: PackedArray::new([self.mu2_grid.len(), self.x1_grid.len(), self.x2_grid.len()]),
            mu2_grid: self.mu2_grid.clone(),
            x1_grid: self.x1_grid.clone(),
            x2_grid: self.x2_grid.clone(),
        }
        .into()
    }

    fn indexed_iter(&self) -> SubgridIndexedIter {
        Box::new(
            self.array
                .indexed_iter()
                .map(|([o, b, c], v)| ((o, b, c), v)),
        )
    }

    fn stats(&self) -> Stats {
        Stats {
            total: self.mu2_grid.len() * self.x1_grid.len() * self.x2_grid.len(),
            allocated: self.array.non_zeros() + self.array.explicit_zeros(),
            zeros: self.array.explicit_zeros(),
            overhead: self.array.overhead(),
            bytes_per_value: mem::size_of::<f64>(),
        }
    }

    fn static_scale(&self) -> Option<Mu2> {
        if let [static_scale] = self.mu2_grid.as_slice() {
            Some(static_scale.clone())
        } else {
            None
        }
    }
}

impl From<&SubgridEnum> for PackedQ1X2SubgridV1 {
    fn from(subgrid: &SubgridEnum) -> Self {
        // find smallest ranges
        let (mu2_range, x1_range, x2_range) = subgrid.indexed_iter().fold(
            (
                subgrid.mu2_grid().len()..0,
                subgrid.x1_grid().len()..0,
                subgrid.x2_grid().len()..0,
            ),
            |prev, ((imu2, ix1, ix2), _)| {
                (
                    prev.0.start.min(imu2)..prev.0.end.max(imu2 + 1),
                    prev.1.start.min(ix1)..prev.1.end.max(ix1 + 1),
                    prev.2.start.min(ix2)..prev.2.end.max(ix2 + 1),
                )
            },
        );

        let (mu2_grid, static_scale) = subgrid.static_scale().map_or_else(
            || (subgrid.mu2_grid()[mu2_range.clone()].to_vec(), false),
            |scale| (vec![scale], true),
        );
        let x1_grid = subgrid.x1_grid()[x1_range.clone()].to_vec();
        let x2_grid = subgrid.x2_grid()[x2_range.clone()].to_vec();

        let mut array = PackedArray::new([mu2_grid.len(), x1_grid.len(), x2_grid.len()]);

        for ((imu2, ix1, ix2), value) in subgrid.indexed_iter() {
            // if there's a static scale we want every value to be added to same grid point
            let index = if static_scale {
                0
            } else {
                imu2 - mu2_range.start
            };

            array[[index, ix1 - x1_range.start, ix2 - x2_range.start]] += value;
        }

        Self {
            array,
            mu2_grid,
            x1_grid,
            x2_grid,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic(expected = "PackedQ1X2SubgridV1 doesn't support the fill operation")]
    fn fill_packed_q1x2_subgrid_v1() {
        let mut subgrid = PackedQ1X2SubgridV1::new(
            PackedArray::new([0, 0, 0]),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );
        subgrid.fill(&[0.0; 4], 0.0);
    }

    #[test]
    fn test_v1() {
        let x = vec![
            0.015625, 0.03125, 0.0625, 0.125, 0.1875, 0.25, 0.375, 0.5, 0.75, 1.0,
        ];
        let mut grid1: SubgridEnum = PackedQ1X2SubgridV1::new(
            PackedArray::new([1, 10, 10]),
            vec![Mu2 {
                ren: 0.0,
                fac: 0.0,
                frg: 0.0,
            }],
            x.clone(),
            x.clone(),
        )
        .into();

        let mu2 = vec![Mu2 {
            ren: 0.0,
            fac: 0.0,
            frg: 0.0,
        }];

        assert_eq!(grid1.mu2_grid().as_ref(), mu2);
        assert_eq!(grid1.x1_grid().as_ref(), x);
        assert_eq!(grid1.x2_grid(), grid1.x1_grid());

        assert!(grid1.is_empty());

        // only use exactly representable numbers here so that we can avoid using approx_eq
        if let SubgridEnum::PackedQ1X2SubgridV1(ref mut x) = grid1 {
            x.array_mut()[[0, 1, 2]] = 1.0;
            x.array_mut()[[0, 1, 3]] = 2.0;
            x.array_mut()[[0, 4, 3]] = 4.0;
            x.array_mut()[[0, 7, 1]] = 8.0;
        }

        assert!(!grid1.is_empty());

        assert_eq!(grid1.indexed_iter().next(), Some(((0, 1, 2), 1.0)));
        assert_eq!(grid1.indexed_iter().nth(1), Some(((0, 1, 3), 2.0)));
        assert_eq!(grid1.indexed_iter().nth(2), Some(((0, 4, 3), 4.0)));
        assert_eq!(grid1.indexed_iter().nth(3), Some(((0, 7, 1), 8.0)));

        // create grid with transposed entries, but different q2
        let mut grid2: SubgridEnum = PackedQ1X2SubgridV1::new(
            PackedArray::new([1, 10, 10]),
            vec![Mu2 {
                ren: 1.0,
                fac: 1.0,
                frg: 1.0,
            }],
            x.clone(),
            x.clone(),
        )
        .into();
        if let SubgridEnum::PackedQ1X2SubgridV1(ref mut x) = grid2 {
            x.array_mut()[[0, 2, 1]] = 1.0;
            x.array_mut()[[0, 3, 1]] = 2.0;
            x.array_mut()[[0, 3, 4]] = 4.0;
            x.array_mut()[[0, 1, 7]] = 8.0;
        }

        assert_eq!(grid2.indexed_iter().next(), Some(((0, 1, 7), 8.0)));
        assert_eq!(grid2.indexed_iter().nth(1), Some(((0, 2, 1), 1.0)));
        assert_eq!(grid2.indexed_iter().nth(2), Some(((0, 3, 1), 2.0)));
        assert_eq!(grid2.indexed_iter().nth(3), Some(((0, 3, 4), 4.0)));

        grid1.merge(&mut grid2, false);

        let mut grid1 = {
            let mut g = grid1.clone_empty();
            g.merge(&mut grid1, false);
            g
        };

        // the luminosity function is symmetric, so after symmetrization the result must be
        // unchanged
        grid1.symmetrize();

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
}
