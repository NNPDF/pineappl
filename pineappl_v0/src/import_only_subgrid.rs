//! TODO

use super::grid::Ntuple;
use super::sparse_array3::SparseArray3;
use super::subgrid::{Mu2, Stats, Subgrid, SubgridEnum, SubgridIndexedIter};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::mem;

/// TODO
#[derive(Clone, Deserialize, Serialize)]
pub struct ImportOnlySubgridV1 {
    array: SparseArray3<f64>,
    q2_grid: Vec<f64>,
    x1_grid: Vec<f64>,
    x2_grid: Vec<f64>,
}

impl ImportOnlySubgridV1 {
    /// Constructor.
    #[must_use]
    pub fn new(
        array: SparseArray3<f64>,
        q2_grid: Vec<f64>,
        x1_grid: Vec<f64>,
        x2_grid: Vec<f64>,
    ) -> Self {
        Self {
            array,
            q2_grid,
            x1_grid,
            x2_grid,
        }
    }

    /// Return the array containing the numerical values of the grid.
    pub fn array_mut(&mut self) -> &mut SparseArray3<f64> {
        &mut self.array
    }
}

impl Subgrid for ImportOnlySubgridV1 {
    fn convolve(
        &self,
        _: &[f64],
        _: &[f64],
        _: &[Mu2],
        lumi: &mut dyn FnMut(usize, usize, usize) -> f64,
    ) -> f64 {
        self.array
            .indexed_iter()
            .map(|((imu2, ix1, ix2), sigma)| sigma * lumi(ix1, ix2, imu2))
            .sum()
    }

    fn fill(&mut self, _: &Ntuple<f64>) {
        panic!("ImportOnlySubgridV1 doesn't support the fill operation");
    }

    fn mu2_grid(&self) -> Cow<'_, [Mu2]> {
        self.q2_grid
            .iter()
            .copied()
            .map(|q2| Mu2 { ren: q2, fac: q2 })
            .collect()
    }

    fn x1_grid(&self) -> Cow<'_, [f64]> {
        Cow::Borrowed(&self.x1_grid)
    }

    fn x2_grid(&self) -> Cow<'_, [f64]> {
        Cow::Borrowed(&self.x2_grid)
    }

    fn is_empty(&self) -> bool {
        self.array.is_empty()
    }

    fn merge(&mut self, other: &mut SubgridEnum, transpose: bool) {
        if let SubgridEnum::ImportOnlySubgridV1(other_grid) = other {
            if self.array.is_empty() && !transpose {
                mem::swap(&mut self.array, &mut other_grid.array);
            } else {
                // TODO: the general case isn't implemented
                assert!(self.x1_grid() == other_grid.x1_grid());
                assert!(self.x2_grid() == other_grid.x2_grid());

                for (other_index, mu2) in other_grid.mu2_grid().iter().enumerate() {
                    // the following should always be the case
                    assert_eq!(mu2.ren, mu2.fac);
                    let q2 = &mu2.ren;

                    let index = match self
                        .q2_grid
                        .binary_search_by(|val| val.partial_cmp(q2).unwrap())
                    {
                        Ok(index) => index,
                        Err(index) => {
                            self.q2_grid.insert(index, *q2);
                            self.array.increase_x_at(index);
                            index
                        }
                    };

                    for ((_, j, k), value) in other_grid
                        .array
                        .indexed_iter()
                        .filter(|&((i, _, _), _)| i == other_index)
                    {
                        let (j, k) = if transpose { (k, j) } else { (j, k) };
                        self.array[[index, j, k]] += value;
                    }
                }
            }
        } else {
            todo!();
        }
    }

    fn scale(&mut self, factor: f64) {
        if factor == 0.0 {
            self.array.clear();
        } else {
            self.array.iter_mut().for_each(|x| *x *= factor);
        }
    }

    fn symmetrize(&mut self) {
        let mut new_array =
            SparseArray3::new(self.q2_grid.len(), self.x1_grid.len(), self.x2_grid.len());

        for ((i, j, k), sigma) in self.array.indexed_iter().filter(|((_, j, k), _)| k >= j) {
            new_array[[i, j, k]] = sigma;
        }
        // do not change the diagonal entries (k==j)
        for ((i, j, k), sigma) in self.array.indexed_iter().filter(|((_, j, k), _)| k < j) {
            new_array[[i, k, j]] += sigma;
        }

        mem::swap(&mut self.array, &mut new_array);
    }

    fn clone_empty(&self) -> SubgridEnum {
        Self {
            array: SparseArray3::new(self.q2_grid.len(), self.x1_grid.len(), self.x2_grid.len()),
            q2_grid: self.q2_grid.clone(),
            x1_grid: self.x1_grid.clone(),
            x2_grid: self.x2_grid.clone(),
        }
        .into()
    }

    fn indexed_iter(&self) -> SubgridIndexedIter<'_> {
        Box::new(self.array.indexed_iter())
    }

    fn stats(&self) -> Stats {
        Stats {
            total: self.q2_grid.len() * self.x1_grid.len() * self.x2_grid.len(),
            allocated: self.array.len() + self.array.zeros(),
            zeros: self.array.zeros(),
            overhead: self.array.overhead(),
            bytes_per_value: mem::size_of::<f64>(),
        }
    }

    fn static_scale(&self) -> Option<Mu2> {
        if let &[static_scale] = self.q2_grid.as_slice() {
            Some(Mu2 {
                ren: static_scale,
                fac: static_scale,
            })
        } else {
            None
        }
    }
}

/// TODO
#[derive(Clone, Deserialize, Serialize)]
pub struct ImportOnlySubgridV2 {
    array: SparseArray3<f64>,
    mu2_grid: Vec<Mu2>,
    x1_grid: Vec<f64>,
    x2_grid: Vec<f64>,
}

impl ImportOnlySubgridV2 {
    /// Constructor.
    #[must_use]
    pub fn new(
        array: SparseArray3<f64>,
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
    pub fn array_mut(&mut self) -> &mut SparseArray3<f64> {
        &mut self.array
    }
}

impl Subgrid for ImportOnlySubgridV2 {
    fn convolve(
        &self,
        _: &[f64],
        _: &[f64],
        _: &[Mu2],
        lumi: &mut dyn FnMut(usize, usize, usize) -> f64,
    ) -> f64 {
        self.array
            .indexed_iter()
            .map(|((imu2, ix1, ix2), sigma)| sigma * lumi(ix1, ix2, imu2))
            .sum()
    }

    fn fill(&mut self, _: &Ntuple<f64>) {
        panic!("ImportOnlySubgridV2 doesn't support the fill operation");
    }

    fn mu2_grid(&self) -> Cow<'_, [Mu2]> {
        Cow::Borrowed(&self.mu2_grid)
    }

    fn x1_grid(&self) -> Cow<'_, [f64]> {
        Cow::Borrowed(&self.x1_grid)
    }

    fn x2_grid(&self) -> Cow<'_, [f64]> {
        Cow::Borrowed(&self.x2_grid)
    }

    fn is_empty(&self) -> bool {
        self.array.is_empty()
    }

    fn merge(&mut self, other: &mut SubgridEnum, transpose: bool) {
        if let SubgridEnum::ImportOnlySubgridV2(other_grid) = other {
            if self.array.is_empty() && !transpose {
                mem::swap(&mut self.array, &mut other_grid.array);
            } else {
                let rhs_x1 = if transpose {
                    other_grid.x2_grid()
                } else {
                    other_grid.x1_grid()
                };
                let rhs_x2 = if transpose {
                    other_grid.x1_grid()
                } else {
                    other_grid.x2_grid()
                };

                if (self.x1_grid() != rhs_x1) || (self.x2_grid() != rhs_x2) {
                    let mut x1_grid = self.x1_grid.clone();
                    let mut x2_grid = self.x2_grid.clone();

                    x1_grid.extend_from_slice(&rhs_x1);
                    x1_grid.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    x1_grid.dedup();
                    x2_grid.extend_from_slice(&rhs_x2);
                    x2_grid.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    x2_grid.dedup();

                    let mut array =
                        SparseArray3::new(self.array.dimensions().0, x1_grid.len(), x2_grid.len());

                    for ((i, j, k), value) in self.array.indexed_iter() {
                        let target_j = x1_grid
                            .iter()
                            .position(|&x| x == self.x1_grid[j])
                            .unwrap_or_else(|| unreachable!());
                        let target_k = x2_grid
                            .iter()
                            .position(|&x| x == self.x2_grid[k])
                            .unwrap_or_else(|| unreachable!());

                        array[[i, target_j, target_k]] = value;
                    }

                    self.array = array;
                    self.x1_grid = x1_grid;
                    self.x2_grid = x2_grid;
                }

                for (other_index, mu2) in other_grid.mu2_grid().iter().enumerate() {
                    let index = match self
                        .mu2_grid
                        .binary_search_by(|val| val.partial_cmp(mu2).unwrap())
                    {
                        Ok(index) => index,
                        Err(index) => {
                            self.mu2_grid.insert(index, mu2.clone());
                            self.array.increase_x_at(index);
                            index
                        }
                    };

                    for ((_, j, k), value) in other_grid
                        .array
                        .indexed_iter()
                        .filter(|&((i, _, _), _)| i == other_index)
                    {
                        let (j, k) = if transpose { (k, j) } else { (j, k) };
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

                        self.array[[index, target_j, target_k]] += value;
                    }
                }
            }
        } else {
            todo!();
        }
    }

    fn scale(&mut self, factor: f64) {
        if factor == 0.0 {
            self.array.clear();
        } else {
            self.array.iter_mut().for_each(|x| *x *= factor);
        }
    }

    fn symmetrize(&mut self) {
        let mut new_array =
            SparseArray3::new(self.mu2_grid.len(), self.x1_grid.len(), self.x2_grid.len());

        for ((i, j, k), sigma) in self.array.indexed_iter().filter(|((_, j, k), _)| k >= j) {
            new_array[[i, j, k]] = sigma;
        }
        // do not change the diagonal entries (k==j)
        for ((i, j, k), sigma) in self.array.indexed_iter().filter(|((_, j, k), _)| k < j) {
            new_array[[i, k, j]] += sigma;
        }

        mem::swap(&mut self.array, &mut new_array);
    }

    fn clone_empty(&self) -> SubgridEnum {
        Self {
            array: SparseArray3::new(self.mu2_grid.len(), self.x1_grid.len(), self.x2_grid.len()),
            mu2_grid: self.mu2_grid.clone(),
            x1_grid: self.x1_grid.clone(),
            x2_grid: self.x2_grid.clone(),
        }
        .into()
    }

    fn indexed_iter(&self) -> SubgridIndexedIter<'_> {
        Box::new(self.array.indexed_iter())
    }

    fn stats(&self) -> Stats {
        Stats {
            total: self.mu2_grid.len() * self.x1_grid.len() * self.x2_grid.len(),
            allocated: self.array.len() + self.array.zeros(),
            zeros: self.array.zeros(),
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

impl From<&SubgridEnum> for ImportOnlySubgridV2 {
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

        let mut array = SparseArray3::new(mu2_grid.len(), x1_grid.len(), x2_grid.len());

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
