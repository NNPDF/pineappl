//! TODO

use super::grid::Ntuple;
use super::lagrange_subgrid::{self, LagrangeSubgridV2};
use super::sparse_array3::SparseArray3;
use super::subgrid::{Subgrid, SubgridEnum};
use either::Either;
use ndarray::Axis;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::mem;
use std::ops::Range;

/// TODO
#[derive(Deserialize, Serialize)]
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
}

impl ImportOnlySubgridV1 {
    /// Return the array containing the numerical values of the grid.
    pub fn array_mut(&mut self) -> &mut SparseArray3<f64> {
        &mut self.array
    }
}

impl Subgrid for ImportOnlySubgridV1 {
    fn convolute(
        &self,
        _: &[f64],
        _: &[f64],
        _: &[f64],
        lumi: Either<&dyn Fn(usize, usize, usize) -> f64, &dyn Fn(f64, f64, f64) -> f64>,
    ) -> f64 {
        let lumi = lumi.left().unwrap();

        self.array
            .indexed_iter()
            .map(|((iq2, ix1, ix2), sigma)| sigma * lumi(ix1, ix2, iq2))
            .sum()
    }

    fn fill(&mut self, _: &Ntuple<f64>) {
        panic!("this grid doesn't support the fill operation");
    }

    fn q2_grid(&self) -> Cow<[f64]> {
        Cow::Borrowed(&self.q2_grid)
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
        if let SubgridEnum::ImportOnlySubgridV1(other_grid) = other {
            if self.array.is_empty() && !transpose {
                mem::swap(&mut self.array, &mut other_grid.array);
            } else {
                // TODO: the general case isn't implemented
                assert!(self.x1_grid() == other_grid.x1_grid());
                assert!(self.x2_grid() == other_grid.x2_grid());

                if self.q2_grid() == other_grid.q2_grid() {
                    if transpose {
                        for ((i, k, j), value) in other_grid.array.indexed_iter() {
                            self.array[[i, j, k]] += value;
                        }
                    } else {
                        for ((i, j, k), value) in other_grid.array.indexed_iter() {
                            self.array[[i, j, k]] += value;
                        }
                    }
                } else {
                    for (other_index, q2) in other_grid.q2_grid().iter().enumerate() {
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

    fn q2_slice(&self) -> Range<usize> {
        self.array.x_range()
    }

    fn fill_q2_slice(&self, q2_slice: usize, grid: &mut [f64]) {
        let x1: Vec<_> = self.x1_grid.iter().map(|&x| 1.0 / x).collect();
        let x2: Vec<_> = self.x2_grid.iter().map(|&x| 1.0 / x).collect();

        for value in grid.iter_mut() {
            *value = 0.0;
        }

        for ((_, ix1, ix2), value) in self
            .array
            .indexed_iter()
            .filter(|((iq2, _, _), _)| *iq2 == q2_slice)
        {
            grid[ix1 * self.x2_grid.len() + ix2] = value * x1[ix1] * x2[ix2];
        }
    }

    fn symmetrize(&mut self) {
        let mut new_array =
            SparseArray3::new(self.q2_grid.len(), self.x1_grid.len(), self.x2_grid.len());

        for ((i, j, k), &sigma) in self.array.indexed_iter().filter(|((_, j, k), _)| k >= j) {
            new_array[[i, j, k]] = sigma;
        }
        // do not change the diagonal entries (k==j)
        for ((i, j, k), &sigma) in self.array.indexed_iter().filter(|((_, j, k), _)| k < j) {
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

    fn iter(&self) -> Box<dyn Iterator<Item = ((usize, usize, usize), &f64)> + '_> {
        Box::new(self.array.indexed_iter())
    }
}

impl From<&LagrangeSubgridV2> for ImportOnlySubgridV1 {
    fn from(subgrid: &LagrangeSubgridV2) -> Self {
        let array = subgrid.grid.as_ref().map_or_else(
            || SparseArray3::new(subgrid.ntau, subgrid.ny1, subgrid.ny2),
            // in the following case we should optimize when ny2 > ny1
            |array| {
                let reweight_x1: Vec<_> = subgrid
                    .x1_grid()
                    .iter()
                    .map(|x| lagrange_subgrid::weightfun(*x))
                    .collect();
                let reweight_x2: Vec<_> = subgrid
                    .x2_grid()
                    .iter()
                    .map(|x| lagrange_subgrid::weightfun(*x))
                    .collect();

                if subgrid.static_q2 > 0.0 {
                    // in this case we've detected a static scale for this bin and we can collapse
                    // the Q^2 axis into a single bin

                    let mut array = array
                        .sum_axis(Axis(0))
                        .into_shape((1, subgrid.ny1, subgrid.ny2))
                        .unwrap();
                    for ((_, ix1, ix2), entry) in array.indexed_iter_mut() {
                        *entry *= reweight_x1[ix1] * reweight_x2[ix2];
                    }
                    SparseArray3::from_ndarray(&array, 0, 1)
                } else {
                    let mut array = array.clone();
                    for ((_, ix1, ix2), entry) in array.indexed_iter_mut() {
                        *entry *= reweight_x1[ix1] * reweight_x2[ix2];
                    }
                    SparseArray3::from_ndarray(&array, subgrid.itaumin, subgrid.ntau)
                }
            },
        );
        let q2_grid = if subgrid.static_q2 > 0.0 {
            vec![subgrid.static_q2]
        } else {
            subgrid.q2_grid().into_owned()
        };
        let x1_grid = subgrid.x1_grid().into_owned();
        let x2_grid = subgrid.x2_grid().into_owned();

        Self {
            array,
            q2_grid,
            x1_grid,
            x2_grid,
        }
    }
}
