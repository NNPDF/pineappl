//! TODO

use super::grid::Ntuple;
use super::lagrange_subgrid::{self, LagrangeSubgridV2};
use super::sparse_array3::SparseArray3;
use super::subgrid::{Mu2, Stats, Subgrid, SubgridEnum, SubgridIter};
use ndarray::Axis;
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
    fn convolute(
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

    fn mu2_grid(&self) -> Cow<[Mu2]> {
        self.q2_grid
            .iter()
            .copied()
            .map(|q2| Mu2 { ren: q2, fac: q2 })
            .collect()
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

    fn iter(&self) -> SubgridIter {
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
                    SparseArray3::from_ndarray(&array, 0, subgrid.itaumax - subgrid.itaumin)
                }
            },
        );
        let q2_grid = if subgrid.static_q2 > 0.0 {
            vec![subgrid.static_q2]
        } else {
            subgrid
                .mu2_grid()
                .iter()
                .skip(subgrid.itaumin)
                .take(subgrid.itaumax - subgrid.itaumin)
                .map(|mu2| mu2.fac)
                .collect()
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
    fn convolute(
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
        if let SubgridEnum::ImportOnlySubgridV2(other_grid) = other {
            if self.array.is_empty() && !transpose {
                mem::swap(&mut self.array, &mut other_grid.array);
            } else {
                // TODO: the general case isn't implemented
                assert!(self.x1_grid() == other_grid.x1_grid());
                assert!(self.x2_grid() == other_grid.x2_grid());

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
            SparseArray3::new(self.mu2_grid.len(), self.x1_grid.len(), self.x2_grid.len());

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
            array: SparseArray3::new(self.mu2_grid.len(), self.x1_grid.len(), self.x2_grid.len()),
            mu2_grid: self.mu2_grid.clone(),
            x1_grid: self.x1_grid.clone(),
            x2_grid: self.x2_grid.clone(),
        }
        .into()
    }

    fn iter(&self) -> SubgridIter {
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::subgrid::{ExtraSubgridParams, SubgridParams};
    use float_cmp::assert_approx_eq;

    #[test]
    fn test_v1() {
        let x = vec![
            0.015625, 0.03125, 0.0625, 0.125, 0.1875, 0.25, 0.375, 0.5, 0.75, 1.0,
        ];
        let mut grid1: SubgridEnum = ImportOnlySubgridV1::new(
            SparseArray3::new(1, 10, 10),
            vec![0.0],
            x.clone(),
            x.clone(),
        )
        .into();

        assert_eq!(
            grid1.stats(),
            Stats {
                total: 100,
                allocated: 0,
                zeros: 0,
                overhead: 2,
                bytes_per_value: 8,
            }
        );

        let mu2 = vec![Mu2 { ren: 0.0, fac: 0.0 }];

        assert_eq!(grid1.mu2_grid().as_ref(), mu2);
        assert_eq!(grid1.x1_grid().as_ref(), x);
        assert_eq!(grid1.x2_grid(), grid1.x1_grid());

        assert!(grid1.is_empty());

        // only use exactly representable numbers here so that we can avoid using approx_eq
        if let SubgridEnum::ImportOnlySubgridV1(ref mut x) = grid1 {
            x.array_mut()[[0, 1, 2]] = 1.0;
            x.array_mut()[[0, 1, 3]] = 2.0;
            x.array_mut()[[0, 4, 3]] = 4.0;
            x.array_mut()[[0, 7, 1]] = 8.0;
        } else {
            unreachable!();
        }

        assert!(!grid1.is_empty());

        assert_eq!(grid1.iter().nth(0), Some(((0, 1, 2), &1.0)));
        assert_eq!(grid1.iter().nth(1), Some(((0, 1, 3), &2.0)));
        assert_eq!(grid1.iter().nth(2), Some(((0, 4, 3), &4.0)));
        assert_eq!(grid1.iter().nth(3), Some(((0, 7, 1), &8.0)));

        // symmetric luminosity function
        let lumi =
            &mut (|ix1, ix2, _| x[ix1] * x[ix2]) as &mut dyn FnMut(usize, usize, usize) -> f64;

        assert_eq!(grid1.convolute(&x, &x, &mu2, lumi), 0.228515625);

        // create grid with transposed entries, but different q2
        let mut grid2: SubgridEnum = ImportOnlySubgridV1::new(
            SparseArray3::new(1, 10, 10),
            vec![1.0],
            x.clone(),
            x.clone(),
        )
        .into();
        if let SubgridEnum::ImportOnlySubgridV1(ref mut x) = grid2 {
            x.array_mut()[[0, 2, 1]] = 1.0;
            x.array_mut()[[0, 3, 1]] = 2.0;
            x.array_mut()[[0, 3, 4]] = 4.0;
            x.array_mut()[[0, 1, 7]] = 8.0;
        } else {
            unreachable!();
        }
        assert_eq!(grid2.convolute(&x, &x, &mu2, lumi), 0.228515625);

        assert_eq!(grid2.iter().nth(0), Some(((0, 1, 7), &8.0)));
        assert_eq!(grid2.iter().nth(1), Some(((0, 2, 1), &1.0)));
        assert_eq!(grid2.iter().nth(2), Some(((0, 3, 1), &2.0)));
        assert_eq!(grid2.iter().nth(3), Some(((0, 3, 4), &4.0)));

        grid1.merge(&mut grid2, false);

        assert_eq!(grid1.convolute(&x, &x, &mu2, lumi), 2.0 * 0.228515625);

        let mut grid1 = {
            let mut g = grid1.clone_empty();
            g.merge(&mut grid1, false);
            g
        };

        // the luminosity function is symmetric, so after symmetrization the result must be
        // unchanged
        grid1.symmetrize();
        assert_eq!(grid1.convolute(&x, &x, &mu2, lumi), 2.0 * 0.228515625);

        grid1.scale(2.0);
        assert_eq!(grid1.convolute(&x, &x, &mu2, lumi), 4.0 * 0.228515625);

        assert_eq!(
            grid1.stats(),
            Stats {
                total: 200,
                allocated: 14,
                zeros: 6,
                overhead: 42,
                bytes_per_value: 8,
            }
        );
    }

    #[test]
    fn test_v2() {
        let x = vec![
            0.015625, 0.03125, 0.0625, 0.125, 0.1875, 0.25, 0.375, 0.5, 0.75, 1.0,
        ];
        let mut grid1: SubgridEnum = ImportOnlySubgridV2::new(
            SparseArray3::new(1, 10, 10),
            vec![Mu2 { ren: 0.0, fac: 0.0 }],
            x.clone(),
            x.clone(),
        )
        .into();

        let mu2 = vec![Mu2 { ren: 0.0, fac: 0.0 }];

        assert_eq!(grid1.mu2_grid().as_ref(), mu2);
        assert_eq!(grid1.x1_grid().as_ref(), x);
        assert_eq!(grid1.x2_grid(), grid1.x1_grid());

        assert!(grid1.is_empty());

        // only use exactly representable numbers here so that we can avoid using approx_eq
        if let SubgridEnum::ImportOnlySubgridV2(ref mut x) = grid1 {
            x.array_mut()[[0, 1, 2]] = 1.0;
            x.array_mut()[[0, 1, 3]] = 2.0;
            x.array_mut()[[0, 4, 3]] = 4.0;
            x.array_mut()[[0, 7, 1]] = 8.0;
        } else {
            unreachable!();
        }

        assert!(!grid1.is_empty());

        assert_eq!(grid1.iter().nth(0), Some(((0, 1, 2), &1.0)));
        assert_eq!(grid1.iter().nth(1), Some(((0, 1, 3), &2.0)));
        assert_eq!(grid1.iter().nth(2), Some(((0, 4, 3), &4.0)));
        assert_eq!(grid1.iter().nth(3), Some(((0, 7, 1), &8.0)));

        // symmetric luminosity function
        let lumi =
            &mut (|ix1, ix2, _| x[ix1] * x[ix2]) as &mut dyn FnMut(usize, usize, usize) -> f64;

        assert_eq!(grid1.convolute(&x, &x, &mu2, lumi), 0.228515625);

        // create grid with transposed entries, but different q2
        let mut grid2: SubgridEnum = ImportOnlySubgridV2::new(
            SparseArray3::new(1, 10, 10),
            vec![Mu2 { ren: 1.0, fac: 1.0 }],
            x.clone(),
            x.clone(),
        )
        .into();
        if let SubgridEnum::ImportOnlySubgridV2(ref mut x) = grid2 {
            x.array_mut()[[0, 2, 1]] = 1.0;
            x.array_mut()[[0, 3, 1]] = 2.0;
            x.array_mut()[[0, 3, 4]] = 4.0;
            x.array_mut()[[0, 1, 7]] = 8.0;
        } else {
            unreachable!();
        }
        assert_eq!(grid2.convolute(&x, &x, &mu2, lumi), 0.228515625);

        assert_eq!(grid2.iter().nth(0), Some(((0, 1, 7), &8.0)));
        assert_eq!(grid2.iter().nth(1), Some(((0, 2, 1), &1.0)));
        assert_eq!(grid2.iter().nth(2), Some(((0, 3, 1), &2.0)));
        assert_eq!(grid2.iter().nth(3), Some(((0, 3, 4), &4.0)));

        grid1.merge(&mut grid2, false);

        assert_eq!(grid1.convolute(&x, &x, &mu2, lumi), 2.0 * 0.228515625);

        let mut grid1 = {
            let mut g = grid1.clone_empty();
            g.merge(&mut grid1, false);
            g
        };

        // the luminosity function is symmetric, so after symmetrization the result must be
        // unchanged
        grid1.symmetrize();
        assert_eq!(grid1.convolute(&x, &x, &mu2, lumi), 2.0 * 0.228515625);

        grid1.scale(2.0);
        assert_eq!(grid1.convolute(&x, &x, &mu2, lumi), 4.0 * 0.228515625);

        assert_eq!(
            grid1.stats(),
            Stats {
                total: 200,
                allocated: 14,
                zeros: 6,
                overhead: 42,
                bytes_per_value: 8,
            }
        );
    }

    #[test]
    #[should_panic(expected = "ImportOnlySubgridV1 doesn't support the fill operation")]
    fn fill_panic_v1() {
        let mut grid =
            ImportOnlySubgridV1::new(SparseArray3::new(1, 1, 1), vec![1.0], vec![1.0], vec![1.0]);

        grid.fill(&Ntuple {
            x1: 0.0,
            x2: 0.0,
            q2: 0.0,
            weight: 1.0,
        });
    }

    #[test]
    #[should_panic(expected = "ImportOnlySubgridV2 doesn't support the fill operation")]
    fn fill_panic_v2() {
        let mut grid = ImportOnlySubgridV2::new(
            SparseArray3::new(1, 1, 1),
            vec![Mu2 { ren: 1.0, fac: 1.0 }],
            vec![1.0],
            vec![1.0],
        );

        grid.fill(&Ntuple {
            x1: 0.0,
            x2: 0.0,
            q2: 0.0,
            weight: 1.0,
        });
    }

    #[test]
    fn from_lagrange_subgrid_v2() {
        let mut lagrange =
            LagrangeSubgridV2::new(&SubgridParams::default(), &ExtraSubgridParams::default());

        // by default this should have 40 grid points
        assert_eq!(lagrange.mu2_grid().len(), 40);

        // only `q2` are important: they're not static and fall between two grid points
        lagrange.fill(&Ntuple {
            x1: 0.25,
            x2: 0.5,
            q2: 10000.0,
            weight: 1.0,
        });
        lagrange.fill(&Ntuple {
            x1: 0.0625,
            x2: 0.125,
            q2: 10001.0,
            weight: 1.0,
        });
        lagrange.fill(&Ntuple {
            x1: 0.5,
            x2: 0.0625,
            q2: 10002.0,
            weight: 1.0,
        });
        lagrange.fill(&Ntuple {
            x1: 0.1,
            x2: 0.2,
            q2: 10003.0,
            weight: 1.0,
        });

        let x1 = lagrange.x1_grid();
        let x2 = lagrange.x2_grid();
        let mu2 = lagrange.mu2_grid();

        let lumi = &mut (|_, _, _| 1.0) as &mut dyn FnMut(usize, usize, usize) -> f64;
        let reference = lagrange.convolute(&x1, &x2, &mu2, lumi);

        let imported = ImportOnlySubgridV1::from(&lagrange);
        let test = imported.convolute(&x1, &x2, &mu2, lumi);

        // make sure the conversion did not change the results
        assert_approx_eq!(f64, reference, test, ulps = 8);

        // all unneccessary grid points should be gone; since we are inserting between two
        // interpolation grid points, the imported grid should have as many interpolation grid
        // points as its interpolation order
        assert_eq!(imported.mu2_grid().len(), 4);
    }
}
