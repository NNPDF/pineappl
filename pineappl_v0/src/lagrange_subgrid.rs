//! Module containing the Lagrange-interpolation subgrid.

use super::convert::f64_from_usize;
use super::sparse_array3::SparseArray3;
use super::subgrid::{
    ExtraSubgridParams, Mu2, Stats, Subgrid, SubgridEnum, SubgridIndexedIter, SubgridParams,
};
use ndarray::Array3;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::iter;
use std::mem;

fn weightfun(x: f64) -> f64 {
    (x.sqrt() / (1.0 - 0.99 * x)).powi(3)
}

fn fx(y: f64) -> f64 {
    let mut yp = y;

    for _ in 0..100 {
        let x = (-yp).exp();
        let delta = y - yp - 5.0 * (1.0 - x);
        if (delta).abs() < 1e-12 {
            return x;
        }
        let deriv = -1.0 - 5.0 * x;
        yp -= delta / deriv;
    }

    unreachable!();
}

fn fy(x: f64) -> f64 {
    (1.0 - x).mul_add(5.0, -x.ln())
}

fn ftau(q2: f64) -> f64 {
    (q2 / 0.0625).ln().ln()
}

fn fq2(tau: f64) -> f64 {
    0.0625 * tau.exp().exp()
}

/// Subgrid which uses Lagrange-interpolation.
#[derive(Clone, Deserialize, Serialize)]
pub struct LagrangeSubgridV1 {
    grid: Option<Array3<f64>>,
    ntau: usize,
    ny: usize,
    yorder: usize,
    tauorder: usize,
    itaumin: usize,
    itaumax: usize,
    reweight: bool,
    ymin: f64,
    ymax: f64,
    taumin: f64,
    taumax: f64,
}

impl LagrangeSubgridV1 {
    /// Constructor.
    #[must_use]
    pub fn new(subgrid_params: &SubgridParams) -> Self {
        Self {
            grid: None,
            ntau: subgrid_params.q2_bins(),
            ny: subgrid_params.x_bins(),
            yorder: subgrid_params.x_order(),
            tauorder: subgrid_params.q2_order(),
            itaumin: 0,
            itaumax: 0,
            reweight: subgrid_params.reweight(),
            ymin: fy(subgrid_params.x_max()),
            ymax: fy(subgrid_params.x_min()),
            taumin: ftau(subgrid_params.q2_min()),
            taumax: ftau(subgrid_params.q2_max()),
        }
    }

    fn deltay(&self) -> f64 {
        (self.ymax - self.ymin) / f64_from_usize(self.ny - 1)
    }

    fn deltatau(&self) -> f64 {
        (self.taumax - self.taumin) / f64_from_usize(self.ntau - 1)
    }

    fn gety(&self, iy: usize) -> f64 {
        f64_from_usize(iy).mul_add(self.deltay(), self.ymin)
    }

    fn gettau(&self, iy: usize) -> f64 {
        f64_from_usize(iy).mul_add(self.deltatau(), self.taumin)
    }

    fn increase_tau(&mut self, new_itaumin: usize, new_itaumax: usize) {
        let min_diff = self.itaumin - new_itaumin;

        let mut new_grid = Array3::zeros((new_itaumax - new_itaumin, self.ny, self.ny));

        for ((i, j, k), value) in self.grid.as_ref().unwrap().indexed_iter() {
            new_grid[[i + min_diff, j, k]] = *value;
        }

        self.itaumin = new_itaumin;
        self.itaumax = new_itaumax;

        mem::swap(&mut self.grid, &mut Some(new_grid));
    }
}

impl Subgrid for LagrangeSubgridV1 {
    fn mu2_grid(&self) -> Cow<'_, [Mu2]> {
        (0..self.ntau)
            .map(|itau| {
                let q2 = fq2(self.gettau(itau));
                Mu2 { ren: q2, fac: q2 }
            })
            .collect()
    }

    fn x1_grid(&self) -> Cow<'_, [f64]> {
        (0..self.ny).map(|iy| fx(self.gety(iy))).collect()
    }

    fn x2_grid(&self) -> Cow<'_, [f64]> {
        self.x1_grid()
    }

    fn is_empty(&self) -> bool {
        self.grid.is_none()
    }

    fn merge(&mut self, other: &mut SubgridEnum, transpose: bool) {
        let x1_equal = self.x1_grid() == other.x1_grid();
        let x2_equal = self.x2_grid() == other.x2_grid();

        if let SubgridEnum::LagrangeSubgridV1(other_grid) = other {
            if let Some(other_grid_grid) = &mut other_grid.grid {
                if self.grid.is_some() {
                    // TODO: the general case isn't implemented
                    assert!(x1_equal);
                    assert!(x2_equal);

                    let new_itaumin = self.itaumin.min(other_grid.itaumin);
                    let new_itaumax = self.itaumax.max(other_grid.itaumax);
                    let offset = other_grid.itaumin.saturating_sub(self.itaumin);

                    // TODO: we need much more checks here if there subgrids are compatible at all

                    if (self.itaumin != new_itaumin) || (self.itaumax != new_itaumax) {
                        self.increase_tau(new_itaumin, new_itaumax);
                    }

                    let self_grid = self.grid.as_mut().unwrap();

                    if transpose {
                        for ((i, k, j), value) in other_grid_grid.indexed_iter() {
                            self_grid[[i + offset, j, k]] += value;
                        }
                    } else {
                        for ((i, j, k), value) in other_grid_grid.indexed_iter() {
                            self_grid[[i + offset, j, k]] += value;
                        }
                    }
                } else {
                    self.grid = other_grid.grid.take();
                    self.itaumin = other_grid.itaumin;
                    self.itaumax = other_grid.itaumax;

                    if transpose {
                        if let Some(grid) = &mut self.grid {
                            grid.swap_axes(1, 2);
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
            self.grid = None;
        } else if let Some(self_grid) = &mut self.grid {
            self_grid.iter_mut().for_each(|x| *x *= factor);
        }
    }

    fn symmetrize(&mut self) {
        if let Some(grid) = self.grid.as_mut() {
            let (i_size, j_size, k_size) = grid.dim();

            for i in 0..i_size {
                for j in 0..j_size {
                    for k in j + 1..k_size {
                        grid[[i, j, k]] += grid[[i, k, j]];
                        grid[[i, k, j]] = 0.0;
                    }
                }
            }
        }
    }

    fn clone_empty(&self) -> SubgridEnum {
        Self {
            grid: None,
            ntau: self.ntau,
            ny: self.ny,
            yorder: self.yorder,
            tauorder: self.tauorder,
            itaumin: 0,
            itaumax: 0,
            reweight: self.reweight,
            ymin: self.ymin,
            ymax: self.ymax,
            taumin: self.taumin,
            taumax: self.taumax,
        }
        .into()
    }

    fn indexed_iter(&self) -> SubgridIndexedIter<'_> {
        self.grid.as_ref().map_or_else(
            || Box::new(iter::empty()) as Box<dyn Iterator<Item = ((usize, usize, usize), f64)>>,
            |grid| {
                Box::new(grid.indexed_iter().filter(|&(_, &value)| value != 0.0).map(
                    |(tuple, &value)| {
                        (
                            (self.itaumin + tuple.0, tuple.1, tuple.2),
                            value
                                * if self.reweight {
                                    weightfun(fx(self.gety(tuple.1)))
                                        * weightfun(fx(self.gety(tuple.2)))
                                } else {
                                    1.0
                                },
                        )
                    },
                ))
            },
        )
    }

    fn stats(&self) -> Stats {
        let (non_zeros, zeros) = self.grid.as_ref().map_or((0, 0), |array| {
            array.iter().fold((0, 0), |mut result, value| {
                if *value == 0.0 {
                    result.0 += 1;
                } else {
                    result.1 += 1;
                }
                result
            })
        });

        Stats {
            total: non_zeros + zeros,
            allocated: non_zeros + zeros,
            zeros,
            overhead: 0,
            bytes_per_value: mem::size_of::<f64>(),
        }
    }

    fn static_scale(&self) -> Option<Mu2> {
        if let [static_scale] = self.mu2_grid().as_ref() {
            Some(static_scale.clone())
        } else {
            None
        }
    }
}

/// Subgrid which uses Lagrange-interpolation.
#[derive(Clone, Deserialize, Serialize)]
pub struct LagrangeSubgridV2 {
    grid: Option<Array3<f64>>,
    ntau: usize,
    ny1: usize,
    ny2: usize,
    y1order: usize,
    y2order: usize,
    tauorder: usize,
    itaumin: usize,
    itaumax: usize,
    reweight1: bool,
    reweight2: bool,
    y1min: f64,
    y1max: f64,
    y2min: f64,
    y2max: f64,
    taumin: f64,
    taumax: f64,
    pub(crate) static_q2: f64,
}

impl LagrangeSubgridV2 {
    /// Constructor.
    #[must_use]
    pub fn new(subgrid_params: &SubgridParams, extra_params: &ExtraSubgridParams) -> Self {
        Self {
            grid: None,
            ntau: subgrid_params.q2_bins(),
            ny1: subgrid_params.x_bins(),
            ny2: extra_params.x2_bins(),
            y1order: subgrid_params.x_order(),
            y2order: extra_params.x2_order(),
            tauorder: subgrid_params.q2_order(),
            itaumin: 0,
            itaumax: 0,
            reweight1: subgrid_params.reweight(),
            reweight2: extra_params.reweight2(),
            y1min: fy(subgrid_params.x_max()),
            y1max: fy(subgrid_params.x_min()),
            y2min: fy(extra_params.x2_max()),
            y2max: fy(extra_params.x2_min()),
            taumin: ftau(subgrid_params.q2_min()),
            taumax: ftau(subgrid_params.q2_max()),
            static_q2: 0.0,
        }
    }

    fn deltay1(&self) -> f64 {
        (self.y1max - self.y1min) / f64_from_usize(self.ny1 - 1)
    }

    fn deltay2(&self) -> f64 {
        (self.y1max - self.y2min) / f64_from_usize(self.ny2 - 1)
    }

    fn deltatau(&self) -> f64 {
        (self.taumax - self.taumin) / f64_from_usize(self.ntau - 1)
    }

    fn gety1(&self, iy: usize) -> f64 {
        if self.y1min == self.y1max {
            debug_assert_eq!(iy, 0);
            self.y1min
        } else {
            f64_from_usize(iy).mul_add(self.deltay1(), self.y1min)
        }
    }

    fn gety2(&self, iy: usize) -> f64 {
        if self.y2min == self.y2max {
            debug_assert_eq!(iy, 0);
            self.y2min
        } else {
            f64_from_usize(iy).mul_add(self.deltay2(), self.y2min)
        }
    }

    fn gettau(&self, iy: usize) -> f64 {
        if self.taumin == self.taumax {
            debug_assert_eq!(iy, 0);
            self.taumin
        } else {
            f64_from_usize(iy).mul_add(self.deltatau(), self.taumin)
        }
    }

    fn increase_tau(&mut self, new_itaumin: usize, new_itaumax: usize) {
        let min_diff = self.itaumin - new_itaumin;

        let mut new_grid = Array3::zeros((new_itaumax - new_itaumin, self.ny1, self.ny2));

        for ((i, j, k), value) in self.grid.as_ref().unwrap().indexed_iter() {
            new_grid[[i + min_diff, j, k]] = *value;
        }

        self.itaumin = new_itaumin;
        self.itaumax = new_itaumax;

        mem::swap(&mut self.grid, &mut Some(new_grid));
    }
}

impl Subgrid for LagrangeSubgridV2 {
    fn mu2_grid(&self) -> Cow<'_, [Mu2]> {
        (0..self.ntau)
            .map(|itau| {
                let q2 = fq2(self.gettau(itau));
                Mu2 { ren: q2, fac: q2 }
            })
            .collect()
    }

    fn x1_grid(&self) -> Cow<'_, [f64]> {
        (0..self.ny1).map(|iy| fx(self.gety1(iy))).collect()
    }

    fn x2_grid(&self) -> Cow<'_, [f64]> {
        (0..self.ny2).map(|iy| fx(self.gety2(iy))).collect()
    }

    fn is_empty(&self) -> bool {
        self.grid.is_none()
    }

    fn merge(&mut self, other: &mut SubgridEnum, transpose: bool) {
        let x1_equal = self.x1_grid() == other.x1_grid();
        let x2_equal = self.x2_grid() == other.x2_grid();

        if let SubgridEnum::LagrangeSubgridV2(other_grid) = other {
            if let Some(other_grid_grid) = &mut other_grid.grid {
                if self.grid.is_some() {
                    // TODO: the general case isn't implemented
                    assert!(x1_equal);
                    assert!(x2_equal);

                    let new_itaumin = self.itaumin.min(other_grid.itaumin);
                    let new_itaumax = self.itaumax.max(other_grid.itaumax);
                    let offset = other_grid.itaumin.saturating_sub(self.itaumin);

                    // TODO: we need much more checks here if there subgrids are compatible at all

                    if (self.itaumin != new_itaumin) || (self.itaumax != new_itaumax) {
                        self.increase_tau(new_itaumin, new_itaumax);
                    }

                    if (other_grid.static_q2 == -1.0) || (self.static_q2 != other_grid.static_q2) {
                        self.static_q2 = -1.0;
                    }

                    let self_grid = self.grid.as_mut().unwrap();

                    if transpose {
                        for ((i, k, j), value) in other_grid_grid.indexed_iter() {
                            self_grid[[i + offset, j, k]] += value;
                        }
                    } else {
                        for ((i, j, k), value) in other_grid_grid.indexed_iter() {
                            self_grid[[i + offset, j, k]] += value;
                        }
                    }
                } else {
                    self.grid = other_grid.grid.take();
                    self.itaumin = other_grid.itaumin;
                    self.itaumax = other_grid.itaumax;
                    self.static_q2 = other_grid.static_q2;

                    if transpose {
                        if let Some(grid) = &mut self.grid {
                            grid.swap_axes(1, 2);
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
            self.grid = None;
        } else if let Some(self_grid) = &mut self.grid {
            self_grid.iter_mut().for_each(|x| *x *= factor);
        }
    }

    fn symmetrize(&mut self) {
        if let Some(grid) = self.grid.as_mut() {
            let (i_size, j_size, k_size) = grid.dim();

            for i in 0..i_size {
                for j in 0..j_size {
                    for k in j + 1..k_size {
                        grid[[i, j, k]] += grid[[i, k, j]];
                        grid[[i, k, j]] = 0.0;
                    }
                }
            }
        }
    }

    fn clone_empty(&self) -> SubgridEnum {
        Self {
            grid: None,
            ntau: self.ntau,
            ny1: self.ny1,
            ny2: self.ny2,
            y1order: self.y1order,
            y2order: self.y2order,
            tauorder: self.tauorder,
            itaumin: 0,
            itaumax: 0,
            reweight1: self.reweight1,
            reweight2: self.reweight2,
            y1min: self.y1min,
            y1max: self.y1max,
            y2min: self.y2min,
            y2max: self.y2max,
            taumin: self.taumin,
            taumax: self.taumax,
            static_q2: 0.0,
        }
        .into()
    }

    fn indexed_iter(&self) -> SubgridIndexedIter<'_> {
        self.grid.as_ref().map_or_else(
            || Box::new(iter::empty()) as Box<dyn Iterator<Item = ((usize, usize, usize), f64)>>,
            |grid| {
                Box::new(grid.indexed_iter().filter(|&(_, &value)| value != 0.0).map(
                    |(tuple, &value)| {
                        (
                            (self.itaumin + tuple.0, tuple.1, tuple.2),
                            value
                                * if self.reweight1 {
                                    weightfun(fx(self.gety1(tuple.1)))
                                } else {
                                    1.0
                                }
                                * if self.reweight2 {
                                    weightfun(fx(self.gety2(tuple.2)))
                                } else {
                                    1.0
                                },
                        )
                    },
                ))
            },
        )
    }

    fn stats(&self) -> Stats {
        let (non_zeros, zeros) = self.grid.as_ref().map_or((0, 0), |array| {
            array.iter().fold((0, 0), |mut result, value| {
                if *value == 0.0 {
                    result.0 += 1;
                } else {
                    result.1 += 1;
                }
                result
            })
        });

        Stats {
            total: non_zeros + zeros,
            allocated: non_zeros + zeros,
            zeros,
            overhead: 0,
            bytes_per_value: mem::size_of::<f64>(),
        }
    }

    fn static_scale(&self) -> Option<Mu2> {
        (self.static_q2 > 0.0).then_some(Mu2 {
            ren: self.static_q2,
            fac: self.static_q2,
        })
    }
}

/// Subgrid which uses Lagrange-interpolation, but also stores its contents in a space-efficient
/// structure.
#[derive(Clone, Deserialize, Serialize)]
pub struct LagrangeSparseSubgridV1 {
    array: SparseArray3<f64>,
    ntau: usize,
    ny: usize,
    yorder: usize,
    tauorder: usize,
    reweight: bool,
    ymin: f64,
    ymax: f64,
    taumin: f64,
    taumax: f64,
}

impl LagrangeSparseSubgridV1 {
    /// Constructor.
    #[must_use]
    pub fn new(subgrid_params: &SubgridParams) -> Self {
        Self {
            array: SparseArray3::new(
                subgrid_params.q2_bins(),
                subgrid_params.x_bins(),
                subgrid_params.x_bins(),
            ),
            ntau: subgrid_params.q2_bins(),
            ny: subgrid_params.x_bins(),
            yorder: subgrid_params.x_order(),
            tauorder: subgrid_params.q2_order(),
            reweight: subgrid_params.reweight(),
            ymin: fy(subgrid_params.x_max()),
            ymax: fy(subgrid_params.x_min()),
            taumin: ftau(subgrid_params.q2_min()),
            taumax: ftau(subgrid_params.q2_max()),
        }
    }

    fn deltay(&self) -> f64 {
        (self.ymax - self.ymin) / f64_from_usize(self.ny - 1)
    }

    fn deltatau(&self) -> f64 {
        (self.taumax - self.taumin) / f64_from_usize(self.ntau - 1)
    }

    fn gety(&self, iy: usize) -> f64 {
        f64_from_usize(iy).mul_add(self.deltay(), self.ymin)
    }

    fn gettau(&self, iy: usize) -> f64 {
        f64_from_usize(iy).mul_add(self.deltatau(), self.taumin)
    }
}

impl Subgrid for LagrangeSparseSubgridV1 {
    fn mu2_grid(&self) -> Cow<'_, [Mu2]> {
        (0..self.ntau)
            .map(|itau| {
                let q2 = fq2(self.gettau(itau));
                Mu2 { ren: q2, fac: q2 }
            })
            .collect()
    }

    fn x1_grid(&self) -> Cow<'_, [f64]> {
        (0..self.ny).map(|iy| fx(self.gety(iy))).collect()
    }

    fn x2_grid(&self) -> Cow<'_, [f64]> {
        self.x1_grid()
    }

    fn is_empty(&self) -> bool {
        self.array.is_empty()
    }

    fn merge(&mut self, other: &mut SubgridEnum, transpose: bool) {
        if let SubgridEnum::LagrangeSparseSubgridV1(other_grid) = other {
            if self.array.is_empty() && !transpose {
                mem::swap(&mut self.array, &mut other_grid.array);
            } else {
                // TODO: the general case isn't implemented
                assert!(self.x1_grid() == other_grid.x1_grid());
                assert!(self.x2_grid() == other_grid.x2_grid());

                // TODO: we need much more checks here if there subgrids are compatible at all

                if transpose {
                    for ((i, k, j), value) in other_grid.array.indexed_iter() {
                        self.array[[i, j, k]] += value;
                    }
                } else {
                    for ((i, j, k), value) in other_grid.array.indexed_iter() {
                        self.array[[i, j, k]] += value;
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
        let mut new_array = SparseArray3::new(self.ntau, self.ny, self.ny);

        for ((i, j, k), sigma) in self.array.indexed_iter().filter(|((_, j, k), _)| k >= j) {
            new_array[[i, j, k]] = sigma;
        }
        for ((i, j, k), sigma) in self.array.indexed_iter().filter(|((_, j, k), _)| k < j) {
            new_array[[i, k, j]] += sigma;
        }

        mem::swap(&mut self.array, &mut new_array);
    }

    fn clone_empty(&self) -> SubgridEnum {
        Self {
            array: SparseArray3::new(self.ntau, self.ny, self.ny),
            ntau: self.ntau,
            ny: self.ny,
            yorder: self.yorder,
            tauorder: self.tauorder,
            reweight: self.reweight,
            ymin: self.ymin,
            ymax: self.ymax,
            taumin: self.taumin,
            taumax: self.taumax,
        }
        .into()
    }

    fn indexed_iter(&self) -> SubgridIndexedIter<'_> {
        Box::new(self.array.indexed_iter().map(|(tuple, value)| {
            (
                tuple,
                value
                    * if self.reweight {
                        weightfun(fx(self.gety(tuple.1))) * weightfun(fx(self.gety(tuple.2)))
                    } else {
                        1.0
                    },
            )
        }))
    }

    fn stats(&self) -> Stats {
        Stats {
            total: self.ntau * self.ny * self.ny,
            allocated: self.array.len() + self.array.zeros(),
            zeros: self.array.zeros(),
            overhead: self.array.overhead(),
            bytes_per_value: mem::size_of::<f64>(),
        }
    }

    fn static_scale(&self) -> Option<Mu2> {
        if let [static_scale] = self.mu2_grid().as_ref() {
            Some(static_scale.clone())
        } else {
            None
        }
    }
}

impl From<&LagrangeSubgridV1> for LagrangeSparseSubgridV1 {
    fn from(subgrid: &LagrangeSubgridV1) -> Self {
        Self {
            array: subgrid.grid.as_ref().map_or_else(
                || SparseArray3::new(subgrid.ntau, subgrid.ny, subgrid.ny),
                |grid| SparseArray3::from_ndarray(grid.view(), subgrid.itaumin, subgrid.ntau),
            ),
            ntau: subgrid.ntau,
            ny: subgrid.ny,
            yorder: subgrid.yorder,
            tauorder: subgrid.tauorder,
            reweight: subgrid.reweight,
            ymin: subgrid.ymin,
            ymax: subgrid.ymax,
            taumin: subgrid.taumin,
            taumax: subgrid.taumax,
        }
    }
}
