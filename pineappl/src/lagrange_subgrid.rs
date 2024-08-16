//! Module containing the Lagrange-interpolation subgrid.

use super::convert::{f64_from_usize, usize_from_f64};
use super::subgrid::{
    ExtraSubgridParams, Mu2, Stats, Subgrid, SubgridEnum, SubgridIndexedIter, SubgridParams,
};
use arrayvec::ArrayVec;
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

fn fi(i: usize, n: usize, u: f64) -> f64 {
    let mut factorials = 1;
    let mut product = 1.0;
    for z in 0..i {
        product *= u - f64_from_usize(z);
        factorials *= i - z;
    }
    for z in i + 1..=n {
        product *= f64_from_usize(z) - u;
        factorials *= z - i;
    }
    product / f64_from_usize(factorials)
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
    fn convolve(
        &self,
        x1: &[f64],
        x2: &[f64],
        _: &[Mu2],
        lumi: &mut dyn FnMut(usize, usize, usize) -> f64,
    ) -> f64 {
        self.grid.as_ref().map_or(0.0, |grid| {
            grid.indexed_iter()
                .map(|((imu2, ix1, ix2), &sigma)| {
                    if sigma == 0.0 {
                        0.0
                    } else {
                        let mut value = sigma * lumi(ix1, ix2, imu2 + self.itaumin);
                        if self.reweight1 {
                            value *= weightfun(x1[ix1]);
                        }
                        if self.reweight2 {
                            value *= weightfun(x2[ix2]);
                        }
                        value
                    }
                })
                .sum()
        })
    }

    fn fill(&mut self, ntuple: &[f64], weight: f64) {
        if weight == 0.0 {
            return;
        }

        let &[x1, x2, q2] = ntuple else {
            unreachable!();
        };
        let y1 = fy(x1);
        let y2 = fy(x2);
        let tau = ftau(q2);

        if self.static_q2 == 0.0 {
            self.static_q2 = q2;
        } else if (self.static_q2 != -1.0) && (self.static_q2 != q2) {
            self.static_q2 = -1.0;
        }

        if (y2 < self.y2min)
            || (y2 > self.y2max)
            || (y1 < self.y1min)
            || (y1 > self.y1max)
            || (tau < self.taumin)
            || (tau > self.taumax)
        {
            return;
        }

        let k1 =
            usize_from_f64((y1 - self.y1min) / self.deltay1() - f64_from_usize(self.y1order / 2))
                .min(self.ny1 - 1 - self.y1order);
        let k2 =
            usize_from_f64((y2 - self.y2min) / self.deltay2() - f64_from_usize(self.y2order / 2))
                .min(self.ny2 - 1 - self.y2order);

        let u_y1 = (y1 - self.gety1(k1)) / self.deltay1();
        let u_y2 = (y2 - self.gety2(k2)) / self.deltay2();

        let fi1: ArrayVec<_, 8> = (0..=self.y1order)
            .map(|i| fi(i, self.y1order, u_y1))
            .collect();
        let fi2: ArrayVec<_, 8> = (0..=self.y2order)
            .map(|i| fi(i, self.y2order, u_y2))
            .collect();

        let k3 = usize_from_f64(
            (tau - self.taumin) / self.deltatau() - f64_from_usize(self.tauorder / 2),
        )
        .min(self.ntau - 1 - self.tauorder);

        let u_tau = (tau - self.gettau(k3)) / self.deltatau();

        let factor = 1.0
            / (if self.reweight1 { weightfun(x1) } else { 1.0 }
                * if self.reweight2 { weightfun(x2) } else { 1.0 });

        let size = self.tauorder + 1;
        let ny1 = self.ny1;
        let ny2 = self.ny2;

        if self.grid.is_none() {
            self.itaumin = k3;
            self.itaumax = k3 + size;
        } else if k3 < self.itaumin || k3 + size > self.itaumax {
            self.increase_tau(self.itaumin.min(k3), self.itaumax.max(k3 + size));
        }

        for i3 in 0..=self.tauorder {
            let fi3i3 = fi(i3, self.tauorder, u_tau);

            for (i1, fi1i1) in fi1.iter().enumerate() {
                for (i2, fi2i2) in fi2.iter().enumerate() {
                    let fillweight = factor * fi1i1 * fi2i2 * fi3i3 * weight;

                    let grid = self
                        .grid
                        .get_or_insert_with(|| Array3::zeros((size, ny1, ny2)));

                    grid[[k3 + i3 - self.itaumin, k1 + i1, k2 + i2]] += fillweight;
                }
            }
        }
    }

    fn mu2_grid(&self) -> Cow<[Mu2]> {
        (0..self.ntau)
            .map(|itau| {
                let q2 = fq2(self.gettau(itau));
                Mu2 {
                    ren: q2,
                    fac: q2,
                    frg: -1.0,
                }
            })
            .collect()
    }

    fn x1_grid(&self) -> Cow<[f64]> {
        (0..self.ny1).map(|iy| fx(self.gety1(iy))).collect()
    }

    fn x2_grid(&self) -> Cow<[f64]> {
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

    fn indexed_iter(&self) -> SubgridIndexedIter {
        self.grid.as_ref().map_or_else(
            || Box::new(iter::empty()) as Box<dyn Iterator<Item = ((usize, usize, usize), f64)>>,
            |grid| {
                Box::new(grid.indexed_iter().filter(|(_, &value)| value != 0.0).map(
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
            frg: -1.0,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::assert_approx_eq;

    fn test_q2_slice_methods<G: Subgrid>(mut grid: G) -> G {
        grid.fill(&[0.1, 0.2, 90.0_f64.powi(2)], 1.0);
        grid.fill(&[0.9, 0.1, 90.0_f64.powi(2)], 1.0);
        grid.fill(&[0.009, 0.01, 90.0_f64.powi(2)], 1.0);
        grid.fill(&[0.009, 0.5, 90.0_f64.powi(2)], 1.0);

        // the grid must not be empty
        assert!(!grid.is_empty());

        let x1 = grid.x1_grid();
        let x2 = grid.x2_grid();
        let mu2 = grid.mu2_grid();

        let reference = grid.convolve(&x1, &x2, &mu2, &mut |ix1, ix2, _| 1.0 / (x1[ix1] * x2[ix2]));

        let mut test = 0.0;

        // check `reference` against manually calculated result from q2 slices
        for ((_, ix1, ix2), value) in grid.indexed_iter() {
            test += value / (x1[ix1] * x2[ix2]);
        }

        assert_approx_eq!(f64, test, reference, ulps = 8);

        grid
    }

    fn test_merge_method<G: Subgrid>(mut grid1: G, mut grid2: G, mut grid3: G)
    where
        SubgridEnum: From<G>,
    {
        grid1.fill(&[0.1, 0.2, 90.0_f64.powi(2)], 1.0);
        grid1.fill(&[0.9, 0.1, 90.0_f64.powi(2)], 1.0);
        grid1.fill(&[0.009, 0.01, 90.0_f64.powi(2)], 1.0);
        grid1.fill(&[0.009, 0.5, 90.0_f64.powi(2)], 1.0);

        assert!(!grid1.is_empty());
        assert!(grid2.is_empty());

        let x1 = grid1.x1_grid().into_owned();
        let x2 = grid1.x2_grid().into_owned();
        let mu2 = grid1.mu2_grid().into_owned();

        let reference =
            grid1.convolve(&x1, &x2, &mu2, &mut |ix1, ix2, _| 1.0 / (x1[ix1] * x2[ix2]));

        // merge filled grid into empty one
        grid2.merge(&mut grid1.into(), false);
        assert!(!grid2.is_empty());

        let merged = grid2.convolve(&x1, &x2, &mu2, &mut |ix1, ix2, _| 1.0 / (x1[ix1] * x2[ix2]));

        assert_approx_eq!(f64, reference, merged, ulps = 8);

        grid3.fill(&[0.1, 0.2, 90.0_f64.powi(2)], 1.0);
        grid3.fill(&[0.9, 0.1, 90.0_f64.powi(2)], 1.0);
        grid3.fill(&[0.009, 0.01, 90.0_f64.powi(2)], 1.0);
        grid3.fill(&[0.009, 0.5, 90.0_f64.powi(2)], 1.0);

        grid2.merge(&mut grid3.into(), false);

        let merged = grid2.convolve(&x1, &x2, &mu2, &mut |ix1, ix2, _| 1.0 / (x1[ix1] * x2[ix2]));

        assert_approx_eq!(f64, 2.0 * reference, merged, ulps = 8);
    }

    fn test_empty_subgrid<G: Subgrid>(mut grid: G) {
        // this following events should be skipped

        // q2 is too large
        grid.fill(&[0.5, 0.5, 2e+8], 1.0);
        // q2 is too small
        grid.fill(&[0.5, 0.5, 5e+1], 1.0);
        // x1 is too large
        grid.fill(&[1.1, 0.5, 1e+3], 1.0);
        // x1 is too small
        grid.fill(&[0.5, 1e-7, 1e+3], 1.0);
        // x1 is too large
        grid.fill(&[0.5, 1.1, 1e+3], 1.0);
        // x1 is too small
        grid.fill(&[1e-7, 0.5, 1e+3], 1.0);

        let x1 = grid.x1_grid();
        let x2 = grid.x2_grid();
        let mu2 = grid.mu2_grid();

        let result = grid.convolve(&x1, &x2, &mu2, &mut |_, _, _| 1.0);

        assert_eq!(result, 0.0);
    }

    #[test]
    fn q2_slice_v2() {
        let subgrid = test_q2_slice_methods(LagrangeSubgridV2::new(
            &SubgridParams::default(),
            &ExtraSubgridParams::default(),
        ));

        assert_eq!(
            subgrid.stats(),
            Stats {
                total: 10000,
                allocated: 10000,
                zeros: 256,
                overhead: 0,
                bytes_per_value: 8
            }
        );
    }

    #[test]
    fn fill_zero_v2() {
        let mut subgrid =
            LagrangeSubgridV2::new(&SubgridParams::default(), &ExtraSubgridParams::default());

        subgrid.fill(&[0.5, 0.5, 1000.0], 0.0);

        assert!(subgrid.is_empty());
        assert_eq!(subgrid.indexed_iter().count(), 0);
    }

    #[test]
    fn merge_dense_v2() {
        test_merge_method(
            LagrangeSubgridV2::new(&SubgridParams::default(), &ExtraSubgridParams::default()),
            LagrangeSubgridV2::new(&SubgridParams::default(), &ExtraSubgridParams::default()),
            LagrangeSubgridV2::new(&SubgridParams::default(), &ExtraSubgridParams::default()),
        );
    }

    #[test]
    fn empty_v2() {
        test_empty_subgrid(LagrangeSubgridV2::new(
            &SubgridParams::default(),
            &ExtraSubgridParams::default(),
        ));
    }
}
