//! Module containing the Lagrange-interpolation subgrid.

use super::grid::{Ntuple, Subgrid, SubgridParams};
use arrayvec::ArrayVec;
use itertools::iproduct;
use ndarray::Array3;
use serde::{Deserialize, Serialize};
use std::any::Any;
use std::mem;

fn weightfun(x: f64) -> f64 {
    (x.sqrt() / (1.0 - 0.99 * x)).powi(3)
}

fn fx(y: f64) -> f64 {
    let eps = 1e-12;
    let imax = 100;

    let mut yp = y;

    for _ in 0..imax {
        let x = (-yp).exp();
        let delta = y - yp - 5.0 * (1.0 - x);
        if (delta).abs() < eps {
            return x;
        }
        let deriv = -1.0 - 5.0 * x;
        yp -= delta / deriv;
    }

    panic!("");
}

fn fy(x: f64) -> f64 {
    -x.ln() + 5.0 * (1.0 - x)
}

fn ftau(q2: f64) -> f64 {
    (q2 / 0.0625).ln().ln()
}

fn fq2(tau: f64) -> f64 {
    0.0625 * tau.exp().exp()
}

fn fi(i: usize, n: usize, u: f64) -> f64 {
    let mut factorials = 1;
    let mut product = if (n - i) % 2 == 0 { 1.0 } else { -1.0 };
    for z in 0..i {
        product *= u - (z as f64);
        factorials *= z + 1;
    }
    for z in i + 1..=n {
        product *= u - (z as f64);
        factorials *= z - i;
    }
    product / (factorials as f64)
}

/// Subgrid which uses Lagrange-interpolation.
#[derive(Deserialize, Serialize)]
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
        (self.ymax - self.ymin) / ((self.ny - 1) as f64)
    }

    fn deltatau(&self) -> f64 {
        (self.taumax - self.taumin) / ((self.ntau - 1) as f64)
    }

    fn gety(&self, iy: usize) -> f64 {
        (iy as f64).mul_add(self.deltay(), self.ymin)
    }

    fn gettau(&self, iy: usize) -> f64 {
        (iy as f64).mul_add(self.deltatau(), self.taumin)
    }
}

impl LagrangeSubgridV1 {
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
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn convolute(&self, lumi: &dyn Fn(f64, f64, f64) -> f64) -> f64 {
        if let Some(self_grid) = &self.grid {
            let qs: Vec<_> = (self.itaumin..self.itaumax)
                .map(|itau| fq2(self.gettau(itau)))
                .collect();
            let x: Vec<_> = (0..self.ny).map(|iy| fx(self.gety(iy))).collect();

            self_grid
                .iter()
                .zip(iproduct!(qs, &x, &x))
                .map(|(&sigma, (q2, &x1, &x2))| {
                    if sigma != 0.0 {
                        let mut value = sigma * lumi(x1, x2, q2);
                        if self.reweight {
                            value *= weightfun(x1) * weightfun(x2);
                        }
                        value
                    } else {
                        0.0
                    }
                })
                .sum()
        } else {
            0.0
        }
    }

    fn fill(&mut self, ntuple: &Ntuple<f64>) {
        let y1 = fy(ntuple.x1);
        let y2 = fy(ntuple.x2);
        let tau = ftau(ntuple.q2);

        if (y2 < self.ymin)
            || (y2 > self.ymax)
            || (y1 < self.ymin)
            || (y1 > self.ymax)
            || (tau < self.taumin)
            || (tau > self.taumax)
        {
            return;
        }

        let k1 = (((y1 - self.ymin) / self.deltay() - ((self.yorder / 2) as f64)).max(0.0)
            as usize)
            .min(self.ny - 1 - self.yorder);

        let k2 = (((y2 - self.ymin) / self.deltay() - ((self.yorder / 2) as f64)).max(0.0)
            as usize)
            .min(self.ny - 1 - self.yorder);

        let u_y1 = (y1 - self.gety(k1)) / self.deltay();
        let u_y2 = (y2 - self.gety(k2)) / self.deltay();

        let fi1: ArrayVec<[_; 8]> = (0..=self.yorder)
            .map(|i| fi(i, self.yorder, u_y1))
            .collect();
        let fi2: ArrayVec<[_; 8]> = (0..=self.yorder)
            .map(|i| fi(i, self.yorder, u_y2))
            .collect();

        let k3 = (((tau - self.taumin) / self.deltatau() - ((self.tauorder / 2) as f64)).max(0.0)
            as usize)
            .min(self.ntau - 1 - self.tauorder);

        let u_tau = (tau - self.gettau(k3)) / self.deltatau();

        let factor = if self.reweight {
            1.0 / (weightfun(ntuple.x1) * weightfun(ntuple.x2))
        } else {
            1.0
        };

        let size = self.tauorder + 1;
        let ny = self.ny;

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
                    let fillweight = factor * fi1i1 * fi2i2 * fi3i3 * ntuple.weight;

                    let grid = self
                        .grid
                        .get_or_insert_with(|| Array3::zeros((size, ny, ny)));

                    grid[[k3 + i3 - self.itaumin, k1 + i1, k2 + i2]] += fillweight;
                }
            }
        }
    }

    fn is_empty(&self) -> bool {
        self.grid.is_none()
    }

    fn merge(&mut self, other: &mut dyn Subgrid) {
        if let Some(other_grid) = other.as_any_mut().downcast_mut::<Self>() {
            if let Some(other_grid_grid) = &mut other_grid.grid {
                if self.grid.is_some() {
                    let new_itaumin = self.itaumin.min(other_grid.itaumin);
                    let new_itaumax = self.itaumax.max(other_grid.itaumax);
                    let offset = other_grid.itaumin.saturating_sub(self.itaumin);

                    // TODO: we need much more checks here if there subgrids are compatible at all

                    if (self.itaumin != new_itaumin) || (self.itaumax != new_itaumax) {
                        self.increase_tau(new_itaumin, new_itaumax);
                    }

                    let self_grid = self.grid.as_mut().unwrap();

                    for ((i, j, k), value) in other_grid_grid.indexed_iter() {
                        self_grid[[i + offset, j, k]] += value;
                    }
                } else {
                    self.grid = other_grid.grid.take();
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
}
