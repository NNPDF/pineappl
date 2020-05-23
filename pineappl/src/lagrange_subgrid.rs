//! Module containing the Lagrange-interpolation subgrid.

use super::grid::{Ntuple, Subgrid, SubgridParams};
use arrayvec::ArrayVec;
use ndarray::Array3;
use serde::{Deserialize, Serialize};
use std::any::Any;

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

fn fac(i: usize) -> f64 {
    if i == 0 {
        1.0
    } else {
        (i as f64) * fac(i - 1)
    }
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
    if n == 0 && i == 0 {
        return 1.0;
    }
    if (u - i as f64).abs() < 1e-8 {
        return 1.0;
    }
    let mut product = (-1.0_f64).powi((n - i) as i32) / (fac(i) * fac(n - i) * (u - i as f64));
    for z in 0..=n {
        product *= u - (z as f64);
    }
    product
}

/// Subgrid which uses Lagrange-interpolation.
#[derive(Deserialize, Serialize)]
pub struct LagrangeSubgrid {
    grid: Option<Array3<f64>>,
    ntau: usize,
    ny1: usize,
    ny2: usize,
    yorder: usize,
    tauorder: usize,
    reweight: bool,
    y1min: f64,
    y1max: f64,
    y2min: f64,
    y2max: f64,
    taumin: f64,
    taumax: f64,
}

impl LagrangeSubgrid {
    /// Constructor.
    pub fn new(subgrid_params: &SubgridParams) -> Self {
        Self {
            grid: None,
            ntau: subgrid_params.q2_bins(),
            ny1: subgrid_params.x_bins(),
            ny2: subgrid_params.x_bins(),
            yorder: subgrid_params.x_order(),
            tauorder: subgrid_params.q2_order(),
            reweight: subgrid_params.reweight(),
            y1min: fy(subgrid_params.x_max()),
            y1max: fy(subgrid_params.x_min()),
            y2min: fy(subgrid_params.x_max()),
            y2max: fy(subgrid_params.x_min()),
            taumin: ftau(subgrid_params.q2_min()),
            taumax: ftau(subgrid_params.q2_max()),
        }
    }

    fn deltay1(&self) -> f64 {
        (self.y1max - self.y1min) / ((self.ny1 - 1) as f64)
    }

    fn deltay2(&self) -> f64 {
        (self.y2max - self.y2min) / ((self.ny2 - 1) as f64)
    }

    fn deltatau(&self) -> f64 {
        (self.taumax - self.taumin) / ((self.ntau - 1) as f64)
    }

    fn gety1(&self, iy: usize) -> f64 {
        (iy as f64).mul_add(self.deltay1(), self.y1min)
    }

    fn gety2(&self, iy: usize) -> f64 {
        (iy as f64).mul_add(self.deltay2(), self.y2min)
    }

    fn gettau(&self, iy: usize) -> f64 {
        (iy as f64).mul_add(self.deltatau(), self.taumin)
    }

    fn fk1(&self, x: f64) -> Option<usize> {
        let y = fy(x);
        if (y < self.y1min) || (y > self.y1max) {
            return None;
        }
        let mut k = ((y - self.y1min) / self.deltay1() - ((self.yorder >> 1) as f64)) as i32;
        if k < 0 {
            k = 0;
        }
        let mut k = k as usize;
        if k + self.yorder >= self.ny1 {
            k = self.ny1 - 1 - self.yorder;
        }

        Some(k)
    }

    fn fk2(&self, x: f64) -> Option<usize> {
        let y = fy(x);
        if (y < self.y2min) || (y > self.y2max) {
            return None;
        }
        let mut k = ((y - self.y2min) / self.deltay2() - ((self.yorder >> 1) as f64)) as i32;
        if k < 0 {
            k = 0;
        }
        let mut k = k as usize;
        if k + self.yorder >= self.ny2 {
            k = self.ny2 - 1 - self.yorder;
        }

        Some(k)
    }

    fn fkappa(&self, q2: f64) -> Option<usize> {
        let tau = ftau(q2);
        if (tau < self.taumin) || (tau > self.taumax) {
            return None;
        }
        let mut kappa =
            ((tau - self.taumin) / self.deltatau() - ((self.tauorder >> 1) as f64)) as i32;
        if kappa + (self.tauorder as i32) >= (self.ntau as i32) {
            kappa = (self.ntau - 1 - self.tauorder) as i32;
        }
        if kappa < 0 {
            kappa = 0;
        }

        Some(kappa as usize)
    }
}

#[typetag::serde]
impl Subgrid for LagrangeSubgrid {
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn convolute(&self, lumi: &dyn Fn(f64, f64, f64) -> f64) -> f64 {
        if let Some(self_grid) = &self.grid {
            let mut dsigma = 0.0;

            for ((itau, iy1, iy2), &sigma) in self_grid.indexed_iter() {
                if sigma != 0.0 {
                    let tau = self.gettau(itau);
                    let q2 = fq2(tau);
                    let y1 = self.gety1(iy1);
                    let x1 = fx(y1);
                    let y2 = self.gety2(iy2);
                    let x2 = fx(y2);
                    let mut lumi_val = lumi(x1, x2, q2);

                    if self.reweight {
                        lumi_val *= weightfun(x1) * weightfun(x2);
                    }

                    dsigma += sigma * lumi_val;
                }
            }

            dsigma
        } else {
            0.0
        }
    }

    fn fill(&mut self, ntuple: &Ntuple<f64>) {
        let k1 = match self.fk1(ntuple.x1) {
            Some(k) => k,
            None => return,
        };
        let k2 = match self.fk2(ntuple.x2) {
            Some(k) => k,
            None => return,
        };
        let k3 = match self.fkappa(ntuple.q2) {
            Some(k) => k,
            None => return,
        };

        let u_y1 = (fy(ntuple.x1) - self.gety1(k1)) / self.deltay1();
        let u_y2 = (fy(ntuple.x2) - self.gety2(k2)) / self.deltay2();

        let fi1: ArrayVec<[_; 8]> = (0..=self.yorder)
            .map(|i| fi(i, self.yorder, u_y1))
            .collect();
        let fi2: ArrayVec<[_; 8]> = (0..=self.yorder)
            .map(|i| fi(i, self.yorder, u_y2))
            .collect();

        let invwfun = 1.0 / (weightfun(ntuple.x1) * weightfun(ntuple.x2));

        let ntau = self.ntau;
        let ny1 = self.ny1;
        let ny2 = self.ny2;

        let u_tau = (ftau(ntuple.q2) - self.gettau(k3)) / self.deltatau();

        for i3 in 0..=self.tauorder {
            let fi3i3 = fi(i3, self.tauorder, u_tau);

            for (i1, fi1i1) in fi1.iter().enumerate() {
                for (i2, fi2i2) in fi2.iter().enumerate() {
                    let mut fi_factor = fi1i1 * fi2i2 * fi3i3;

                    if self.reweight {
                        fi_factor *= invwfun;
                    }

                    let fillweight = ntuple.weight * fi_factor;
                    let grid = self
                        .grid
                        .get_or_insert_with(|| Array3::zeros((ntau, ny1, ny2)));
                    grid[[k3 + i3, k1 + i1, k2 + i2]] += fillweight;
                }
            }
        }
    }

    fn is_empty(&self) -> bool {
        self.grid.is_none()
    }

    fn merge(&mut self, other: &mut dyn Subgrid) {
        if let Some(other_grid) = other.as_any_mut().downcast_mut::<Self>() {
            if let Some(self_grid) = &mut self.grid {
                if let Some(other_grid_grid) = &mut other_grid.grid {
                    self_grid.scaled_add(1.0, other_grid_grid);
                }
            } else {
                self.grid = other_grid.grid.take();
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
