//! Module containing the Lagrange-interpolation subgrid.

use super::grid::{Ntuple, Subgrid};
use ndarray::{Array3, Dimension};
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
    grid: Array3<f64>,
    m_yorder: usize,
    m_tauorder: usize,
    m_reweight: bool,
    xmin: f64,
    xmax: f64,
    q2min: f64,
    q2max: f64,
}

impl Default for LagrangeSubgrid {
    fn default() -> Self {
        let nq2 = 30;
        let nx = 50;
        let xorder = 3;
        let q2order = 3;

        Self {
            grid: Array3::zeros((nq2, nx, nx)),
            m_yorder: xorder,
            m_tauorder: q2order,
            m_reweight: false,
            xmin: 2e-7,
            xmax: 1.0,
            q2min: 100.0,
            q2max: 1_000_000.0,
        }
    }
}

impl LagrangeSubgrid {
    fn ny1(&self) -> usize {
        self.grid.raw_dim().into_pattern().1
    }

    fn y1min(&self) -> f64 {
        fy(self.xmax)
    }

    fn y1max(&self) -> f64 {
        fy(self.xmin)
    }

    fn deltay1(&self) -> f64 {
        (self.y1max() - self.y1min()) / ((self.ny1() - 1) as f64)
    }

    fn ny2(&self) -> usize {
        self.grid.raw_dim().into_pattern().2
    }

    fn y2min(&self) -> f64 {
        fy(self.xmax)
    }

    fn y2max(&self) -> f64 {
        fy(self.xmin)
    }

    fn deltay2(&self) -> f64 {
        (self.y2max() - self.y2min()) / ((self.ny2() - 1) as f64)
    }

    fn ntau(&self) -> usize {
        self.grid.raw_dim().into_pattern().0
    }

    fn taumin(&self) -> f64 {
        ftau(self.q2min)
    }

    fn taumax(&self) -> f64 {
        ftau(self.q2max)
    }

    fn deltatau(&self) -> f64 {
        (self.taumax() - self.taumin()) / ((self.ntau() - 1) as f64)
    }

    fn gety1(&self, iy: usize) -> f64 {
        (iy as f64).mul_add(self.deltay1(), self.y1min())
    }

    fn gety2(&self, iy: usize) -> f64 {
        (iy as f64).mul_add(self.deltay2(), self.y2min())
    }

    fn gettau(&self, iy: usize) -> f64 {
        (iy as f64).mul_add(self.deltatau(), self.taumin())
    }

    fn fk1(&self, x: f64) -> usize {
        let y = fy(x);
        assert!((y > self.y1min()) && (y < self.y1max()));
        let mut k = ((y - self.y1min()) / self.deltay1() - ((self.m_yorder >> 1) as f64)) as i32;
        if k < 0 {
            k = 0;
        }
        let mut k = k as usize;
        if k + self.m_yorder >= self.ny1() {
            k = self.ny1() - 1 - self.m_yorder;
        }

        k
    }

    fn fk2(&self, x: f64) -> usize {
        let y = fy(x);
        assert!((y > self.y2min()) && (y < self.y2max()));
        let mut k = ((y - self.y2min()) / self.deltay2() - ((self.m_yorder >> 1) as f64)) as i32;
        if k < 0 {
            k = 0;
        }
        let mut k = k as usize;
        if k + self.m_yorder >= self.ny2() {
            k = self.ny2() - 1 - self.m_yorder;
        }

        k
    }

    fn fkappa(&self, q2: f64) -> usize {
        let tau = ftau(q2);
        assert!((tau > self.taumin()) && (tau < self.taumax()));
        let mut kappa =
            ((tau - self.taumin()) / self.deltatau() - ((self.m_tauorder >> 1) as f64)) as i32;
        if kappa + (self.m_tauorder as i32) >= (self.ntau() as i32) {
            kappa = (self.ntau() - 1 - self.m_tauorder) as i32;
        }
        if kappa < 0 {
            kappa = 0;
        }

        kappa as usize
    }
}

#[typetag::serde]
impl Subgrid for LagrangeSubgrid {
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn convolute(&self, lumi: &dyn Fn(f64, f64, f64) -> f64) -> f64 {
        let mut dsigma = 0.0;

        for itau in 0..self.ntau() {
            for iy1 in (0..self.ny1()).rev() {
                for iy2 in (0..self.ny2()).rev() {
                    let sig = self.grid[[itau, iy1, iy2]];

                    if sig != 0.0 {
                        let tau = self.gettau(itau);
                        let q2 = fq2(tau);
                        let y1 = self.gety1(iy1);
                        let x1 = fx(y1);
                        let y2 = self.gety2(iy2);
                        let x2 = fx(y2);
                        let mut lumi_val = lumi(x1, x2, q2);

                        if self.m_reweight {
                            lumi_val *= weightfun(x1) * weightfun(x2);
                        }

                        dsigma += sig * lumi_val;
                    }
                }
            }
        }

        dsigma
    }

    fn fill(&mut self, ntuple: &Ntuple<f64>) {
        let k1 = self.fk1(ntuple.x1);
        let k2 = self.fk2(ntuple.x2);
        let k3 = self.fkappa(ntuple.q2);

        let u_y1 = (fy(ntuple.x1) - self.gety1(k1)) / self.deltay1();
        let u_y2 = (fy(ntuple.x2) - self.gety2(k2)) / self.deltay2();
        let u_tau = (ftau(ntuple.q2) - self.gettau(k3)) / self.deltatau();

        let fi1: Vec<f64> = (0..=self.m_yorder)
            .map(|i| fi(i, self.m_yorder, u_y1))
            .collect();
        let fi2: Vec<f64> = (0..=self.m_yorder)
            .map(|i| fi(i, self.m_yorder, u_y2))
            .collect();
        let fi3: Vec<f64> = (0..=self.m_tauorder)
            .map(|i| fi(i, self.m_tauorder, u_tau))
            .collect();

        let invwfun = 1.0 / (weightfun(ntuple.x1) * weightfun(ntuple.x2));

        for (i3, fi3i3) in fi3.iter().enumerate() {
            for (i1, fi1i1) in fi1.iter().enumerate() {
                for (i2, fi2i2) in fi2.iter().enumerate() {
                    let mut fi_factor = fi1i1 * fi2i2 * fi3i3;

                    if self.m_reweight {
                        fi_factor *= invwfun;
                    }

                    let fillweight = ntuple.weight * fi_factor;
                    self.grid[[k3 + i3, k1 + i1, k2 + i2]] += fillweight;
                }
            }
        }
    }

    fn is_empty(&self) -> bool {
        !self.grid.iter().any(|x| *x != 0.0)
    }

    fn merge(&mut self, other: &mut dyn Subgrid) {
        if let Some(other_grid) = other.as_any_mut().downcast_mut::<Self>() {
            self.grid.scaled_add(1.0, &other_grid.grid);
        } else {
            todo!();
        }
    }

    fn scale(&mut self, factor: f64) {
        self.grid.iter_mut().for_each(|x| *x *= factor);
    }
}
