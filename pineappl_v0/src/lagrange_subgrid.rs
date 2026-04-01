//! Module containing the Lagrange-interpolation subgrid.

use super::convert::f64_from_usize;
use super::sparse_array3::SparseArray3;
use super::subgrid::{Mu2, Subgrid, SubgridIndexedIter};
use ndarray::Array3;
use serde::Deserialize;
use std::borrow::Cow;
use std::iter;

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

fn fq2(tau: f64) -> f64 {
    0.0625 * tau.exp().exp()
}

/// Subgrid which uses Lagrange-interpolation.
#[derive(Deserialize)]
pub struct LagrangeSubgridV1 {
    grid: Option<Array3<f64>>,
    ntau: usize,
    ny: usize,
    _yorder: usize,
    _tauorder: usize,
    itaumin: usize,
    _itaumax: usize,
    reweight: bool,
    ymin: f64,
    ymax: f64,
    taumin: f64,
    taumax: f64,
}

impl LagrangeSubgridV1 {
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
}

/// Subgrid which uses Lagrange-interpolation.
#[derive(Deserialize)]
pub struct LagrangeSubgridV2 {
    grid: Option<Array3<f64>>,
    ntau: usize,
    ny1: usize,
    ny2: usize,
    _y1order: usize,
    _y2order: usize,
    _tauorder: usize,
    itaumin: usize,
    _itaumax: usize,
    reweight1: bool,
    reweight2: bool,
    y1min: f64,
    y1max: f64,
    y2min: f64,
    y2max: f64,
    taumin: f64,
    taumax: f64,
    pub(crate) _static_q2: f64,
}

impl LagrangeSubgridV2 {
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
}

/// Subgrid which uses Lagrange-interpolation, but also stores its contents in a space-efficient
/// structure.
#[derive(Deserialize)]
pub struct LagrangeSparseSubgridV1 {
    array: SparseArray3<f64>,
    ntau: usize,
    ny: usize,
    _yorder: usize,
    _tauorder: usize,
    reweight: bool,
    ymin: f64,
    ymax: f64,
    taumin: f64,
    taumax: f64,
}

impl LagrangeSparseSubgridV1 {
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
}
