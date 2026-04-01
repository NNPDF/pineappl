//! Module containing the Lagrange-interpolation subgrid.

use super::convert::f64_from_usize;
use super::sparse_array3::SparseArray3;
use super::subgrid::{ExtraSubgridParams, Mu2, Subgrid, SubgridIndexedIter, SubgridParams};
use ndarray::Array3;
use serde::{Deserialize, Serialize};
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
