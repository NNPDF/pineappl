//! TODO

use super::convert;
use arrayvec::ArrayVec;
use std::ops::IndexMut;

const INTERP_ORDER_MAX_PLUS_ONE: usize = 8;

mod applgrid {
    pub fn reweight_x(x: f64) -> f64 {
        (x.sqrt() / (1.0 - 0.99 * x)).powi(3)
    }

    pub fn fx2(y: f64) -> f64 {
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

    pub fn fy2(x: f64) -> f64 {
        (1.0 - x).mul_add(5.0, -x.ln())
    }

    pub fn ftau0(q2: f64) -> f64 {
        (q2 / 0.0625).ln().ln()
    }

    pub fn fq20(tau: f64) -> f64 {
        0.0625 * tau.exp().exp()
    }
}

fn no_reweight(_: f64) -> f64 {
    1.0
}

fn lagrange_weights(i: usize, n: usize, u: f64) -> f64 {
    let mut factorials = 1;
    let mut product = 1.0;
    for z in 0..i {
        product *= u - convert::f64_from_usize(z);
        factorials *= i - z;
    }
    for z in i + 1..=n {
        product *= convert::f64_from_usize(z) - u;
        factorials *= z - i;
    }
    product / convert::f64_from_usize(factorials)
}

/// TODO
pub enum ReweightMeth {
    /// TODO
    ApplGridX,
    /// TODO
    NoReweight,
}

/// TODO
pub enum Map {
    /// TODO
    ApplGridF2,
    /// TODO
    ApplGridH0,
}

/// TODO
pub enum InterpMeth {
    /// TODO
    Lagrange,
}

/// TODO
pub struct Interp {
    min: f64,
    max: f64,
    nodes: usize,
    order: usize,
    reweight_x: fn(f64) -> f64,
    map_x_to_y: fn(f64) -> f64,
    _map_y_to_x: fn(f64) -> f64,
    node_weights: fn(usize, usize, f64) -> f64,
}

impl Interp {
    /// TODO
    ///
    /// # Panics
    ///
    /// Panics if `nodes` is `0`, or if `nodes` is smaller or equal to `order` or if `order` is
    /// larger than some internally specified maximum.
    pub fn new(
        min: f64,
        max: f64,
        nodes: usize,
        order: usize,
        reweight: ReweightMeth,
        map: Map,
        interp_meth: InterpMeth,
    ) -> Self {
        // for interpolation to work `nodes` has to be at least `1`
        assert!(nodes > 0);
        // for each interpolated point `order + 1` nodes are updated
        assert!(nodes > order);
        // for `order`
        assert!(order < INTERP_ORDER_MAX_PLUS_ONE);

        let mut result = Self {
            min: 0.0,
            max: 0.0,
            nodes,
            order,
            reweight_x: match reweight {
                ReweightMeth::ApplGridX => applgrid::reweight_x,
                ReweightMeth::NoReweight => no_reweight,
            },
            map_x_to_y: match map {
                Map::ApplGridF2 => applgrid::fx2,
                Map::ApplGridH0 => applgrid::fq20,
            },
            _map_y_to_x: match map {
                Map::ApplGridF2 => applgrid::fy2,
                Map::ApplGridH0 => applgrid::ftau0,
            },
            node_weights: match interp_meth {
                InterpMeth::Lagrange => lagrange_weights,
            },
        };

        result.min = (result.map_x_to_y)(min);
        result.max = (result.map_x_to_y)(max);

        result
    }

    fn deltay(&self) -> f64 {
        (self.max - self.min) / convert::f64_from_usize(self.nodes - 1)
    }

    fn gety(&self, index: usize) -> f64 {
        convert::f64_from_usize(index).mul_add(self.deltay(), self.min)
    }

    /// TODO
    pub fn reweight(&self, x: f64) -> f64 {
        (self.reweight_x)(x)
    }

    /// TODO
    pub fn interpolate(&self, x: f64) -> Option<(usize, f64)> {
        let y = (self.map_x_to_y)(x);

        // points falling outside the interpolation range shouldn't happen very often, because when
        // it does it degrades the interpolation quality
        if (self.min > y) || (y > self.max) {
            return None;
        }

        if self.nodes == 1 {
            // TODO: is the `y_fraction` correct?
            Some((0, 0.0))
        } else {
            let index = convert::usize_from_f64(
                (y - self.min) / self.deltay() - convert::f64_from_usize(self.order / 2),
            )
            .min(self.nodes - self.order - 1);
            let y_fraction = (y - self.gety(index)) / self.deltay();

            Some((index, y_fraction))
        }
    }

    /// TODO
    pub fn node_weights(&self, fraction: f64) -> ArrayVec<f64, INTERP_ORDER_MAX_PLUS_ONE> {
        (0..=self.order)
            .map(|i| (self.node_weights)(i, self.order, fraction))
            .collect()
    }
}

/// TODO
pub fn interpolate<const D: usize>(
    interps: &[Interp],
    ntuple: &[f64],
    weight: f64,
    _array: &impl IndexMut<usize>,
) {
    if weight == 0.0 {
        return;
    }

    // we must have as many variables as we want to interpolate
    debug_assert_eq!(interps.len(), ntuple.len());
    debug_assert_eq!(interps.len(), D);

    let Some(result): Option<ArrayVec<_, D>> = interps
        .iter()
        .zip(ntuple)
        .map(|(interp, &x)| interp.interpolate(x))
        .collect()
    else {
        return;
    };

    //if self.static_q2 == 0.0 {
    //    self.static_q2 = q2;
    //} else if (self.static_q2 != -1.0) && (self.static_q2 != q2) {
    //    self.static_q2 = -1.0;
    //}

    let _weight = weight
        / interps
            .iter()
            .zip(ntuple)
            .map(|(interp, &x)| interp.reweight(x))
            .product::<f64>();

    let _node_weights: ArrayVec<_, D> = interps
        .iter()
        .zip(result)
        .map(|(interp, (_, fraction))| interp.node_weights(fraction))
        .collect();

    //for i3 in 0..=self.tauorder {
    //    let fi3i3 = fi(i3, self.tauorder, u_tau);

    //    for (i1, fi1i1) in fi1.iter().enumerate() {
    //        for (i2, fi2i2) in fi2.iter().enumerate() {
    //            let fillweight = factor * fi1i1 * fi2i2 * fi3i3 * weight;

    //            let grid = self
    //                .grid
    //                .get_or_insert_with(|| Array3::zeros((size, ny1, ny2)));

    //            grid[[k3 + i3 - self.itaumin, k1 + i1, k2 + i2]] += fillweight;
    //        }
    //    }
    //}
}
