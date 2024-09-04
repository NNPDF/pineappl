//! TODO

use super::convert;
use arrayvec::ArrayVec;
use std::mem;
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
                Map::ApplGridF2 => applgrid::fy2,
                Map::ApplGridH0 => applgrid::ftau0,
            },
            _map_y_to_x: match map {
                Map::ApplGridF2 => applgrid::fx2,
                Map::ApplGridH0 => applgrid::fq20,
            },
            node_weights: match interp_meth {
                InterpMeth::Lagrange => lagrange_weights,
            },
        };

        result.min = (result.map_x_to_y)(min);
        result.max = (result.map_x_to_y)(max);

        // for some maps the minimum in x is mapped to the maximum in y
        if result.min > result.max {
            mem::swap(&mut result.min, &mut result.max);
        }

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

    /// TODO
    pub fn order(&self) -> usize {
        self.order
    }
}

/// TODO
pub fn interpolate<const D: usize>(
    interps: &[Interp],
    ntuple: &[f64],
    weight: f64,
    array: &mut impl IndexMut<[usize; D], Output = f64>,
) {
    use super::packed_array;
    use itertools::Itertools;

    if weight == 0.0 {
        return;
    }

    // we must have as many variables as we want to interpolate
    debug_assert_eq!(interps.len(), ntuple.len());
    debug_assert_eq!(interps.len(), D);

    let Some((indices, fractions)): Option<(ArrayVec<_, D>, ArrayVec<_, D>)> = interps
        .iter()
        .zip(ntuple)
        .map(|(interp, &x)| interp.interpolate(x))
        .collect()
    else {
        return;
    };

    // TODO: add static value detection

    let weight = weight
        / interps
            .iter()
            .zip(ntuple)
            .map(|(interp, &x)| interp.reweight(x))
            .product::<f64>();

    let node_weights: ArrayVec<_, D> = interps
        .iter()
        .zip(fractions)
        .map(|(interp, fraction)| interp.node_weights(fraction))
        .collect();

    let shape: ArrayVec<_, D> = interps.iter().map(|interp| interp.order() + 1).collect();

    for (i, node_weights) in node_weights
        .into_iter()
        .multi_cartesian_product()
        .enumerate()
    {
        let mut index = packed_array::unravel_index::<D>(i, &shape);
        for (entry, start_index) in index.iter_mut().zip(&indices) {
            *entry += start_index;
        }
        array[index] += weight * node_weights.iter().product::<f64>();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::assert_approx_eq;

    #[test]
    fn test_interpolation() {
        let interps = vec![
            Interp::new(
                1e2,
                1e8,
                40,
                3,
                ReweightMeth::NoReweight,
                Map::ApplGridH0,
                InterpMeth::Lagrange,
            ),
            Interp::new(
                2e-7,
                1.0,
                50,
                3,
                ReweightMeth::ApplGridX,
                Map::ApplGridF2,
                InterpMeth::Lagrange,
            ),
            Interp::new(
                2e-7,
                1.0,
                50,
                3,
                ReweightMeth::ApplGridX,
                Map::ApplGridF2,
                InterpMeth::Lagrange,
            ),
        ];
        let mut array = crate::packed_array::PackedArray::<f64, 3>::new([40, 50, 50]);
        let ntuples = [[100000.0, 0.25, 0.5], [1000.0, 0.5, 0.5]];
        let weight = 1.0;

        for ntuple in &ntuples {
            interpolate(&interps, ntuple, weight, &mut array);
        }

        let reference = [
            ([9, 6, 6], -4.0913584971505212e-6),
            ([9, 6, 7], 3.0858594463668783e-5),
            ([9, 6, 8], 6.0021251939206686e-5),
            ([9, 6, 9], -5.0714506160633226e-6),
            ([9, 7, 6], 3.0858594463668783e-5),
            ([9, 7, 7], -2.3274735101712016e-4),
            ([9, 7, 8], -4.5270329502624643e-4),
            ([9, 7, 9], 3.8250825004119329e-5),
            ([9, 8, 6], 6.0021251939206680e-5),
            ([9, 8, 7], -4.5270329502624637e-4),
            ([9, 8, 8], -8.8052677047459023e-4),
            ([9, 8, 9], 7.4399448333843429e-5),
            ([9, 9, 6], -5.0714506160633217e-6),
            ([9, 9, 7], 3.8250825004119329e-5),
            ([9, 9, 8], 7.4399448333843443e-5),
            ([9, 9, 9], -6.2863255246593026e-6),
            ([10, 6, 6], 3.2560454032038003e-4),
            ([10, 6, 7], -2.4558342839606324e-3),
            ([10, 6, 8], -4.7767000033681270e-3),
            ([10, 6, 9], 4.0360368023258439e-4),
            ([10, 7, 6], -2.4558342839606324e-3),
            ([10, 7, 7], 1.8522843767295388e-2),
            ([10, 7, 8], 3.6027702872090658e-2),
            ([10, 7, 9], -3.0441337030269453e-3),
            ([10, 8, 6], -4.7767000033681270e-3),
            ([10, 8, 7], 3.6027702872090658e-2),
            ([10, 8, 8], 7.0075383161814372e-2),
            ([10, 8, 9], -5.9209668846429931e-3),
            ([10, 9, 6], 4.0360368023258439e-4),
            ([10, 9, 7], -3.0441337030269453e-3),
            ([10, 9, 8], -5.9209668846429931e-3),
            ([10, 9, 9], 5.0028765120106755e-4),
            ([11, 6, 6], 1.3274904136884986e-5),
            ([11, 6, 7], -1.0012441676511976e-4),
            ([11, 6, 8], -1.9474616224017421e-4),
            ([11, 6, 9], 1.6454930754680843e-5),
            ([11, 7, 6], -1.0012441676511976e-4),
            ([11, 7, 7], 7.5517674019963063e-4),
            ([11, 7, 8], 1.4688502237364042e-3),
            ([11, 7, 9], -1.2410939677862364e-4),
            ([11, 8, 6], -1.9474616224017418e-4),
            ([11, 8, 7], 1.4688502237364042e-3),
            ([11, 8, 8], 2.8569748840518382e-3),
            ([11, 8, 9], -2.4139794768822075e-4),
            ([11, 9, 6], 1.6454930754680843e-5),
            ([11, 9, 7], -1.2410939677862364e-4),
            ([11, 9, 8], -2.4139794768822075e-4),
            ([11, 9, 9], 2.0396738337944602e-5),
            ([12, 6, 6], -2.1682835394615433e-6),
            ([12, 6, 7], 1.6354025801721504e-5),
            ([12, 6, 8], 3.1809261566371142e-5),
            ([12, 6, 9], -2.6876996722875166e-6),
            ([12, 7, 6], 1.6354025801721504e-5),
            ([12, 7, 7], -1.2334833293517984e-4),
            ([12, 7, 8], -2.3991764680339134e-4),
            ([12, 7, 9], 2.0271661426154572e-5),
            ([12, 8, 6], 3.1809261566371142e-5),
            ([12, 8, 7], -2.3991764680339134e-4),
            ([12, 8, 8], -4.6664981907720756e-4),
            ([12, 8, 9], 3.9429226082154630e-5),
            ([12, 9, 6], -2.6876996722875166e-6),
            ([12, 9, 7], 2.0271661426154572e-5),
            ([12, 9, 8], 3.9429226082154623e-5),
            ([12, 9, 9], -3.3315428526512343e-6),
            ([23, 11, 6], -2.4353100307613186e-4),
            ([23, 11, 7], 1.8368041980410083e-3),
            ([23, 11, 8], 3.5726606946862392e-3),
            ([23, 11, 9], -3.0186928289005667e-4),
            ([23, 12, 6], 2.9987494527093064e-3),
            ([23, 12, 7], -2.2617718130482554e-2),
            ([23, 12, 8], -4.3992404119311192e-2),
            ([23, 12, 9], 3.7171051546702580e-3),
            ([23, 13, 6], 1.4248943085993610e-3),
            ([23, 13, 7], -1.0747099197804599e-2),
            ([23, 13, 8], -2.0903555712057060e-2),
            ([23, 13, 9], 1.7662302446007081e-3),
            ([23, 14, 6], -1.9189233197773798e-4),
            ([23, 14, 7], 1.4473255417027965e-3),
            ([23, 14, 8], 2.8151084806817320e-3),
            ([23, 14, 9], -2.3786047737056168e-4),
            ([24, 11, 6], 2.4624842908465045e-3),
            ([24, 11, 7], -1.8573000668924675e-2),
            ([24, 11, 8], -3.6125260135520983e-2),
            ([24, 11, 9], 3.0523767307502974e-3),
            ([24, 12, 6], -3.0322108175987520e-2),
            ([24, 12, 7], 2.2870096573985707e-1),
            ([24, 12, 8], 4.4483290707142076e-1),
            ([24, 12, 9], -3.7585822483302452e-2),
            ([24, 13, 6], -1.4407939057950724e-2),
            ([24, 13, 7], 1.0867020055959625e-1),
            ([24, 13, 8], 2.1136806777609088e-1),
            ([24, 13, 9], -1.7859386182495846e-2),
            ([24, 14, 6], 1.9403355098954720e-3),
            ([24, 14, 7], -1.4634754364600849e-2),
            ([24, 14, 8], -2.8465206988616668e-2),
            ([24, 14, 9], 2.4051462916002946e-3),
            ([25, 11, 6], 1.7967411488022474e-3),
            ([25, 11, 7], -1.3551710637356816e-2),
            ([25, 11, 8], -2.6358641814670535e-2),
            ([25, 11, 9], 2.2271536489275397e-3),
            ([25, 12, 6], -2.2124396764984615e-2),
            ([25, 12, 7], 1.6687068317270662e-1),
            ([25, 12, 8], 3.2457043135158342e-1),
            ([25, 12, 9], -2.7424334895599013e-2),
            ([25, 13, 6], -1.0512691216379747e-2),
            ([25, 13, 7], 7.9290747851593249e-2),
            ([25, 13, 8], 1.5422380817933001e-1),
            ([25, 13, 9], -1.3031024874237778e-2),
            ([25, 14, 6], 1.4157575201882570e-3),
            ([25, 14, 7], -1.0678186036448784e-2),
            ([25, 14, 8], -2.0769516741988788e-2),
            ([25, 14, 9], 1.7549047224670190e-3),
            ([26, 11, 6], -2.1941078639412583e-4),
            ([26, 11, 7], 1.6548802758317395e-3),
            ([26, 11, 8], 3.2188110862231270e-3),
            ([26, 11, 9], -2.7197102590848567e-4),
            ([26, 12, 6], 2.7017421490774826e-3),
            ([26, 12, 7], -2.0377575170166206e-2),
            ([26, 12, 8], -3.9635232727098596e-2),
            ([26, 12, 9], 3.3489492294371172e-3),
            ([26, 13, 6], 1.2837674744868705e-3),
            ([26, 13, 7], -9.6826665051300553e-3),
            ([26, 13, 8], -1.8833189775767707e-2),
            ([26, 13, 9], 1.5912962293338163e-3),
            ([26, 14, 6], -1.7288660142000884e-4),
            ([26, 14, 7], 1.3039770348009471e-3),
            ([26, 14, 8], 2.5362896622162690e-3),
            ([26, 14, 9], -2.1430189065349477e-4),
        ];

        for ((index, value), (ref_index, ref_value)) in array.indexed_iter().zip(reference) {
            assert_eq!(index, ref_index);
            assert_approx_eq!(f64, value, ref_value, ulps = 4);
        }
    }
}
