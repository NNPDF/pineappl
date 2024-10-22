//! Interpolation module.

use super::convert;
use super::packed_array::PackedArray;
use arrayvec::ArrayVec;
use serde::{Deserialize, Serialize};
use std::mem;

const MAX_INTERP_ORDER_PLUS_ONE: usize = 8;
const MAX_DIMENSIONS: usize = 8;

mod applgrid {
    pub fn reweight_x(x: f64) -> f64 {
        (x.sqrt() / (1.0 - 0.99 * x)).powi(3)
    }

    pub fn fx2(y: f64) -> f64 {
        let mut yp = y;

        for _ in 0..100 {
            let x = (-yp).exp();
            let delta = y - yp - 5.0 * (1.0 - x);
            if delta.abs() < 1e-12 {
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
#[derive(Clone, Copy, Deserialize, Serialize)]
pub enum ReweightMeth {
    /// TODO
    ApplGridX,
    /// TODO
    NoReweight,
}

/// TODO
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub enum Map {
    /// TODO
    ApplGridF2,
    /// TODO
    ApplGridH0,
}

/// TODO
#[derive(Clone, Deserialize, Serialize)]
pub enum InterpMeth {
    /// TODO
    Lagrange,
}

/// TODO
#[derive(Clone, Deserialize, Serialize)]
pub struct Interp {
    min: f64,
    max: f64,
    nodes: usize,
    order: usize,
    reweight: ReweightMeth,
    map: Map,
    interp_meth: InterpMeth,
}

impl Interp {
    /// TODO
    ///
    /// # Panics
    ///
    /// Panics if `nodes` is `0`, or if `nodes` is smaller or equal to `order` or if `order` is
    /// larger than some internally specified maximum.
    #[must_use]
    pub fn new(
        min: f64,
        max: f64,
        nodes: usize,
        order: usize,
        reweight: ReweightMeth,
        map: Map,
        interp_meth: InterpMeth,
    ) -> Self {
        // minimum must be larger or equal to the maximum
        assert!(min <= max);
        // for interpolation to work `nodes` has to be at least `1`
        assert!(nodes > 0);
        // for each interpolated point `order + 1` nodes are updated
        assert!(nodes > order);
        // using arrays with fixed size limit the possible max value of `order`
        assert!(order < MAX_INTERP_ORDER_PLUS_ONE);

        let mut result = Self {
            min: 0.0,
            max: 0.0,
            nodes,
            order,
            reweight,
            map,
            interp_meth,
        };

        result.min = result.map_x_to_y(min);
        result.max = result.map_x_to_y(max);

        // for some maps the minimum in x is mapped to the maximum in y
        if result.min > result.max {
            // TODO: alternatively we have to modify our range check in `Self::interpolate`, which
            // has the advantage that we don't swap min and max in `x` space
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
    #[must_use]
    pub fn reweight(&self, x: f64) -> f64 {
        match self.reweight {
            ReweightMeth::ApplGridX => applgrid::reweight_x(x),
            ReweightMeth::NoReweight => 1.0,
        }
    }

    /// TODO
    #[must_use]
    pub fn interpolate(&self, x: f64) -> Option<(usize, f64)> {
        let y = self.map_x_to_y(x);

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
    #[must_use]
    pub fn node_weights(&self, fraction: f64) -> ArrayVec<f64, MAX_INTERP_ORDER_PLUS_ONE> {
        (0..=self.order)
            .map(|i| match self.interp_meth {
                InterpMeth::Lagrange => lagrange_weights(i, self.order, fraction),
            })
            .collect()
    }

    /// TODO
    #[must_use]
    pub const fn order(&self) -> usize {
        self.order
    }

    /// TODO
    #[must_use]
    pub fn node_values(&self) -> Vec<f64> {
        (0..self.nodes)
            .map(|node| self.map_y_to_x(self.gety(node)))
            .collect()
    }

    fn map_y_to_x(&self, y: f64) -> f64 {
        match self.map {
            Map::ApplGridF2 => applgrid::fx2(y),
            Map::ApplGridH0 => applgrid::fq20(y),
        }
    }

    fn map_x_to_y(&self, x: f64) -> f64 {
        match self.map {
            Map::ApplGridF2 => applgrid::fy2(x),
            Map::ApplGridH0 => applgrid::ftau0(x),
        }
    }

    /// TODO
    #[must_use]
    pub const fn nodes(&self) -> usize {
        self.nodes
    }

    /// TODO
    #[must_use]
    pub fn min(&self) -> f64 {
        self.map_y_to_x(self.min).min(self.map_y_to_x(self.max))
    }

    /// TODO
    #[must_use]
    pub fn max(&self) -> f64 {
        self.map_y_to_x(self.min).max(self.map_y_to_x(self.max))
    }

    /// TODO
    #[must_use]
    pub const fn map(&self) -> Map {
        self.map
    }
}

/// TODO
pub fn interpolate(
    interps: &[Interp],
    ntuple: &[f64],
    weight: f64,
    array: &mut PackedArray<f64>,
) -> bool {
    use itertools::Itertools;

    if weight == 0.0 {
        return false;
    }

    // we must have as many variables as we want to interpolate
    debug_assert_eq!(interps.len(), ntuple.len());
    debug_assert!(interps.len() <= MAX_DIMENSIONS);

    let Some((indices, fractions)): Option<(
        ArrayVec<_, MAX_DIMENSIONS>,
        ArrayVec<_, MAX_DIMENSIONS>,
    )> = interps
        .iter()
        .zip(ntuple)
        .map(|(interp, &x)| interp.interpolate(x))
        .collect()
    else {
        return false;
    };

    // TODO: add static value detection

    let weight = weight
        / interps
            .iter()
            .zip(ntuple)
            .map(|(interp, &x)| interp.reweight(x))
            .product::<f64>();

    let node_weights: ArrayVec<_, MAX_DIMENSIONS> = interps
        .iter()
        .zip(fractions)
        .map(|(interp, fraction)| interp.node_weights(fraction))
        .collect();

    let shape: ArrayVec<_, MAX_DIMENSIONS> =
        interps.iter().map(|interp| interp.order() + 1).collect();

    for (i, node_weight) in node_weights
        .into_iter()
        // TODO: replace this with something else to avoid allocating memory
        .multi_cartesian_product()
        .map(|weights| weights.iter().product::<f64>())
        .enumerate()
    {
        let idx = array.sub_block_idx(&indices, i, &shape);
        array[idx] += weight * node_weight;
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::assert_approx_eq;

    #[test]
    fn interpolate_two_points() {
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

        let node_values: Vec<_> = interps.iter().map(Interp::node_values).collect();

        let q2_reference = [
            9.999999999999999e1,
            1.2242682307575689e2,
            1.507173582975839e2,
            1.8660624792652183e2,
            2.3239844323901826e2,
            2.911750445478316e2,
            3.670799619445291e2,
            4.657216764869711e2,
            5.947399998930223e2,
            7.646109579666331e2,
            9.897977073478313e2,
            1.2904078604330668e3,
            1.694597307328949e3,
            2.2420826491130997e3,
            2.989312590729525e3,
            4.017141299790263e3,
            5.442305429193529e3,
            7.434731381687921e3,
            1.024385467001917e4,
            1.4238990475802799e4,
            1.9971806922234402e4,
            2.8273883344269376e4,
            4.041048232844362e4,
            5.832525318921733e4,
            8.503347534094655e4,
            1.2526040013230646e5,
            1.864882133214792e5,
            2.806914902174795e5,
            4.272453808062111e5,
            6.578537431299294e5,
            1.0249965523865514e6,
            1.6165812577807596e6,
            2.581663421106388e6,
            4.1761634755570055e6,
            6.845167341538921e6,
            1.1373037585359517e7,
            1.916090997202005e7,
            3.2746801715531096e7,
            5.679435282347418e7,
            9.99999999999995e7,
        ];

        assert_eq!(node_values[0].len(), interps[0].nodes());

        for (&node, ref_node) in node_values[0].iter().zip(q2_reference) {
            assert_approx_eq!(f64, node, ref_node, ulps = 4);
        }

        let x_reference = [
            1.0,
            9.309440808717544e-1,
            8.627839323906108e-1,
            7.956242522922756e-1,
            7.295868442414312e-1,
            6.648139482473823e-1,
            6.01472197967335e-1,
            5.397572337880445e-1,
            4.798989029610255e-1,
            4.221667753589648e-1,
            3.668753186482242e-1,
            3.1438740076927585e-1,
            2.651137041582823e-1,
            2.195041265003886e-1,
            1.7802566042569432e-1,
            1.4112080644440345e-1,
            1.0914375746330703e-1,
            8.228122126204893e-2,
            6.0480028754447364e-2,
            4.341491741702269e-2,
            3.0521584007828916e-2,
            2.108918668378717e-2,
            1.4375068581090129e-2,
            9.699159574043399e-3,
            6.496206194633799e-3,
            4.328500638820811e-3,
            2.8738675812817515e-3,
            1.9034634022867384e-3,
            1.2586797144272762e-3,
            8.314068836488144e-4,
            5.487795323670796e-4,
            3.6205449638139736e-4,
            2.3878782918561914e-4,
            1.5745605600841445e-4,
            1.0381172986576898e-4,
            6.843744918967897e-5,
            4.511438394964044e-5,
            2.97384953722449e-5,
            1.9602505002391748e-5,
            1.292101569074731e-5,
            8.516806677573355e-6,
            5.613757716930151e-6,
            3.7002272069854957e-6,
            2.438943292891682e-6,
            1.607585498470808e-6,
            1.0596094959101024e-6,
            6.984208530700364e-7,
            4.6035014748963906e-7,
            3.034304765867952e-7,
            1.9999999999999954e-7,
        ];

        assert_eq!(node_values[1].len(), interps[1].nodes());
        assert_eq!(node_values[2].len(), interps[2].nodes());

        for (&node, ref_node) in node_values[1].iter().zip(x_reference) {
            assert_approx_eq!(f64, node, ref_node, ulps = 4);
        }

        for (&node, ref_node) in node_values[2].iter().zip(x_reference) {
            assert_approx_eq!(f64, node, ref_node, ulps = 4);
        }

        let mut array = crate::packed_array::PackedArray::<f64>::new(vec![40, 50, 50]);
        let ntuples = [[100000.0, 0.25, 0.5], [1000.0, 0.5, 0.5]];
        let weight = 1.0;

        for ntuple in &ntuples {
            interpolate(&interps, ntuple, weight, &mut array);
        }

        let reference = [
            ([9, 6, 6], -4.091358497150521e-6),
            ([9, 6, 7], 3.085859446366878e-5),
            ([9, 6, 8], 6.0021251939206686e-5),
            ([9, 6, 9], -5.0714506160633226e-6),
            ([9, 7, 6], 3.085859446366878e-5),
            ([9, 7, 7], -2.3274735101712016e-4),
            ([9, 7, 8], -4.527032950262464e-4),
            ([9, 7, 9], 3.825082500411933e-5),
            ([9, 8, 6], 6.002125193920668e-5),
            ([9, 8, 7], -4.5270329502624637e-4),
            ([9, 8, 8], -8.805267704745902e-4),
            ([9, 8, 9], 7.439944833384343e-5),
            ([9, 9, 6], -5.071450616063322e-6),
            ([9, 9, 7], 3.825082500411933e-5),
            ([9, 9, 8], 7.439944833384344e-5),
            ([9, 9, 9], -6.286325524659303e-6),
            ([10, 6, 6], 3.25604540320380e-4),
            ([10, 6, 7], -2.4558342839606324e-3),
            ([10, 6, 8], -4.776700003368127e-3),
            ([10, 6, 9], 4.036036802325844e-4),
            ([10, 7, 6], -2.4558342839606324e-3),
            ([10, 7, 7], 1.8522843767295388e-2),
            ([10, 7, 8], 3.602770287209066e-2),
            ([10, 7, 9], -3.0441337030269453e-3),
            ([10, 8, 6], -4.776700003368127e-3),
            ([10, 8, 7], 3.602770287209066e-2),
            ([10, 8, 8], 7.007538316181437e-2),
            ([10, 8, 9], -5.920966884642993e-3),
            ([10, 9, 6], 4.036036802325844e-4),
            ([10, 9, 7], -3.0441337030269453e-3),
            ([10, 9, 8], -5.920966884642993e-3),
            ([10, 9, 9], 5.002876512010676e-4),
            ([11, 6, 6], 1.3274904136884986e-5),
            ([11, 6, 7], -1.0012441676511976e-4),
            ([11, 6, 8], -1.947461622401742e-4),
            ([11, 6, 9], 1.6454930754680843e-5),
            ([11, 7, 6], -1.0012441676511976e-4),
            ([11, 7, 7], 7.551767401996306e-4),
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
            ([12, 6, 8], 3.180926156637114e-5),
            ([12, 6, 9], -2.6876996722875166e-6),
            ([12, 7, 6], 1.6354025801721504e-5),
            ([12, 7, 7], -1.2334833293517984e-4),
            ([12, 7, 8], -2.3991764680339134e-4),
            ([12, 7, 9], 2.0271661426154572e-5),
            ([12, 8, 6], 3.180926156637114e-5),
            ([12, 8, 7], -2.3991764680339134e-4),
            ([12, 8, 8], -4.6664981907720756e-4),
            ([12, 8, 9], 3.942922608215463e-5),
            ([12, 9, 6], -2.6876996722875166e-6),
            ([12, 9, 7], 2.0271661426154572e-5),
            ([12, 9, 8], 3.942922608215462e-5),
            ([12, 9, 9], -3.3315428526512343e-6),
            ([23, 11, 6], -2.4353100307613186e-4),
            ([23, 11, 7], 1.8368041980410083e-3),
            ([23, 11, 8], 3.572660694686239e-3),
            ([23, 11, 9], -3.0186928289005667e-4),
            ([23, 12, 6], 2.9987494527093064e-3),
            ([23, 12, 7], -2.2617718130482554e-2),
            ([23, 12, 8], -4.399240411931119e-2),
            ([23, 12, 9], 3.717105154670258e-3),
            ([23, 13, 6], 1.424894308599361e-3),
            ([23, 13, 7], -1.0747099197804599e-2),
            ([23, 13, 8], -2.090355571205706e-2),
            ([23, 13, 9], 1.766230244600708e-3),
            ([23, 14, 6], -1.9189233197773798e-4),
            ([23, 14, 7], 1.4473255417027965e-3),
            ([23, 14, 8], 2.815108480681732e-3),
            ([23, 14, 9], -2.3786047737056168e-4),
            ([24, 11, 6], 2.4624842908465045e-3),
            ([24, 11, 7], -1.8573000668924675e-2),
            ([24, 11, 8], -3.612526013552098e-2),
            ([24, 11, 9], 3.0523767307502974e-3),
            ([24, 12, 6], -3.032210817598752e-2),
            ([24, 12, 7], 2.2870096573985707e-1),
            ([24, 12, 8], 4.4483290707142076e-1),
            ([24, 12, 9], -3.758582248330245e-2),
            ([24, 13, 6], -1.4407939057950724e-2),
            ([24, 13, 7], 1.0867020055959625e-1),
            ([24, 13, 8], 2.1136806777609088e-1),
            ([24, 13, 9], -1.7859386182495846e-2),
            ([24, 14, 6], 1.940335509895472e-3),
            ([24, 14, 7], -1.463475436460085e-2),
            ([24, 14, 8], -2.8465206988616668e-2),
            ([24, 14, 9], 2.4051462916002946e-3),
            ([25, 11, 6], 1.7967411488022474e-3),
            ([25, 11, 7], -1.3551710637356816e-2),
            ([25, 11, 8], -2.6358641814670535e-2),
            ([25, 11, 9], 2.2271536489275397e-3),
            ([25, 12, 6], -2.2124396764984615e-2),
            ([25, 12, 7], 1.6687068317270662e-1),
            ([25, 12, 8], 3.245704313515834e-1),
            ([25, 12, 9], -2.7424334895599013e-2),
            ([25, 13, 6], -1.0512691216379747e-2),
            ([25, 13, 7], 7.929074785159325e-2),
            ([25, 13, 8], 1.54223808179330e-1),
            ([25, 13, 9], -1.3031024874237778e-2),
            ([25, 14, 6], 1.415757520188257e-3),
            ([25, 14, 7], -1.0678186036448784e-2),
            ([25, 14, 8], -2.0769516741988788e-2),
            ([25, 14, 9], 1.754904722467019e-3),
            ([26, 11, 6], -2.1941078639412583e-4),
            ([26, 11, 7], 1.6548802758317395e-3),
            ([26, 11, 8], 3.218811086223127e-3),
            ([26, 11, 9], -2.7197102590848567e-4),
            ([26, 12, 6], 2.7017421490774826e-3),
            ([26, 12, 7], -2.0377575170166206e-2),
            ([26, 12, 8], -3.9635232727098596e-2),
            ([26, 12, 9], 3.3489492294371172e-3),
            ([26, 13, 6], 1.2837674744868705e-3),
            ([26, 13, 7], -9.682666505130055e-3),
            ([26, 13, 8], -1.8833189775767707e-2),
            ([26, 13, 9], 1.5912962293338163e-3),
            ([26, 14, 6], -1.7288660142000884e-4),
            ([26, 14, 7], 1.303977034800947e-3),
            ([26, 14, 8], 2.536289662216269e-3),
            ([26, 14, 9], -2.1430189065349477e-4),
        ];

        for ((index, value), (ref_index, ref_value)) in array.indexed_iter().zip(reference) {
            assert_eq!(index, ref_index);
            assert_approx_eq!(f64, value, ref_value, ulps = 4);
        }
    }

    #[test]
    fn interpolate_zero_and_outside() {
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
        let mut array = crate::packed_array::PackedArray::<f64>::new(vec![40, 50, 50]);

        let ntuple = [1000.0, 0.5, 0.5];
        let weight = 0.0;
        interpolate(&interps, &ntuple, weight, &mut array);

        assert_eq!(array.non_zeros(), 0);
        assert_eq!(array.explicit_zeros(), 0);

        let ntuple = [10.0, 0.5, 0.5];
        let weight = 1.0;
        interpolate(&interps, &ntuple, weight, &mut array);

        assert_eq!(array.non_zeros(), 0);
        assert_eq!(array.explicit_zeros(), 0);
    }

    #[test]
    fn interpolate_with_one_node() {
        // TODO: does it make sense for an interpolation to have `min = max`? There will be
        // numerical problems if the `x` value doesn't exactly hit the limits
        let interps = vec![Interp::new(
            90.0_f64.powi(2),
            90.0_f64.powi(2),
            1,
            0,
            ReweightMeth::NoReweight,
            Map::ApplGridH0,
            InterpMeth::Lagrange,
        )];
        let mut array = crate::packed_array::PackedArray::<f64>::new(vec![1]);

        let ntuple = [90.0_f64.powi(2)];
        let weight = 1.0;
        interpolate(&interps, &ntuple, weight, &mut array);

        assert_eq!(array[[0]], 1.0);

        let node_values = interps[0].node_values();

        assert_eq!(node_values.len(), 1);

        // TODO: the return value is not the one expected (90^2), because `deltay` is zero
        assert!(node_values[0].is_nan());
    }
}
