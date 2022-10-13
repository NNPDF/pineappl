//! Supporting classes and functions for [`Grid::evolve`].

use super::grid::{GridError, Order};
use super::subgrid::{Mu2, Subgrid, SubgridEnum};
use float_cmp::approx_eq;
use itertools::Itertools;
use ndarray::{s, Array3, Array5, ArrayView1, Axis};
use std::iter;

/// Information about the evolution kernel operator (EKO) passed to [`Grid::evolve`] as `operator`,
/// which is used to convert a [`Grid`] into an [`FkTable`]. The dimensions of the EKO must
/// correspond to the values given in [`fac1`], [`pids0`], [`x0`], [`pids1`] and [`x1`], exactly in
/// this order. Members with a `1` are defined at the squared factorization scales given in
/// [`fac1`] (often called process scales) and are found in the [`Grid`] that [`Grid::evolve`] is
/// called with. Members with a `0` are defined at the squared factorization scale [`fac0`] (often
/// called fitting scale or starting scale) and are found in the [`FkTable`] resulting from
/// [`Grid::evolve`].
///
/// The EKO may convert a `Grid` from a basis given by the particle identifiers [`pids1`] to a
/// possibly different basis given by [`pids0`]. This basis must also be identified using
/// [`lumi_id_types`], which tells [`FkTable::convolute`] how to perform a convolution. The members
/// [`ren1`] and [`alphas`] must be the strong couplings given at the respective renormalization
/// scales. Finally, [`xir`] and [`xif`] can be used to vary the renormalization and factorization
/// scales, respectively, around their central values.
pub struct OperatorInfo {
    /// Squared factorization scales of the `Grid`.
    pub fac1: Vec<f64>,
    /// Particle identifiers of the `FkTable`.
    pub pids0: Vec<i32>,
    /// `x`-grid coordinates of the `FkTable`
    pub x0: Vec<f64>,
    /// Particle identifiers of the `Grid`. If the `Grid` contains more particle identifiers than
    /// given here, the contributions of them are silently ignored.
    pub pids1: Vec<i32>,
    /// `x`-grid coordinates of the `Grid`.
    pub x1: Vec<f64>,

    /// Squared factorization scale of the `FkTable`.
    pub fac0: f64,
    /// Renormalization scales of the `Grid`.
    pub ren1: Vec<f64>,
    /// Strong couplings corresponding to the order given in [`ren1`].
    pub alphas: Vec<f64>,
    /// Multiplicative factor for the central renormalization scale.
    pub xir: f64,
    /// Multiplicative factor for the central factorization scale.
    pub xif: f64,
    /// Identifier of the particle basis for the `FkTable`.
    pub lumi_id_types: String,
}

pub(crate) fn pids(
    operator: &Array5<f64>,
    info: &OperatorInfo,
    pid1_nonzero: &dyn Fn(i32) -> bool,
) -> (Vec<(usize, usize)>, Vec<(i32, i32)>) {
    // list of all non-zero PID indices
    let pid_indices: Vec<_> = (0..operator.dim().3)
        .cartesian_product(0..operator.dim().1)
        .filter(|&(pid0_idx, pid1_idx)| {
            // 1) at least one element of the operator must be non-zero, and 2) the pid must be
            // contained in the lumi somewhere
            operator
                .slice(s![.., pid1_idx, .., pid0_idx, ..])
                .iter()
                .any(|&value| value != 0.0)
                && pid1_nonzero(info.pids1[pid1_idx])
        })
        .collect();

    // list of all non-zero (pid0, pid1) combinations
    let pids = pid_indices
        .iter()
        .map(|&(pid0_idx, pid1_idx)| (info.pids0[pid0_idx], info.pids1[pid1_idx]))
        .collect();

    (pid_indices, pids)
}

pub(crate) fn lumi0_with_two(pids_a: &[(i32, i32)], pids_b: &[(i32, i32)]) -> Vec<(i32, i32)> {
    let mut pids0_a: Vec<_> = pids_a.iter().map(|&(pid0, _)| pid0).collect();
    pids0_a.sort_unstable();
    pids0_a.dedup();
    let mut pids0_b: Vec<_> = pids_b.iter().map(|&(pid0, _)| pid0).collect();
    pids0_b.sort_unstable();
    pids0_b.dedup();

    pids0_a
        .iter()
        .copied()
        .cartesian_product(pids0_b.iter().copied())
        .collect()
}

pub(crate) fn operators(
    operator: &Array5<f64>,
    info: &OperatorInfo,
    pid_indices: &[(usize, usize)],
    x1: &[f64],
) -> Result<Vec<Array3<f64>>, GridError> {
    // permutation between the grid x values and the operator x1 values
    let x1_indices: Vec<_> = if let Some(x1_indices) = x1
        .iter()
        .map(|&x1p| {
            info.x1
                .iter()
                .position(|&x1| approx_eq!(f64, x1p, x1, ulps = 64))
        })
        .collect()
    {
        x1_indices
    } else {
        return Err(GridError::EvolutionFailure(
            "operator information does not match grid's x-grid values".to_string(),
        ));
    };

    // create the corresponding operators accessible in the form [muf2, x0, x1]
    let operators: Vec<_> = pid_indices
        .iter()
        .map(|&(pid0_idx, pid1_idx)| {
            operator
                .slice(s![.., pid1_idx, .., pid0_idx, ..])
                .select(Axis(1), &x1_indices)
                .permuted_axes([0, 2, 1])
                .as_standard_layout()
                .into_owned()
        })
        .collect();

    Ok(operators)
}

pub(crate) fn ndarray_from_subgrid_orders(
    info: &OperatorInfo,
    subgrids: &ArrayView1<SubgridEnum>,
    orders: &[Order],
    order_mask: &[bool],
) -> Result<(Vec<f64>, Vec<f64>, Array3<f64>), GridError> {
    let mut x1_a: Vec<_> = subgrids
        .iter()
        .flat_map(|subgrid| subgrid.x1_grid().into_owned())
        .collect();
    let mut x1_b: Vec<_> = subgrids
        .iter()
        .flat_map(|subgrid| subgrid.x2_grid().into_owned())
        .collect();

    x1_a.sort_by(|a, b| a.partial_cmp(b).unwrap());
    x1_a.dedup_by(|a, b| approx_eq!(f64, *a, *b, ulps = 64));
    x1_b.sort_by(|a, b| a.partial_cmp(b).unwrap());
    x1_b.dedup_by(|a, b| approx_eq!(f64, *a, *b, ulps = 64));

    let mut array = Array3::<f64>::zeros((info.fac1.len(), x1_a.len(), x1_b.len()));

    // add subgrids for different orders, but the same bin and lumi, using the right
    // couplings
    for (subgrid, order) in subgrids
        .iter()
        .zip(orders.iter())
        .zip(order_mask.iter().chain(iter::repeat(&true)))
        .filter_map(|((subgrid, order), &enabled)| {
            (enabled && !subgrid.is_empty()).then(|| (subgrid, order))
        })
    {
        let mut logs = 1.0;

        if order.logxir > 0 {
            if approx_eq!(f64, info.xir, 1.0, ulps = 4) {
                continue;
            }

            logs *= (info.xir * info.xir).ln();
        }

        if order.logxif > 0 {
            if approx_eq!(f64, info.xif, 1.0, ulps = 4) {
                continue;
            }

            logs *= (info.xif * info.xif).ln();
        }

        let xa_indices: Vec<_> = subgrid
            .x1_grid()
            .iter()
            .map(|&xa| {
                x1_a.iter()
                    .position(|&x1a| approx_eq!(f64, x1a, xa, ulps = 64))
                    .unwrap()
            })
            .collect();
        let xb_indices: Vec<_> = subgrid
            .x2_grid()
            .iter()
            .map(|&xb| {
                x1_b.iter()
                    .position(|&x1b| approx_eq!(f64, x1b, xb, ulps = 64))
                    .unwrap()
            })
            .collect();

        for ((imu2, ix1, ix2), value) in subgrid.iter() {
            let Mu2 {
                ren: mur2,
                fac: muf2,
            } = subgrid.mu2_grid()[imu2];

            let als = if let Some(alphas) = info
                .ren1
                .iter()
                .zip(info.alphas.iter())
                .find_map(|(&ren1, &alphas)| approx_eq!(f64, ren1, mur2, ulps = 64).then(|| alphas))
            {
                alphas.powi(order.alphas.try_into().unwrap())
            } else {
                return Err(GridError::EvolutionFailure(format!(
                    "could not find alphas for mur2 = {}",
                    mur2
                )));
            };

            // TODO: get rid of the `unwrap`
            let mu2_index = info
                .fac1
                .iter()
                .position(|&fac| approx_eq!(f64, fac, muf2, ulps = 64))
                .unwrap();

            array[[mu2_index, xa_indices[ix1], xb_indices[ix2]]] += als * logs * value;
        }
    }

    Ok((x1_a, x1_b, array))
}
