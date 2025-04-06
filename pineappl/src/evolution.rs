//! Supporting classes and functions for [`Grid::evolve`].

use super::boc::{Channel, Kinematics, Order};
use super::convolutions::ConvType;
use super::error::{Error, Result};
use super::grid::Grid;
use super::packed_array::PackedArray;
use super::pids::PidBasis;
use super::subgrid::{self, ImportSubgridV1, Subgrid, SubgridEnum};
use float_cmp::approx_eq;
use itertools::izip;
use itertools::Itertools;
use ndarray::linalg;
use ndarray::{
    s, Array1, Array2, Array3, ArrayD, ArrayView1, ArrayView4, ArrayViewD, ArrayViewMutD, Axis,
    Ix1, Ix2,
};
use std::iter;

/// This structure captures the information needed to create an evolution kernel operator (EKO) for
/// a specific [`Grid`].
pub struct EvolveInfo {
    /// Squared factorization scales of the `Grid`.
    pub fac1: Vec<f64>,
    /// Squared fragmentation scales of the `Grid`.
    pub frg1: Vec<f64>,
    /// Particle identifiers of the `Grid`.
    pub pids1: Vec<i32>,
    /// `x`-grid coordinates of the `Grid`.
    pub x1: Vec<f64>,
    /// Renormalization scales of the `Grid`.
    pub ren1: Vec<f64>,
}

/// Information about the evolution kernel operator slice (EKO) passed to [`Grid::evolve`] as
/// `operator`, which is used to convert a [`Grid`] into an [`FkTable`](super::fk_table::FkTable).
/// The dimensions of the EKO must correspond to the values given in [`fac1`](Self::fac1),
/// [`pids0`](Self::pids0), [`x0`](Self::x0), [`pids1`](Self::pids1) and [`x1`](Self::x1), exactly
/// in this order. Members with a `1` are defined at the squared factorization scale given as
/// `fac1` (often called process scale) and are found in the [`Grid`] that [`Grid::evolve`] is
/// called with. Members with a `0` are defined at the squared factorization scale
/// [`fac0`](Self::fac0) (often called fitting scale or starting scale) and are found in the
/// `FkTable` resulting from [`Grid::evolve`].
///
/// The EKO slice may convert a `Grid` from a basis given by the particle identifiers `pids1` to a
/// possibly different basis given by `pids0`. This basis must also be identified using
/// [`pid_basis`](Self::pid_basis), which tells
/// [`FkTable::convolve`](super::fk_table::FkTable::convolve) how to perform a convolution.
///
/// [`Grid::evolve`]: super::grid::Grid::evolve
#[derive(Clone)]
pub struct OperatorSliceInfo {
    /// Squared factorization/fragmentation scale of the `FkTable`.
    pub fac0: f64,
    /// Particle identifiers of the `FkTable`.
    pub pids0: Vec<i32>,
    /// `x`-grid coordinates of the `FkTable`
    pub x0: Vec<f64>,
    /// Squared factorization/fragmentation scale of the slice of `Grid` that should be evolved.
    pub fac1: f64,
    /// Particle identifiers of the `Grid`. If the `Grid` contains more particle identifiers than
    /// given here, the contributions of them are silently ignored.
    pub pids1: Vec<i32>,
    /// `x`-grid coordinates of the `Grid`.
    pub x1: Vec<f64>,

    /// Particle ID basis for `FkTable`.
    pub pid_basis: PidBasis,
    /// Type of convolution which this operator evolves. This also determines whether
    /// [`Self::fac0`] and [`Self::fac1`] is a factorization ([`ConvType::UnpolPDF`] or
    /// [`ConvType::PolPDF`]) or a fragmentation ([`ConvType::UnpolFF`] or [`ConvType::PolFF`])
    /// scale.
    pub conv_type: ConvType,
}

/// A mapping of squared renormalization scales in `ren1` to strong couplings in `alphas`. The
/// ordering of both members defines the mapping.
pub struct AlphasTable {
    /// Renormalization scales of the `Grid`.
    pub ren1: Vec<f64>,
    /// Strong couplings corresponding to the order given in [`ren1`](Self::ren1).
    pub alphas: Vec<f64>,
}

impl AlphasTable {
    /// Create an `AlphasTable` for `grid`, varying the renormalization scale by `xir` for the
    /// strong couplings given by `alphas`. The only argument of `alphas` must be the squared
    /// renormalization scale.
    pub fn from_grid(grid: &Grid, xir: f64, alphas: &dyn Fn(f64) -> f64) -> Self {
        let mut ren1: Vec<_> = grid
            .subgrids()
            .iter()
            .flat_map(|subgrid| {
                grid.scales()
                    .ren
                    .calc(&subgrid.node_values(), grid.kinematics())
                    .into_owned()
                    .into_iter()
                    .map(|ren| xir * xir * ren)
            })
            .collect();
        ren1.sort_by(f64::total_cmp);
        ren1.dedup_by(subgrid::node_value_eq_ref_mut);
        let ren1 = ren1;
        let alphas: Vec<_> = ren1.iter().map(|&mur2| alphas(mur2)).collect();

        Self { ren1, alphas }
    }
}

fn gluon_has_pid_zero(grid: &Grid) -> bool {
    // if there are any PID zero particles ...
    grid.channels()
        .iter()
        .any(|entry| entry.entry().iter().any(|(pids, _)| pids.contains(&0)))
        // and if the particle IDs are encoded using PDG MC IDs
        && *grid.pid_basis() == PidBasis::Pdg
}

type Pid01IndexTuples = Vec<(usize, usize)>;
type Pid01Tuples = Vec<(i32, i32)>;

fn pid_slices(
    operator: &ArrayView4<f64>,
    info: &OperatorSliceInfo,
    gluon_has_pid_zero: bool,
    pid1_nonzero: &dyn Fn(i32) -> bool,
) -> Result<(Pid01IndexTuples, Pid01Tuples)> {
    // list of all non-zero PID indices
    let pid_indices: Vec<_> = (0..operator.dim().2)
        .cartesian_product(0..operator.dim().0)
        .filter(|&(pid0_idx, pid1_idx)| {
            // 1) at least one element of the operator must be non-zero, and 2) the pid must be
            // contained in some channel
            operator
                .slice(s![pid1_idx, .., pid0_idx, ..])
                .iter()
                .any(|&value| value != 0.0)
                && pid1_nonzero(if gluon_has_pid_zero && info.pids1[pid1_idx] == 21 {
                    0
                } else {
                    info.pids1[pid1_idx]
                })
        })
        .collect();

    if pid_indices.is_empty() {
        return Err(Error::General(
            "no non-zero operator found; result would be an empty FkTable".to_owned(),
        ));
    }

    // list of all non-zero (pid0, pid1) combinations
    let pids = pid_indices
        .iter()
        .map(|&(pid0_idx, pid1_idx)| {
            (
                info.pids0[pid0_idx],
                if gluon_has_pid_zero && info.pids1[pid1_idx] == 21 {
                    0
                } else {
                    info.pids1[pid1_idx]
                },
            )
        })
        .collect();

    Ok((pid_indices, pids))
}

fn operator_slices(
    operator: &ArrayView4<f64>,
    info: &OperatorSliceInfo,
    pid_indices: &[(usize, usize)],
    x1: &[f64],
) -> Result<Vec<Array2<f64>>> {
    // permutation between the grid x values and the operator x1 values
    let x1_indices: Vec<_> = x1
        .iter()
        .map(|&x1p| {
            info.x1
                .iter()
                .position(|&x1| subgrid::node_value_eq(x1p, x1))
                .ok_or_else(|| Error::General(format!("no operator for x = {x1p} found")))
        })
        // TODO: use `try_collect` once stabilized
        .collect::<Result<_>>()?;

    // create the corresponding operators accessible in the form [muf2, x0, x1]
    let operators: Vec<_> = pid_indices
        .iter()
        .map(|&(pid0_idx, pid1_idx)| {
            operator
                .slice(s![pid1_idx, .., pid0_idx, ..])
                .select(Axis(0), &x1_indices)
                .reversed_axes()
                .as_standard_layout()
                .into_owned()
        })
        .collect();

    Ok(operators)
}

type X1aX1bOpDTuple = (Vec<Vec<f64>>, Option<ArrayD<f64>>);

fn ndarray_from_subgrid_orders_slice(
    grid: &Grid,
    fac1: f64,
    kinematics: &[Kinematics],
    subgrids: &ArrayView1<SubgridEnum>,
    orders: &[Order],
    order_mask: &[bool],
    (xir, xif, xia): (f64, f64, f64),
    alphas_table: &AlphasTable,
) -> Result<X1aX1bOpDTuple> {
    // TODO: remove these assumptions from the following code
    assert_eq!(grid.kinematics()[0], Kinematics::Scale(0));
    let x_start = grid
        .kinematics()
        .iter()
        .position(|kin| matches!(kin, Kinematics::X(_)))
        // UNWRAP: there is always at least one X
        .unwrap();

    // create a Vec of all x values for each dimension
    let mut x1n: Vec<_> = kinematics
        .iter()
        .enumerate()
        .filter_map(|(idx, kin)| matches!(kin, Kinematics::X(_)).then_some(idx))
        .map(|kin_idx| {
            subgrids
                .iter()
                .enumerate()
                .filter(|&(ord_idx, subgrid)| {
                    order_mask.get(ord_idx).copied().unwrap_or(true)
                        // TODO: empty subgrids don't have node values
                        && !subgrid.is_empty()
                })
                .flat_map(|(_, subgrid)| subgrid.node_values()[kin_idx].clone())
                .collect::<Vec<_>>()
        })
        .collect();

    for x1 in &mut x1n {
        x1.sort_by(f64::total_cmp);
        x1.dedup_by(subgrid::node_value_eq_ref_mut);
    }

    let dim: Vec<_> = x1n.iter().map(Vec::len).collect();
    let mut array = ArrayD::<f64>::zeros(dim);
    let mut zero = true;
    let mut x1_idx = vec![0; grid.convolutions().len()];

    // for the same bin and channel, sum subgrids of different orders, using the right couplings
    for (subgrid, order) in subgrids
        .iter()
        .zip(orders.iter())
        .zip(order_mask.iter().chain(iter::repeat(&true)))
        .filter_map(|((subgrid, order), &enabled)| {
            (enabled && !subgrid.is_empty()).then_some((subgrid, order))
        })
    {
        let mut logs = 1.0;

        if order.logxir > 0 {
            if approx_eq!(f64, xir, 1.0, ulps = 4) {
                continue;
            }

            logs *= (xir * xir).ln();
        }

        if order.logxif > 0 {
            if approx_eq!(f64, xif, 1.0, ulps = 4) {
                continue;
            }

            logs *= (xif * xif).ln();
        }

        if order.logxia > 0 {
            if approx_eq!(f64, xia, 1.0, ulps = 4) {
                continue;
            }

            logs *= (xia * xia).ln();
        }

        let node_values = subgrid.node_values();
        let scale_dims: Vec<_> = grid
            .kinematics()
            .iter()
            .zip(&node_values)
            .filter_map(|(kin, node_values)| {
                matches!(kin, Kinematics::Scale(_)).then_some(node_values.len())
            })
            .collect();

        let x1_indices: Vec<Vec<_>> = kinematics
            .iter()
            .enumerate()
            .filter_map(|(idx, kin)| matches!(kin, Kinematics::X(_)).then_some(idx))
            .zip(&x1n)
            .map(|(kin_idx, x1)| {
                node_values[kin_idx]
                    .iter()
                    .map(|&xs| {
                        x1.iter()
                            .position(|&x| subgrid::node_value_eq(x, xs))
                            // UNWRAP: `x1n` contains all x-values, so we must find each `x`
                            .unwrap()
                    })
                    .collect()
            })
            .collect();

        let rens = grid.scales().ren.calc(&node_values, grid.kinematics());
        let facs = grid.scales().fac.calc(&node_values, grid.kinematics());

        for (indices, value) in subgrid.indexed_iter() {
            // TODO: implement evolution for non-zero fragmentation scales
            let ren = rens[grid.scales().ren.idx(&indices, &scale_dims)];
            let fac = facs[grid.scales().fac.idx(&indices, &scale_dims)];

            if !subgrid::node_value_eq(xif * xif * fac, fac1) {
                continue;
            }

            let mur2 = xir * xir * ren;

            let als = if order.alphas == 0 {
                1.0
            } else if let Some(alphas) = alphas_table
                .ren1
                .iter()
                .zip(alphas_table.alphas.iter())
                .find_map(|(&ren1, &alphas)| subgrid::node_value_eq(ren1, mur2).then_some(alphas))
            {
                alphas.powi(order.alphas.into())
            } else {
                return Err(Error::General(format!("no alphas for mur2 = {mur2} found")));
            };

            zero = false;

            for (i, &index) in indices.iter().skip(x_start).enumerate() {
                x1_idx[i] = x1_indices[i][index];
            }

            array[x1_idx.as_slice()] += als * logs * value;
        }
    }

    Ok((x1n, (!zero).then_some(array)))
}

pub(crate) fn evolve_slice(
    grid: &Grid,
    operators: &[ArrayView4<f64>],
    infos: &[OperatorSliceInfo],
    scale_values: &[f64],
    order_mask: &[bool],
    xi: (f64, f64, f64),
    alphas_table: &AlphasTable,
) -> Result<(Array3<SubgridEnum>, Vec<Channel>)> {
    let gluon_has_pid_zero = gluon_has_pid_zero(grid);

    // TODO: implement matching of different scales for different EKOs
    let mut fac1_scales: Vec<_> = infos.iter().map(|info| info.fac1).collect();
    fac1_scales.sort_by(f64::total_cmp);
    assert!(fac1_scales
        .windows(2)
        .all(|scales| subgrid::node_value_eq(scales[0], scales[1])));
    let fac1 = fac1_scales[0];

    assert_eq!(operators.len(), infos.len());
    assert_eq!(operators.len(), grid.convolutions().len());

    let (pid_indices, pids01): (Vec<_>, Vec<_>) = izip!(0..infos.len(), operators, infos)
        .map(|(d, operator, info)| {
            pid_slices(operator, info, gluon_has_pid_zero, &|pid1| {
                grid.channels()
                    .iter()
                    .flat_map(Channel::entry)
                    .any(|(pids, _)| pids[d] == pid1)
            })
        })
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .unzip();

    let mut channels0: Vec<_> = pids01
        .iter()
        .map(|pids| pids.iter().map(|&(pid0, _)| pid0))
        .multi_cartesian_product()
        .collect();
    channels0.sort_unstable();
    channels0.dedup();
    let channels0 = channels0;

    let mut sub_fk_tables = Vec::with_capacity(grid.bwfl().len() * channels0.len());

    // TODO: generalize to `n`
    let mut last_x1 = vec![Vec::new(); infos.len()];
    let mut eko_slices = vec![Vec::new(); infos.len()];
    let dim: Vec<_> = infos.iter().map(|info| info.x0.len()).collect();

    for subgrids_oc in grid.subgrids().axis_iter(Axis(1)) {
        let mut tables = vec![ArrayD::zeros(dim.clone()); channels0.len()];

        for (subgrids_o, channel1) in subgrids_oc.axis_iter(Axis(1)).zip(grid.channels()) {
            let (x1, array) = ndarray_from_subgrid_orders_slice(
                grid,
                fac1,
                grid.kinematics(),
                &subgrids_o,
                grid.orders(),
                order_mask,
                xi,
                alphas_table,
            )?;

            // skip over zero arrays to speed up evolution and avoid problems with NaNs
            let Some(array) = array else {
                continue;
            };

            for (last_x1, x1, pid_indices, slices, operator, info) in izip!(
                &mut last_x1,
                x1,
                &pid_indices,
                &mut eko_slices,
                operators,
                infos
            ) {
                if (last_x1.len() != x1.len())
                    || last_x1
                        .iter()
                        .zip(x1.iter())
                        .any(|(&lhs, &rhs)| !subgrid::node_value_eq(lhs, rhs))
                {
                    *slices = operator_slices(operator, info, pid_indices, &x1)?;
                    *last_x1 = x1;
                }
            }

            for (pids1, factor) in channel1.entry() {
                for (fk_table, ops) in
                    channels0
                        .iter()
                        .zip(tables.iter_mut())
                        .filter_map(|(pids0, fk_table)| {
                            izip!(pids0, pids1, &pids01, &eko_slices)
                                .map(|(&pid0, &pid1, pids, slices)| {
                                    pids.iter().zip(slices).find_map(|(&(p0, p1), op)| {
                                        ((p0 == pid0) && (p1 == pid1)).then_some(op)
                                    })
                                })
                                // TODO: avoid using `collect`
                                .collect::<Option<Vec<_>>>()
                                .map(|ops| (fk_table, ops))
                        })
                {
                    general_tensor_mul(*factor, array.view(), &ops, fk_table.view_mut());
                }
            }
        }

        let mut node_values = vec![scale_values.to_vec()];

        for info in infos {
            node_values.push(info.x0.clone());
        }

        sub_fk_tables.extend(tables.into_iter().map(|table| {
            ImportSubgridV1::new(
                // TODO: insert one more axis if initial frg scale unequals fac
                PackedArray::from(table.insert_axis(Axis(0)).view()),
                node_values.clone(),
            )
            .into()
        }));
    }

    Ok((
        Array1::from_iter(sub_fk_tables)
            .into_shape((1, grid.bwfl().len(), channels0.len()))
            .unwrap(),
        channels0
            .into_iter()
            .map(|c| Channel::new(vec![(c, 1.0)]))
            .collect(),
    ))
}

fn general_tensor_mul(
    factor: f64,
    array: ArrayViewD<f64>,
    ops: &[&Array2<f64>],
    mut fk_table: ArrayViewMutD<f64>,
) {
    let convolutions = array.shape().len();

    match convolutions {
        // the case `0` is not possible
        1 => {
            let array = array
                .into_dimensionality::<Ix1>()
                // UNWRAP: cannot fail due to previous `match`
                .unwrap();
            let mut fk_table = fk_table
                .into_dimensionality::<Ix1>()
                // UNWRAP: cannot fail due to previous `match`
                .unwrap();

            // fk_table += factor * ops[0] * array
            linalg::general_mat_vec_mul(factor, ops[0], &array, 1.0, &mut fk_table);
        }
        2 => {
            let array = array
                .into_dimensionality::<Ix2>()
                // UNWRAP: cannot fail due to previous `match`
                .unwrap();
            let mut fk_table = fk_table
                .into_dimensionality::<Ix2>()
                // UNWRAP: cannot fail due to previous `match`
                .unwrap();

            let mut tmp = Array2::zeros((array.shape()[0], ops[1].shape()[0]));
            // tmp = array * ops[1]^T
            linalg::general_mat_mul(1.0, &array, &ops[1].t(), 0.0, &mut tmp);
            // fk_table += factor * ops[0] * tmp
            linalg::general_mat_mul(factor, ops[0], &tmp, 1.0, &mut fk_table);
        }
        _ => {
            let (ops_0, ops_dm1) = ops.split_first().unwrap();

            for (mut fk_table_i, ops_0_i) in fk_table
                .axis_iter_mut(Axis(0))
                .zip(ops_0.axis_iter(Axis(0)))
            {
                for (array_j, ops_0_ij) in array.axis_iter(Axis(0)).zip(ops_0_i.iter()) {
                    general_tensor_mul(factor * ops_0_ij, array_j, ops_dm1, fk_table_i.view_mut());
                }
            }
        }
    }
}
