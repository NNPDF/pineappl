use anyhow::{bail, Result};
use cxx::{let_cxx_string, UniquePtr};
use float_cmp::assert_approx_eq;
use lhapdf::Pdf;
use ndarray::{s, Axis};
use pineappl::boc::{Kinematics, Order};
use pineappl::grid::Grid;
use pineappl::interpolation::{Interp, InterpMeth, Map, ReweightMeth};
use pineappl::subgrid::{self, Subgrid};
use pineappl_applgrid::ffi::{self, grid};
use std::f64::consts::TAU;
use std::iter;
use std::path::Path;
use std::pin::Pin;

fn reconstruct_subgrid_params(grid: &Grid, order: usize, bin: usize) -> Result<Vec<Interp>> {
    if grid
        .kinematics()
        .iter()
        .filter(|kin| matches!(kin, Kinematics::Scale(_)))
        .count()
        > 1
    {
        bail!("APPLgrid does not support grids with more than one scale");
    }

    let mut mu2_grid: Vec<_> = grid
        .subgrids()
        .slice(s![order, bin, ..])
        .iter()
        .filter(|subgrid| !subgrid.is_empty())
        .flat_map(|subgrid| {
            grid.scales()
                .fac
                .calc(&subgrid.node_values(), grid.kinematics())
                .into_owned()
        })
        .collect();
    mu2_grid.dedup_by(subgrid::node_value_eq_ref_mut);

    // TODO: implement the general case
    let mut result = vec![
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

    if let &[fac] = mu2_grid.as_slice() {
        result.insert(
            0,
            Interp::new(
                fac,
                fac,
                1,
                0,
                ReweightMeth::NoReweight,
                Map::ApplGridH0,
                InterpMeth::Lagrange,
            ),
        );
    } else {
        result.insert(
            0,
            Interp::new(
                1e2,
                1e8,
                40,
                3,
                ReweightMeth::NoReweight,
                Map::ApplGridH0,
                InterpMeth::Lagrange,
            ),
        );
    }

    Ok(result)
}

pub fn convert_into_applgrid(
    grid: &Grid,
    output: &Path,
    discard_non_matching_scales: bool,
) -> Result<(UniquePtr<grid>, Vec<bool>)> {
    let dim = grid.bwfl().dimensions();

    if dim > 1 {
        bail!(
            "grid has {} dimensions, but APPLgrid only supports one-dimensional distributions",
            dim
        );
    }

    // APPLgrid can only be used with one-dimensional consecutive bin limits
    if grid.bwfl().slices().len() != 1 {
        bail!("grid has non-consecutive bin limits, which APPLgrid does not support");
    }

    if grid.convolutions().len() > 2 {
        bail!("APPLgrid does not support grids with more than two convolutions");
    }

    let lumis = grid.channels().len();
    let has_pdf1 = !grid.convolutions().is_empty();
    let has_pdf2 = grid.convolutions().get(1).is_some();

    // TODO: check that PDG MC IDs are used

    let combinations: Vec<_> =
        iter::once(lumis.try_into().unwrap())
            .chain(
                grid.channels()
                    .iter()
                    .enumerate()
                    .flat_map(|(index, entry)| {
                        [
                            index.try_into().unwrap(),
                            entry.entry().len().try_into().unwrap(),
                        ]
                        .into_iter()
                        .chain(entry.entry().iter().flat_map(|&(ref pids, factor)| {
                            // TODO: if the factors aren't trivial, we have to find some other way
                            // to propagate them
                            assert_approx_eq!(f64, factor, 1.0, ulps = 4);

                            pids.iter()
                                .copied()
                                .chain(iter::repeat(0))
                                .take(2)
                                .collect::<Vec<_>>()
                        }))
                    }),
            )
            .collect();

    // `id` must end with '.config' for APPLgrid to know its type is `lumi_pdf`
    let id = "PineAPPL-Lumi.config";
    // this object is managed by APPLgrid internally
    ffi::make_lumi_pdf(id, &combinations).into_raw();

    let limits: Vec<_> = grid
        .bwfl()
        .bins()
        .iter()
        .map(|bin| {
            // TODO: instead of `bin.limits()[0]` we should use `bin.fill_limits()`, but this
            // requires changing the normalization
            bin.limits()[0].0
        })
        .chain(Some(
            grid.bwfl()
                .bins()
                .last()
                // UNWRAP: every `grid` should have at least one bin
                .unwrap()
                .limits()[0]
                .1,
        ))
        .collect();

    let order_mask = Order::create_mask(grid.orders(), 3, 0, false);
    let orders_with_mask: Vec<_> = grid
        .orders()
        .iter()
        .cloned()
        .zip(order_mask.iter().copied())
        .collect();
    let lo_alphas = orders_with_mask
        .iter()
        .filter_map(|&(Order { alphas, .. }, keep)| keep.then_some(alphas))
        .min()
        .unwrap();
    let loops = orders_with_mask
        .iter()
        .filter_map(|&(Order { alphas, .. }, keep)| keep.then_some(alphas))
        .max()
        .unwrap()
        - lo_alphas;

    let mut applgrid =
        ffi::make_empty_grid(&limits, id, lo_alphas.into(), loops.into(), "f2", "h0");

    for (appl_order, order) in order_mask
        .iter()
        .enumerate()
        .filter_map(|(index, keep)| keep.then_some(index))
        .enumerate()
    {
        let factor = TAU.powi(grid.orders()[order].alphas.into());

        for (bin, subgrids) in grid
            .subgrids()
            .index_axis(Axis(0), order)
            .axis_iter(Axis(0))
            .enumerate()
        {
            let interps = reconstruct_subgrid_params(grid, order, bin)?;
            // TODO: support DIS case
            assert_eq!(interps.len(), 3);

            // TODO: make sure interps[1] is the same as interps[2]

            let mut igrid = ffi::make_igrid(
                interps[0].nodes().try_into().unwrap(),
                interps[0].min(),
                interps[0].max(),
                interps[0].order().try_into().unwrap(),
                interps[1].nodes().try_into().unwrap(),
                interps[1].min(),
                interps[1].max(),
                interps[1].order().try_into().unwrap(),
                match interps[1].map() {
                    Map::ApplGridF2 => "f2",
                    map @ Map::ApplGridH0 => panic!("export does not support {map:?}"),
                },
                match interps[0].map() {
                    Map::ApplGridH0 => "h0",
                    map @ Map::ApplGridF2 => panic!("export does not support {map:?}"),
                },
                grid.channels().len().try_into().unwrap(),
                grid.convolutions().len() == 1,
            );
            let appl_q2: Vec<_> = (0..igrid.Ntau()).map(|i| igrid.getQ2(i)).collect();
            let appl_x1: Vec<_> = (0..igrid.Ny1()).map(|i| igrid.getx1(i)).collect();
            let appl_x2: Vec<_> = (0..igrid.Ny2()).map(|i| igrid.getx2(i)).collect();

            for (channel, subgrid) in subgrids
                .iter()
                .enumerate()
                .filter(|(_, subgrid)| !subgrid.is_empty())
            {
                let appl_q2_idx: Vec<_> = grid.scales().fac.calc(&subgrid.node_values(), grid.kinematics())
                    .iter()
                    .map(|&fac| {
                        appl_q2
                            .iter()
                            .position(|&x| subgrid::node_value_eq(x, fac))
                            .map_or_else(
                                || {
                                    if discard_non_matching_scales {
                                        Ok(-1)
                                    } else {
                                        bail!(
                                            "factorization scale muf2 = {fac} not found in APPLgrid",
                                        )
                                    }
                                },
                                |idx| Ok(idx.try_into().unwrap()),
                            )
                    })
                    .collect::<Result<_>>()?;

                // in the DIS case APPLgrid always uses the first x dimension

                let (x1_grid, x2_grid) = if has_pdf1 && has_pdf2 {
                    (
                        grid.kinematics()
                            .iter()
                            .zip(subgrid.node_values())
                            .find_map(|(kin, node_values)| {
                                matches!(kin, &Kinematics::X(idx) if idx == 0)
                                    .then_some(node_values)
                            })
                            // TODO: convert this into an error
                            .unwrap(),
                        grid.kinematics()
                            .iter()
                            .zip(subgrid.node_values())
                            .find_map(|(kin, node_values)| {
                                matches!(kin, &Kinematics::X(idx) if idx == 1)
                                    .then_some(node_values)
                            })
                            // TODO: convert this into an error
                            .unwrap(),
                    )
                } else if has_pdf1 {
                    (
                        grid.kinematics()
                            .iter()
                            .zip(subgrid.node_values())
                            .find_map(|(kin, node_values)| {
                                matches!(kin, &Kinematics::X(idx) if idx == 0)
                                    .then_some(node_values)
                            })
                            // TODO: convert this into an error
                            .unwrap(),
                        Vec::new(),
                    )
                } else {
                    (
                        grid.kinematics()
                            .iter()
                            .zip(subgrid.node_values())
                            .find_map(|(kin, node_values)| {
                                matches!(kin, &Kinematics::X(idx) if idx == 1)
                                    .then_some(node_values)
                            })
                            // TODO: convert this into an error
                            .unwrap(),
                        Vec::new(),
                    )
                };

                let appl_x1_idx: Vec<_> = x1_grid
                    .iter()
                    .map(|&x1| {
                        appl_x1
                            .iter()
                            .position(|&x| subgrid::node_value_eq(x, x1))
                            .map_or_else(
                                || bail!("momentum fraction x1 = {x1} not found in APPLgrid"),
                                |idx| Ok(idx.try_into().unwrap()),
                            )
                    })
                    .collect::<Result<_>>()?;
                let appl_x2_idx: Vec<_> = x2_grid
                    .iter()
                    .map(|&x2| {
                        appl_x2
                            .iter()
                            .position(|&x| subgrid::node_value_eq(x, x2))
                            .map_or_else(
                                || bail!("momentum fraction x2 = {x2} not found in APPLgrid"),
                                |idx| Ok(idx.try_into().unwrap()),
                            )
                    })
                    .collect::<Result<_>>()?;

                let mut weightgrid = ffi::igrid_weightgrid(igrid.pin_mut(), channel);

                for (indices, value) in subgrid.indexed_iter() {
                    // TODO: here we assume that all X are consecutive starting from the second
                    // element and are in ascending order
                    let iq2 = indices[0];
                    let appl_q2_idx = appl_q2_idx[iq2];

                    if appl_q2_idx == -1 {
                        if value != 0.0 {
                            println!(
                                "WARNING: discarding non-matching scale muf2 = {}",
                                grid.scales()
                                    .fac
                                    .calc(&subgrid.node_values(), grid.kinematics())[iq2]
                            );
                        }

                        continue;
                    }

                    ffi::sparse_matrix_set(
                        weightgrid.as_mut(),
                        appl_q2_idx,
                        appl_x1_idx[indices[1]],
                        if has_pdf1 && has_pdf2 {
                            appl_x2_idx[indices[2]]
                        } else {
                            0
                        },
                        factor * value,
                    );
                }

                // TODO: is this call needed?
                weightgrid.trim();
            }

            igrid.pin_mut().setlimits();

            unsafe {
                applgrid.pin_mut().add_igrid(
                    bin.try_into().unwrap(),
                    appl_order.try_into().unwrap(),
                    igrid.into_raw(),
                );
            }
        }
    }

    applgrid.pin_mut().include_photon(true);

    let_cxx_string!(filename = output.to_str().unwrap());
    let_cxx_string!(empty = "");

    applgrid.pin_mut().Write(&filename, &empty, &empty);

    Ok((applgrid, order_mask))
}

// TODO: deduplicate this function from import
pub fn convolve_applgrid(grid: Pin<&mut grid>, conv_funs: &mut [Pdf]) -> Vec<f64> {
    let nloops = grid.nloops();

    // TODO: add support for convolving an APPLgrid with two functions
    assert_eq!(conv_funs.len(), 1);

    pineappl_applgrid::grid_convolve_with_one(grid, &mut conv_funs[0], nloops, 1.0, 1.0, 1.0)
}
