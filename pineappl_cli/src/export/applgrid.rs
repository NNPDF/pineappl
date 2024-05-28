use anyhow::{anyhow, bail, Result};
use cxx::{let_cxx_string, UniquePtr};
use float_cmp::approx_eq;
use ndarray::{s, Axis};
use pineappl::grid::{Convolution, Grid};
use pineappl::order::Order;
use pineappl::subgrid::{Mu2, Subgrid, SubgridParams};
use pineappl_applgrid::ffi::{self, grid};
use std::borrow::Cow;
use std::f64::consts::TAU;
use std::iter;
use std::path::Path;
use std::pin::Pin;

fn reconstruct_subgrid_params(grid: &Grid, order: usize, bin: usize) -> Result<SubgridParams> {
    let mut result = SubgridParams::default();

    let mu2_grid: Vec<_> = grid
        .subgrids()
        .slice(s![order, bin, ..])
        .iter()
        .map(|subgrid| {
            subgrid
                .mu2_grid()
                .iter()
                .map(|&Mu2 { ren, fac }| {
                    if !approx_eq!(f64, ren, fac, ulps = 128) {
                        bail!("subgrid has mur2 != muf2, which APPLgrid does not support");
                    }

                    Ok(fac)
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<_>>()?;
    let mut mu2_grid: Vec<_> = mu2_grid.into_iter().flatten().collect();
    mu2_grid.dedup_by(|a, b| approx_eq!(f64, *a, *b, ulps = 128));
    let mu2_grid = mu2_grid.as_slice();

    if let &[fac] = mu2_grid {
        result.set_q2_bins(1);
        result.set_q2_max(fac);
        result.set_q2_min(fac);
        result.set_q2_order(0);
    }

    // TODO: implement the general case
    Ok(result)
}

pub fn convert_into_applgrid(
    grid: &Grid,
    output: &Path,
    discard_non_matching_scales: bool,
) -> Result<(UniquePtr<grid>, Vec<bool>)> {
    let bin_info = grid.bin_info();
    let dim = bin_info.dimensions();

    if dim > 1 {
        bail!(
            "grid has {} dimensions, but APPLgrid only supports one-dimensional distributions",
            dim
        );
    }

    if bin_info.slices().len() != 1 {
        bail!("grid has non-consecutive bin limits, which APPLgrid does not support");
    }

    let lumis = grid.lumi().len();
    let has_pdf1 = grid.convolutions()[0] != Convolution::None;
    let has_pdf2 = grid.convolutions()[1] != Convolution::None;

    // TODO: check that PDG MC IDs are used

    let combinations: Vec<_> = iter::once(lumis.try_into().unwrap())
        .chain(grid.lumi().iter().enumerate().flat_map(|(index, entry)| {
            [
                index.try_into().unwrap(),
                entry.entry().len().try_into().unwrap(),
            ]
            .into_iter()
            .chain(entry.entry().iter().flat_map(|&(a, b, factor)| {
                // TODO: if the factors aren't trivial, we have to find some other way to
                // propagate them
                assert_eq!(factor, 1.0);

                match (has_pdf1, has_pdf2) {
                    (true, true) => [a, b],
                    (true, false) => [a, 0],
                    (false, true) => [b, 0],
                    (false, false) => unreachable!(),
                }
            }))
        }))
        .collect();

    // `id` must end with '.config' for APPLgrid to know its type is `lumi_pdf`
    let id = "PineAPPL-Lumi.config";
    // this object is managed by APPLgrid internally
    ffi::make_lumi_pdf(id, &combinations).into_raw();

    let limits = &bin_info.limits();
    let limits: Vec<_> = limits
        .iter()
        .map(|vec| vec[0].0)
        .chain(limits.last().map(|vec| vec[0].1))
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

    let mut applgrid = ffi::make_empty_grid(
        &limits,
        id,
        lo_alphas.try_into().unwrap(),
        loops.try_into().unwrap(),
        "f2",
        "h0",
    );

    for (appl_order, order) in order_mask
        .iter()
        .enumerate()
        .filter_map(|(index, keep)| keep.then_some(index))
        .enumerate()
    {
        let factor = TAU.powi(grid.orders()[order].alphas.try_into().unwrap());

        for (bin, subgrids) in grid
            .subgrids()
            .index_axis(Axis(0), order)
            .axis_iter(Axis(0))
            .enumerate()
        {
            let p = reconstruct_subgrid_params(grid, order, bin)?;

            let mut igrid = ffi::make_igrid(
                p.q2_bins().try_into().unwrap(),
                p.q2_min(),
                p.q2_max(),
                p.q2_order().try_into().unwrap(),
                p.x_bins().try_into().unwrap(),
                p.x_min(),
                p.x_max(),
                p.x_order().try_into().unwrap(),
                "f2",
                "h0",
                grid.lumi().len().try_into().unwrap(),
                has_pdf1 != has_pdf2,
            );
            let appl_q2: Vec<_> = (0..igrid.Ntau()).map(|i| igrid.getQ2(i)).collect();
            let appl_x1: Vec<_> = (0..igrid.Ny1()).map(|i| igrid.getx1(i)).collect();
            let appl_x2: Vec<_> = (0..igrid.Ny2()).map(|i| igrid.getx2(i)).collect();

            for (lumi, subgrid) in subgrids.iter().enumerate() {
                let appl_q2_idx: Vec<_> = subgrid
                    .mu2_grid()
                    .iter()
                    .map(|&Mu2 { ren, fac }| {
                        if !approx_eq!(f64, ren, fac, ulps = 128) {
                            bail!("subgrid has mur2 != muf2, which APPLgrid does not support");
                        }
                        appl_q2
                            .iter()
                            .position(|&x| approx_eq!(f64, x, fac, ulps = 128))
                            .map_or_else(
                                || {
                                    if discard_non_matching_scales {
                                        Ok(-1)
                                    } else {
                                        Err(anyhow!(
                                            "factorization scale muf2 = {} not found in APPLgrid",
                                            fac
                                        ))
                                    }
                                },
                                |idx| Ok(idx.try_into().unwrap()),
                            )
                    })
                    .collect::<Result<_>>()?;

                // in the DIS case APPLgrid always uses the first x dimension
                let (x1_grid, x2_grid) = if has_pdf1 && has_pdf2 {
                    (subgrid.x1_grid(), subgrid.x2_grid())
                } else if has_pdf1 {
                    (subgrid.x1_grid(), Cow::Owned(vec![]))
                } else {
                    (subgrid.x2_grid(), Cow::Owned(vec![]))
                };

                let appl_x1_idx: Vec<_> = x1_grid
                    .iter()
                    .map(|&x1| {
                        appl_x1
                            .iter()
                            .position(|&x| approx_eq!(f64, x, x1, ulps = 128))
                            .map_or_else(
                                || {
                                    Err(anyhow!(
                                        "momentum fraction x1 = {} not found in APPLgrid",
                                        x1
                                    ))
                                },
                                |idx| Ok(idx.try_into().unwrap()),
                            )
                    })
                    .collect::<Result<_>>()?;
                let appl_x2_idx: Vec<_> = x2_grid
                    .iter()
                    .map(|&x2| {
                        appl_x2
                            .iter()
                            .position(|&x| approx_eq!(f64, x, x2, ulps = 128))
                            .map_or_else(
                                || {
                                    Err(anyhow!(
                                        "momentum fraction x2 = {} not found in APPLgrid",
                                        x2
                                    ))
                                },
                                |idx| Ok(idx.try_into().unwrap()),
                            )
                    })
                    .collect::<Result<_>>()?;

                let mut weightgrid = ffi::igrid_weightgrid(igrid.pin_mut(), lumi);

                for ((iq2, ix1, ix2), value) in subgrid.indexed_iter() {
                    let appl_q2_idx = appl_q2_idx[iq2];

                    if appl_q2_idx == -1 {
                        if value != 0.0 {
                            println!(
                                "WARNING: discarding non-matching scale muf2 = {}",
                                subgrid.mu2_grid()[iq2].fac
                            );
                        }

                        continue;
                    }

                    ffi::sparse_matrix_set(
                        weightgrid.as_mut(),
                        appl_q2_idx,
                        appl_x1_idx[ix1],
                        if has_pdf1 && has_pdf2 {
                            appl_x2_idx[ix2]
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
pub fn convolve_applgrid(grid: Pin<&mut grid>, pdfset: &str, member: usize) -> Vec<f64> {
    let nloops = grid.nloops();

    ffi::grid_convolve(
        grid,
        pdfset,
        member.try_into().unwrap(),
        nloops,
        1.0,
        1.0,
        1.0,
    )
}
