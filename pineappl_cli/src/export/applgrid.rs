use anyhow::{bail, Result};
use cxx::{let_cxx_string, UniquePtr};
use ndarray::s;
use pineappl::grid::{Grid, Order};
use pineappl::subgrid::SubgridParams;
use pineappl_applgrid::ffi::{self, grid};
use std::iter;
use std::mem;
use std::path::Path;
use std::pin::Pin;

fn reconstruct_subgrid_params(_: &Grid) -> Result<(&str, &str, SubgridParams)> {
    // TODO: implement the general case
    Ok(("f2", "h0", SubgridParams::default()))
}

pub fn convert_into_applgrid(grid: &Grid, output: &Path) -> Result<UniquePtr<grid>> {
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

    // TODO: check that PDG MC IDs are used

    let combinations: Vec<_> = iter::once(i32::try_from(lumis).unwrap())
        .chain(grid.lumi().iter().enumerate().flat_map(|(index, entry)| {
            [
                i32::try_from(index).unwrap(),
                i32::try_from(entry.entry().len()).unwrap(),
            ]
            .into_iter()
            .chain(entry.entry().iter().flat_map(|&(a, b, factor)| {
                // TODO: if the factors aren't trivial, we have to find some other way to
                // propagate them
                assert_eq!(factor, 1.0);

                [a, b]
            }))
        }))
        .collect();

    // TODO: `id` should be a unique string for each call of this function, because the lumis are
    // registered via their pointers in APPLgrid. When this function goes out of scope, the
    // pointers will be dangling. However, this function should, in practice, only be called once
    let id = "PineAPPL-Lumi";
    let lumi_pdf = ffi::make_lumi_pdf(id, &combinations);

    let limits = &bin_info.limits()[0];
    let limits: Vec<_> = limits
        .iter()
        .map(|(left, _)| left)
        .chain(limits.last().map(|(_, right)| right))
        .copied()
        .collect();

    let (xtrans, qtrans, p) = reconstruct_subgrid_params(grid)?;

    let mut applgrid = ffi::make_empty_grid(&limits, id, 0, 0, xtrans, qtrans);

    for (appl_order, order) in grid
        .orders()
        .iter()
        .enumerate()
        .filter_map(|(index, &Order { logxir, logxif, .. })| {
            ((logxir != 0) || (logxir != 0)).then_some(index)
        })
        .enumerate()
    {
        for bin in 0..bin_info.bins() {
            // we need to convert this slice of subgrids ...
            let subgrids = grid.subgrids().slice(s![order, bin, ..]);

            // into this this igrid
            let mut igrid = ffi::make_igrid(
                p.q2_bins().try_into().unwrap(),
                p.q2_min(),
                p.q2_max(),
                p.q2_order().try_into().unwrap(),
                p.x_bins().try_into().unwrap(),
                p.x_min(),
                p.x_max(),
                p.x_order().try_into().unwrap(),
                xtrans,
                qtrans,
                lumis.try_into().unwrap(),
                false, // TODO: implement the DIS case
            )?;

            let x1 = 0;
            let x2 = 0;
            let q2 = 0;
            let lumi_results = vec![0.0; lumis];
            unsafe {
                igrid
                    .pin_mut()
                    .fill_index(x1, x2, q2, lumi_results.as_ptr());
            }

            ffi::grid_add_igrid(
                applgrid.pin_mut(),
                bin.try_into().unwrap(),
                order.try_into().unwrap(),
                igrid,
            );
        }
    }

    let_cxx_string!(filename = output.to_str().unwrap());
    let_cxx_string!(dirname = "");
    let_cxx_string!(pdfname = id);

    applgrid.pin_mut().Write(&filename, &dirname, &pdfname);

    mem::drop(lumi_pdf);

    Ok(applgrid)
}

// TODO: deduplicate this function from import
pub fn convolute_applgrid(grid: Pin<&mut grid>, pdfset: &str, member: usize) -> Vec<f64> {
    let nloops = grid.nloops();

    ffi::grid_convolute(
        grid,
        pdfset,
        member.try_into().unwrap(),
        nloops,
        1.0,
        1.0,
        1.0,
    )
}
