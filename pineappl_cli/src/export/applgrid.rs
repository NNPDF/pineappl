use anyhow::{bail, Result};
use cxx::{let_cxx_string, UniquePtr};
use ndarray::Axis;
use pineappl::grid::{Grid, Order};
use pineappl::subgrid::{Subgrid, SubgridParams};
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

    // `id` must end with '.config' for APPLgrid to know its type is `lumi_pdf`
    let id = "PineAPPL-Lumi.config";
    let lumi_pdf = ffi::make_lumi_pdf(id, &combinations).into_raw();

    let limits = &bin_info.limits();
    let limits: Vec<_> = limits
        .iter()
        .map(|vec| vec[0].0)
        .chain(limits.last().map(|vec| vec[0].1))
        .collect();

    let (xtrans, qtrans, p) = reconstruct_subgrid_params(grid)?;

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

    let mut applgrid = ffi::make_new_grid(
        &limits,
        p.q2_bins().try_into().unwrap(),
        p.q2_min(),
        p.q2_max(),
        p.q2_order().try_into().unwrap(),
        p.x_bins().try_into().unwrap(),
        p.x_min(),
        p.x_max(),
        p.x_order().try_into().unwrap(),
        id,
        lo_alphas.try_into().unwrap(),
        loops.try_into().unwrap(),
        xtrans,
        qtrans,
        false, // TODO: implement the DIS case
    );

    for Order {
        alphas,
        alpha,
        logxir,
        logxif,
    } in orders_with_mask
        .iter()
        .filter_map(|(order, keep)| (!keep).then_some(order.clone()))
    {
        println!("WARNING: the order O(as^{alphas} a^{alpha} lr^{logxir} lf^{logxif}) isn't supported by APPLgrid and will be skipped.");
    }

    for (appl_order, order) in order_mask
        .into_iter()
        .enumerate()
        .filter_map(|(index, keep)| keep.then_some(index))
        .enumerate()
    {
        for ((bin, lumi), subgrid) in grid.subgrids().index_axis(Axis(0), order).indexed_iter() {
            let mut igrid = ffi::grid_get_igrid(applgrid.pin_mut(), appl_order, bin);
            let mut weightgrid = ffi::igrid_weightgrid(igrid.as_mut(), lumi);

            for ((x1, x2, q2), value) in subgrid.indexed_iter() {
                ffi::sparse_matrix_set(
                    weightgrid.as_mut(),
                    x1.try_into().unwrap(),
                    x2.try_into().unwrap(),
                    q2.try_into().unwrap(),
                    value,
                );
            }

            // TODO: is this call needed?
            weightgrid.as_mut().trim();
        }

        for bin in 0..bin_info.bins() {
            let mut igrid = ffi::grid_get_igrid(applgrid.pin_mut(), appl_order, bin);
            // TODO: is this call needed?
            igrid.as_mut().setlimits();
        }
    }

    let_cxx_string!(filename = output.to_str().unwrap());
    let_cxx_string!(empty = "");

    applgrid.pin_mut().Write(&filename, &empty, &empty);

    // silence the warning that `lumi_pdf` isn't used
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
