use anyhow::{bail, Result};
use cxx::UniquePtr;
use float_cmp::assert_approx_eq;
use lhapdf::Pdf;
use pineappl::boc::Order;
use pineappl::grid::Grid;
use pineappl_fastnlo::ffi::{self, fastNLOLHAPDF};
use std::path::Path;
use std::pin::Pin;

pub fn convert_into_fastnlo(
    grid: &Grid,
    _output: &Path,
    _discard_non_matching_scales: bool,
) -> Result<(UniquePtr<fastNLOLHAPDF>, Vec<bool>)> {
    let bin_info = grid.bin_info();
    let dim = bin_info.dimensions();

    if dim > 3 {
        bail!(
            "grid has {} dimensions, but fastNLO only supports up to three-dimensional distributions",
            dim
        );
    }

    let bin_limits = bin_info.limits();
    let left_bin_limits: Vec<Vec<_>> = bin_limits
        .iter()
        .map(|limits| limits.iter().map(|&(left, _)| left).collect())
        .collect();
    let right_bin_limits: Vec<Vec<_>> = bin_limits
        .iter()
        .map(|limits| limits.iter().map(|&(_, right)| right).collect())
        .collect();
    let normalizations = bin_info.normalizations();

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
        // UNWRAP: this will fail for `Grid` with no orders, but this shouldn't happen
        .unwrap();
    //let loops = orders_with_mask
    //    .iter()
    //    .filter_map(|&(Order { alphas, .. }, keep)| keep.then_some(alphas))
    //    .max()
    //    .unwrap()
    //    - lo_alphas;

    let convolutions: Vec<i32> = grid
        .convolutions()
        .iter()
        .filter_map(|conv| conv.pid())
        .collect();

    let channels: Vec<Vec<_>> = grid
        .channels()
        .iter()
        .map(|channel| {
            channel
                .entry()
                .iter()
                .map(|&(a, b, factor)| {
                    assert_approx_eq!(f64, factor, 1.0, ulps = 4);
                    ffi::pair_int_int {
                        first: a,
                        second: b,
                    }
                })
                .collect()
        })
        .collect();

    //for (fnlo_order, order) in order_mask
    //    .iter()
    //    .enumerate()
    //    .filter_map(|(index, keep)| keep.then_some(index))
    //    .enumerate()
    //{}

    let _fastnlo = ffi::make_fastnlo_create(
        // UNWRAP: negative numbers and overflow should not happen
        lo_alphas.try_into().unwrap(),
        &left_bin_limits,
        &right_bin_limits,
        &normalizations,
        // TODO: calculate channels for each order separately
        // UNWRAP: negative numbers and overflow should not happen
        channels.len().try_into().unwrap(),
        // UNWRAP: negative numbers and overflow should not happen
        channels.len().try_into().unwrap(),
        // UNWRAP: negative numbers and overflow should not happen
        channels.len().try_into().unwrap(),
        &convolutions,
        &channels,
    );
}

pub fn convolve_fastnlo(_grid: Pin<&mut fastNLOLHAPDF>, conv_funs: &mut [Pdf]) -> Vec<f64> {
    // TODO: add support for convolving an fastNLO table with two functions
    assert_eq!(conv_funs.len(), 1);

    todo!()
}
