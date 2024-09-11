//! The C-language interface for `PineAPPL`.
//!
//! The `PineAPPL` Application Programming Interface for the C language (CAPI) defines types and
//! functions that allow `PineAPPL` to be used without having to write Rust code, and instead
//! offering
//!
//! * C or C++, or
//! * Fortran as programming languages. Fortran is supported using the `iso_c_binding` module to
//!   wrap CAPI's functions (see the [Fortran example] in the repository).
//!
//! Note that the CAPI only defines a subset of functionality available in the Rust module, which
//! is extended on an as-needed basis. If you happen to miss a function please open a new [issue]
//! on the github page.
//!
//! [Fortran example]: https://github.com/NNPDF/pineappl/tree/master/examples/fortran
//! [issue]: https://github.com/NNPDF/pineappl/issues
//!
//! # Using and linking the CAPI
//!
//! To use `PineAPPL`'s CAPI you must include its header `pineappl_capi.h`. To successfully build
//! your program, you should rely on [`pkgconf`] or [`pkg-conf`], as follows:
//!
//! ```shell
//! pkg-config --cflags pineappl_capi
//! pkg-config --libs pineappl_capi
//! ```
//!
//! This will read `PineAPPL`'s `.pc` file and print the neccessary `CFLAGS` (`--cflags`) and
//! linker flags (`--libs`). This procedure is supported by many build systems, such as
//!
//! * [Autoconf]/[Automake], using the `PKG_CHECK_MODULES` macro, see the [Autotools mythbuster]
//!   page for correct usage,
//! * [CMake], using [`FindPkgConfig`] and
//! * [meson], using [`dependency`]; it usually suffices to write `dependency('pineappl_capi')`.
//!
//! [Autoconf]: https://www.gnu.org/software/autoconf/
//! [Automake]: https://www.gnu.org/software/automake/
//! [Autotools mythbuster]: https://autotools.io/pkgconfig/pkg_check_modules.html
//! [CMake]: https://cmake.org/
//! [`dependency`]: https://mesonbuild.com/Reference-manual.html#dependency
//! [`FindPkgConfig`]: https://cmake.org/cmake/help/latest/module/FindPkgConfig.html
//! [meson]: https://mesonbuild.com/
//! [`pkgconf`]: https://github.com/pkgconf/pkgconf
//! [`pkg-conf`]: https://www.freedesktop.org/wiki/Software/pkg-config
//!
//! # Strings and other types
//!
//! Strings used in this library are have the Rust type `*const c_char` or `*mut c_char` and
//! correspond to the C types `const char*` and `char*`, respectively. The strings are assumed to
//! be encoded using UTF-8, which Rust uses internally.
//!
//! All other built-in types in Rust (`f64`, `usize`, `i32`, `u32`, ...) correspond to types in C
//! as defined by the [translation tables] of cbindgen. If in doubt, consult the generated header
//! `pineappl_capi.h`.
//!
//! [translation tables]: https://github.com/eqrion/cbindgen/blob/master/docs.md#std-types

use itertools::izip;
use pineappl::bin::BinRemapper;
use pineappl::boc::{Channel, Kinematics, Order};
use pineappl::convolutions::{Convolution, LumiCache};
use pineappl::grid::{Grid, GridOptFlags};
use pineappl::interpolation::{Interp, InterpMeth, Map, ReweightMeth};
use std::collections::HashMap;
use std::ffi::{CStr, CString};
use std::fs::File;
use std::mem;
use std::os::raw::{c_char, c_void};
use std::path::Path;
use std::slice;

// TODO: make sure no `panic` calls leave functions marked as `extern "C"`

fn grid_interpolation_params(key_vals: Option<&KeyVal>) -> Vec<Interp> {
    let mut q2_min = 1e2;
    let mut q2_max = 1e8;
    let mut q2_nodes = 40;
    let mut q2_order = 3;
    let mut x1_min = 2e-7;
    let mut x1_max = 1.0;
    let mut x1_nodes = 50;
    let mut x1_order = 3;
    let mut x2_min = 2e-7;
    let mut x2_max = 1.0;
    let mut x2_nodes = 50;
    let mut x2_order = 3;
    let mut reweight = ReweightMeth::ApplGridX;

    if let Some(keyval) = key_vals {
        if let Some(value) = keyval
            .ints
            .get("q2_bins")
            .or_else(|| keyval.ints.get("nq2"))
        {
            q2_nodes = usize::try_from(*value).unwrap();
        }

        if let Some(value) = keyval
            .doubles
            .get("q2_max")
            .or_else(|| keyval.doubles.get("q2max"))
        {
            q2_max = *value;
        }

        if let Some(value) = keyval
            .doubles
            .get("q2_min")
            .or_else(|| keyval.doubles.get("q2min"))
        {
            q2_min = *value;
        }

        if let Some(value) = keyval
            .ints
            .get("q2_order")
            .or_else(|| keyval.ints.get("q2order"))
        {
            q2_order = usize::try_from(*value).unwrap();
        }

        if let Some(value) = keyval.bools.get("reweight") {
            if !value {
                reweight = ReweightMeth::NoReweight;
            }
        }

        if let Some(value) = keyval.ints.get("x_bins").or_else(|| keyval.ints.get("nx")) {
            let value = usize::try_from(*value).unwrap();
            x1_nodes = value;
            x2_nodes = value;
        }

        if let Some(value) = keyval
            .doubles
            .get("x_max")
            .or_else(|| keyval.doubles.get("xmax"))
        {
            x1_max = *value;
            x2_max = *value;
        }

        if let Some(value) = keyval
            .doubles
            .get("x_min")
            .or_else(|| keyval.doubles.get("xmin"))
        {
            x1_min = *value;
            x2_min = *value;
        }

        if let Some(value) = keyval
            .ints
            .get("x_order")
            .or_else(|| keyval.ints.get("xorder"))
        {
            let value = usize::try_from(*value).unwrap();
            x1_order = value;
            x2_order = value;
        }

        if let Some(value) = keyval.ints.get("x1_bins") {
            x1_nodes = usize::try_from(*value).unwrap();
        }

        if let Some(value) = keyval.doubles.get("x1_max") {
            x1_max = *value;
        }

        if let Some(value) = keyval.doubles.get("x1_min") {
            x1_min = *value;
        }

        if let Some(value) = keyval.ints.get("x1_order") {
            x1_order = usize::try_from(*value).unwrap();
        }

        if let Some(value) = keyval.ints.get("x2_bins") {
            x2_nodes = usize::try_from(*value).unwrap();
        }

        if let Some(value) = keyval.doubles.get("x2_max") {
            x2_max = *value;
        }

        if let Some(value) = keyval.doubles.get("x2_min") {
            x2_min = *value;
        }

        if let Some(value) = keyval.ints.get("x2_order") {
            x2_order = usize::try_from(*value).unwrap();
        }
    }

    vec![
        Interp::new(
            q2_min,
            q2_max,
            q2_nodes,
            q2_order,
            ReweightMeth::NoReweight,
            Map::ApplGridH0,
            InterpMeth::Lagrange,
        ),
        Interp::new(
            x1_min,
            x1_max,
            x1_nodes,
            x1_order,
            reweight,
            Map::ApplGridF2,
            InterpMeth::Lagrange,
        ),
        Interp::new(
            x2_min,
            x2_max,
            x2_nodes,
            x2_order,
            reweight,
            Map::ApplGridF2,
            InterpMeth::Lagrange,
        ),
    ]
}

/// Type for defining a luminosity function.
#[derive(Default)]
pub struct Lumi(Vec<Channel>);

/// Returns the number of bins in `grid`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
#[must_use]
pub unsafe extern "C" fn pineappl_grid_bin_count(grid: *const Grid) -> usize {
    let grid = unsafe { &*grid };

    grid.bin_info().bins()
}

/// Returns the number of dimensions of the bins this grid has.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_bin_dimensions(grid: *const Grid) -> usize {
    let grid = unsafe { &*grid };

    grid.bin_info().dimensions()
}

/// Stores the bin sizes of `grid` in `bin_sizes`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. The parameter `bin_sizes` must point to an array that is as
/// long as `grid` has bins.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_bin_normalizations(grid: *const Grid, bin_sizes: *mut f64) {
    let grid = unsafe { &*grid };
    let sizes = grid.bin_info().normalizations();
    let bin_sizes = unsafe { slice::from_raw_parts_mut(bin_sizes, sizes.len()) };

    for (i, size) in sizes.iter().copied().enumerate() {
        bin_sizes[i] = size;
    }
}

/// Write the left limits for the specified dimension into `left`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. The parameter `left` must point to an array that is as large
/// as `grid` has bins. If `dimension` is larger or equal the number of dimensions for this grid,
/// nothing is written into `left`, the result is undefined.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_bin_limits_left(
    grid: *const Grid,
    dimension: usize,
    left: *mut f64,
) {
    let grid = unsafe { &*grid };
    let limits = grid.bin_info().left(dimension);
    let left_limits = unsafe { slice::from_raw_parts_mut(left, limits.len()) };

    left_limits.copy_from_slice(&limits);
}

/// Write the right limits for the specified dimension into `right`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. The parameter `right` must point to an array that is as
/// large as `grid` has bins. If `dimension` is larger or equal the number of dimensions for this
/// grid, nothing is written into `right`, the result is undefined.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_bin_limits_right(
    grid: *const Grid,
    dimension: usize,
    right: *mut f64,
) {
    let grid = unsafe { &*grid };
    let limits = grid.bin_info().right(dimension);
    let right_limits = unsafe { slice::from_raw_parts_mut(right, limits.len()) };

    right_limits.copy_from_slice(&limits);
}

/// Returns a cloned object of `grid`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_clone(grid: *const Grid) -> Box<Grid> {
    let grid = unsafe { &*grid };

    Box::new(grid.clone())
}

/// Wrapper for [`pineappl_grid_convolve_with_one`].
///
/// # Safety
///
/// See [`pineappl_grid_convolve_with_one`].
#[deprecated(
    since = "0.8.0",
    note = "please use `pineappl_grid_convolve_with_one` instead"
)]
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_convolute_with_one(
    grid: *const Grid,
    pdg_id: i32,
    xfx: extern "C" fn(pdg_id: i32, x: f64, q2: f64, state: *mut c_void) -> f64,
    alphas: extern "C" fn(q2: f64, state: *mut c_void) -> f64,
    state: *mut c_void,
    order_mask: *const bool,
    channel_mask: *const bool,
    xi_ren: f64,
    xi_fac: f64,
    results: *mut f64,
) {
    unsafe {
        pineappl_grid_convolve_with_one(
            grid,
            pdg_id,
            xfx,
            alphas,
            state,
            order_mask,
            channel_mask,
            xi_ren,
            xi_fac,
            results,
        );
    }
}

/// Wrapper for [`pineappl_grid_convolve_with_two`].
///
/// # Safety
///
/// See [`pineappl_grid_convolve_with_two`].
#[deprecated(
    since = "0.8.0",
    note = "please use `pineappl_grid_convolve_with_two` instead"
)]
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_convolute_with_two(
    grid: *const Grid,
    pdg_id1: i32,
    xfx1: extern "C" fn(pdg_id: i32, x: f64, q2: f64, state: *mut c_void) -> f64,
    pdg_id2: i32,
    xfx2: extern "C" fn(pdg_id: i32, x: f64, q2: f64, state: *mut c_void) -> f64,
    alphas: extern "C" fn(q2: f64, state: *mut c_void) -> f64,
    state: *mut c_void,
    order_mask: *const bool,
    channel_mask: *const bool,
    xi_ren: f64,
    xi_fac: f64,
    results: *mut f64,
) {
    unsafe {
        pineappl_grid_convolve_with_two(
            grid,
            pdg_id1,
            xfx1,
            pdg_id2,
            xfx2,
            alphas,
            state,
            order_mask,
            channel_mask,
            xi_ren,
            xi_fac,
            results,
        );
    }
}

/// Convolutes the specified grid with the PDF `xfx`, which is the PDF for a hadron with the PDG id
/// `pdg_id`, and strong coupling `alphas`. These functions must evaluate the PDFs for the given
/// `x` and `q2` for the parton with the given PDG id, `pdg_id`, and return the result. Note that
/// the value must be the PDF multiplied with its argument `x`. The value of the pointer `state`
/// provided to these functions is the same one given to this function. The parameter `order_mask`
/// must be as long as there are perturbative orders contained in `grid` and is used to selectively
/// disable (`false`) or enable (`true`) individual orders. If `order_mask` is set to `NULL`, all
/// orders are active. The parameter `channel_mask` can be used similarly, but must be as long as
/// the channels `grid` was created with has entries, or `NULL` to enable all channels. The values
/// `xi_ren` and `xi_fac` can be used to vary the renormalization and factorization from its
/// central value, which corresponds to `1.0`. After convolution of the grid with the PDFs the
/// differential cross section for each bin is written into `results`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. The function pointers `xfx1`, `xfx2`, and `alphas` must not
/// be null pointers and point to valid functions. The parameters `order_mask` and `channel_mask`
/// must either be null pointers or point to arrays that are as long as `grid` has orders and
/// channels, respectively. Finally, `results` must be as long as `grid` has bins.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_convolve_with_one(
    grid: *const Grid,
    pdg_id: i32,
    xfx: extern "C" fn(pdg_id: i32, x: f64, q2: f64, state: *mut c_void) -> f64,
    alphas: extern "C" fn(q2: f64, state: *mut c_void) -> f64,
    state: *mut c_void,
    order_mask: *const bool,
    channel_mask: *const bool,
    xi_ren: f64,
    xi_fac: f64,
    results: *mut f64,
) {
    let grid = unsafe { &*grid };
    let mut pdf = |id, x, q2| xfx(id, x, q2, state);
    let mut als = |q2| alphas(q2, state);
    let order_mask = if order_mask.is_null() {
        vec![]
    } else {
        unsafe { slice::from_raw_parts(order_mask, grid.orders().len()) }.to_owned()
    };
    let channel_mask = if channel_mask.is_null() {
        vec![]
    } else {
        unsafe { slice::from_raw_parts(channel_mask, grid.channels().len()) }.to_vec()
    };
    let results = unsafe { slice::from_raw_parts_mut(results, grid.bin_info().bins()) };
    let mut lumi_cache = LumiCache::with_one(pdg_id, &mut pdf, &mut als);

    results.copy_from_slice(&grid.convolve(
        &mut lumi_cache,
        &order_mask,
        &[],
        &channel_mask,
        &[(xi_ren, xi_fac, 1.0)],
    ));
}

/// Convolutes the specified grid with the PDFs `xfx1` and `xfx2`, which are the PDFs of hadrons
/// with PDG ids `pdg_id1` and `pdg_id2`, respectively, and strong coupling `alphas`. These
/// functions must evaluate the PDFs for the given `x` and `q2` for the parton with the given PDG
/// id, `pdg_id`, and return the result. Note that the value must be the PDF multiplied with its
/// argument `x`. The value of the pointer `state` provided to these functions is the same one
/// given to this function. The parameter `order_mask` must be as long as there are perturbative
/// orders contained in `grid` and is used to selectively disable (`false`) or enable (`true`)
/// individual orders. If `order_mask` is set to `NULL`, all orders are active. The parameter
/// `channel_mask` can be used similarly, but must be as long as the channels `grid` was created
/// with has entries, or `NULL` to enable all channels. The values `xi_ren` and `xi_fac` can be
/// used to vary the renormalization and factorization from its central value, which corresponds to
/// `1.0`. After convolution of the grid with the PDFs the differential cross section for each bin
/// is written into `results`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. The function pointers `xfx1`, `xfx2`, and `alphas` must not
/// be null pointers and point to valid functions. The parameters `order_mask` and `channel_mask`
/// must either be null pointers or point to arrays that are as long as `grid` has orders and
/// channels, respectively. Finally, `results` must be as long as `grid` has bins.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_convolve_with_two(
    grid: *const Grid,
    pdg_id1: i32,
    xfx1: extern "C" fn(pdg_id: i32, x: f64, q2: f64, state: *mut c_void) -> f64,
    pdg_id2: i32,
    xfx2: extern "C" fn(pdg_id: i32, x: f64, q2: f64, state: *mut c_void) -> f64,
    alphas: extern "C" fn(q2: f64, state: *mut c_void) -> f64,
    state: *mut c_void,
    order_mask: *const bool,
    channel_mask: *const bool,
    xi_ren: f64,
    xi_fac: f64,
    results: *mut f64,
) {
    let grid = unsafe { &*grid };
    let mut pdf1 = |id, x, q2| xfx1(id, x, q2, state);
    let mut pdf2 = |id, x, q2| xfx2(id, x, q2, state);
    let mut als = |q2| alphas(q2, state);
    let order_mask = if order_mask.is_null() {
        vec![]
    } else {
        unsafe { slice::from_raw_parts(order_mask, grid.orders().len()) }.to_vec()
    };
    let channel_mask = if channel_mask.is_null() {
        vec![]
    } else {
        unsafe { slice::from_raw_parts(channel_mask, grid.channels().len()) }.to_vec()
    };
    let results = unsafe { slice::from_raw_parts_mut(results, grid.bin_info().bins()) };
    let mut lumi_cache = LumiCache::with_two(pdg_id1, &mut pdf1, pdg_id2, &mut pdf2, &mut als);

    results.copy_from_slice(&grid.convolve(
        &mut lumi_cache,
        &order_mask,
        &[],
        &channel_mask,
        &[(xi_ren, xi_fac, 1.0)],
    ));
}

/// Try to deduplicate channels of `grid` by detecting pairs of them that contain the same
/// subgrids. The numerical equality is tested using a tolerance of `ulps`, given in [units of
/// least precision](https://docs.rs/float-cmp/latest/float_cmp/index.html#some-explanation).
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_dedup_channels(grid: *mut Grid, ulps: i64) {
    let grid = unsafe { &mut *grid };

    grid.dedup_channels(ulps);
}

/// Delete a grid previously created with `pineappl_grid_new`.
#[no_mangle]
#[allow(unused_variables)]
pub extern "C" fn pineappl_grid_delete(grid: Option<Box<Grid>>) {}

/// Fill `grid` for the given momentum fractions `x1` and `x2`, at the scale `q2` for the given
/// value of the `order`, `observable`, and `lumi` with `weight`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_fill(
    grid: *mut Grid,
    x1: f64,
    x2: f64,
    q2: f64,
    order: usize,
    observable: f64,
    lumi: usize,
    weight: f64,
) {
    let grid = unsafe { &mut *grid };

    grid.fill(order, observable, lumi, &[q2, x1, x2], weight);
}

/// Fill `grid` for the given momentum fractions `x1` and `x2`, at the scale `q2` for the given
/// value of the `order` and `observable` with `weights`. The parameter of weight must contain a
/// result for entry of the luminosity function the grid was created with.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_fill_all(
    grid: *mut Grid,
    x1: f64,
    x2: f64,
    q2: f64,
    order: usize,
    observable: f64,
    weights: *const f64,
) {
    let grid = unsafe { &mut *grid };
    let weights = unsafe { slice::from_raw_parts(weights, grid.channels().len()) };

    for (channel, &weight) in weights.iter().enumerate() {
        grid.fill(order, observable, channel, &[q2, x1, x2], weight);
    }
}

/// Fill `grid` with as many points as indicated by `size`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. Additionally, all remaining pointer parameters must be
/// arrays as long as specified by `size`.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_fill_array(
    grid: *mut Grid,
    x1: *const f64,
    x2: *const f64,
    q2: *const f64,
    orders: *const usize,
    observables: *const f64,
    lumis: *const usize,
    weights: *const f64,
    size: usize,
) {
    let grid = unsafe { &mut *grid };
    let x1 = unsafe { slice::from_raw_parts(x1, size) };
    let x2 = unsafe { slice::from_raw_parts(x2, size) };
    let q2 = unsafe { slice::from_raw_parts(q2, size) };
    let orders = unsafe { slice::from_raw_parts(orders, size) };
    let observables = unsafe { slice::from_raw_parts(observables, size) };
    let lumis = unsafe { slice::from_raw_parts(lumis, size) };
    let weights = unsafe { slice::from_raw_parts(weights, size) };

    for (&x1, &x2, &q2, &order, &observable, &lumi, &weight) in
        izip!(x1, x2, q2, orders, observables, lumis, weights)
    {
        grid.fill(order, observable, lumi, &[q2, x1, x2], weight);
    }
}

/// Return the luminosity function of `grid`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_lumi(grid: *const Grid) -> Box<Lumi> {
    let grid = unsafe { &*grid };

    Box::new(Lumi(grid.channels().to_vec()))
}

/// Write the order parameters of `grid` into `order_params`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. The pointer `order_params` must point to an array as large
/// as four times the number of orders in `grid`.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_order_params(grid: *const Grid, order_params: *mut u32) {
    let grid = unsafe { &*grid };
    let orders = grid.orders();
    let order_params = unsafe { slice::from_raw_parts_mut(order_params, 4 * orders.len()) };

    for (i, order) in orders.iter().enumerate() {
        order_params[4 * i] = order.alphas;
        order_params[4 * i + 1] = order.alpha;
        order_params[4 * i + 2] = order.logxir;
        order_params[4 * i + 3] = order.logxif;
    }
}

/// Return the number of orders in `grid`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
#[must_use]
pub unsafe extern "C" fn pineappl_grid_order_count(grid: *const Grid) -> usize {
    let grid = unsafe { &*grid };

    grid.orders().len()
}

// TODO: add a function that supports fragmentation scale logs

/// Creates a new and empty grid. The creation requires four different sets of parameters:
/// - The luminosity function `lumi`: A pointer to the luminosity function that specifies how the
///   cross section should be reconstructed.
/// - Order specification `orders` and `order_params`. Each `PineAPPL` grid contains a number of
///   different perturbative orders, specified by `orders`. The array `order_params` stores the
///   exponent of each perturbative order and must contain 4 integers denoting the exponent of the
///   string coupling, of the electromagnetic coupling, of the logarithm of the renormalization
///   scale, and finally of the logarithm of the factorization scale.
/// - The observable definition `bins` and `bin_limits`. Each `PineAPPL` grid can store observables
///   from a one-dimensional distribution. To this end `bins` specifies how many observables are
///   stored and `bin_limits` must contain `bins + 1` entries denoting the left and right limit for
///   each bin.
/// - More (optional) information can be given in a key-value storage `key_vals`, which might be
///   a null pointer, to signal there are no further parameters that need to be set.
///
/// # Safety
///
/// The parameter `lumi` must point a valid luminosity function created by `pineappl_lumi_new`.
/// `order_params` must be an array with a length of `4 * orders`, and `bin_limits` an array with
/// length `bins + 1`. `key_vals` must be a valid `KeyVal` object created by `pineappl_keyval_new`.
///
/// # Panics
///
/// TODO
#[no_mangle]
#[must_use]
pub unsafe extern "C" fn pineappl_grid_new(
    lumi: *const Lumi,
    orders: usize,
    order_params: *const u32,
    bins: usize,
    bin_limits: *const f64,
    key_vals: *const KeyVal,
) -> Box<Grid> {
    let order_params = unsafe { slice::from_raw_parts(order_params, 4 * orders) };
    let orders: Vec<_> = order_params
        .chunks(4)
        .map(|s| Order {
            alphas: s[0],
            alpha: s[1],
            logxir: s[2],
            logxif: s[3],
            // this function doesn't support fragmentation scale logs
            logxia: 0,
        })
        .collect();

    let key_vals = unsafe { key_vals.as_ref() };
    let interps = grid_interpolation_params(key_vals);

    let lumi = unsafe { &*lumi };

    let mut convolutions = vec![Convolution::UnpolPDF(2212); 2];

    if let Some(keyval) = key_vals {
        if let Some(value) = keyval.strings.get("initial_state_1") {
            convolutions[0] = Convolution::UnpolPDF(value.to_string_lossy().parse().unwrap());
        }

        if let Some(value) = keyval.strings.get("initial_state_2") {
            convolutions[1] = Convolution::UnpolPDF(value.to_string_lossy().parse().unwrap());
        }
    }

    Box::new(Grid::new(
        lumi.0.clone(),
        orders,
        unsafe { slice::from_raw_parts(bin_limits, bins + 1) }.to_vec(),
        convolutions,
        interps,
        vec![Kinematics::MU2_RF, Kinematics::X1, Kinematics::X2],
    ))
}

/// Read a `PineAPPL` grid from a file with name `filename`.
///
/// # Safety
///
/// The parameter `filename` must be a C string pointing to an existing grid file.
///
/// # Panics
///
/// TODO
#[no_mangle]
#[must_use]
pub unsafe extern "C" fn pineappl_grid_read(filename: *const c_char) -> Box<Grid> {
    let filename = unsafe { CStr::from_ptr(filename) };
    let filename = filename.to_string_lossy();
    let reader = File::open(filename.as_ref()).unwrap();

    Box::new(Grid::read(reader).unwrap())
}

/// Merges the bins of corresponding to the indices from the half-open interval `[from, to]` into a
/// single bin.
///
/// # Safety
///
/// The parameter `grid` must be valid `Grid` object created by either `pineappl_grid_new` or
/// `pineappl_grid_read`.
///
/// # Panics
///
/// TODO
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_merge_bins(grid: *mut Grid, from: usize, to: usize) {
    let grid = unsafe { &mut *grid };

    grid.merge_bins(from..to).unwrap();
}

/// Merges `other` into `grid` and subsequently deletes `other`.
///
/// # Safety
///
/// Both `grid` and `other` must be valid `Grid` objects created by either `pineappl_grid_new` or
/// `pineappl_grid_read`. If `other` is a `NULL` pointer, this function does not do anything.
///
/// # Panics
///
/// TODO
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_merge_and_delete(grid: *mut Grid, other: Option<Box<Grid>>) {
    if let Some(other) = other {
        let grid = unsafe { &mut *grid };

        grid.merge(*other).unwrap();
    }
}

/// Scale all grids in `grid` by `factor`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_scale(grid: *mut Grid, factor: f64) {
    let grid = unsafe { &mut *grid };

    grid.scale(factor);
}

/// Splits the grid such that the luminosity function contains only a single combination per
/// channel.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_split_lumi(grid: *mut Grid) {
    let grid = unsafe { &mut *grid };

    grid.split_channels();
}

/// Optimizes the grid representation for space efficiency.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_optimize(grid: *mut Grid) {
    let grid = unsafe { &mut *grid };

    grid.optimize();
}

/// Optimizes the grid representation for space efficiency. The parameter `flags` determines which
/// optimizations are applied, see [`GridOptFlags`].
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_optimize_using(grid: *mut Grid, flags: GridOptFlags) {
    let grid = unsafe { &mut *grid };

    grid.optimize_using(flags);
}

/// Scales each subgrid by a bin-dependent factor given in `factors`. If a bin does not have a
/// corresponding entry in `factors` it is not rescaled. If `factors` has more entries than there
/// are bins the superfluous entries do not have an effect.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. The pointer `factors` must be an array of at least the size
/// given by `count`.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_scale_by_bin(
    grid: *mut Grid,
    count: usize,
    factors: *const f64,
) {
    let grid = unsafe { &mut *grid };
    let factors = unsafe { slice::from_raw_parts(factors, count) };

    grid.scale_by_bin(factors);
}

/// Scales each subgrid by a factor which is the product of the given values `alphas`, `alpha`,
/// `logxir`, and `logxif`, each raised to the corresponding powers for each subgrid. In addition,
/// every subgrid is scaled by a factor `global` independently of its order.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_scale_by_order(
    grid: *mut Grid,
    alphas: f64,
    alpha: f64,
    logxir: f64,
    logxif: f64,
    global: f64,
) {
    let grid = unsafe { &mut *grid };

    grid.scale_by_order(alphas, alpha, logxir, logxif, global);
}

/// Return the value for `key` stored in `grid`. If `key` isn't found, `NULL` will be returned.
/// After usage the string must be deallocated using [`pineappl_string_delete`].
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the `NULL`
/// pointer, this function is not safe to call. The parameter `key` must be non-`NULL` and a valid
/// C string.
///
/// # Panics
///
/// TODO
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_key_value(
    grid: *const Grid,
    key: *const c_char,
) -> *mut c_char {
    let grid = unsafe { &*grid };
    let key = unsafe { CStr::from_ptr(key) };
    let key = key.to_string_lossy();

    // backwards compatibility
    let index = match key.as_ref() {
        "initial_state_1" => Some(0),
        "initial_state_2" => Some(1),
        _ => None,
    };

    if let Some(index) = index {
        return CString::new(grid.convolutions()[index].pid().unwrap().to_string())
            .unwrap()
            .into_raw();
    }

    CString::new(grid.metadata().get(key.as_ref()).map_or("", String::as_str))
        .unwrap()
        .into_raw()
}

/// Sets an internal key-value pair for the grid.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. The parameters `key` and `value` must be non-`NULL` and
/// valid C strings.
///
/// # Panics
///
/// TODO
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_set_key_value(
    grid: *mut Grid,
    key: *const c_char,
    value: *const c_char,
) {
    let grid = unsafe { &mut *grid };
    let key = unsafe { CStr::from_ptr(key) }
        .to_string_lossy()
        .into_owned();
    let value = unsafe { CStr::from_ptr(value) }
        .to_string_lossy()
        .into_owned();

    // backwards compatibility
    let index = match key.as_str() {
        "initial_state_1" => Some(0),
        "initial_state_2" => Some(1),
        _ => None,
    };

    if let Some(index) = index {
        grid.convolutions_mut()[index] = Convolution::UnpolPDF(value.parse().unwrap());
    }

    grid.metadata_mut().insert(key, value);
}

/// Sets a remapper for the grid. This can be used to 'upgrade' one-dimensional bin limits to
/// N-dimensional ones. The new bin limits must be given in the form of tuples giving the left and
/// right limits, and a tuple for each dimension.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. The arrays `normalizations` and `limits` must be at least as
/// long as the number of bins of the grid and `2 * dimensions * bins`, respectively.
///
/// # Panics
///
/// TODO
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_set_remapper(
    grid: *mut Grid,
    dimensions: usize,
    normalizations: *const f64,
    limits: *const f64,
) {
    let grid = unsafe { &mut *grid };
    let bins = grid.bin_info().bins();
    let normalizations = unsafe { slice::from_raw_parts(normalizations, bins) };
    let limits = unsafe { slice::from_raw_parts(limits, 2 * dimensions * bins) };

    grid.set_remapper(
        BinRemapper::new(
            normalizations.to_vec(),
            limits
                .chunks_exact(2)
                .map(|chunk| (chunk[0], chunk[1]))
                .collect(),
        )
        .unwrap(),
    )
    .unwrap();
}

/// Write `grid` to a file with name `filename`. If `filename` ends in `.lz4` the grid is
/// automatically LZ4 compressed.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. The parameter `filename` must be a non-`NULL`, non-empty,
/// and valid C string pointing to a non-existing, but writable file.
///
/// # Panics
///
/// TODO
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_write(grid: *const Grid, filename: *const c_char) {
    let grid = unsafe { &*grid };
    let filename = unsafe { CStr::from_ptr(filename) };
    let filename = filename.to_string_lossy();
    let path = Path::new(filename.as_ref());
    let writer = File::create(path).unwrap();

    if path.extension().map_or(false, |ext| ext == "lz4") {
        grid.write_lz4(writer).unwrap();
    } else {
        grid.write(writer).unwrap();
    }
}

/// Adds a linear combination of initial states to the luminosity function `lumi`.
///
/// # Safety
///
/// The parameter `lumi` must point to a valid `Lumi` object created by `pineappl_lumi_new`.
/// `pdg_id_pairs` must be an array with length `2 * combinations`, and `factors` with length of
/// `combinations`.
#[no_mangle]
pub unsafe extern "C" fn pineappl_lumi_add(
    lumi: *mut Lumi,
    combinations: usize,
    pdg_id_pairs: *const i32,
    factors: *const f64,
) {
    let lumi = unsafe { &mut *lumi };
    let pdg_id_pairs = unsafe { slice::from_raw_parts(pdg_id_pairs, 2 * combinations) };
    let factors = if factors.is_null() {
        vec![1.0; combinations]
    } else {
        unsafe { slice::from_raw_parts(factors, combinations) }.to_vec()
    };

    lumi.0.push(Channel::new(
        pdg_id_pairs
            .chunks(2)
            .zip(factors)
            .map(|x| (vec![x.0[0], x.0[1]], x.1))
            .collect(),
    ));
}

/// Returns the number of combinations of the luminosity function `lumi` for the specified entry.
///
/// # Safety
///
/// The parameter `lumi` must point to a valid `Lumi` object created by `pineappl_lumi_new` or
/// `pineappl_grid_lumi`.
#[no_mangle]
pub unsafe extern "C" fn pineappl_lumi_combinations(lumi: *const Lumi, entry: usize) -> usize {
    let lumi = unsafe { &*lumi };

    lumi.0[entry].entry().len()
}

/// Returns the number of channel for the luminosity function `lumi`.
///
/// # Safety
///
/// The parameter `lumi` must point to a valid `Lumi` object created by `pineappl_lumi_new` or
/// `pineappl_grid_lumi`.
#[no_mangle]
pub unsafe extern "C" fn pineappl_lumi_count(lumi: *const Lumi) -> usize {
    let lumi = unsafe { &*lumi };

    lumi.0.len()
}

/// Delete luminosity function previously created with `pineappl_lumi_new`.
#[no_mangle]
#[allow(unused_variables)]
pub extern "C" fn pineappl_lumi_delete(lumi: Option<Box<Lumi>>) {}

/// Read out the channel with index `entry` of the luminosity function `lumi`. The PDG ids and
/// factors will be copied into `pdg_ids` and `factors`.
///
/// # Safety
///
/// The parameter `lumi` must point to a valid `Lumi` object created by `pineappl_lumi_new` or
/// `pineappl_grid_lumi`. The parameter `factors` must point to an array as long as the size
/// returned by `pineappl_lumi_combinations` and `pdg_ids` must point to an array that is twice as
/// long.
#[no_mangle]
pub unsafe extern "C" fn pineappl_lumi_entry(
    lumi: *const Lumi,
    entry: usize,
    pdg_ids: *mut i32,
    factors: *mut f64,
) {
    let lumi = unsafe { &*lumi };
    let entry = lumi.0[entry].entry();
    let pdg_ids = unsafe { slice::from_raw_parts_mut(pdg_ids, 2 * entry.len()) };
    let factors = unsafe { slice::from_raw_parts_mut(factors, entry.len()) };

    entry
        .iter()
        .flat_map(|(pids, _)| pids)
        .zip(pdg_ids.iter_mut())
        .for_each(|(from, to)| *to = *from);
    entry
        .iter()
        .map(|(_, factor)| factor)
        .zip(factors.iter_mut())
        .for_each(|(from, to)| *to = *from);
}

/// Creates a new luminosity function and returns a pointer to it. If no longer needed, the object
/// should be deleted using `pineappl_lumi_delete`.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_lumi_new() -> Box<Lumi> {
    Box::default()
}

/// Key-value storage for passing optional information during grid creation with
/// `pineappl_grid_new`.
#[derive(Default)]
pub struct KeyVal {
    bools: HashMap<String, bool>,
    doubles: HashMap<String, f64>,
    ints: HashMap<String, i32>,
    strings: HashMap<String, CString>,
}

/// Delete the previously created object pointed to by `key_vals`.
#[no_mangle]
#[allow(unused_variables)]
pub extern "C" fn pineappl_keyval_delete(key_vals: Option<Box<KeyVal>>) {}

/// Get the boolean-valued parameter with name `key` stored in `key_vals`.
///
/// # Safety
///
/// The parameter `key_vals` must point to a valid `KeyVal` object created by
/// `pineappl_keyval_new`. `key` must be a valid C string.
#[no_mangle]
#[must_use]
pub unsafe extern "C" fn pineappl_keyval_bool(key_vals: *const KeyVal, key: *const c_char) -> bool {
    let key_vals = unsafe { &*key_vals };
    let key = unsafe { CStr::from_ptr(key) };

    key_vals.bools[key.to_string_lossy().as_ref()]
}

/// Get the double-valued parameter with name `key` stored in `key_vals`.
///
/// # Safety
///
/// The parameter `key_vals` must point to a valid `KeyVal` object created by
/// `pineappl_keyval_new`. `key` must be a valid C string.
#[no_mangle]
#[must_use]
pub unsafe extern "C" fn pineappl_keyval_double(
    key_vals: *const KeyVal,
    key: *const c_char,
) -> f64 {
    let key_vals = unsafe { &*key_vals };
    let key = unsafe { CStr::from_ptr(key) };

    key_vals.doubles[key.to_string_lossy().as_ref()]
}

/// Get the string-valued parameter with name `key` stored in `key_vals`.
///
/// # Safety
///
/// The parameter `key_vals` must point to a valid `KeyVal` object created by
/// `pineappl_keyval_new`. `key` must be a valid C string.
#[no_mangle]
#[must_use]
pub unsafe extern "C" fn pineappl_keyval_int(key_vals: *const KeyVal, key: *const c_char) -> i32 {
    let key_vals = unsafe { &*key_vals };
    let key = unsafe { CStr::from_ptr(key) };

    key_vals.ints[key.to_string_lossy().as_ref()]
}

/// Get the int-valued parameter with name `key` stored in `key_vals`.
///
/// # Safety
///
/// The parameter `key_vals` must point to a valid `KeyVal` object created by
/// `pineappl_keyval_new`. `key` must be a valid C string.
#[no_mangle]
#[must_use]
pub unsafe extern "C" fn pineappl_keyval_string(
    key_vals: *const KeyVal,
    key: *const c_char,
) -> *const c_char {
    let key_vals = unsafe { &*key_vals };
    let key = unsafe { CStr::from_ptr(key) };

    key_vals.strings[key.to_string_lossy().as_ref()].as_ptr()
}

/// Return a pointer to newly-created `pineappl_keyval` object.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_keyval_new() -> Box<KeyVal> {
    Box::default()
}

/// Set the double-valued parameter with name `key` to `value` in `key_vals`.
///
/// # Safety
///
/// The parameter `key_vals` must point to a valid `KeyVal` object created by
/// `pineappl_keyval_new`. `key` must be a valid C string.
#[no_mangle]
pub unsafe extern "C" fn pineappl_keyval_set_bool(
    key_vals: *mut KeyVal,
    key: *const c_char,
    value: bool,
) {
    let key_vals = unsafe { &mut *key_vals };
    let key = unsafe { CStr::from_ptr(key) };

    key_vals
        .bools
        .insert(key.to_string_lossy().into_owned(), value);
}

/// Set the double-valued parameter with name `key` to `value` in `key_vals`.
///
/// # Safety
///
/// The parameter `key_vals` must point to a valid `KeyVal` object created by
/// `pineappl_keyval_new`. `key` must be a valid C string.
#[no_mangle]
pub unsafe extern "C" fn pineappl_keyval_set_double(
    key_vals: *mut KeyVal,
    key: *const c_char,
    value: f64,
) {
    let key_vals = unsafe { &mut *key_vals };
    let key = unsafe { CStr::from_ptr(key) };

    key_vals
        .doubles
        .insert(key.to_string_lossy().into_owned(), value);
}

/// Set the int-valued parameter with name `key` to `value` in `key_vals`.
///
/// # Safety
///
/// The parameter `key_vals` must point to a valid `KeyVal` object created by
/// `pineappl_keyval_new`. `key` must be a valid C string.
#[no_mangle]
pub unsafe extern "C" fn pineappl_keyval_set_int(
    key_vals: *mut KeyVal,
    key: *const c_char,
    value: i32,
) {
    let key_vals = unsafe { &mut *key_vals };
    let key = unsafe { CStr::from_ptr(key) };

    key_vals
        .ints
        .insert(key.to_string_lossy().into_owned(), value);
}

/// Set the string-valued parameter with name `key` to `value` in `key_vals`.
///
/// # Safety
///
/// The parameter `key_vals` must point to a valid `KeyVal` object created by
/// `pineappl_keyval_new`. `key` must be a valid C string.
#[no_mangle]
pub unsafe extern "C" fn pineappl_keyval_set_string(
    key_vals: *mut KeyVal,
    key: *const c_char,
    value: *const c_char,
) {
    let key_vals = unsafe { &mut *key_vals };
    let key = unsafe { CStr::from_ptr(key) };
    let value = unsafe { CStr::from_ptr(value) };

    key_vals
        .strings
        .insert(key.to_string_lossy().into_owned(), CString::from(value));
}

/// Deletes a string previously allocated by [`pineappl_grid_key_value`]. If `string` is a `NULL`
/// pointer this function does nothing.
///
/// # Safety
///
/// The parameter `string` must be a pointer to string created by [`pineappl_grid_key_value`] or
/// `NULL`, otherwise this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_string_delete(string: *mut c_char) {
    if !string.is_null() {
        mem::drop(unsafe { CString::from_raw(string) });
    }
}
