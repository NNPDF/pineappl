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
use ndarray::{Array4, CowArray, Ix4};
use pineappl::boc::{Bin, BinsWithFillLimits, Channel, Kinematics, Order, ScaleFuncForm, Scales};
use pineappl::convolutions::{Conv, ConvType, ConvolutionCache};
use pineappl::evolution::{AlphasTable, OperatorSliceInfo};
use pineappl::fk_table::FkTable;
use pineappl::grid::{Grid, GridOptFlags};
use pineappl::interpolation::{Interp, InterpMeth, Map, ReweightMeth};
use pineappl::packed_array::ravel_multi_index;
use pineappl::pids::PidBasis;
use pineappl::subgrid::Subgrid;
use std::collections::HashMap;
use std::collections::HashSet;
use std::ffi::{CStr, CString};
use std::fs::File;
use std::mem;
use std::os::raw::{c_char, c_int, c_void};
use std::path::Path;
use std::slice;

/// TODO
pub const PINEAPPL_GOF_OPTIMIZE_SUBGRID_TYPE: GridOptFlags = GridOptFlags::OPTIMIZE_SUBGRID_TYPE;

/// TODO
pub const PINEAPPL_GOF_OPTIMIZE_NODES: GridOptFlags = GridOptFlags::OPTIMIZE_NODES;

/// TODO
pub const PINEAPPL_GOF_SYMMETRIZE_CHANNELS: GridOptFlags = GridOptFlags::SYMMETRIZE_CHANNELS;

/// TODO
pub const PINEAPPL_GOF_STRIP_EMPTY_ORDERS: GridOptFlags = GridOptFlags::STRIP_EMPTY_ORDERS;

/// TODO
pub const PINEAPPL_GOF_MERGE_SAME_CHANNELS: GridOptFlags = GridOptFlags::MERGE_SAME_CHANNELS;

/// TODO
pub const PINEAPPL_GOF_STRIP_EMPTY_CHANNELS: GridOptFlags = GridOptFlags::STRIP_EMPTY_CHANNELS;

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

    grid.bwfl().len()
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

    grid.bwfl().dimensions()
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
    let bins = grid.bwfl().len();
    let bin_sizes = unsafe { slice::from_raw_parts_mut(bin_sizes, bins) };

    for (i, size) in grid.bwfl().normalizations().into_iter().enumerate() {
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
    let bins = grid.bwfl().len();
    let result = unsafe { slice::from_raw_parts_mut(left, bins) };

    for (lhs, rhs) in result.iter_mut().zip(
        grid.bwfl()
            .bins()
            .iter()
            .map(|bin| bin.limits()[dimension].0),
    ) {
        *lhs = rhs;
    }
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
    let bins = grid.bwfl().len();
    let result = unsafe { slice::from_raw_parts_mut(right, bins) };

    for (lhs, rhs) in result.iter_mut().zip(
        grid.bwfl()
            .bins()
            .iter()
            .map(|bin| bin.limits()[dimension].1),
    ) {
        *lhs = rhs;
    }
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
/// this function is not safe to call. The function pointers `xfx` and `alphas` must not
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
    let mut xfx = |id, x, q2| xfx(id, x, q2, state);
    let mut als = |q2| alphas(q2, state);
    let order_mask = if order_mask.is_null() {
        &[]
    } else {
        unsafe { slice::from_raw_parts(order_mask, grid.orders().len()) }
    };
    let channel_mask = if channel_mask.is_null() {
        &[]
    } else {
        unsafe { slice::from_raw_parts(channel_mask, grid.channels().len()) }
    };
    let bins = grid.bwfl().len();
    let results = unsafe { slice::from_raw_parts_mut(results, bins) };
    let mut convolution_cache = ConvolutionCache::new(
        vec![Conv::new(ConvType::UnpolPDF, pdg_id)],
        vec![&mut xfx],
        &mut als,
    );

    results.copy_from_slice(&grid.convolve(
        &mut convolution_cache,
        order_mask,
        &[],
        channel_mask,
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
    let mut xfx1 = |id, x, q2| xfx1(id, x, q2, state);
    let mut xfx2 = |id, x, q2| xfx2(id, x, q2, state);
    let mut als = |q2| alphas(q2, state);
    let order_mask = if order_mask.is_null() {
        &[]
    } else {
        unsafe { slice::from_raw_parts(order_mask, grid.orders().len()) }
    };
    let channel_mask = if channel_mask.is_null() {
        &[]
    } else {
        unsafe { slice::from_raw_parts(channel_mask, grid.channels().len()) }
    };
    let bins = grid.bwfl().len();
    let results = unsafe { slice::from_raw_parts_mut(results, bins) };
    let mut convolution_cache = ConvolutionCache::new(
        vec![
            Conv::new(ConvType::UnpolPDF, pdg_id1),
            Conv::new(ConvType::UnpolPDF, pdg_id2),
        ],
        vec![&mut xfx1, &mut xfx2],
        &mut als,
    );

    results.copy_from_slice(&grid.convolve(
        &mut convolution_cache,
        order_mask,
        &[],
        channel_mask,
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
        order_params[4 * i] = order.alphas.into();
        order_params[4 * i + 1] = order.alpha.into();
        order_params[4 * i + 2] = order.logxir.into();
        order_params[4 * i + 3] = order.logxif.into();
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
            // UNWRAP: there shouldn't be orders with exponents larger than 255
            alphas: s[0].try_into().unwrap(),
            alpha: s[1].try_into().unwrap(),
            logxir: s[2].try_into().unwrap(),
            logxif: s[3].try_into().unwrap(),
            // this function doesn't support fragmentation scale logs
            logxia: 0,
        })
        .collect();

    let key_vals = unsafe { key_vals.as_ref() };
    let interps = grid_interpolation_params(key_vals);

    let lumi = unsafe { &*lumi };

    let mut convolutions = vec![Conv::new(ConvType::UnpolPDF, 2212); 2];

    if let Some(keyval) = key_vals {
        if let Some(value) = keyval.strings.get("initial_state_1") {
            convolutions[0] =
                Conv::new(ConvType::UnpolPDF, value.to_string_lossy().parse().unwrap());
        }

        if let Some(value) = keyval.strings.get("initial_state_2") {
            convolutions[1] =
                Conv::new(ConvType::UnpolPDF, value.to_string_lossy().parse().unwrap());
        }
    }

    let bins = BinsWithFillLimits::from_fill_limits(
        unsafe { slice::from_raw_parts(bin_limits, bins + 1) }.to_vec(),
    )
    .unwrap();

    Box::new(Grid::new(
        bins,
        orders,
        lumi.0.clone(),
        PidBasis::Pdg,
        convolutions,
        interps,
        vec![Kinematics::Scale(0), Kinematics::X(0), Kinematics::X(1)],
        Scales {
            ren: ScaleFuncForm::Scale(0),
            fac: ScaleFuncForm::Scale(0),
            frg: ScaleFuncForm::NoScale,
        },
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

/// Change the particle ID basis of a given Grid.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object (for example when `grid` is the null pointer)
/// or if `pid_basis` does not refer to a correct basis, then this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_rotate_pid_basis(grid: *mut Grid, pid_basis: PidBasis) {
    let grid = unsafe { &mut *grid };

    grid.rotate_pid_basis(pid_basis);
}

/// Get the particle ID basis of a Grid.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the `NULL`
/// pointer, this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_pid_basis(grid: *mut Grid) -> PidBasis {
    let grid = unsafe { &mut *grid };

    *grid.pid_basis()
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

    grid.scale_by_order(alphas, alpha, logxir, logxif, 1.0, global);
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
        return CString::new(grid.convolutions()[index].pid().to_string())
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
        grid.convolutions_mut()[index] = Conv::new(ConvType::UnpolPDF, value.parse().unwrap());
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
    let bins = grid.bwfl().len();
    let normalizations = unsafe { slice::from_raw_parts(normalizations, bins) };
    let limits = unsafe { slice::from_raw_parts(limits, 2 * dimensions * bins) };

    let new_bins: Vec<_> = limits
        .chunks_exact(2 * dimensions)
        .zip(normalizations)
        .map(|(limits, &normalization)| {
            Bin::new(
                limits
                    .chunks_exact(2)
                    .map(|limits| (limits[0], limits[1]))
                    .collect(),
                normalization,
            )
        })
        .collect();

    grid.set_bwfl(
        BinsWithFillLimits::new(new_bins, grid.bwfl().fill_limits().to_vec())
            // UNWRAP: error handling in the CAPI is to abort
            .unwrap(),
    )
    // UNWRAP: error handling in the CAPI is to abort
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

    if path.extension().is_some_and(|ext| ext == "lz4") {
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

// Here starts the generalized C-API interface.

/// Type for defining a Channel function.
#[derive(Default)]
pub struct Channels(Vec<Channel>);

/// Type for defining the interpolation object
#[repr(C)]
pub struct InterpTuples {
    node_min: f64,
    node_max: f64,
    nb_nodes: usize,
    interp_degree: usize,
    reweighting_method: ReweightMeth,
    mapping: Map,
    interpolation_method: InterpMeth,
}

#[must_use]
fn construct_interpolation(interp: &InterpTuples) -> Interp {
    Interp::new(
        interp.node_min,
        interp.node_max,
        interp.nb_nodes,
        interp.interp_degree,
        interp.reweighting_method,
        interp.mapping,
        interp.interpolation_method,
    )
}

/// Type for defining the Operator info.
#[repr(C)]
#[derive(Debug)]
pub struct OperatorInfo {
    fac0: f64,
    fac1: f64,
    pid_basis: PidBasis,
    conv_type: ConvType,
}

/// An exact duplicate of `pineappl_lumi_new` to make naming (lumi -> channel) consistent.
/// should be deleted using `pineappl_channels_delete`.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_channels_new() -> Box<Channels> {
    Box::default()
}

/// Adds a generalized linear combination of initial states to the Luminosity.
///
/// # Safety
///
/// The parameter `channels` must point to a valid `Channels` object created by `pineappl_channels_new`.
/// `pdg_id_combinations` must be an array with length `nb_combinations * combinations`, and
/// `factors` with length of `combinations`. The `nb_convolutions` describe the number of
/// parton distributions involved, while `combinations` represent the number of different
/// channel combinations.
#[no_mangle]
pub unsafe extern "C" fn pineappl_channels_add(
    channels: *mut Channels,
    combinations: usize,
    nb_convolutions: usize,
    pdg_id_combinations: *const i32,
    factors: *const f64,
) {
    let channels = unsafe { &mut *channels };
    let pdg_id_pairs =
        unsafe { slice::from_raw_parts(pdg_id_combinations, nb_convolutions * combinations) };
    let factors = if factors.is_null() {
        vec![1.0; combinations]
    } else {
        unsafe { slice::from_raw_parts(factors, combinations) }.to_vec()
    };

    channels.0.push(Channel::new(
        pdg_id_pairs
            .chunks(nb_convolutions)
            .zip(factors)
            .map(|x| ((0..nb_convolutions).map(|i| x.0[i]).collect(), x.1))
            .collect(),
    ));
}

/// An exact duplicate of `pineappl_grid_lumi` to make naming (lumi -> channel) consistent.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_channels(grid: *const Grid) -> Box<Channels> {
    let grid = unsafe { &*grid };

    Box::new(Channels(grid.channels().to_vec()))
}

/// An exact duplicate of `pineappl_lumi_count` to make naming (lumi -> channel) consistent.
///
/// # Safety
///
/// The parameter `channels` must point to a valid `Lumi` object created by `pineappl_channels_new` or
/// `pineappl_grid_channels`.
#[no_mangle]
pub unsafe extern "C" fn pineappl_channels_count(channels: *const Channels) -> usize {
    let channels = unsafe { &*channels };

    channels.0.len()
}

/// An exact duplicate of `pineappl_lumi_combinations` to make naming (lumi -> channel) consistent.
///
/// # Safety
///
/// The parameter `channels` must point to a valid `Channels` object created by `pineappl_channels_new` or
/// `pineappl_grid_channels`.
#[no_mangle]
pub unsafe extern "C" fn pineappl_channels_combinations(
    channels: *const Channels,
    entry: usize,
) -> usize {
    let channels = unsafe { &*channels };

    channels.0[entry].entry().len()
}

/// An exact duplicate of `pineappl_lumi_delete` to make naming (lumi -> channel) consistent.
#[no_mangle]
#[allow(unused_variables)]
pub extern "C" fn pineappl_channels_delete(channels: Option<Box<Channels>>) {}

/// Creates a new and empty grid that can accept any number of convolutions. The creation requires
/// the following different sets of parameters:
/// - The PID basis `pid_basis`: The basis onto which the partons are mapped, can be `Evol` or `Pdg`.
/// - The channel function `channels`: A pointer to the channels function that specifies how the
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
/// - The types of convolutions `convolution_types` and their numbers `nb_convolutions`: specify how
///   how many different convolutions are involved and their types - which are a cross product of the
///   the following combination: (unpolarized, polarized) ⊗ (PDF, Fragmentation Function).
/// - The PDG IDs of the involved initial- or final-state hadrons `pdg_ids`.
/// - The types of kinematics `kinematics`: specify the various kinematics required to construct the
///   Grid. These can be the energy scales and the various momentum fractions.
/// - The specifications of the interpolation methods `interpolations`: provide the specifications on
///   how each of the kinematics should be interpolated.
/// - The unphysical renormalization, factorization, and fragmentation scales: `mu_scales`. Its entries
///   have to be ordered following {ren, fac, frg}. The mapping is as follows:
///   `0` -> `ScaleFuncForm::NoScale`, ..., `n` -> `ScaleFuncForm::Scale(n - 1)`.
///
/// # Safety
/// TODO
///
/// # Panics
/// TODO
#[no_mangle]
#[must_use]
pub unsafe extern "C" fn pineappl_grid_new2(
    pid_basis: PidBasis,
    channels: *const Channels,
    orders: usize,
    order_params: *const u8,
    bins: usize,
    bin_limits: *const f64,
    nb_convolutions: usize,
    convolution_types: *const ConvType,
    pdg_ids: *const c_int,
    kinematics: *const Kinematics,
    interpolations: *const InterpTuples,
    mu_scales: *const usize,
) -> Box<Grid> {
    // Luminosity channels
    let channels = unsafe { &*channels };

    // Perturbative orders
    let order_params = unsafe { slice::from_raw_parts(order_params, 5 * orders) };
    let orders: Vec<_> = order_params
        .chunks(5)
        .map(|s| Order {
            alphas: s[0],
            alpha: s[1],
            logxir: s[2],
            logxif: s[3],
            logxia: s[4],
        })
        .collect();

    let bins = BinsWithFillLimits::from_fill_limits(
        unsafe { slice::from_raw_parts(bin_limits, bins + 1) }.to_vec(),
    )
    .unwrap();

    // Construct the convolution objects
    let convolution_types =
        unsafe { slice::from_raw_parts(convolution_types, nb_convolutions).to_vec() };
    let pdg_ids = unsafe { slice::from_raw_parts(pdg_ids, nb_convolutions).to_vec() };
    let convolutions = izip!(convolution_types.iter(), pdg_ids.iter())
        .map(|(&conv, &pdg_value)| Conv::new(conv, pdg_value))
        .collect();

    // Grid interpolations
    let interp_slices = unsafe { std::slice::from_raw_parts(interpolations, nb_convolutions + 1) };
    let interp_vecs: Vec<Interp> = interp_slices.iter().map(construct_interpolation).collect();

    // Construct the kinematic variables
    let kinematics = unsafe { slice::from_raw_parts(kinematics, interp_vecs.len()).to_vec() };

    // Scales. An array containing the values of {ren, fac, frg}
    let mu_scales = unsafe { std::slice::from_raw_parts(mu_scales, 3) };
    let mu_scales_vec: Vec<ScaleFuncForm> = mu_scales
        .iter()
        .map(|&scale| {
            if scale == 0 {
                ScaleFuncForm::NoScale
            } else {
                ScaleFuncForm::Scale(scale - 1)
            }
        })
        .collect();

    Box::new(Grid::new(
        bins,
        orders,
        channels.0.clone(),
        pid_basis,
        convolutions,
        interp_vecs,
        kinematics,
        Scales {
            ren: mu_scales_vec[0].clone(),
            fac: mu_scales_vec[1].clone(),
            frg: mu_scales_vec[2].clone(),
        },
    ))
}

/// Similar to  `pineappl_grid_fill` but accepts any given momentum fractions {`x1`, ...,`xn`} at
/// various energy scalesfor the given value of the `order`, `observable`, and `lumi` with `weight`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_fill2(
    grid: *mut Grid,
    order: usize,
    observable: f64,
    channel: usize,
    ntuple: *const f64,
    weight: f64,
) {
    let grid = unsafe { &mut *grid };
    let ntuple = unsafe { slice::from_raw_parts(ntuple, grid.kinematics().len()) };

    grid.fill(order, observable, channel, ntuple, weight);
}

/// Similar to  `pineappl_grid_fill_all` but accepts any given momentum fractions {`x1`, ...,`xn`} at
/// various energy scalesfor the given value of the `order`, `observable`, and `lumi` with `weight`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_fill_all2(
    grid: *mut Grid,
    order: usize,
    observable: f64,
    ntuple: *const f64,
    weights: *const f64,
) {
    let grid = unsafe { &mut *grid };
    let ntuple = unsafe { slice::from_raw_parts(ntuple, grid.kinematics().len()) };
    let weights = unsafe { slice::from_raw_parts(weights, grid.channels().len()) };

    for (channel, &weight) in weights.iter().enumerate() {
        grid.fill(order, observable, channel, ntuple, weight);
    }
}

/// Similar to  `pineappl_grid_fill_array` but accepts any given momentum fractions
/// {`x1`, ...,`xn`} at various energy scalesfor the given value of the `order`, `observable`,
/// and `lumi` with `weight`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. Additionally, all remaining pointer parameters must be
/// arrays as long as specified by `size`.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_fill_array2(
    grid: *mut Grid,
    orders: *const usize,
    observables: *const f64,
    ntuples: *const f64,
    channels: *const usize,
    weights: *const f64,
    size: usize,
) {
    let grid = unsafe { &mut *grid };
    let orders = unsafe { slice::from_raw_parts(orders, size) };
    let observables = unsafe { slice::from_raw_parts(observables, size) };
    let channels = unsafe { slice::from_raw_parts(channels, size) };
    let weights = unsafe { slice::from_raw_parts(weights, size) };

    // Convert the 1D slice into a 2D array
    let ntuples = unsafe { slice::from_raw_parts(ntuples, size * grid.kinematics().len()) };
    let ntuples_2d: Vec<&[f64]> = ntuples.chunks(grid.kinematics().len()).collect();

    for (ntuple, &order, &observable, &channel, &weight) in
        izip!(ntuples_2d, orders, observables, channels, weights)
    {
        grid.fill(order, observable, channel, ntuple, weight);
    }
}

/// Similar to `pineappl_grid_split_channels`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_split_channels(grid: *mut Grid) {
    let grid = unsafe { &mut *grid };

    grid.split_channels();
}

/// Read out the channel with index `entry` of the given `channels`. The PDG ids and factors will
/// be copied into `pdg_ids` and `factors`.
///
/// # Safety
///
/// The parameter `channels` must point to a valid [`Channels`] object created by
/// [`pineappl_channels_new`] or [`pineappl_grid_channels`]. The parameter `factors` must point to
/// an array as long as the result of [`pineappl_channels_combinations`] and `pdg_ids` must
/// point to an array as long as the result multiplied with the number of convolutions.
#[no_mangle]
pub unsafe extern "C" fn pineappl_channels_entry(
    channels: *const Channels,
    entry: usize,
    pdg_ids: *mut i32,
    factors: *mut f64,
) {
    let channels = unsafe { &*channels };
    let entry = channels.0[entry].entry();
    // if the channel has no entries we assume no convolutions, which is OK we don't copy anything
    // in this case
    let convolutions = entry.first().map_or(0, |x| x.0.len());
    let pdg_ids = unsafe { slice::from_raw_parts_mut(pdg_ids, convolutions * entry.len()) };
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

/// An extension of `pineappl_grid_order_params` that accounts for the order of the fragmentation
/// logs.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. The pointer `order_params` must point to an array as large
/// as four times the number of orders in `grid`.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_order_params2(grid: *const Grid, order_params: *mut u32) {
    let grid = unsafe { &*grid };
    let orders = grid.orders();
    let order_params = unsafe { slice::from_raw_parts_mut(order_params, 5 * orders.len()) };

    for (i, order) in orders.iter().enumerate() {
        order_params[5 * i] = order.alphas.into();
        order_params[5 * i + 1] = order.alpha.into();
        order_params[5 * i + 2] = order.logxir.into();
        order_params[5 * i + 3] = order.logxif.into();
        order_params[5 * i + 4] = order.logxia.into();
    }
}

/// A generalization of the convolution function.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. The function pointers `xfx` and `alphas` must not
/// be null pointers and point to valid functions. The parameters `order_mask` and `channel_mask`
/// must either be null pointers or point to arrays that are as long as `grid` has orders and
/// channels, respectively. Finally, `results` must be as long as `grid` has bins.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_convolve(
    grid: *const Grid,
    xfx: extern "C" fn(pdg_id: i32, x: f64, q2: f64, state: *mut c_void) -> f64,
    alphas: extern "C" fn(q2: f64, state: *mut c_void) -> f64,
    pdfs_state: *mut *mut c_void,
    alphas_state: *mut c_void,
    order_mask: *const bool,
    channel_mask: *const bool,
    bin_indices: *const usize,
    nb_scales: usize,
    mu_scales: *const f64,
    results: *mut f64,
) {
    let grid = unsafe { &*grid };

    let order_mask = if order_mask.is_null() {
        &[]
    } else {
        unsafe { slice::from_raw_parts(order_mask, grid.orders().len()) }
    };

    let channel_mask = if channel_mask.is_null() {
        &[]
    } else {
        unsafe { slice::from_raw_parts(channel_mask, grid.channels().len()) }
    };

    let bin_indices = if bin_indices.is_null() {
        &[]
    } else {
        unsafe { slice::from_raw_parts(bin_indices, grid.bwfl().len()) }
    };

    // Construct the alphas and PDFs functions
    let mut als = |q2| alphas(q2, alphas_state);

    let pdfs_slices = unsafe { slice::from_raw_parts(pdfs_state, grid.convolutions().len()) };
    let mut xfx_funcs: Vec<_> = pdfs_slices
        .iter()
        .map(|&state| move |id, x, q2| xfx(id, x, q2, state))
        .collect();

    // Construct the Convolution cache
    let mut convolution_cache = ConvolutionCache::new(
        grid.convolutions().to_vec(),
        xfx_funcs
            .iter_mut()
            .map(|fx| fx as &mut dyn FnMut(i32, f64, f64) -> f64)
            .collect(),
        &mut als,
    );

    // The factorization, renormalization, and fragmentation scale factors
    let mu_scales = if mu_scales.is_null() {
        &[(1.0, 1.0, 1.0)]
    } else {
        unsafe { slice::from_raw_parts(mu_scales.cast::<(f64, f64, f64)>(), nb_scales) }
    };

    let results = unsafe { slice::from_raw_parts_mut(results, grid.bwfl().len()) };

    results.copy_from_slice(&grid.convolve(
        &mut convolution_cache,
        order_mask,
        bin_indices,
        channel_mask,
        mu_scales,
    ));
}

/// Get the type of convolutions for this Grid given an index.
///
/// # Safety
///
/// TODO
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_conv_type(grid: *const Grid, index: usize) -> ConvType {
    let grid = unsafe { &*grid };

    grid.convolutions()[index].conv_type()
}

/// Get the number of convolutions for a given Grid.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_convolutions_len(grid: *mut Grid) -> usize {
    let grid = unsafe { &mut *grid };

    grid.convolutions().len()
}

/// Get the number of different kinematics for a given Grid.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_kinematics_len(grid: *mut Grid) -> usize {
    let grid = unsafe { &mut *grid };

    grid.kinematics().len()
}

/// Get the shape of a subgrid for a given bin, channel, and order.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. Additionally, the pointer that specifies the shape of the
/// subgrid has to be an array whose size must be as given by `pineappl_grid_kinematics_len`.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_subgrid_shape(
    grid: *const Grid,
    bin: usize,
    order: usize,
    channel: usize,
    shape: *mut usize,
) {
    let grid = unsafe { &*grid };
    let subgrid = &grid.subgrids()[[order, bin, channel]];
    let subgrid_shape = if subgrid.is_empty() {
        // avoid calling `Subgrid::shape()` for empty grids, which may panic
        let subgrid_dim = grid.kinematics().len();
        &vec![0; subgrid_dim]
    } else {
        subgrid.shape()
    };
    let shape = unsafe { slice::from_raw_parts_mut(shape, grid.kinematics().len()) };

    shape.copy_from_slice(subgrid_shape);
}

/// Get the nodes of the subgrid.
///
/// # Safety
///
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_subgrid_node_values(
    grid: *const Grid,
    bin: usize,
    order: usize,
    channel: usize,
    node_values: *mut f64,
) {
    let grid = unsafe { &*grid };
    let subgrid = &grid.subgrids()[[order, bin, channel]];

    let flatten_subgrid: Vec<f64> = subgrid.node_values().into_iter().flatten().collect();
    let node_values = unsafe { slice::from_raw_parts_mut(node_values, flatten_subgrid.len()) };

    node_values.copy_from_slice(flatten_subgrid.as_slice());
}

/// Get the subgrid for a given bin, channel, and order.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. Additionally, the pointer that specifies the size of the subgrid
/// when flattened must be an array; its size must be computed by multiplying the shape dimension as
/// given by `pineappl_grid_subgrid_shape`.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_subgrid_array(
    grid: *const Grid,
    bin: usize,
    order: usize,
    channel: usize,
    subgrid_array: *mut f64,
) {
    let grid = unsafe { &*grid };
    let subgrid = &grid.subgrids()[[order, bin, channel]];

    // avoid calling `Subgrid::shape()` for empty grids, which may panic
    if !subgrid.is_empty() {
        let shape = subgrid.shape();
        let subgrid_array =
            unsafe { slice::from_raw_parts_mut(subgrid_array, shape.iter().product()) };

        for (index, value) in subgrid.indexed_iter() {
            let ravel_index = ravel_multi_index(index.as_slice(), shape);
            subgrid_array[ravel_index] = value;
        }
    }
}

/// Get the shape of the objects represented in the evolve info.
///
/// # Safety
///
/// TODO
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_evolve_info_shape(
    grid: *const Grid,
    orders: *const u8,
    shape_info: *mut usize,
) {
    let grid = unsafe { &*grid };

    let orders = unsafe { slice::from_raw_parts(orders, 2) };
    let order_mask = Order::create_mask(grid.orders(), orders[0], orders[1], true);
    let shape_info = unsafe { slice::from_raw_parts_mut(shape_info, 5) };

    let evolve_info = grid.evolve_info(&order_mask);
    let evinfo_shapes = [
        evolve_info.fac1.len(),
        evolve_info.frg1.len(),
        evolve_info.pids1.len(),
        evolve_info.x1.len(),
        evolve_info.ren1.len(),
    ];

    shape_info.copy_from_slice(&evinfo_shapes);
}

/// Get the values of the parameters contained in evolve info.
///
/// # Safety
///
/// TODO
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_evolve_info(
    grid: *const Grid,
    orders: *const u8,
    fac1: *mut f64,
    frg1: *mut f64,
    pids1: *mut i32,
    x1: *mut f64,
    ren1: *mut f64,
) {
    let grid = unsafe { &*grid };

    let orders = unsafe { slice::from_raw_parts(orders, 2) };
    let order_mask = Order::create_mask(grid.orders(), orders[0], orders[1], true);

    let evolve_info = grid.evolve_info(&order_mask);
    let fac1 = unsafe { slice::from_raw_parts_mut(fac1, evolve_info.fac1.len()) };
    let frg1 = unsafe { slice::from_raw_parts_mut(frg1, evolve_info.frg1.len()) };
    let pids1 = unsafe { slice::from_raw_parts_mut(pids1, evolve_info.pids1.len()) };
    let x1 = unsafe { slice::from_raw_parts_mut(x1, evolve_info.x1.len()) };
    let ren1 = unsafe { slice::from_raw_parts_mut(ren1, evolve_info.ren1.len()) };

    fac1.copy_from_slice(&grid.evolve_info(&order_mask).fac1);
    frg1.copy_from_slice(&grid.evolve_info(&order_mask).frg1);
    pids1.copy_from_slice(&grid.evolve_info(&order_mask).pids1);
    x1.copy_from_slice(&grid.evolve_info(&order_mask).x1);
    ren1.copy_from_slice(&grid.evolve_info(&order_mask).ren1);
}

/// Evolve a grid and dump the resulting FK table.
///
/// # Safety
///
/// TODO
///
/// # Panics
///
/// TODO
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_evolve(
    grid: *mut Grid,
    op_info: *mut OperatorInfo,
    orders: *const u8,
    operators: *mut f64,
    x_grid: *mut f64,
    x_fktable: *mut f64,
    pids_grid: *mut i32,
    pids_fktable: *mut i32,
    eko_shape: *mut usize,
    xi: *mut f64,
    ren1: *mut f64,
    alphas: *mut f64,
) -> Box<FkTable> {
    let grid = unsafe { &mut *grid };

    let orders = unsafe { slice::from_raw_parts(orders, 2) };
    let eko_shape = unsafe { slice::from_raw_parts(eko_shape, 4) };
    let pids_fktable = unsafe { slice::from_raw_parts(pids_fktable, eko_shape[0]) };
    let x_fktable = unsafe { slice::from_raw_parts(x_fktable, eko_shape[1]) };
    let pids_grid = unsafe { slice::from_raw_parts(pids_grid, eko_shape[2]) };
    let x_grid = unsafe { slice::from_raw_parts(x_grid, eko_shape[3]) };

    let order_mask = Order::create_mask(grid.orders(), orders[0], orders[1], true);
    let evolve_info = grid.evolve_info(&order_mask);

    let ren1_len = evolve_info.ren1.len();
    let ren1 = unsafe { Vec::from_raw_parts(ren1, ren1_len, ren1_len) };
    let alphas = unsafe { Vec::from_raw_parts(alphas, ren1_len, ren1_len) };
    let xi = unsafe { Vec::from_raw_parts(xi, 3, 3) };

    let conv_types: HashSet<_> = grid
        .convolutions()
        .iter()
        .map(pineappl::convolutions::Conv::conv_type)
        .collect();

    let op_info =
        unsafe { slice::from_raw_parts(op_info, conv_types.len() * evolve_info.fac1.len()) };
    let opinfo_chunk = op_info.chunks_exact(conv_types.len());

    let flattened_shapes: usize = op_info
        .iter()
        .map(|_| eko_shape.iter().product::<usize>())
        .sum();
    let operators =
        unsafe { slice::from_raw_parts(operators, conv_types.len() * flattened_shapes) };

    let mut start_idx = 0;
    let op_split: Vec<_> = op_info
        .iter()
        .map(|_| {
            let end_idx = start_idx + eko_shape.iter().product::<usize>();
            let op_range = operators[start_idx..end_idx].to_vec();
            start_idx = end_idx;
            op_range
        })
        .collect();
    let ops_chunk = op_split.chunks_exact(conv_types.len());

    let slices = opinfo_chunk
        .into_iter()
        .zip(ops_chunk)
        .map(|(op_subinfo, op_range)| {
            op_subinfo.iter().zip(op_range).map(|(opinfo, op_values)| {
                let operator_slice_info = OperatorSliceInfo {
                    pid_basis: opinfo.pid_basis,
                    fac0: opinfo.fac0,
                    pids0: pids_fktable.to_vec(),
                    x0: x_fktable.to_vec(),
                    fac1: opinfo.fac1,
                    pids1: pids_grid.to_vec(),
                    x1: x_grid.to_vec(),
                    conv_type: opinfo.conv_type,
                };

                let array = Array4::from_shape_vec(
                    Ix4(eko_shape[0], eko_shape[1], eko_shape[2], eko_shape[3]),
                    op_values.clone(),
                )
                .expect("Shape mismatch or invalid input.");

                Ok::<_, std::io::Error>((operator_slice_info, CowArray::from(array)))
            })
        })
        .collect();

    let fk_table = grid.evolve(
        slices,
        &order_mask,
        (xi[0], xi[1], xi[2]),
        &AlphasTable { ren1, alphas },
    );

    Box::new(fk_table.unwrap())
}

/// Delete an FK table.
#[no_mangle]
#[allow(unused_variables)]
pub extern "C" fn pineappl_fk_table_delete(fktable: Option<Box<FkTable>>) {}

/// Write `fktable` to a file with name `filename`. If `filename` ends in `.lz4` the Fk table is
/// automatically LZ4 compressed.
///
/// # Safety
///
/// If `fktable` does not point to a valid `FkTable` object, for example when `fktable` is a null
/// pointer, this function is not safe to call. The parameter `filename` must be a non-`NULL`,
/// non-empty, and valid C string pointing to a non-existing, but writable file.
///
/// # Panics
///
/// TODO
#[no_mangle]
pub unsafe extern "C" fn pineappl_fktable_write(fktable: *const FkTable, filename: *const c_char) {
    let fktable = unsafe { &*fktable };
    let filename = unsafe { CStr::from_ptr(filename) };
    let filename = filename.to_string_lossy();
    let path = Path::new(filename.as_ref());
    let writer = File::create(path).unwrap();

    if path.extension().is_some_and(|ext| ext == "lz4") {
        fktable.grid().write_lz4(writer).unwrap();
    } else {
        fktable.grid().write(writer).unwrap();
    }
}
