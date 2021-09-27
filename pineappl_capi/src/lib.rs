#![warn(clippy::all, clippy::cargo, clippy::nursery, clippy::pedantic)]
#![warn(missing_docs)]

//! The C-language interface for `PineAPPL`.
//!
//! The `PineAPPL` Application Programming Interface for the C language (CAPI) defines types and
//! functions that allow one to use [`PineAPPL`] without having to write Rust code, and instead
//! using
//!
//! * C or C++ (see the [C++ example] in the repository).
//! * Fortran can be used as well, using the `iso_c_binding` module to wrap CAPI's functions (see
//!   the [Fortran example] in the repository).
//!
//! Note that the CAPI only defines a subset of functionality available in the Rust module, which
//! is extended on an as-needed basis. If you happen to miss a function please open a new [issue]
//! on the github page.
//!
//! [`PineAPPL`]: https://docs.rs/pineappl
//! [C++ example]: https://github.com/N3PDF/pineappl/tree/master/examples/cpp
//! [Fortran example]: https://github.com/N3PDF/pineappl/tree/master/examples/fortran
//! [issue]: https://github.com/N3PDF/pineappl/issues
//!
//! # Using and linking the CAPI
//!
//! To use `PineAPPL`'s CAPI you have to include its header `pineappl_capi.h`. To successfully
//! build your program, you should rely on [`pkgconf`] or [`pkg-conf`], as follows:
//!
//! ```shell
//! pkg-config --cflags pineappl_capi
//! pkg-config --libs pineappl_capi
//! ```
//!
//! This will read `PineAPPL`'s `.pc` file and print the neccessary `CFLAGS` (`--cflags`) and
//! linker flags (`--libs`). This procedure is supported by many build systems, such as
//!
//! * Autotools (using the `PKG_CHECK_MODULES` macro, see the [Autotools mythbuster] page for
//!   correct usage)
//! * `CMake`, using [`FindPkgConfig`], and
//! * `meson`, using [`dependency`]; it usually suffices to write `dependency('pineappl_capi')`.
//!
//! [`pkgconf`]: https://github.com/pkgconf/pkgconf
//! [`pkg-conf`]: https://www.freedesktop.org/wiki/Software/pkg-config
//! [Autotools mythbuster]: https://autotools.io/pkgconfig/pkg_check_modules.html
//! [`FindPkgConfig`]: https://cmake.org/cmake/help/latest/module/FindPkgConfig.html
//! [`dependency`]: https://mesonbuild.com/Reference-manual.html#dependency
//!
//! # Strings and other types
//!
//! All strings used in this library, denoted with the Rust type `*const c_char` or `*mut c_char`,
//! which correspond to the C types `const char*` and `char*`, respectively, are assumed to be
//! encoded using UTF-8, which Rust uses internally. The conversions are performed using
//! [`CStr::to_string_lossy`], which you can consult for how the strings are converted.
//!
//! All other built-in types in Rust (`f64`, `usize`, `i32`, `u32`, ...) correspond to types in C
//! as defined by the [translation tables] of cbindgen. If in doubt, consult the generated header
//! `pineappl_capi.h`.
//!
//! [translation tables]: https://github.com/eqrion/cbindgen/blob/master/docs.md#std-types

use itertools::izip;
use pineappl::bin::BinRemapper;
use pineappl::empty_subgrid::EmptySubgridV1;
use pineappl::grid::{Grid, Ntuple, Order};
use pineappl::import_only_subgrid::ImportOnlySubgridV2;
use pineappl::lumi::LumiEntry;
use pineappl::sparse_array3::SparseArray3;
use pineappl::subgrid::{ExtraSubgridParams, Mu2, Subgrid, SubgridParams};
use std::collections::HashMap;
use std::convert::TryFrom;
use std::ffi::{CStr, CString};
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::os::raw::{c_char, c_void};
use std::slice;

// TODO: make sure no `panic` calls leave functions marked as `extern "C"`

unsafe fn grid_params(key_vals: *const KeyVal) -> (String, SubgridParams, ExtraSubgridParams) {
    let mut subgrid_type = "LagrangeSubgrid".to_string();
    let mut subgrid_params = SubgridParams::default();
    let mut extra = ExtraSubgridParams::default();

    if !key_vals.is_null() {
        let keyval = &*key_vals;

        if let Some(value) = keyval.ints.get("q2_bins") {
            subgrid_params.set_q2_bins(usize::try_from(*value).unwrap());
        }

        if let Some(value) = keyval.doubles.get("q2_max") {
            subgrid_params.set_q2_max(*value);
        }

        if let Some(value) = keyval.doubles.get("q2_min") {
            subgrid_params.set_q2_min(*value);
        }

        if let Some(value) = keyval.ints.get("q2_order") {
            subgrid_params.set_q2_order(usize::try_from(*value).unwrap());
        }

        if let Some(value) = keyval.bools.get("reweight") {
            subgrid_params.set_reweight(*value);
        }

        if let Some(value) = keyval.ints.get("x_bins") {
            let value = usize::try_from(*value).unwrap();
            subgrid_params.set_x_bins(value);
            extra.set_x2_bins(value);
        }

        if let Some(value) = keyval.doubles.get("x_max") {
            subgrid_params.set_x_max(*value);
            extra.set_x2_max(*value);
        }

        if let Some(value) = keyval.doubles.get("x_min") {
            subgrid_params.set_x_min(*value);
            extra.set_x2_min(*value);
        }

        if let Some(value) = keyval.ints.get("x_order") {
            let value = usize::try_from(*value).unwrap();
            subgrid_params.set_x_order(value);
            extra.set_x2_order(value);
        }

        if let Some(value) = keyval.ints.get("x1_bins") {
            subgrid_params.set_x_bins(usize::try_from(*value).unwrap());
        }

        if let Some(value) = keyval.doubles.get("x1_max") {
            subgrid_params.set_x_max(*value);
        }

        if let Some(value) = keyval.doubles.get("x1_min") {
            subgrid_params.set_x_min(*value);
        }

        if let Some(value) = keyval.ints.get("x1_order") {
            subgrid_params.set_x_order(usize::try_from(*value).unwrap());
        }

        if let Some(value) = keyval.ints.get("x2_bins") {
            extra.set_x2_bins(usize::try_from(*value).unwrap());
        }

        if let Some(value) = keyval.doubles.get("x2_max") {
            extra.set_x2_max(*value);
        }

        if let Some(value) = keyval.doubles.get("x2_min") {
            extra.set_x2_min(*value);
        }

        if let Some(value) = keyval.ints.get("x2_order") {
            extra.set_x2_order(usize::try_from(*value).unwrap());
        }

        if let Some(value) = keyval.strings.get("subgrid_type") {
            subgrid_type = value.to_str().unwrap().to_string();
        }
    }

    (subgrid_type, subgrid_params, extra)
}

/// Type for defining a luminosity function.
#[derive(Default)]
pub struct Lumi(Vec<LumiEntry>);

/// Type for reading and accessing subgrids.
pub struct SubGrid(ImportOnlySubgridV2);

/// Returns the number of bins in `grid`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
#[must_use]
pub unsafe extern "C" fn pineappl_grid_bin_count(grid: *const Grid) -> usize {
    (*grid).bin_info().bins()
}

/// Returns the number of dimensions of the bins this grid has.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_bin_dimensions(grid: *const Grid) -> usize {
    (*grid).bin_info().dimensions()
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
    let sizes = (*grid).bin_info().normalizations();
    let bin_sizes = slice::from_raw_parts_mut(bin_sizes, sizes.len());

    for (i, size) in sizes.iter().enumerate() {
        bin_sizes[i] = *size;
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
    let limits = (*grid).bin_info().left(dimension);
    slice::from_raw_parts_mut(left, limits.len()).copy_from_slice(&limits);
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
    let limits = (*grid).bin_info().right(dimension);
    slice::from_raw_parts_mut(right, limits.len()).copy_from_slice(&limits);
}

/// Convolutes the specified grid with the PDFs `xfx1` and `xfx2` and strong coupling `alphas`.
/// These functions must evaluate the PDFs for the given `x` and `q2` for the parton with the given
/// PDG id, `pdg_id`, and return the result. Note that the value must be the PDF multiplied with
/// its argument `x`. The value of the pointer `state` provided to these functions is the same one
/// given to this function. The parameter `order_mask` must be as long as there are perturbative
/// orders contained in `grid` and is used to selectively disable (`false`) or enable (`true`)
/// individual orders. If `order_mask` is set to `NULL`, all orders are active. The parameter
/// `lumi_mask` can be used similarly, but must be as long as the luminosity function `grid` was
/// created with has entries, or `NULL` to enable all luminosities. The values `xi_ren` and
/// `xi_fac` can be used to vary the renormalization and factorization from its central value,
/// which corresponds to `1.0`. After convolution of the grid with the PDFs the differential cross
/// section for each bin is written into `results`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. The function pointers `xfx1`, `xfx2`, and `alphas` must be
/// null pointers and point to valid functions. The parameters `order_mask` and `lumi_mask` must
/// either be null pointers or point to arrays that are as long as `grid` has orders and lumi
/// entries, respectively. Finally, `results` must be as long as `grid` has bins.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_convolute(
    grid: *const Grid,
    xfx1: extern "C" fn(pdg_id: i32, x: f64, q2: f64, state: *mut c_void) -> f64,
    xfx2: extern "C" fn(pdg_id: i32, x: f64, q2: f64, state: *mut c_void) -> f64,
    alphas: extern "C" fn(q2: f64, state: *mut c_void) -> f64,
    state: *mut c_void,
    order_mask: *const bool,
    lumi_mask: *const bool,
    xi_ren: f64,
    xi_fac: f64,
    results: *mut f64,
) {
    let grid = &*grid;
    let pdf1 = |id, x, q2| xfx1(id, x, q2, state);
    let pdf2 = |id, x, q2| xfx2(id, x, q2, state);
    let als = |q2| alphas(q2, state);
    let order_mask = if order_mask.is_null() {
        vec![]
    } else {
        slice::from_raw_parts(order_mask, grid.orders().len()).to_vec()
    };
    let lumi_mask = if lumi_mask.is_null() {
        vec![]
    } else {
        slice::from_raw_parts(lumi_mask, grid.lumi().len()).to_vec()
    };
    let results = slice::from_raw_parts_mut(results, grid.bin_info().bins());

    results.copy_from_slice(&grid.convolute(
        &pdf1,
        &pdf2,
        &als,
        &order_mask,
        &[],
        &lumi_mask,
        &[(xi_ren, xi_fac)],
    ));
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
    (*grid).fill(order, observable, lumi, &Ntuple { x1, x2, q2, weight });
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
    let grid = &mut *grid;

    grid.fill_all(
        order,
        observable,
        &Ntuple {
            x1,
            x2,
            q2,
            weight: (),
        },
        slice::from_raw_parts(weights, grid.lumi().len()),
    );
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
    let x1 = slice::from_raw_parts(x1, size);
    let x2 = slice::from_raw_parts(x2, size);
    let q2 = slice::from_raw_parts(q2, size);
    let orders = slice::from_raw_parts(orders, size);
    let observables = slice::from_raw_parts(observables, size);
    let lumis = slice::from_raw_parts(lumis, size);
    let weights = slice::from_raw_parts(weights, size);

    for (&x1, &x2, &q2, &order, &observable, &lumi, &weight) in
        izip!(x1, x2, q2, orders, observables, lumis, weights)
    {
        (*grid).fill(order, observable, lumi, &Ntuple { x1, x2, q2, weight });
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
    Box::new(Lumi((*grid).lumi().to_vec()))
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
    let orders = (*grid).orders();
    let order_params = slice::from_raw_parts_mut(order_params, 4 * orders.len());

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
    (*grid).orders().len()
}

/// Creates a new and empty grid. The creation requires four different sets of parameters:
/// - The luminosity function `lumi`: A pointer to the luminosity function that specifies how the
/// cross section should be reconstructed.
/// - Order specification `orders` and `order_params`. Each `PineAPPL` grid contains a number of
/// different perturbative orders, specified by `orders`. The array `order_params` stores the
/// exponent of each perturbative order and must contain 4 integers denoting the exponent of the
/// string coupling, of the electromagnetic coupling, of the logarithm of the renormalization
/// scale, and finally of the logarithm of the factorization scale.
/// - The observable definition `bins` and `bin_limits`. Each `PineAPPL` grid can store observables
/// from a one-dimensional distribution. To this end `bins` specifies how many observables are
/// stored and `bin_limits` must contain `bins + 1` entries denoting the left and right limit for
/// each bin.
/// - More (optional) information can be given in a key-value storage `key_vals`, which might be
/// a null pointer, to signal there are no further parameters that need to be set.
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
    let order_params = slice::from_raw_parts(order_params, 4 * orders);
    let orders: Vec<_> = order_params
        .chunks(4)
        .map(|s| Order {
            alphas: s[0],
            alpha: s[1],
            logxir: s[2],
            logxif: s[3],
        })
        .collect();

    let (subgrid_type, subgrid_params, extra) = grid_params(key_vals);

    let mut grid = Box::new(
        Grid::with_subgrid_type(
            (*lumi).0.clone(),
            orders,
            slice::from_raw_parts(bin_limits, bins + 1).to_vec(),
            subgrid_params,
            extra,
            &subgrid_type,
        )
        .unwrap(),
    );

    if !key_vals.is_null() {
        let keyval = &*key_vals;

        if let Some(value) = keyval.strings.get("initial_state_1") {
            grid.set_key_value("initial_state_1", value.to_str().unwrap());
        }

        if let Some(value) = keyval.strings.get("initial_state_2") {
            grid.set_key_value("initial_state_2", value.to_str().unwrap());
        }
    }

    grid
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
    let filename = CStr::from_ptr(filename).to_string_lossy();
    let reader = BufReader::new(File::open(filename.as_ref()).unwrap());

    Box::new(Grid::read(reader).unwrap())
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
        (*grid).merge(*other).unwrap();
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
    (&mut *grid).scale(factor);
}

/// Optimizes the grid representation for space efficiency.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_optimize(grid: *mut Grid) {
    (&mut *grid).optimize();
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
    (&mut *grid).scale_by_order(alphas, alpha, logxir, logxif, global);
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
    (*grid).set_key_value(
        CStr::from_ptr(key).to_string_lossy().as_ref(),
        CStr::from_ptr(value).to_string_lossy().as_ref(),
    );
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
    let grid = &mut *grid;
    let bins = grid.bin_info().bins();

    grid.set_remapper(
        BinRemapper::new(
            slice::from_raw_parts(normalizations, bins).to_vec(),
            slice::from_raw_parts(limits, 2 * dimensions * bins)
                .chunks_exact(2)
                .map(|chunk| (chunk[0], chunk[1]))
                .collect(),
        )
        .unwrap(),
    )
    .unwrap();
}

/// Write `grid` to a file with name `filename`.
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
    let filename = CStr::from_ptr(filename).to_string_lossy();
    let writer = BufWriter::new(File::create(filename.as_ref()).unwrap());

    (*grid).write(writer).unwrap();
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
    let factors = if factors.is_null() {
        vec![1.0; combinations]
    } else {
        slice::from_raw_parts(factors, combinations).to_vec()
    };

    (*lumi).0.push(LumiEntry::new(
        slice::from_raw_parts(pdg_id_pairs, 2 * combinations)
            .chunks(2)
            .zip(factors)
            .map(|x| (x.0[0], x.0[1], x.1))
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
    (*lumi).0[entry].entry().len()
}

/// Returns the number of channel for the luminosity function `lumi`.
///
/// # Safety
///
/// The parameter `lumi` must point to a valid `Lumi` object created by `pineappl_lumi_new` or
/// `pineappl_grid_lumi`.
#[no_mangle]
pub unsafe extern "C" fn pineappl_lumi_count(lumi: *mut Lumi) -> usize {
    (*lumi).0.len()
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
    lumi: *mut Lumi,
    entry: usize,
    pdg_ids: *mut i32,
    factors: *mut f64,
) {
    let entry = (*lumi).0[entry].entry();
    let pdg_ids = slice::from_raw_parts_mut(pdg_ids, 2 * entry.len());
    let factors = slice::from_raw_parts_mut(factors, entry.len());

    entry
        .iter()
        .flat_map(|(id1, id2, _)| vec![id1, id2])
        .zip(pdg_ids.iter_mut())
        .for_each(|(from, to)| *to = *from);
    entry
        .iter()
        .map(|(_, _, factor)| factor)
        .zip(factors.iter_mut())
        .for_each(|(from, to)| *to = *from);
}

/// Creates a new luminosity function and returns a pointer to it. If no longer needed, the object
/// should be deleted using `pineappl_lumi_delete`.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_lumi_new() -> Box<Lumi> {
    Box::new(Lumi::default())
}

/// Fills `buffer` with the q2-slice for the index `q2_slice` of the grid for the specified
/// `order`, `bin`, and `lumi`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. `buffer` must be as large as the square of the return value
/// of `pineappl_subgrid_x_grid_count`.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_export_mu2_slice(
    grid: *const Grid,
    order: usize,
    bin: usize,
    lumi: usize,
    q2_slice: usize,
    buffer: *mut f64,
) {
    let subgrid = (*grid).subgrid(order, bin, lumi);
    let x1_len = subgrid.x1_grid().len();
    let slice = slice::from_raw_parts_mut(buffer, x1_len * subgrid.x2_grid().len());
    subgrid
        .iter()
        .filter(|((imu2, _, _), _)| *imu2 == q2_slice)
        .for_each(|((_, ix1, ix2), &value)| slice[ix1 + x1_len * ix2] = value);
}

/// Write into `tuple` the lower and upper limit of filled q2 slices for the grid with the
/// specified indices. The last slice that is filled is one minus the upper limit.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. `tuple` must point to an array with two elements.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_nonzero_mu2_slices(
    grid: *const Grid,
    order: usize,
    bin: usize,
    lumi: usize,
    tuple: *mut usize,
) {
    let tuple = slice::from_raw_parts_mut(tuple, 2);
    let mut iter = (*grid).subgrid(order, bin, lumi).iter();

    if let Some(((first, _, _), _)) = iter.next() {
        tuple[0] = first;

        if let Some(((last, _, _), _)) = iter.last() {
            tuple[1] = last + 1;
        } else {
            tuple[0] = first + 1;
        }
    } else {
        tuple[0] = 0;
        tuple[1] = 0;
    }
}

/// Deletes a subgrid created with [`pineappl_subgrid_new2`]. If `subgrid` is the null pointer,
/// nothing is done.
#[no_mangle]
#[allow(unused_variables)]
pub extern "C" fn pineappl_subgrid_delete(subgrid: Option<Box<SubGrid>>) {}

/// This function takes replaces the subgrid in `grid` with the corresponding indices `order`,
/// `bin` and `lumi` with the one given in `subgrid`. If `subgrid` is the null pointer, the specied
/// subgrid is replaced with an empty one.
///
/// # Safety
///
/// Both `grid` and `subgrid` must point to valid objects. The parameter `subgrid` can be the null
/// pointer.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_replace_and_delete(
    grid: *mut Grid,
    subgrid: Option<Box<SubGrid>>,
    order: usize,
    bin: usize,
    lumi: usize,
) {
    (*grid).set_subgrid(
        order,
        bin,
        lumi,
        subgrid.map_or_else(
            || EmptySubgridV1::default().into(),
            |subgrid| subgrid.0.into(),
        ),
    );
}

/// Creates a new subgrid, which can be filled with [`pineappl_subgrid_import_mu2_slice`]. The
/// array `mu2_grid` must contain the (squared) values of the renormalization and then the
/// factorization scale, such that twice the value of `mu2_grid_len` gives the length of the array.
///
/// # Safety
///
/// The pointers `mu2_grid`, `x1_grid` and `x2_grid` must be non-`NULL` and array. Furthermore, the
/// array `mu2_grid` must be *twice* as long as given in `mu2_grid_len`, and the arrays `x1_grid`
/// and `x2_grid` as long as specified by `x1_grid_len` and `x2_grid_len`, respectively.
#[no_mangle]
pub unsafe extern "C" fn pineappl_subgrid_new2(
    mu2_grid_len: usize,
    mu2_grid: *const f64,
    x1_grid_len: usize,
    x1_grid: *const f64,
    x2_grid_len: usize,
    x2_grid: *const f64,
) -> Box<SubGrid> {
    let mu2: Vec<_> = slice::from_raw_parts(mu2_grid, 2 * mu2_grid_len)
        .chunks_exact(2)
        .map(|mu2| Mu2 {
            ren: mu2[0],
            fac: mu2[1],
        })
        .collect();
    let x1 = slice::from_raw_parts(x1_grid, x1_grid_len).to_vec();
    let x2 = slice::from_raw_parts(x2_grid, x2_grid_len).to_vec();

    Box::new(SubGrid(ImportOnlySubgridV2::new(
        SparseArray3::new(mu2.len(), x1.len(), x2.len()),
        mu2,
        x1,
        x2,
    )))
}

/// Imports `slice` for the given index into `subgrid`.
///
/// # Safety
///
/// The parameter `subgrid` and the array `slice` must be non-`NULL` and `slice` must be at least
/// as long as the product `x1_grid_len * x2_grid_len` that were used to create the subgrid with.
/// The index `mu2_slice` must be smaller than `mu2_grid_len` that was used in
/// [`pineappl_subgrid_new2`] to create `subgrid`.
#[no_mangle]
pub unsafe extern "C" fn pineappl_subgrid_import_mu2_slice(
    subgrid: *mut SubGrid,
    mu2_slice: usize,
    slice: *const f64,
) {
    let array = (*subgrid).0.array_mut();
    let (_, nx1, nx2) = array.dimensions();
    let slice = slice::from_raw_parts(slice, nx1 * nx2);

    for (index, &value) in slice.iter().enumerate().filter(|(_, &value)| value != 0.0) {
        let ix1 = index / nx2;
        let ix2 = index % nx2;

        array[[mu2_slice, ix1, ix2]] = value;
    }
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
    (*key_vals).bools[CStr::from_ptr(key).to_string_lossy().as_ref()]
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
    (*key_vals).doubles[CStr::from_ptr(key).to_string_lossy().as_ref()]
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
    (*key_vals).ints[CStr::from_ptr(key).to_string_lossy().as_ref()]
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
    (*key_vals).strings[CStr::from_ptr(key).to_string_lossy().as_ref()].as_ptr()
}

/// Return a pointer to newly-created `pineappl_keyval` object.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_keyval_new() -> Box<KeyVal> {
    Box::new(KeyVal::default())
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
    (*key_vals)
        .bools
        .insert(CStr::from_ptr(key).to_string_lossy().into_owned(), value);
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
    (*key_vals)
        .doubles
        .insert(CStr::from_ptr(key).to_string_lossy().into_owned(), value);
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
    (*key_vals)
        .ints
        .insert(CStr::from_ptr(key).to_string_lossy().into_owned(), value);
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
    (*key_vals).strings.insert(
        CStr::from_ptr(key).to_string_lossy().into_owned(),
        CString::from(CStr::from_ptr(value)),
    );
}
