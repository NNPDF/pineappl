#![warn(clippy::all, clippy::cargo, clippy::nursery, clippy::pedantic)]
#![warn(missing_docs)]

//! C-language interface for `PineAPPL`.

use itertools::izip;
use pineappl::grid::{Grid, Ntuple, Order};
use pineappl::lagrange_subgrid::{LagrangeSparseSubgridV1, LagrangeSubgridV1, LagrangeSubgridV2};
use pineappl::lumi::LumiEntry;
use pineappl::subgrid::{ExtraSubgridParams, Subgrid, SubgridEnum, SubgridParams};
use std::collections::HashMap;
use std::convert::TryFrom;
use std::ffi::{CStr, CString};
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::os::raw::{c_char, c_void};
use std::slice;

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
pub struct SubGrid(SubgridEnum);

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
pub unsafe extern "C" fn pineappl_grid_bin_sizes(grid: *const Grid, bin_sizes: *mut f64) {
    let sizes = (*grid).bin_info().normalizations();
    let bin_sizes = slice::from_raw_parts_mut(bin_sizes, sizes.len());

    for (i, size) in sizes.iter().enumerate() {
        bin_sizes[i] = *size;
    }
}

/// Stores the bin limits of `grid` in `bin_limits`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. The parameter `bin_limits` must point to an array that is
/// one element longer as `grid` has bins.
#[deprecated(
    since = "0.4.0",
    note = "this function only returns the limits for the zeroth dimension"
)]
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_bin_limits(grid: *const Grid, bin_limits: *mut f64) {
    let mut limits = (*grid).bin_info().left(0);
    limits.push(*(*grid).bin_info().right(0).last().unwrap());
    slice::from_raw_parts_mut(bin_limits, limits.len()).copy_from_slice(&limits);
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
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, this function is not safe to call. If `grid`
/// is a `NULL` pointer, this function does not do anything.
#[no_mangle]
#[allow(unused_variables)]
pub unsafe extern "C" fn pineappl_grid_delete(grid: Option<Box<Grid>>) {}

/// Performs an operation `name` on `grid` using as input or output parameters `key_vals`. This is
/// used to get access to functions that are otherwise not available through other functions. If
/// the operation was successful, returns `true`. Otherwise, or if the `name` wasn't recognized
/// `false` is returned. The parameter `_key_vals` is currently ignored.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. The parameter `name` must be a valid and non-`NULL` C
/// string.
#[no_mangle]
#[must_use]
pub unsafe extern "C" fn pineappl_grid_ext(
    grid: *mut Grid,
    name: *const c_char,
    _key_vals: *mut KeyVal,
) -> bool {
    let grid = &mut *grid;
    let name = CStr::from_ptr(name).to_str().unwrap();

    if name == "optimise" {
        grid.scale(0.0);

        true
    } else {
        false
    }
}

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
#[no_mangle]
#[must_use]
pub unsafe extern "C" fn pineappl_grid_read(filename: *const c_char) -> Box<Grid> {
    let filename = CStr::from_ptr(filename).to_str().unwrap();
    let reader = BufReader::new(File::open(filename).unwrap());

    Box::new(Grid::read(reader).unwrap())
}

/// Merges `other` into `grid` and subsequently deletes `other`.
///
/// # Safety
///
/// Both `grid` and `other` must be valid `Grid` objects created by either `pineappl_grid_new` or
/// `pineappl_grid_read`. If `other` is a `NULL` pointer, this function does not do anything.
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

/// Write `grid` to a file with name `filename`.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. The parameter must be a non-`NULL`, non-empty, and valid C
/// string pointing to a non-existing, but writable file.
#[no_mangle]
pub unsafe extern "C" fn pineappl_grid_write(grid: *const Grid, filename: *const c_char) {
    let filename = CStr::from_ptr(filename).to_str().unwrap();
    let writer = BufWriter::new(File::create(filename).unwrap());

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
///
/// # Safety
///
/// The parameter `lumi` must point to a valid `Lumi` object created by `pineappl_lumi_new`. If
/// `lumi` is a `NULL` pointer, this function does not do anything.
#[no_mangle]
#[allow(unused_variables)]
pub unsafe extern "C" fn pineappl_lumi_delete(lumi: Option<Box<Lumi>>) {}

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
pub unsafe extern "C" fn pineappl_subgrid_q2_slice(
    grid: *const Grid,
    order: usize,
    bin: usize,
    lumi: usize,
    q2_slice: usize,
    buffer: *mut f64,
) {
    let size = pineappl_subgrid_x_grid_count(grid).pow(2);
    (*grid)
        .subgrid(order, bin, lumi)
        .fill_q2_slice(q2_slice, slice::from_raw_parts_mut(buffer, size));
}

/// Write into `tuple` the lower and upper limit of filled q2 slices for the grid with the
/// specified indices. The last slice that is filled is one minus the upper limit.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. `tuple` must point to an array with two elements.
#[no_mangle]
pub unsafe extern "C" fn pineappl_subgrid_filled_q2_slices(
    grid: *const Grid,
    order: usize,
    bin: usize,
    lumi: usize,
    tuple: *mut usize,
) {
    let tuple = slice::from_raw_parts_mut(tuple, 2);
    let slice = (*grid).subgrid(order, bin, lumi).q2_slice();
    tuple[0] = slice.0;
    tuple[1] = slice.1;
}

/// Writes into `buffer` the q2 values that the grid spans.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. The parameter `buffer` must point to an array which is as
/// large as the value returned by `pineappl_subgrid_q2_grid_count`.
#[no_mangle]
pub unsafe extern "C" fn pineappl_subgrid_q2_grid(grid: *const Grid, buffer: *mut f64) {
    let size = pineappl_subgrid_q2_grid_count(grid);
    slice::from_raw_parts_mut(buffer, size).copy_from_slice(&(*grid).subgrid(0, 0, 0).q2_grid());
}

/// Returns the number grid values in q2 direction.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_subgrid_q2_grid_count(grid: *const Grid) -> usize {
    (*grid).subgrid(0, 0, 0).q2_grid().len()
}

/// Writes into `buffer` the x values that the grid spans.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call. The parameter `buffer` must point to an array which is as
/// large as the value returned by `pineappl_subgrid_x_grid_count`.
#[no_mangle]
pub unsafe extern "C" fn pineappl_subgrid_x_grid(grid: *const Grid, buffer: *mut f64) {
    let size = pineappl_subgrid_x_grid_count(grid);
    slice::from_raw_parts_mut(buffer, size).copy_from_slice(&(*grid).subgrid(0, 0, 0).x1_grid());
}

/// Returns the number grid values in x direction.
///
/// # Safety
///
/// If `grid` does not point to a valid `Grid` object, for example when `grid` is the null pointer,
/// this function is not safe to call.
#[no_mangle]
pub unsafe extern "C" fn pineappl_subgrid_x_grid_count(grid: *const Grid) -> usize {
    (*grid).subgrid(0, 0, 0).x1_grid().len()
}

/// Creates a new subgrid, using the paramters stored in `keyvals`. If `key_vals` is the null
/// pointer, default parameters are chosen.
///
/// # Safety
///
/// The parameter `key_vals` must be either the null pointer, or point to a valid `KeyVal` object.
#[no_mangle]
pub unsafe extern "C" fn pineappl_subgrid_new(key_vals: *const KeyVal) -> Box<SubGrid> {
    let (subgrid_type, subgrid_params, extra) = grid_params(key_vals);

    match subgrid_type.as_str() {
        "LagrangeSubgrid" | "LagrangeSubgridV1" => {
            Box::new(SubGrid(LagrangeSubgridV1::new(&subgrid_params).into()))
        }
        "LagrangeSubgridV2" => Box::new(SubGrid(
            LagrangeSubgridV2::new(&subgrid_params, &extra).into(),
        )),
        "LagrangeSparseSubgrid" => Box::new(SubGrid(
            LagrangeSparseSubgridV1::new(&subgrid_params).into(),
        )),
        _ => {
            panic!();
        }
    }
}

/// Deletes a subgrid created with `pineappl_subgrid_new`. If `subgrid` is the null pointer,
/// nothing is done.
///
/// # Safety
///
/// The parameter `subgrid` must point to a valid `Subgrid` object created by
/// `pineappl_subgrid_new` or be the null pointer.
#[no_mangle]
#[allow(unused_variables)]
pub unsafe extern "C" fn pineappl_subgrid_delete(subgrid: Option<Box<SubGrid>>) {}

/// This function takes replaces the subgrid in `grid` with the corresponding indices `order`,
/// `bin` and `lumi` with the one given in `subgrid`. If `subgrid` is the null pointer, nothing is
/// done.
///
/// # Safety
///
/// Both `grid` and `subgrid` must point to valid objects. The parameter `subgrid` can be the null
/// pointer.
#[no_mangle]
pub unsafe extern "C" fn pineappl_subgrid_replace_and_delete(
    grid: *mut Grid,
    subgrid: Option<Box<SubGrid>>,
    order: usize,
    bin: usize,
    lumi: usize,
) {
    if let Some(subgrid) = subgrid {
        (*grid).set_subgrid(order, bin, lumi, (*subgrid).0);
    }
}

/// Fills the slice with the given `index` of `subgrid` with the contents of `slice`.
///
/// # Safety
///
/// The parameter `subgrid` must point to a valid `SubGrid` object and be non-null, `index` must be
/// an index smaller than the parameter `q2_bins` the `SubGrid` was created with, and `slice` must
/// be an array at least as large as the square of the parameter `x_bins` the `SubGrid` was
/// created with.
#[no_mangle]
pub unsafe extern "C" fn pineappl_subgrid_fill_q2_slice(
    subgrid: *mut SubGrid,
    index: usize,
    slice: *const f64,
) {
    (*subgrid).0.write_q2_slice(
        index,
        slice::from_raw_parts(slice, (*subgrid).0.x1_grid().len().pow(2)),
    );
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
///
/// # Safety
///
/// The parameter `key_vals` must point to a valid `KeyVal` object created by
/// `pineappl_keyval_new`. If `key_vals` is a `NULL` pointer, this function does not do anything.
#[no_mangle]
#[allow(unused_variables)]
pub unsafe extern "C" fn pineappl_keyval_delete(key_vals: Option<Box<KeyVal>>) {}

/// Get the boolean-valued parameter with name `key` stored in `key_vals`.
///
/// # Safety
///
/// The parameter `key_vals` must point to a valid `KeyVal` object created by
/// `pineappl_keyval_new`. `key` must be a valid C string.
#[no_mangle]
#[must_use]
pub unsafe extern "C" fn pineappl_keyval_bool(key_vals: *const KeyVal, key: *const c_char) -> bool {
    (*key_vals).bools[CStr::from_ptr(key).to_str().unwrap()]
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
    (*key_vals).doubles[CStr::from_ptr(key).to_str().unwrap()]
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
    (*key_vals).ints[CStr::from_ptr(key).to_str().unwrap()]
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
    (*key_vals).strings[CStr::from_ptr(key).to_str().unwrap()].as_ptr()
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
        .insert(CStr::from_ptr(key).to_str().unwrap().to_string(), value);
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
        .insert(CStr::from_ptr(key).to_str().unwrap().to_string(), value);
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
        .insert(CStr::from_ptr(key).to_str().unwrap().to_string(), value);
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
        CStr::from_ptr(key).to_str().unwrap().to_string(),
        CString::from(CStr::from_ptr(value)),
    );
}
