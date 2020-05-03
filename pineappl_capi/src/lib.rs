#![warn(clippy::all, clippy::cargo, clippy::nursery, clippy::pedantic)]
#![warn(missing_docs)]

//! C-language interface for `PineAPPL`.

use pineappl_core::grid::{Grid, Ntuple, Order};
use pineappl_core::lumi::{Lumi, LumiEntry};
use std::collections::HashMap;
use std::ffi::{CStr, CString};
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::os::raw::{c_char, c_void};
use std::slice;

// TODO: change raw pointer and Option of pointer to Box if possible, as soon as cbindgen supports
// this; see github issue: https://github.com/eqrion/cbindgen/issues/474

/// Convolutes the specified grid with the PDFs `xfx1` and `xfx2` and strong coupling `alphas`.
/// These functions must evaluate the PDFs for the given `x` and `q2` and write the results for all
/// partons into `pdfs`. Note that the value must be the PDF for the given `pdg_id` multiplied with
/// `x`. The value of the pointer `state` provided to these functions is the same one given to this
/// function. The parameter `order_mask` must be as long as there are perturbative orders contained
/// in `grid` and is used to selectively disable (`false`) or enable (`true`) individual orders. If
/// `order_mask` is set to `NULL`, all orders are active. The parameter `lumi_mask` can be used
/// similarly, but must be as long as the luminosity function `grid` was created with has entries,
/// or `NULL`. The values `xi_ren` and `xi_fac` can be used to vary the renormalization and
/// factorization from its central value, which corresponds to `1.0`. After convolution of the grid
/// with the PDFS value of the observable for each bin is written into `results`.
#[no_mangle]
pub extern "C" fn pineappl_grid_convolute(
    grid: Option<*const Grid>,
    xfx1: extern "C" fn(x: f64, q2: f64, pdg_id: i32, state: *mut c_void) -> f64,
    xfx2: extern "C" fn(x: f64, q2: f64, pdg_id: i32, state: *mut c_void) -> f64,
    alphas: extern "C" fn(q2: f64, state: *mut c_void) -> f64,
    state: *mut c_void,
    order_mask: Option<*const bool>,
    lumi_mask: Option<*const bool>,
    xi_ren: f64,
    xi_fac: f64,
    results: Option<*mut f64>,
) {
    let grid = unsafe { &*grid.unwrap() };
    let xfx1 = |x, q2, id| xfx1(x, q2, id, state);
    let xfx2 = |x, q2, id| xfx2(x, q2, id, state);
    let alphas = |q2| alphas(q2, state);
    let order_mask = if let Some(order_mask) = order_mask {
        unsafe { slice::from_raw_parts(order_mask, grid.orders().len()) }.to_vec()
    } else {
        vec![]
    };
    let lumi_mask = if let Some(lumi_mask) = lumi_mask {
        unsafe { slice::from_raw_parts(lumi_mask, grid.lumi().len()) }.to_vec()
    } else {
        vec![]
    };
    let results = unsafe { slice::from_raw_parts_mut(results.unwrap(), grid.bin_limits().bins()) };

    results.copy_from_slice(&grid.convolute(
        &xfx1,
        &xfx2,
        &alphas,
        &order_mask,
        &lumi_mask,
        &(xi_ren, xi_fac),
    ));
}

/// Delete a grid previously created with `pineappl_grid_new`.
#[no_mangle]
pub extern "C" fn pineappl_grid_delete(grid: Option<*mut Grid>) {
    unsafe {
        Box::from_raw(grid.unwrap());
    }
}

/// Fill `grid` for the given momentum fractions `x1` and `x2`, at the scale `q2` for the given
/// value of the `order`, `observable`, and `lumi` with `weight`.
#[no_mangle]
pub extern "C" fn pineappl_grid_fill(
    grid: Option<*mut Grid>,
    x1: f64,
    x2: f64,
    q2: f64,
    order: usize,
    observable: f64,
    lumi: usize,
    weight: f64,
) {
    let grid = unsafe { &mut *grid.unwrap() };

    grid.fill(order, observable, lumi, Ntuple { x1, x2, q2, weight });
}

/// Fill `grid` for the given momentum fractions `x1` and `x2`, at the scale `q2` for the given
/// value of the `order` and `observable` with `weights`. The parameter of weight must contain a
/// result for entry of the luminosity function the grid was created with.
#[no_mangle]
pub extern "C" fn pineappl_grid_fill_all(
    grid: Option<*mut Grid>,
    x1: f64,
    x2: f64,
    q2: f64,
    order: usize,
    observable: f64,
    weights: Option<*const f64>,
) {
    let grid = unsafe { &mut *grid.unwrap() };
    let weights = unsafe { slice::from_raw_parts(weights.unwrap(), grid.lumi().len()) };

    grid.fill_all(
        order,
        observable,
        Ntuple {
            x1,
            x2,
            q2,
            weight: (),
        },
        weights,
    );
}

/// Write the order parameters of `grid` into `order_params`.
#[no_mangle]
pub extern "C" fn pineappl_grid_get_order_params(
    grid: Option<*const Grid>,
    order_params: Option<*mut u32>,
) {
    let orders = unsafe { &*grid.unwrap() }.orders();
    let order_params =
        unsafe { slice::from_raw_parts_mut(order_params.unwrap(), 4 * orders.len()) };

    for (i, order) in orders.iter().enumerate() {
        order_params[4 * i] = order.alphas;
        order_params[4 * i + 1] = order.alpha;
        order_params[4 * i + 2] = order.logxir;
        order_params[4 * i + 3] = order.logxif;
    }
}

/// Return the number of orders in `grid`.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_grid_get_order_count(grid: Option<*const Grid>) -> usize {
    let subgrids = unsafe { &*grid.unwrap() }.orders();

    subgrids.len()
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
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_grid_new(
    lumi: Option<*const Lumi>,
    orders: usize,
    order_params: Option<*const u32>,
    bins: usize,
    bin_limits: Option<*const f64>,
    key_vals: Option<*const KeyVal>,
) -> *mut Grid {
    let lumi = unsafe { &*lumi.unwrap() };
    let order_params = unsafe { slice::from_raw_parts(order_params.unwrap(), 4 * orders) };
    let orders = order_params
        .chunks(4)
        .map(|s| Order {
            alphas: s[0],
            alpha: s[1],
            logxir: s[2],
            logxif: s[3],
        })
        .collect::<Vec<_>>();
    let bin_limits = unsafe { slice::from_raw_parts(bin_limits.unwrap(), bins + 1) };

    if let Some(key_vals) = key_vals {
        // TODO: do something with the contents
        let _keyval = unsafe { &*key_vals };
    }

    Box::into_raw(Box::new(Grid::new(
        lumi.clone(),
        orders,
        bin_limits.to_vec(),
    )))
}

/// Read a `PineAPPL` grid from a file with name `filename`.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_grid_read(filename: Option<*const c_char>) -> *mut Grid {
    let filename = String::from(
        unsafe { CStr::from_ptr(filename.unwrap()) }
            .to_str()
            .unwrap(),
    );
    let reader = BufReader::new(File::open(filename).unwrap());
    let grid = Box::new(bincode::deserialize_from(reader).unwrap());

    Box::into_raw(grid)
}

/// Scale all grids in `grid` by `factor`.
#[no_mangle]
pub extern "C" fn pineappl_grid_scale(grid: Option<*mut Grid>, factor: f64) {
    let grid = unsafe { &mut *grid.unwrap() };

    grid.scale(factor);
}

/// Write `grid` to a file with name `filename`.
#[no_mangle]
pub extern "C" fn pineappl_grid_write(grid: Option<*const Grid>, filename: Option<*const c_char>) {
    let grid = unsafe { &*grid.unwrap() };
    let filename = String::from(
        unsafe { CStr::from_ptr(filename.unwrap()) }
            .to_str()
            .unwrap(),
    );
    let writer = BufWriter::new(File::create(filename).unwrap());

    bincode::serialize_into(writer, &grid).unwrap();
}

/// Adds a linear combination of initial states to the luminosity function `lumi`.
#[no_mangle]
pub extern "C" fn pineappl_lumi_add(
    lumi: Option<*mut Lumi>,
    combinations: usize,
    pdg_id_pairs: Option<*const i32>,
    factors: Option<*const f64>,
) {
    let lumi = unsafe { &mut *lumi.unwrap() };
    let ids = unsafe { slice::from_raw_parts(pdg_id_pairs.unwrap(), 2 * combinations) };

    let factors = if let Some(factors) = factors {
        unsafe { slice::from_raw_parts(factors, combinations) }.to_vec()
    } else {
        vec![1.0; combinations]
    };

    lumi.add(LumiEntry::new(
        ids.chunks(2)
            .zip(factors)
            .map(|x| (x.0[0], x.0[1], x.1))
            .collect(),
    ));
}

/// Delete luminosity function previously created with `pineappl_lumi_new`.
#[no_mangle]
pub extern "C" fn pineappl_lumi_delete(lumi: Option<*mut Lumi>) {
    unsafe {
        Box::from_raw(lumi.unwrap());
    }
}

/// Creates a new luminosity function and returns a pointer to it. If no longer needed, the object
/// should be deleted using `pineappl_lumi_delete`.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_lumi_new() -> *mut Lumi {
    Box::into_raw(Box::new(Lumi::default()))
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
pub extern "C" fn pineappl_keyval_delete(key_vals: Option<*mut KeyVal>) {
    unsafe {
        Box::from_raw(key_vals.unwrap());
    }
}

/// Get the boolean-valued parameter with name `key` stored in `key_vals`.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_keyval_get_bool(
    key_vals: Option<*const KeyVal>,
    key: Option<*const c_char>,
) -> bool {
    let key_vals = unsafe { &*key_vals.unwrap() };
    let key = String::from(unsafe { CStr::from_ptr(key.unwrap()) }.to_str().unwrap());

    key_vals.bools[&key]
}

/// Get the double-valued parameter with name `key` stored in `key_vals`.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_keyval_get_double(
    key_vals: Option<*const KeyVal>,
    key: Option<*const c_char>,
) -> f64 {
    let key_vals = unsafe { &*key_vals.unwrap() };
    let key = String::from(unsafe { CStr::from_ptr(key.unwrap()) }.to_str().unwrap());

    key_vals.doubles[&key]
}

/// Get the string-valued parameter with name `key` stored in `key_vals`.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_keyval_get_int(
    key_vals: Option<*const KeyVal>,
    key: Option<*const c_char>,
) -> i32 {
    let key_vals = unsafe { &*key_vals.unwrap() };
    let key = String::from(unsafe { CStr::from_ptr(key.unwrap()) }.to_str().unwrap());

    key_vals.ints[&key]
}

/// Get the int-valued parameter with name `key` stored in `key_vals`.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_keyval_get_string(
    key_vals: Option<*const KeyVal>,
    key: Option<*const c_char>,
) -> *const c_char {
    let key_vals = unsafe { &*key_vals.unwrap() };
    let key = String::from(unsafe { CStr::from_ptr(key.unwrap()) }.to_str().unwrap());

    key_vals.strings[&key].as_ptr()
}

/// Return a pointer to newly-created `pineappl_keyval` object.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_keyval_new() -> *mut KeyVal {
    Box::into_raw(Box::new(KeyVal::default()))
}

/// Set the double-valued parameter with name `key` to `value` in `key_vals`.
#[no_mangle]
pub extern "C" fn pineappl_keyval_set_bool(
    key_vals: Option<*mut KeyVal>,
    key: Option<*const c_char>,
    value: bool,
) {
    let key_vals = unsafe { &mut *key_vals.unwrap() };
    let key = String::from(unsafe { CStr::from_ptr(key.unwrap()) }.to_str().unwrap());

    key_vals.bools.insert(key, value);
}

/// Set the double-valued parameter with name `key` to `value` in `key_vals`.
#[no_mangle]
pub extern "C" fn pineappl_keyval_set_double(
    key_vals: Option<*mut KeyVal>,
    key: Option<*const c_char>,
    value: f64,
) {
    let key_vals = unsafe { &mut *key_vals.unwrap() };
    let key = String::from(unsafe { CStr::from_ptr(key.unwrap()) }.to_str().unwrap());

    key_vals.doubles.insert(key, value);
}

/// Set the int-valued parameter with name `key` to `value` in `key_vals`.
#[no_mangle]
pub extern "C" fn pineappl_keyval_set_int(
    key_vals: Option<*mut KeyVal>,
    key: Option<*const c_char>,
    value: i32,
) {
    let key_vals = unsafe { &mut *key_vals.unwrap() };
    let key = String::from(unsafe { CStr::from_ptr(key.unwrap()) }.to_str().unwrap());

    key_vals.ints.insert(key, value);
}

/// Set the string-valued parameter with name `key` to `value` in `key_vals`.
#[no_mangle]
pub extern "C" fn pineappl_keyval_set_string(
    key_vals: Option<*mut KeyVal>,
    key: Option<*const c_char>,
    value: Option<*const c_char>,
) {
    let key_vals = unsafe { &mut *key_vals.unwrap() };
    let key = String::from(unsafe { CStr::from_ptr(key.unwrap()) }.to_str().unwrap());
    let value = CString::from(unsafe { CStr::from_ptr(value.unwrap()) });

    key_vals.strings.insert(key, value);
}
