#![warn(clippy::all, clippy::cargo, clippy::nursery, clippy::pedantic)]
#![warn(missing_docs)]

//! C-language interface for `PineAPPL`.

use pineappl_core::grid::{Grid, Order};
use pineappl_core::lumi::{Lumi, LumiEntry};
use std::collections::HashMap;
use std::ffi::{CStr, CString};
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::os::raw::c_char;
use std::slice;

// TODO: change raw pointer and Option of pointer to Box if possible, as soon as cbindgen supports
// this; see github issue: https://github.com/eqrion/cbindgen/issues/474

/// Delete a grid previously created with `pineappl_grid_new`.
#[no_mangle]
pub extern "C" fn pineappl_grid_delete(grid: Option<*mut Grid>) {
    unsafe {
        Box::from_raw(grid.unwrap());
    }
}

/// long the corresponding luminosity function the grid was created with and contain the
/// corresponding weights at each index. The value `grid_index` selects one of the subgrids whose
/// meaning was specified during creation with `grid_parameters` in `pineappl_grid_new`.
#[no_mangle]
pub extern "C" fn pineappl_grid_fill(
    grid: Option<*mut Grid>,
    x1: f64,
    x2: f64,
    q2: f64,
    observable: f64,
    weights: Option<*const f64>,
    subgrid: usize,
) {
    let grid = unsafe { &mut *grid.unwrap() };
    let weights = unsafe { slice::from_raw_parts(weights.unwrap(), grid.lumi().len()) };

    grid.fill(x1, x2, q2, observable, weights, subgrid);
}

/// Write the subgrid parameters of `grid` into `subgrid_params`.
#[no_mangle]
pub extern "C" fn pineappl_grid_get_subgrid_params(
    grid: Option<*const Grid>,
    subgrid_params: Option<*mut u32>,
) {
    let orders = unsafe { &*grid.unwrap() }.orders();
    let subgrid_params =
        unsafe { slice::from_raw_parts_mut(subgrid_params.unwrap(), 4 * orders.len()) };

    for (i, order) in orders.iter().enumerate() {
        subgrid_params[4 * i] = order.alphas;
        subgrid_params[4 * i + 1] = order.alpha;
        subgrid_params[4 * i + 2] = order.logxir;
        subgrid_params[4 * i + 3] = order.logxif;
    }
}

/// Return the number of subgrids in `grid`.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_grid_get_subgrids(grid: Option<*const Grid>) -> usize {
    let subgrids = unsafe { &*grid.unwrap() }.orders();

    subgrids.len()
}

/// Create a new `pineappl_grid`. The creation requires four different sets of parameters:
/// see `pineappl_lumi_new`.
/// - Subgrid specification: `format`, `subgrids`, and `subgrids_params`. Each `PineAPPL` grid
/// contains a number of subgrids, given as `subgrids`, which usually store the grids for each
/// perturbative order separately. The concrete meaning of the subgrids is determined by `format`
/// and `subgrids_params`.
/// - Observable definition: `bins` and `bin_limits`. Each subgrid can store observables from a
/// one-dimensional distribution. To this end `bins` specifies how many observables are stored and
/// with `bins + 1` entries.
/// - More complex information can be given in a key-value storage `key_vals`.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_grid_new(
    lumi: Option<*const Lumi>,
    subgrids: usize,
    subgrid_params: Option<*const u32>,
    bins: usize,
    bin_limits: Option<*const f64>,
    key_vals: Option<*const KeyVal>,
) -> *mut Grid {
    let lumi = unsafe { &*lumi.unwrap() };
    let subgrid_params = unsafe { slice::from_raw_parts(subgrid_params.unwrap(), 4 * subgrids) };
    let subgrid_data = subgrid_params
        .chunks(4)
        .map(|s| Order {
            alphas: s[0],
            alpha: s[1],
            logxir: s[2],
            logxif: s[3],
        })
        .collect::<Vec<_>>();
    let bin_limits = unsafe { slice::from_raw_parts(bin_limits.unwrap(), bins + 1) };

    // TODO: do something with the contents
    let _keyval = unsafe { &*key_vals.unwrap() };

    Box::into_raw(Box::new(Grid::new(
        lumi.clone(),
        subgrid_data,
        bin_limits.to_vec(),
    )))
}

/// Read a `pineappl_grid` from a file with name `filename`.
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

/// Key-value storage.
#[derive(Default)]
pub struct KeyVal {
    _method: String,
    bools: HashMap<String, bool>,
    doubles: HashMap<String, f64>,
    ints: HashMap<String, i32>,
    strings: HashMap<String, CString>,
}

/// Delete the previously created object pointed to by `storage`.
#[no_mangle]
pub extern "C" fn pineappl_keyval_delete(storage: Option<*mut KeyVal>) {
    unsafe {
        Box::from_raw(storage.unwrap());
    }
}

/// Get the boolean-valued parameter with name `key` for `storage`.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_keyval_get_bool(
    storage: Option<*const KeyVal>,
    key: Option<*const c_char>,
) -> bool {
    let storage = unsafe { &*storage.unwrap() };
    let key = String::from(unsafe { CStr::from_ptr(key.unwrap()) }.to_str().unwrap());

    storage.bools[&key]
}

/// Get the double-valued parameter with name `key` for `storage`.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_keyval_get_double(
    storage: Option<*const KeyVal>,
    key: Option<*const c_char>,
) -> f64 {
    let storage = unsafe { &*storage.unwrap() };
    let key = String::from(unsafe { CStr::from_ptr(key.unwrap()) }.to_str().unwrap());

    storage.doubles[&key]
}

/// Get the string-valued parameter with name `key` for `storage`.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_keyval_get_int(
    storage: Option<*const KeyVal>,
    key: Option<*const c_char>,
) -> i32 {
    let storage = unsafe { &*storage.unwrap() };
    let key = String::from(unsafe { CStr::from_ptr(key.unwrap()) }.to_str().unwrap());

    storage.ints[&key]
}

/// Get the int-valued parameter with name `key` for `storage`.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_keyval_get_string(
    storage: Option<*const KeyVal>,
    key: Option<*const c_char>,
) -> *const c_char {
    let storage = unsafe { &*storage.unwrap() };
    let key = String::from(unsafe { CStr::from_ptr(key.unwrap()) }.to_str().unwrap());

    storage.strings[&key].as_ptr()
}

/// Return a pointer to newly-created `pineappl_storage` object. The parameter method determines the
/// storage layout.
#[no_mangle]
#[must_use]
pub extern "C" fn pineappl_keyval_new(method: Option<*const c_char>) -> *mut KeyVal {
    Box::into_raw(Box::new(KeyVal {
        _method: String::from(unsafe { CStr::from_ptr(method.unwrap()) }.to_str().unwrap()),
        ..KeyVal::default()
    }))
}

/// Set the double-valued parameter with name `key` to `value` for `storage`.
#[no_mangle]
pub extern "C" fn pineappl_keyval_set_bool(
    storage: Option<*mut KeyVal>,
    key: Option<*const c_char>,
    value: bool,
) {
    let storage = unsafe { &mut *storage.unwrap() };
    let key = String::from(unsafe { CStr::from_ptr(key.unwrap()) }.to_str().unwrap());

    storage.bools.insert(key, value);
}

/// Set the double-valued parameter with name `key` to `value` for `storage`.
#[no_mangle]
pub extern "C" fn pineappl_keyval_set_double(
    storage: Option<*mut KeyVal>,
    key: Option<*const c_char>,
    value: f64,
) {
    let storage = unsafe { &mut *storage.unwrap() };
    let key = String::from(unsafe { CStr::from_ptr(key.unwrap()) }.to_str().unwrap());

    storage.doubles.insert(key, value);
}

/// Set the int-valued parameter with name `key` to `value` for `storage`.
#[no_mangle]
pub extern "C" fn pineappl_keyval_set_int(
    storage: Option<*mut KeyVal>,
    key: Option<*const c_char>,
    value: i32,
) {
    let storage = unsafe { &mut *storage.unwrap() };
    let key = String::from(unsafe { CStr::from_ptr(key.unwrap()) }.to_str().unwrap());

    storage.ints.insert(key, value);
}

/// Set the string-valued parameter with name `key` to `value` for `storage`.
#[no_mangle]
pub extern "C" fn pineappl_keyval_set_string(
    storage: Option<*mut KeyVal>,
    key: Option<*const c_char>,
    value: Option<*const c_char>,
) {
    let storage = unsafe { &mut *storage.unwrap() };
    let key = String::from(unsafe { CStr::from_ptr(key.unwrap()) }.to_str().unwrap());
    let value = CString::from(unsafe { CStr::from_ptr(value.unwrap()) });

    storage.strings.insert(key, value);
}
