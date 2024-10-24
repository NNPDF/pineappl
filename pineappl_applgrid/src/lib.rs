//! TODO

#![allow(clippy::missing_errors_doc)]
#![allow(clippy::missing_safety_doc)]
// doesn't work with cxx
#![allow(clippy::multiple_inherent_impl)]
// doesn't work with cxx
#![allow(clippy::multiple_unsafe_ops_per_block)]
#![allow(clippy::must_use_candidate)]
#![allow(clippy::too_many_arguments)]
#![allow(missing_docs)]

use lhapdf::Pdf;
use std::mem;
use std::pin::Pin;
use std::ptr;
use std::slice;
use std::sync::{Mutex, OnceLock};

static MUTEX: OnceLock<Mutex<()>> = OnceLock::new();

#[cxx::bridge]
pub mod ffi {
    #[repr(u32)]
    enum grid_CALCULATION {
        STANDARD = 0,
        AMCATNLO = 1,
        SHERPA = 2,
    }

    unsafe extern "C++" {
        // this header is needed to make the enum a member of the class `appl::grid`
        include!("pineappl_applgrid/src/helpers.hpp");

        type grid_CALCULATION;

        // TODO: `std::ffi::c_void` isn't support by cxx, see:
        // https://github.com/dtolnay/cxx/issues/1049, though there is a PR:
        // https://github.com/dtolnay/cxx/pull/1204
        type c_void;
    }

    #[namespace = "appl"]
    unsafe extern "C++" {
        include!("appl_grid/appl_grid.h");

        type grid;

        unsafe fn add_igrid(self: Pin<&mut Self>, _: i32, _: i32, _: *mut igrid);
        fn calculation(&self) -> grid_CALCULATION;
        fn genpdf(&self, _: i32, _: bool) -> *const appl_pdf;
        fn getApplyCorrection(&self, _: u32) -> bool;
        fn getApplyCorrections(&self) -> bool;
        fn getDynamicScale(&self) -> f64;
        fn getNormalised(&self) -> bool;
        fn include_photon(self: Pin<&mut Self>, _: bool);
        fn isDIS(&self) -> bool;
        fn leadingOrder(&self) -> i32;
        fn nloops(&self) -> i32;
        fn Nobs_internal(&self) -> i32;
        fn obslow_internal(&self, _: i32) -> f64;
        fn run(self: Pin<&mut Self>) -> &mut f64;
        fn weightgrid(&self, _: i32, _: i32) -> *const igrid;
        fn Write(self: Pin<&mut Self>, _: &CxxString, _: &CxxString, _: &CxxString);
    }

    #[namespace = "appl"]
    unsafe extern "C++" {
        include!("appl_igrid.h");

        type igrid;

        fn getQ2(&self, _: i32) -> f64;
        fn getx1(&self, _: i32) -> f64;
        fn getx2(&self, _: i32) -> f64;
        fn Ntau(&self) -> i32;
        fn Ny1(&self) -> i32;
        fn Ny2(&self) -> i32;
        fn setlimits(self: Pin<&mut Self>);
        fn weightgrid(&self, _: i32) -> *const SparseMatrix3d;
    }

    #[namespace = "appl"]
    unsafe extern "C++" {
        include!("appl_grid/appl_pdf.h");

        type appl_pdf;

        fn Nproc(&self) -> i32;
        unsafe fn evaluate(&self, _: *const f64, _: *const f64, _: *mut f64);
    }

    unsafe extern "C++" {
        type SparseMatrix3d;

        fn trim(self: Pin<&mut Self>);
    }

    unsafe extern "C++" {
        include!("appl_grid/lumi_pdf.h");

        type lumi_pdf;
    }

    unsafe extern "C++" {
        include!("pineappl_applgrid/src/applgrid.hpp");

        fn make_grid(_: &str) -> Result<UniquePtr<grid>>;
        fn make_igrid(
            _: i32,
            _: f64,
            _: f64,
            _: i32,
            _: i32,
            _: f64,
            _: f64,
            _: i32,
            _: &str,
            _: &str,
            _: i32,
            _: bool,
        ) -> UniquePtr<igrid>;
        fn make_empty_grid(_: &[f64], _: &str, _: i32, _: i32, _: &str, _: &str)
            -> UniquePtr<grid>;
        fn make_lumi_pdf(_: &str, _: &[i32]) -> UniquePtr<lumi_pdf>;

        fn grid_combine(_: &grid) -> Vec<i32>;

        unsafe fn grid_convolve_with_one(
            _: Pin<&mut grid>,
            _: unsafe fn(&f64, &f64, *mut f64, *mut c_void),
            _: unsafe fn(&f64, *mut c_void) -> f64,
            _: *mut c_void,
            _: i32,
            _: f64,
            _: f64,
            _: f64,
        ) -> Vec<f64>;

        fn sparse_matrix_get(_: &SparseMatrix3d, _: i32, _: i32, _: i32) -> f64;
        fn sparse_matrix_set(_: Pin<&mut SparseMatrix3d>, _: i32, _: i32, _: i32, _: f64);

        // TODO: class member functions aren't supported yet by cxx, see
        // https://github.com/dtolnay/cxx/issues/447
        fn weightfun(_: f64) -> f64;

        fn igrid_m_reweight(_: &igrid) -> bool;
        fn igrid_weightgrid(_: Pin<&mut igrid>, _: usize) -> Pin<&mut SparseMatrix3d>;
    }
}

pub fn grid_convolve_with_one(
    grid: Pin<&mut ffi::grid>,
    pdf: &mut Pdf,
    nloops: i32,
    rscale: f64,
    fscale: f64,
    escale: f64,
) -> Vec<f64> {
    let xfx = |x: &f64, q: &f64, results: *mut f64, pdf: *mut ffi::c_void| {
        let pdf = unsafe { &mut *pdf.cast::<Pdf>() };
        let results = unsafe { slice::from_raw_parts_mut(results, 14) };
        for (pid, result) in [-6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6, 22]
            .into_iter()
            .zip(results.iter_mut())
        {
            // some grids have x-nodes at slightly larger values than `1.0`; in that cases these
            // are numerical problems which we 'fix' by evaluating at exactly `1.0` instead
            *result = pdf.xfx_q2(pid, x.min(1.0), *q * *q);
        }
    };

    let alphas = |q: &f64, pdf: *mut ffi::c_void| -> f64 {
        let pdf = unsafe { &mut *pdf.cast::<Pdf>() };
        pdf.alphas_q2(*q * *q)
    };

    let lock = MUTEX
        .get_or_init(|| Mutex::new(()))
        .lock()
        // UNWRAP: if this fails there's an unexpected bug somewhere
        .unwrap_or_else(|_| unreachable!());

    let results = unsafe {
        ffi::grid_convolve_with_one(
            grid,
            xfx,
            alphas,
            ptr::from_mut(pdf).cast::<ffi::c_void>(),
            nloops,
            rscale,
            fscale,
            escale,
        )
    };

    mem::drop(lock);

    results
}
