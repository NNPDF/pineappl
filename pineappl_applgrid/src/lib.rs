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
        include!("pineappl_applgrid/src/calculation.hpp");

        type grid_CALCULATION;
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
        fn grid_convolute(
            _: Pin<&mut grid>,
            _: &str,
            _: i32,
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
