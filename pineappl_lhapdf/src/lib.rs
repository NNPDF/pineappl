#[macro_use]
extern crate cpp;

use std::ffi::{c_void, CString};

cpp! {{
    #include <LHAPDF/LHAPDF.h>
}}

pub struct Pdf {
    ptr: *mut c_void,
}

impl Pdf {
    pub fn new(setname: &str, member: i32) -> Self {
        let setname = CString::new(setname).unwrap();
        let setname_ptr = setname.as_ptr();

        Self {
            ptr: unsafe {
                cpp!([setname_ptr as "const char *",
                      member as "int"] -> *mut c_void as "LHAPDF::PDF*" {
                    return LHAPDF::mkPDF(setname_ptr, member);
                })
            },
        }
    }

    pub fn xfx_q2(&self, id: i32, x: f64, q2: f64) -> f64 {
        let self_ptr = self.ptr;

        unsafe {
            cpp!([self_ptr as "LHAPDF::PDF *",
                           id as "int",
                           x as "double",
                           q2 as "double"] -> f64 as "double" {
                return self_ptr->xfxQ2(id, x, q2);
            })
        }
    }

    pub fn alphas_q2(&self, q2: f64) -> f64 {
        let self_ptr = self.ptr;

        unsafe {
            cpp!([self_ptr as "LHAPDF::PDF *", q2 as "double"] -> f64 as "double" {
                return self_ptr->alphasQ2(q2);
            })
        }
    }
}

impl Drop for Pdf {
    fn drop(&mut self) {
        let self_ptr = self.ptr;

        unsafe {
            cpp!([self_ptr as "LHAPDF::PDF *"] -> () as "void" {
                delete self_ptr;
            });
        }
    }
}
