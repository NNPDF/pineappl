#[macro_use]
extern crate cpp;

use std::ffi::{c_void, CStr, CString};
use std::os::raw::c_char;

cpp! {{
    #include <LHAPDF/LHAPDF.h>
}}

pub fn available_pdf_sets() -> Vec<String> {
    let pdfs = unsafe {
        cpp!([] -> usize as "size_t" {
            return static_cast<unsigned> (LHAPDF::availablePDFSets().size());
        })
    };

    let mut pdf_sets: Vec<String> = Vec::with_capacity(pdfs);

    for i in 0..pdfs {
        let cstr_ptr = unsafe {
            cpp!([i as "size_t"] -> *const c_char as "const char *" {
                return LHAPDF::availablePDFSets().at(i).c_str();
            })
        };
        pdf_sets.push(
            unsafe { CStr::from_ptr(cstr_ptr) }
                .to_str()
                .unwrap()
                .to_string(),
        );
    }

    pdf_sets
}

pub struct Pdf {
    ptr: *mut c_void,
}

impl Pdf {
    pub fn with_setname_and_member(setname: &str, member: i32) -> Self {
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

    pub fn with_lhaid(lhaid: i32) -> Self {
        Self {
            ptr: unsafe {
                cpp!([lhaid as "int"] -> *mut c_void as "LHAPDF::PDF*" {
                    return LHAPDF::mkPDF(lhaid);
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
