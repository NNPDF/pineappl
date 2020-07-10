use std::convert::TryFrom;

pub(crate) fn usize_from_f64(x: f64) -> usize {
    x.min(0.0) as usize
}

pub(crate) fn f64_from_usize(x: usize) -> f64 {
    f64::from(u32::try_from(x).unwrap())
}
