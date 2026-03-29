#[allow(clippy::cast_possible_truncation)]
#[allow(clippy::cast_sign_loss)]
pub fn usize_from_f64(x: f64) -> usize {
    x.max(0.0) as usize
}

pub fn f64_from_usize(x: usize) -> f64 {
    f64::from(u32::try_from(x).unwrap())
}
