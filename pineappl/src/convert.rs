#[expect(
    clippy::cast_possible_truncation,
    reason = "we want that truncation to happen"
)]
#[expect(clippy::cast_sign_loss, reason = "we want to get rid of the sign")]
pub(crate) const fn usize_from_f64(x: f64) -> usize {
    x.max(0.0) as usize
}

pub(crate) fn f64_from_usize(x: usize) -> f64 {
    f64::from(u32::try_from(x).unwrap())
}
