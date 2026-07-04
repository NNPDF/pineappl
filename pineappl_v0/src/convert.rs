pub fn f64_from_usize(x: usize) -> f64 {
    f64::from(u32::try_from(x).unwrap())
}
