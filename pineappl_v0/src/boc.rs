//! Module containing structures for the 3 dimensions of a [`Grid`]: bins, [`Order`] and channels
//! (`boc`).

use serde::Deserialize;

/// Coupling powers for each grid.
#[derive(Deserialize)]
pub struct Order {
    /// Exponent of the strong coupling.
    pub alphas: u32,
    /// Exponent of the electromagnetic coupling.
    pub alpha: u32,
    /// Exponent of the logarithm of the scale factor of the renomalization scale.
    pub logxir: u32,
    /// Exponent of the logarithm of the scale factor of the factorization scale.
    pub logxif: u32,
}

/// This structure represents a channel. Each channel consists of a tuple containing in the
/// following order, the particle ID of the first incoming parton, then the particle ID of the
/// second parton, and finally a numerical factor that will multiply the result for this specific
/// combination.
#[derive(Deserialize)]
pub struct Channel {
    entry: Vec<(i32, i32, f64)>,
}

impl Channel {
    /// Returns a tuple representation of this entry.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use pineappl_v0::channel;
    /// use pineappl_v0::boc::Channel;
    ///
    /// let entry = channel![4, 4, 1.0; 2, 2, 1.0];
    ///
    /// assert_eq!(entry.entry(), [(2, 2, 1.0), (4, 4, 1.0)]);
    /// ```
    #[must_use]
    pub fn entry(&self) -> &[(i32, i32, f64)] {
        &self.entry
    }
}
