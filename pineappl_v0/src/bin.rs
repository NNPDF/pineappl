//! Module that contains helpers for binning observables

use super::convert::f64_from_usize;
use serde::Deserialize;
use std::f64;

#[derive(Deserialize)]
enum Limits {
    Equal { left: f64, right: f64, bins: usize },
    Unequal { limits: Vec<f64> },
}

/// Structure representing bin limits.
#[derive(Deserialize)]
pub struct BinLimits(Limits);

/// Structure for remapping bin limits.
#[derive(Deserialize)]
pub struct BinRemapper {
    normalizations: Vec<f64>,
    limits: Vec<(f64, f64)>,
}

/// Captures all information about the bins in a grid.
pub struct BinInfo<'a> {
    limits: &'a BinLimits,
    remapper: Option<&'a BinRemapper>,
}

impl<'a> BinInfo<'a> {
    /// Constructor.
    #[must_use]
    pub const fn new(limits: &'a BinLimits, remapper: Option<&'a BinRemapper>) -> Self {
        Self { limits, remapper }
    }

    /// Returns the number of bins.
    #[must_use]
    pub const fn bins(&self) -> usize {
        self.limits.bins()
    }

    /// Returns the number of dimensions.
    #[must_use]
    pub fn dimensions(&self) -> usize {
        self.remapper.map_or(1, BinRemapper::dimensions)
    }

    /// For each bin return a vector of `(left, right)` limits for each dimension.
    #[must_use]
    pub fn limits(&self) -> Vec<Vec<(f64, f64)>> {
        self.remapper.map_or_else(
            || {
                self.limits
                    .limits()
                    .windows(2)
                    .map(|window| vec![(window[0], window[1])])
                    .collect()
            },
            |remapper| {
                remapper
                    .limits()
                    .to_vec()
                    .chunks_exact(self.dimensions())
                    .map(<[(f64, f64)]>::to_vec)
                    .collect()
            },
        )
    }

    /// Returns all normalization factors.
    #[must_use]
    pub fn normalizations(&self) -> Vec<f64> {
        self.remapper.map_or_else(
            || self.limits.bin_sizes(),
            |remapper| remapper.normalizations().to_vec(),
        )
    }
}

impl BinRemapper {
    /// Return the number of dimensions.
    #[must_use]
    pub const fn dimensions(&self) -> usize {
        self.limits.len() / self.normalizations.len()
    }

    /// Return tuples of left and right bin limits for all dimensions and all bins.
    #[must_use]
    pub fn limits(&self) -> &[(f64, f64)] {
        &self.limits
    }

    /// Return the normalization factors for all bins.
    #[must_use]
    pub fn normalizations(&self) -> &[f64] {
        &self.normalizations
    }
}

impl BinLimits {
    /// Returns the number of bins.
    #[must_use]
    pub const fn bins(&self) -> usize {
        match &self.0 {
            Limits::Equal { bins, .. } => *bins,
            Limits::Unequal { limits } => limits.len() - 1,
        }
    }

    /// Returns the limits in a `Vec`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use pineappl_v0::bin::BinLimits;
    ///
    /// // example with equally sized bins
    /// let equal_bins = BinLimits::new(vec![0.25, 0.5, 0.75, 1.0]);
    /// assert_eq!(equal_bins.limits(), vec![0.25, 0.5, 0.75, 1.0]);
    ///
    /// // example with unequally sized bins
    /// let unequal_bins = BinLimits::new(vec![0.125, 0.25, 1.0, 1.5]);
    /// assert_eq!(unequal_bins.limits(), vec![0.125, 0.25, 1.0, 1.5]);
    /// ```
    #[must_use]
    pub fn limits(&self) -> Vec<f64> {
        match &self.0 {
            Limits::Equal { left, right, bins } => (0..=*bins)
                .map(|b| (*right - *left).mul_add(f64_from_usize(b) / f64_from_usize(*bins), *left))
                .collect(),
            Limits::Unequal { limits } => limits.clone(),
        }
    }

    /// Returns the size for each bin.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use pineappl_v0::bin::BinLimits;
    ///
    /// // example with equally sized bins
    /// let equal_bins = BinLimits::new(vec![0.25, 0.5, 0.75, 1.0]);
    /// assert_eq!(equal_bins.bin_sizes(), vec![0.25, 0.25, 0.25]);
    ///
    /// // example with unequally sized bins
    /// let unequal_bins = BinLimits::new(vec![0.125, 0.25, 1.0, 1.5]);
    /// assert_eq!(unequal_bins.bin_sizes(), vec![0.125, 0.75, 0.5]);
    /// ```
    #[must_use]
    pub fn bin_sizes(&self) -> Vec<f64> {
        match &self.0 {
            Limits::Equal { left, right, bins } => {
                vec![(*right - *left) / f64_from_usize(*bins); *bins]
            }
            Limits::Unequal { limits } => limits.windows(2).map(|x| x[1] - x[0]).collect(),
        }
    }
}
