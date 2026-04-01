//! Module that contains helpers for binning observables

use super::convert::f64_from_usize;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::f64;
use std::ops::Range;
use thiserror::Error;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
enum Limits {
    Equal { left: f64, right: f64, bins: usize },
    Unequal { limits: Vec<f64> },
}

/// Error type which is returned when two `BinLimits` objects are merged which are not
/// connected/non-consecutive.
#[derive(Debug, Error)]
pub enum MergeBinError {
    /// Returned when two `BinLimits` objects `a` and `b` were tried to be merged using
    /// `a.merge(b)`, but when the right-most limit of `a` does not match the left-most limit of
    /// `b`.
    #[error("can not merge bins which end at {lhs} with bins that start at {rhs}")]
    NonConsecutiveBins {
        /// right-most limit of the `BinLimits` object that is being merged into.
        lhs: f64,
        /// left-most limit of the `BinLimits` object that is being merged.
        rhs: f64,
    },

    /// Returned by [`BinRemapper::merge_bins`] whenever it can not merge bins.
    #[error("can not merge bins with indices {0:?}")]
    NonConsecutiveRange(Range<usize>),

    /// Returned by [`BinLimits::merge_bins`] whenever the range is outside the available bins.
    #[error("tried to merge bins with indices {range:?}, but there are only {bins} bins")]
    InvalidRange {
        /// Range given to [`BinLimits::merge_bins`].
        range: Range<usize>,
        /// Number of bins.
        bins: usize,
    },

    /// Returned by [`BinRemapper::merge`] whenever the dimensions of two `BinRemapper` are not the
    /// same.
    #[error("tried to merge bins with different dimensions {lhs} and {rhs}")]
    IncompatibleDimensions {
        /// Dimension of the bins of the first `BinRemapper`.
        lhs: usize,
        /// Dimension of the bins of the second `BinRemapper`.
        rhs: usize,
    },
}

/// Structure representing bin limits.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct BinLimits(Limits);

/// Error type that is returned by the constructor of `BinRemapper`.
#[derive(Debug, Error)]
pub enum BinRemapperNewError {
    /// Returned if the lengths of the normalization and limits vectors do not allow to determine a
    /// well-defined number of dimensions.
    #[error("could not determine the dimensions from a normalization vector with length {normalizations_len} and limits vector with length {limits_len}")]
    DimensionUnknown {
        /// Length of the normalization vector.
        normalizations_len: usize,
        /// Length of the limits vector.
        limits_len: usize,
    },
    /// Returned if bins overlap.
    #[error("the bin limits for the bins with indices {} overlap with other bins", overlaps.iter().map(ToString::to_string).join(","))]
    OverlappingBins {
        /// Indices of the bins that overlap with other bins.
        overlaps: Vec<usize>,
    },
}

/// Structure for remapping bin limits.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct BinRemapper {
    normalizations: Vec<f64>,
    limits: Vec<(f64, f64)>,
}

/// Captures all information about the bins in a grid.
#[derive(Debug)]
pub struct BinInfo<'a> {
    limits: &'a BinLimits,
    remapper: Option<&'a BinRemapper>,
}

/// Error type returned by [`BinRemapper::from_str`]
#[derive(Debug, Error)]
pub enum ParseBinRemapperError {
    /// An error that occured while parsing the string in [`BinRemapper::from_str`].
    #[error("{0}")]
    Error(String),
    /// An error that occured while constructing the remapper with [`BinRemapper::new`].
    #[error("{source}")]
    BinRemapperNewError {
        // TODO: enable #[backtrace] whenever the feature is stable
        /// The error returned by [`BinRemapper::new`].
        source: BinRemapperNewError,
    },
}

impl<'a> BinInfo<'a> {
    /// Constructor.
    #[must_use]
    pub const fn new(limits: &'a BinLimits, remapper: Option<&'a BinRemapper>) -> Self {
        Self { limits, remapper }
    }

    /// Returns the number of bins.
    #[must_use]
    pub fn bins(&self) -> usize {
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

impl PartialEq<BinInfo<'_>> for BinInfo<'_> {
    fn eq(&self, other: &BinInfo) -> bool {
        (self.limits() == other.limits()) && (self.normalizations() == other.normalizations())
    }
}

impl BinRemapper {
    /// Return the number of dimensions.
    #[must_use]
    pub fn dimensions(&self) -> usize {
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

impl PartialEq<Self> for BinRemapper {
    fn eq(&self, other: &Self) -> bool {
        (self.limits == other.limits) && (self.normalizations == other.normalizations)
    }
}

impl BinLimits {
    /// Returns the number of bins.
    #[must_use]
    pub fn bins(&self) -> usize {
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
