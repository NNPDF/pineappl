//! Module that contains helpers for binning observables

use super::convert::{f64_from_usize, usize_from_f64};
use float_cmp::approx_eq;
use itertools::izip;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::f64;
use std::ops::Range;
use std::str::FromStr;
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

impl FromStr for BinRemapper {
    type Err = ParseBinRemapperError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let remaps: Result<Vec<Vec<Vec<_>>>, Self::Err> = s
            .split(';')
            .map(|string| {
                string
                    .split('|')
                    .map(|string| {
                        string
                            .split_once(':')
                            .map_or(Ok(string), |(lhs, rhs)| {
                                match (lhs.trim().parse::<usize>(), rhs.trim().parse::<usize>()) {
                                    (Err(lhs), Err(rhs)) => Err(ParseBinRemapperError::Error(format!(
                                        "unable to parse 'N:M' syntax from: '{string}' (N: '{lhs}', M: '{rhs}')"
                                    ))),
                                    // skip :N specification
                                    (Err(_), Ok(_)) => Ok(lhs),
                                    // skip N: specification
                                    (Ok(_), Err(_)) => Ok(rhs),
                                    // skip N:M specification
                                    (Ok(_), Ok(_)) => Ok(""),
                                }
                            })?
                            .split(',')
                            .filter_map(|string| {
                                let string = string.trim();
                                if string.is_empty() {
                                    None
                                } else {
                                    Some(string.parse::<f64>().map_err(|err| {
                                        ParseBinRemapperError::Error(format!(
                                            "unable to parse limit '{string}': '{err}')"
                                        ))
                                    }))
                                }
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect();
        let mut remaps = remaps?;

        if let Some(first) = remaps.first() {
            if first.len() != 1 {
                return Err(ParseBinRemapperError::Error(
                    "'|' syntax not meaningful for first dimension".to_owned(),
                ));
            }
        }

        // go over `remaps` again, and repeat previous entries as requested with the `|` syntax
        for vec in &mut remaps {
            for i in 1..vec.len() {
                if vec[i].is_empty() {
                    if vec[i - 1].is_empty() {
                        return Err(ParseBinRemapperError::Error(
                            "empty repetition with '|'".to_owned(),
                        ));
                    }

                    vec[i] = vec[i - 1].clone();
                }
            }
        }

        // go over `remaps` again, this time remove bin as requested with the `:N` or `N:` syntax
        for (vec, string) in remaps.iter_mut().zip(s.split(';')) {
            for (vec, string) in vec.iter_mut().zip(string.split('|')) {
                let (lhs, rhs) = {
                    if let Some((lhs, rhs)) = string.split_once(':') {
                        (lhs.parse::<usize>(), rhs.parse::<usize>())
                    } else {
                        // there's no colon
                        continue;
                    }
                };

                if let Ok(num) = rhs {
                    vec.truncate(vec.len() - num);
                }

                if let Ok(num) = lhs {
                    vec.drain(0..num);
                }

                if vec.len() <= 1 {
                    return Err(ParseBinRemapperError::Error(
                        "no limits due to ':' syntax".to_owned(),
                    ));
                }
            }
        }

        let dimensions = remaps.len();
        let mut normalizations = Vec::new();
        let mut limits = Vec::new();
        let mut buffer = Vec::with_capacity(dimensions);
        let mut pipe_indices = vec![0; dimensions];
        let mut last_indices = vec![0; dimensions];

        'looop: for indices in remaps
            .iter()
            .map(|vec| 0..vec.iter().map(|vec| vec.len() - 1).max().unwrap())
            .multi_cartesian_product()
        {
            // calculate `pipe_indices`, which stores the indices for the second dimension of `remaps`
            for d in 0..dimensions - 1 {
                if indices[d] > last_indices[d] {
                    for dp in d + 1..dimensions {
                        if remaps[dp].len() != 1 {
                            pipe_indices[dp] += 1;
                        }
                    }
                }
            }

            last_indices.clone_from(&indices);

            let mut normalization = 1.0;

            for (remap, &pipe_index, &i) in izip!(&remaps, &pipe_indices, &indices) {
                if let Some(r) = remap.get(pipe_index) {
                    if r.len() <= (i + 1) {
                        buffer.clear();

                        // this index doesn't exist
                        continue 'looop;
                    }

                    let left = r[i];
                    let right = r[i + 1];

                    buffer.push((left, right));
                    normalization *= right - left;
                } else {
                    return Err(ParseBinRemapperError::Error(
                        "missing '|' specification: number of variants too small".to_owned(),
                    ));
                }
            }

            limits.append(&mut buffer);
            normalizations.push(normalization);
        }

        Self::new(normalizations, limits)
            .map_err(|err| ParseBinRemapperError::BinRemapperNewError { source: err })
    }
}

impl<'a> BinInfo<'a> {
    /// Constructor.
    #[must_use]
    pub const fn new(limits: &'a BinLimits, remapper: Option<&'a BinRemapper>) -> Self {
        Self { limits, remapper }
    }

    /// Return the bin limits for the bin with index `bin`.
    #[must_use]
    pub fn bin_limits(&self, bin: usize) -> Vec<(f64, f64)> {
        // TODO: make return type a Cow
        self.remapper.map_or_else(
            || {
                let limits = &self.limits.limits()[bin..=bin + 1];
                vec![(limits[0], limits[1])]
            },
            |remapper| {
                let dim = remapper.dimensions();
                remapper.limits()[bin * dim..(bin + 1) * dim].to_vec()
            },
        )
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

    /// Return the index of the bin corresponding to `limits`. If no bin is found `None` is
    /// returned.
    #[must_use]
    pub fn find_bin(&self, limits: &[(f64, f64)]) -> Option<usize> {
        (0..self.bins())
            .map(|bin| self.bin_limits(bin))
            .position(|lim| lim == limits)
    }

    /// Returns all left-limits for the specified dimension. If the dimension does not exist, an
    /// empty vector is returned.
    #[must_use]
    pub fn left(&self, dimension: usize) -> Vec<f64> {
        if dimension >= self.dimensions() {
            vec![]
        } else {
            self.remapper.map_or_else(
                || {
                    self.limits
                        .limits()
                        .iter()
                        .take(self.bins())
                        .copied()
                        .collect()
                },
                |remapper| {
                    remapper
                        .limits()
                        .iter()
                        .skip(dimension)
                        .step_by(self.dimensions())
                        .take(self.bins())
                        .map(|tuple| tuple.0)
                        .collect()
                },
            )
        }
    }

    /// Returns all right-limits for the specified dimension. If the dimension does not exist, an
    /// empty vector is returned.
    #[must_use]
    pub fn right(&self, dimension: usize) -> Vec<f64> {
        if dimension >= self.dimensions() {
            vec![]
        } else {
            self.remapper.map_or_else(
                || {
                    self.limits
                        .limits()
                        .iter()
                        .skip(1)
                        .take(self.bins())
                        .copied()
                        .collect()
                },
                |remapper| {
                    remapper
                        .limits()
                        .iter()
                        .skip(dimension)
                        .step_by(self.dimensions())
                        .take(self.bins())
                        .map(|tuple| tuple.1)
                        .collect()
                },
            )
        }
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

    /// Returns a vector of half-open intervals that show how multi-dimensional bins can be
    /// efficiently sliced into one-dimensional histograms.
    #[must_use]
    pub fn slices(&self) -> Vec<(usize, usize)> {
        // TODO: convert this to Vec<Range<usize>>
        self.remapper
            .map_or_else(|| vec![(0, self.limits.bins())], BinRemapper::slices)
    }
}

impl PartialEq<BinInfo<'_>> for BinInfo<'_> {
    fn eq(&self, other: &BinInfo) -> bool {
        (self.limits() == other.limits()) && (self.normalizations() == other.normalizations())
    }
}

impl BinRemapper {
    /// Create a new `BinRemapper` object with the specified number of bins and dimensions and
    /// limits.
    ///
    /// # Errors
    ///
    /// Returns an error if the length of `limits` is not a multiple of the length of
    /// `normalizations`, or if the limits of at least two bins overlap.
    pub fn new(
        normalizations: Vec<f64>,
        limits: Vec<(f64, f64)>,
    ) -> Result<Self, BinRemapperNewError> {
        if limits.len() % normalizations.len() == 0 {
            let dimensions = limits.len() / normalizations.len();
            let mut overlaps = Vec::new();

            for (i, bin_i) in limits.chunks_exact(dimensions).enumerate() {
                for (j, bin_j) in limits.chunks_exact(dimensions).enumerate().skip(i + 1) {
                    if bin_i.iter().zip(bin_j).all(|((l1, r1), (l2, r2))| {
                        ((l2 >= l1) && (l2 < r1)) || ((l1 >= l2) && (l1 < r2))
                    }) {
                        overlaps.push(j);
                    }
                }
            }

            overlaps.sort_unstable();
            overlaps.dedup();

            if overlaps.is_empty() {
                Ok(Self {
                    normalizations,
                    limits,
                })
            } else {
                Err(BinRemapperNewError::OverlappingBins { overlaps })
            }
        } else {
            Err(BinRemapperNewError::DimensionUnknown {
                normalizations_len: normalizations.len(),
                limits_len: limits.len(),
            })
        }
    }

    /// Return the number of bins.
    #[must_use]
    pub fn bins(&self) -> usize {
        self.normalizations.len()
    }

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

    /// Merges the bins for the corresponding range together in a single one.
    ///
    /// # Errors
    ///
    /// When `range` refers to non-consecutive bins, an error is returned.
    pub fn merge_bins(&mut self, range: Range<usize>) -> Result<(), MergeBinError> {
        if self
            .slices()
            .iter()
            .any(|&(start, end)| (start <= range.start) && (range.end <= end))
        {
            for bin in range.start + 1..range.end {
                self.normalizations[range.start] += self.normalizations[bin];
            }

            let dim = self.dimensions();

            self.normalizations.drain(range.start + 1..range.end);
            self.limits[dim * (range.start + 1) - 1].1 = self.limits[dim * range.end - 1].1;
            self.limits.drain(dim * (range.start + 1)..dim * range.end);

            Ok(())
        } else {
            Err(MergeBinError::NonConsecutiveRange(range))
        }
    }

    /// Merge the `BinRemapper` of `other` into `self` on the right-hand-side.
    ///
    /// # Errors
    ///
    /// If the dimensions of both remappers are not the same an error is returned.
    pub fn merge(&mut self, other: &Self) -> Result<(), MergeBinError> {
        let lhs_dim = self.dimensions();
        let rhs_dim = other.dimensions();

        if lhs_dim != rhs_dim {
            return Err(MergeBinError::IncompatibleDimensions {
                lhs: lhs_dim,
                rhs: rhs_dim,
            });
        }

        // TODO: we shouldn't allow overlapping bins
        self.normalizations.extend_from_slice(&other.normalizations);
        self.limits.extend_from_slice(&other.limits);

        Ok(())
    }

    /// Return the normalization factors for all bins.
    #[must_use]
    pub fn normalizations(&self) -> &[f64] {
        &self.normalizations
    }

    /// Returns a vector of half-open intervals that show how multi-dimensional bins can be
    /// efficiently sliced into one-dimensional histograms.
    #[must_use]
    pub fn slices(&self) -> Vec<(usize, usize)> {
        if self.dimensions() == 1 {
            vec![(0, self.bins())]
        } else {
            self.limits()
                .iter()
                .enumerate()
                .filter_map(|(index, x)| {
                    ((index % self.dimensions()) != (self.dimensions() - 1)).then_some(x)
                })
                .collect::<Vec<_>>()
                .chunks_exact(self.dimensions() - 1)
                .enumerate()
                .dedup_by_with_count(|(_, x), (_, y)| x == y)
                .map(|(count, (index, _))| (index, index + count))
                .collect()
        }
    }

    /// Deletes all bins whose corresponding indices are in one of the ranges of `bins`.
    pub fn delete_bins(&mut self, bins: &[Range<usize>]) {
        let dim = self.dimensions();

        for range in bins.iter().cloned().rev() {
            self.normalizations.drain(range);
        }

        for range in bins.iter().rev() {
            self.limits.drain((range.start * dim)..(range.end * dim));
        }
    }
}

impl PartialEq<Self> for BinRemapper {
    fn eq(&self, other: &Self) -> bool {
        (self.limits == other.limits) && (self.normalizations == other.normalizations)
    }
}

impl BinLimits {
    /// Constructor for `BinLimits`.
    ///
    /// # Panics
    ///
    /// TODO
    #[must_use]
    pub fn new(mut limits: Vec<f64>) -> Self {
        limits.sort_by(|left, right| left.partial_cmp(right).unwrap());

        if limits
            .iter()
            .zip(limits.iter().skip(1))
            .map(|(current, next)| next - current)
            .collect::<Vec<f64>>()
            .windows(2)
            .all(|val| approx_eq!(f64, val[0], val[1], ulps = 8))
        {
            Self(Limits::Equal {
                left: *limits.first().unwrap(),
                right: *limits.last().unwrap(),
                bins: limits.len() - 1,
            })
        } else {
            Self(Limits::Unequal { limits })
        }
    }

    /// Returns the number of bins.
    #[must_use]
    pub fn bins(&self) -> usize {
        match &self.0 {
            Limits::Equal { bins, .. } => *bins,
            Limits::Unequal { limits } => limits.len() - 1,
        }
    }

    /// Returns the bin index for observable `value`. If the value over- or underflows, the return
    /// value is `None`.
    ///
    /// # Panics
    ///
    /// TODO
    #[must_use]
    pub fn index(&self, value: f64) -> Option<usize> {
        match &self.0 {
            Limits::Equal { left, right, bins } => {
                if value < *left || value >= *right {
                    None
                } else {
                    Some(usize_from_f64(
                        (value - left) / (right - left) * f64_from_usize(*bins),
                    ))
                }
            }
            Limits::Unequal { limits } => {
                match limits.binary_search_by(|left| left.partial_cmp(&value).unwrap()) {
                    Err(0) => None,
                    Err(index) if index == limits.len() => None,
                    Ok(index) if index == (limits.len() - 1) => None,
                    Ok(index) => Some(index),
                    Err(index) => Some(index - 1),
                }
            }
        }
    }

    /// Returns the left-most bin limit
    ///
    /// # Panics
    ///
    /// TODO
    #[must_use]
    pub fn left(&self) -> f64 {
        match &self.0 {
            Limits::Unequal { limits } => *limits.first().unwrap(),
            Limits::Equal { left, .. } => *left,
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

    /// Merges the bins for the corresponding range together in a single one.
    ///
    /// # Errors
    ///
    /// When `bins` contains any indices that do not correspond to bins this method returns an
    /// error.
    pub fn merge_bins(&mut self, range: Range<usize>) -> Result<(), MergeBinError> {
        if range.end > self.bins() {
            return Err(MergeBinError::InvalidRange {
                range,
                bins: self.bins(),
            });
        }

        let mut new_limits = self.limits();
        new_limits.drain(range.start + 1..range.end);
        *self = Self::new(new_limits);

        Ok(())
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

    /// Merge the limits of `other` into `self` on the right-hand-side. If both limits are
    /// non-consecutive, an error is returned.
    ///
    /// # Errors
    ///
    /// If the right-most limit of `self` is different from the left-most limit of `other`, the
    /// bins are non-consecutive and an error is returned.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn merge(&mut self, other: &Self) -> Result<(), MergeBinError> {
        if !approx_eq!(f64, self.right(), other.left(), ulps = 8) {
            return Err(MergeBinError::NonConsecutiveBins {
                lhs: self.right(),
                rhs: other.left(),
            });
        }

        let mut limits = self.limits();
        let add_limits = other.limits();

        // average over the shared limit
        *limits.last_mut().unwrap() =
            0.5 * (*limits.last().unwrap() + *add_limits.first().unwrap());
        // add the new limits
        limits.extend_from_slice(&add_limits[1..]);

        // use the constructor to get a valid state
        *self = Self::new(limits);

        Ok(())
    }

    /// Returns the right-most bin limit
    ///
    /// # Panics
    ///
    /// TODO
    #[must_use]
    pub fn right(&self) -> f64 {
        match &self.0 {
            Limits::Unequal { limits } => *limits.last().unwrap(),
            Limits::Equal { right, .. } => *right,
        }
    }

    /// Delete `bins` number of bins from the start.
    pub fn delete_bins_left(&mut self, bins: usize) {
        let mut limits = self.limits();
        limits.drain(..bins);
        *self = Self::new(limits);
    }

    /// Delete `bins` number of bins from the end.
    pub fn delete_bins_right(&mut self, bins: usize) {
        let mut limits = self.limits();
        limits.drain((limits.len() - bins)..);
        *self = Self::new(limits);
    }
}
