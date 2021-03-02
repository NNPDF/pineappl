//! Module that contains helpers for binning observables

use super::convert::{f64_from_usize, usize_from_f64};
use float_cmp::approx_eq;
use serde::{Deserialize, Serialize};
use std::f64;
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

    /// Returns all left-limits for the specified dimension. If the dimension does not exist, an
    /// empty vector is returned.
    #[must_use]
    pub fn left(&self, dimension: usize) -> Vec<f64> {
        if dimension >= self.dimensions() {
            vec![]
        } else if let Some(remapper) = self.remapper {
            remapper
                .limits()
                .iter()
                .skip(dimension)
                .step_by(self.dimensions())
                .take(self.bins())
                .map(|tuple| tuple.0)
                .collect()
        } else {
            self.limits
                .limits()
                .iter()
                .skip(0)
                .take(self.bins())
                .copied()
                .collect()
        }
    }

    /// Returns all right-limits for the specified dimension. If the dimension does not exist, an
    /// empty vector is returned.
    #[must_use]
    pub fn right(&self, dimension: usize) -> Vec<f64> {
        if dimension >= self.dimensions() {
            vec![]
        } else if let Some(remapper) = self.remapper {
            remapper
                .limits()
                .iter()
                .skip(dimension)
                .step_by(self.dimensions())
                .take(self.bins())
                .map(|tuple| tuple.1)
                .collect()
        } else {
            self.limits
                .limits()
                .iter()
                .skip(1)
                .take(self.bins())
                .copied()
                .collect()
        }
    }

    /// Returns all normalization factors.
    #[must_use]
    pub fn normalizations(&self) -> Vec<f64> {
        if let Some(remapper) = self.remapper {
            remapper.normalizations().to_vec()
        } else {
            self.limits.bin_sizes()
        }
    }
}

impl PartialEq<BinInfo<'_>> for BinInfo<'_> {
    fn eq(&self, other: &BinInfo) -> bool {
        (self.limits == other.limits) && (self.remapper == other.remapper)
    }
}

impl BinRemapper {
    /// Create a new `BinRemapper` object with the specified number of bins and dimensions and
    /// limits.
    ///
    /// # Errors
    ///
    /// Returns an error if the length of `limits` is not a multiple of the length of
    /// `normalizations`.
    pub fn new(
        normalizations: Vec<f64>,
        limits: Vec<(f64, f64)>,
    ) -> Result<Self, BinRemapperNewError> {
        if limits.len() % normalizations.len() == 0 {
            Ok(Self {
                normalizations,
                limits,
            })
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

    /// Return the normalization factors for all bins.
    #[must_use]
    pub fn normalizations(&self) -> &[f64] {
        &self.normalizations
    }
}

impl PartialEq<BinRemapper> for BinRemapper {
    fn eq(&self, other: &Self) -> bool {
        (self.limits == other.limits) && (self.normalizations == other.normalizations)
    }
}

impl BinLimits {
    /// Constructor for `BinLimits`.
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
    /// use pineappl::bin::BinLimits;
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
    /// use pineappl::bin::BinLimits;
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
    #[must_use]
    pub fn right(&self) -> f64 {
        match &self.0 {
            Limits::Unequal { limits } => *limits.last().unwrap(),
            Limits::Equal { right, .. } => *right,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn bin_limits_merge() {
        let mut limits = BinLimits::new(vec![0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0]);

        // right merge
        assert!(limits
            .merge(&BinLimits::new(vec![
                1.0,
                1.0 + 1.0 / 3.0,
                1.0 + 2.0 / 3.0,
                2.0
            ]))
            .is_ok());

        assert_eq!(limits.left(), 0.0);
        assert_eq!(limits.right(), 2.0);
        assert_eq!(limits.bins(), 6);

        let non_consecutive_bins = BinLimits::new(vec![3.0, 4.0]);

        assert!(limits.merge(&non_consecutive_bins).is_err());

        assert_eq!(limits.left(), 0.0);
        assert_eq!(limits.right(), 2.0);
        assert_eq!(limits.bins(), 6);

        // left merge
        assert!(limits
            .merge(&BinLimits::new(vec![
                -1.0,
                -1.0 + 1.0 / 3.0,
                -1.0 + 2.0 / 3.0,
                0.0
            ]))
            .is_err());

        assert_eq!(limits.left(), 0.0);
        assert_eq!(limits.right(), 2.0);
        assert_eq!(limits.bins(), 6);
    }

    #[test]
    fn bin_info_without_remapper() {
        let limits = BinLimits::new(vec![0.0, 0.125, 0.25, 0.375, 0.5]);
        let info = BinInfo::new(&limits, None);

        assert_eq!(info.bins(), 4);
        assert_eq!(info.dimensions(), 1);
        assert_eq!(info.left(0), vec![0.0, 0.125, 0.25, 0.375]);
        assert_eq!(info.right(0), vec![0.125, 0.25, 0.375, 0.5]);
        assert_eq!(info.normalizations(), vec![0.125; 4]);

        assert_eq!(info.left(1), vec![]);
        assert_eq!(info.right(1), vec![]);
    }

    #[test]
    fn bin_info_with_remapper() {
        let limits = BinLimits::new(vec![0.0, 0.125, 0.25, 0.375, 0.5]);
        let remapper = BinRemapper::new(
            vec![1.0; 4],
            vec![
                (0.0, 0.5),
                (0.25, 0.75),
                (1.0, 2.0),
                (0.5, 1.0),
                (0.75, 1.0),
                (2.0, 5.0),
                (1.0, 2.0),
                (1.75, 2.0),
                (5.0, 5.5),
                (2.5, 3.0),
                (2.0, 2.5),
                (6.0, 8.0),
            ],
        )
        .unwrap();
        let info = BinInfo::new(&limits, Some(&remapper));

        assert_ne!(info, BinInfo::new(&limits, None));
        assert_eq!(info, BinInfo::new(&limits, Some(&remapper)));

        assert_eq!(info.bins(), 4);
        assert_eq!(info.dimensions(), 3);
        assert_eq!(info.left(0), vec![0.0, 0.5, 1.0, 2.5]);
        assert_eq!(info.left(1), vec![0.25, 0.75, 1.75, 2.0]);
        assert_eq!(info.left(2), vec![1.0, 2.0, 5.0, 6.0]);
        assert_eq!(info.right(0), vec![0.5, 1.0, 2.0, 3.0]);
        assert_eq!(info.right(1), vec![0.75, 1.0, 2.0, 2.5]);
        assert_eq!(info.right(2), vec![2.0, 5.0, 5.5, 8.0]);
        assert_eq!(info.normalizations(), vec![1.0; 4]);

        assert_eq!(info.left(3), vec![]);
        assert_eq!(info.right(3), vec![]);
    }

    #[test]
    fn bin_limits() {
        // first check BinLimits with exactly representable bin sizes
        let limits = BinLimits::new(vec![0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0]);

        assert_eq!(limits.bins(), 8);
        assert_eq!(limits.index(-0.1), None);
        assert_eq!(limits.index(0.1), Some(0));
        assert_eq!(limits.index(0.2), Some(1));
        assert_eq!(limits.index(0.3), Some(2));
        assert_eq!(limits.index(0.4), Some(3));
        assert_eq!(limits.index(0.55), Some(4));
        assert_eq!(limits.index(0.65), Some(5));
        assert_eq!(limits.index(0.8), Some(6));
        assert_eq!(limits.index(0.9), Some(7));
        assert_eq!(limits.index(1.1), None);

        // check bin limits that are equally sized, with values on the limits
        let limits = BinLimits::new(vec![0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0]);
        assert_eq!(limits.index(0.0), Some(0));
        assert_eq!(limits.index(0.125), Some(1));
        assert_eq!(limits.index(0.25), Some(2));
        assert_eq!(limits.index(0.375), Some(3));
        assert_eq!(limits.index(0.5), Some(4));
        assert_eq!(limits.index(0.625), Some(5));
        assert_eq!(limits.index(0.75), Some(6));
        assert_eq!(limits.index(0.875), Some(7));
        assert_eq!(limits.index(1.0), None);

        // now, check with bin sizes that are not exactly representable
        let limits = BinLimits::new(vec![0.0, 0.1, 0.2, 0.3, 0.4, 0.5]);

        assert_eq!(limits.bins(), 5);
        assert_eq!(limits.index(-1.0), None);
        assert_eq!(limits.index(0.05), Some(0));
        assert_eq!(limits.index(0.15), Some(1));
        assert_eq!(limits.index(0.25), Some(2));
        assert_eq!(limits.index(0.35), Some(3));
        assert_eq!(limits.index(0.45), Some(4));
        assert_eq!(limits.index(1.1), None);

        // check the special case of one bin
        let limits = BinLimits::new(vec![0.0, 1.0]);
        assert_eq!(limits.bins(), 1);
        assert_eq!(limits.index(-0.1), None);
        assert_eq!(limits.index(0.5), Some(0));
        assert_eq!(limits.index(1.1), None);

        // check bin limits that are unequally sized, with ascending bin sizes
        let limits = BinLimits::new(vec![0.0, 0.1, 0.3, 0.6, 1.0]);
        assert_eq!(limits.bins(), 4);
        assert_eq!(limits.index(-1.0), None);
        assert_eq!(limits.index(0.05), Some(0));
        assert_eq!(limits.index(0.2), Some(1));
        assert_eq!(limits.index(0.4), Some(2));
        assert_eq!(limits.index(0.9), Some(3));
        assert_eq!(limits.index(1.3), None);

        // check bin limits that are unequally sized, with values on the limits
        let limits = BinLimits::new(vec![0.0, 0.25, 0.75, 0.875, 1.0]);
        assert_eq!(limits.index(0.0), Some(0));
        assert_eq!(limits.index(0.25), Some(1));
        assert_eq!(limits.index(0.75), Some(2));
        assert_eq!(limits.index(0.875), Some(3));
        assert_eq!(limits.index(1.0), None);

        // check bin limits that are unequally sized, with descending bin sizes
        let limits = BinLimits::new(vec![0.0, 0.4, 0.7, 0.9, 1.0]);
        assert_eq!(limits.bins(), 4);
        assert_eq!(limits.index(-1.0), None);
        assert_eq!(limits.index(0.2), Some(0));
        assert_eq!(limits.index(0.5), Some(1));
        assert_eq!(limits.index(0.8), Some(2));
        assert_eq!(limits.index(0.95), Some(3));
        assert_eq!(limits.index(1.3), None);
    }

    #[test]
    fn bin_remapper() {
        let remapper = BinRemapper::new(
            vec![1.0; 4],
            vec![
                (0.0, 0.5),
                (0.25, 0.75),
                (0.5, 1.0),
                (0.75, 1.0),
                (1.0, 2.0),
                (1.75, 2.0),
                (2.5, 3.0),
                (2.0, 2.5),
            ],
        )
        .unwrap();

        assert_ne!(
            remapper,
            BinRemapper::new(
                vec![1.0; 4],
                vec![(0.0, 1.0), (1.0, 2.0), (2.0, 3.0), (4.0, 5.0)]
            )
            .unwrap()
        );

        assert!(matches!(
            BinRemapper::new(vec![1.0; 8], vec![(0.0, 1.0); 2]),
            Err(BinRemapperNewError::DimensionUnknown{normalizations_len, limits_len})
                if (normalizations_len == 8) && (limits_len == 2)
        ));

        assert_eq!(remapper.bins(), 4);
        assert_eq!(remapper.dimensions(), 2);
        assert_eq!(
            remapper.limits(),
            &[
                (0.0, 0.5),
                (0.25, 0.75),
                (0.5, 1.0),
                (0.75, 1.0),
                (1.0, 2.0),
                (1.75, 2.0),
                (2.5, 3.0),
                (2.0, 2.5)
            ]
        );
        assert_eq!(remapper.normalizations(), vec![1.0; 4]);
    }
}
