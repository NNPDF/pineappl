//! Module that contains helpers for binning observables

use serde::{Deserialize, Serialize};
use std::f64;
use thiserror::Error;

fn float_eq_within(lhs: f64, rhs: f64, ulps: usize) -> bool {
    // TODO: only works well enough if the numbers are far enough from zero or exacly zero
    if (lhs != 0.0) && (rhs != 0.0) {
        (lhs / rhs).abs().max((rhs / lhs).abs()) < f64::EPSILON.mul_add(ulps as f64, 1.0)
    } else {
        lhs == rhs
    }
}

#[derive(Deserialize, PartialEq, Serialize)]
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
#[derive(Deserialize, PartialEq, Serialize)]
pub struct BinLimits(Limits);

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
            .all(|val| float_eq_within(val[0], val[1], 8))
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
                    Some(((value - left) / (right - left) * (*bins as f64)) as usize)
                }
            }
            Limits::Unequal { limits } => {
                let index = limits
                    .binary_search_by(|left| left.partial_cmp(&value).unwrap())
                    .unwrap_or_else(|e| e);

                if index == 0 || index == limits.len() {
                    None
                } else {
                    Some(index - 1)
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
                .map(|b| (*right - *left).mul_add((b as f64) / (*bins as f64), *left))
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
            Limits::Equal { left, right, bins } => vec![(*right - *left) / (*bins as f64); *bins],
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
        if !float_eq_within(self.right(), other.left(), 8) {
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
    fn test() {
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
}
