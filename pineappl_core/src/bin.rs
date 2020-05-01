//! Module that contains helpers for binning observables

use std::cmp::Ordering;
use std::f64;

enum Limits {
    Equal { left: f64, right: f64, bins: usize },
    Unequal { limits: Vec<f64> },
}

/// Structure representing bin limits.
pub struct BinLimits(Limits);

impl BinLimits {
    /// Constructor for `BinLimits`.
    #[must_use]
    pub fn new(mut limits: Vec<f64>) -> Self {
        limits.sort_by(|left, right| {
            if left < right {
                Ordering::Less
            } else {
                Ordering::Greater
            }
        });

        if limits
            .iter()
            .zip(limits.iter().skip(1))
            .map(|(current, next)| next - current)
            .collect::<Vec<f64>>()
            .windows(2)
            .all(|val| (val[0] / val[1]).max(val[1] / val[0]) <= f64::EPSILON.mul_add(8.0, 1.0))
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

    /// Returns the bin index for observable `value`. If the value over- or underflows, the return
    /// value is `Option::None`.
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
                    .binary_search_by(|left| {
                        if left < &value {
                            Ordering::Less
                        } else {
                            Ordering::Greater
                        }
                    })
                    .unwrap_or_else(|e| e);

                if index == 0 || index == limits.len() {
                    None
                } else {
                    Some(index - 1)
                }
            }
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
}

#[cfg(test)]
mod test {
    use super::*;

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
