//! Module containing structures for the 3 dimensions of a [`Grid`]: bins, [`Order`] and channels
//! (`boc`).
//!
//! [`Grid`]: super::grid::Grid

use super::subgrid::NodeValues;
use float_cmp::approx_eq;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::str::FromStr;
use thiserror::Error;

/// TODO
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub enum Kinematics {
    /// TODO
    Scale(usize),
    /// TODO
    X(usize),
}

impl Kinematics {
    /// TODO
    pub const X1: Self = Self::X(0);

    /// TODO
    pub const X2: Self = Self::X(1);
}

/// TODO
#[derive(Clone, Deserialize, Serialize)]
pub enum ScaleFuncForm {
    /// TODO
    NoScale,
    /// TODO
    Scale(usize),
    /// TODO
    QuadraticSum(usize, usize),
}

impl ScaleFuncForm {
    /// TODO
    pub fn calc(&self, node_values: &[NodeValues], kinematics: &[Kinematics]) -> Option<Vec<f64>> {
        match self {
            ScaleFuncForm::NoScale => None,
            &ScaleFuncForm::Scale(index) => Some(if node_values.is_empty() {
                // TODO: empty subgrid should have as many node values as dimensions
                Vec::new()
            } else {
                node_values[kinematics
                    .iter()
                    .position(|&kin| kin == Kinematics::Scale(index))
                    // UNWRAP: this should be guaranteed by `Grid::new`
                    .unwrap()]
                .values()
            }),
            ScaleFuncForm::QuadraticSum(_, _) => todo!(),
        }
    }
}

/// TODO
#[derive(Clone, Deserialize, Serialize)]
pub struct Scales {
    /// TODO
    pub ren: ScaleFuncForm,
    /// TODO
    pub fac: ScaleFuncForm,
    /// TODO
    pub frg: ScaleFuncForm,
}

impl Scales {
    /// TODO
    pub fn compatible_with(&self, kinematics: &[Kinematics]) -> bool {
        for scale in [&self.ren, &self.fac, &self.frg].map(Clone::clone) {
            match scale {
                ScaleFuncForm::NoScale => {}
                ScaleFuncForm::Scale(index)
                    if kinematics
                        .iter()
                        .any(|&kin| kin == Kinematics::Scale(index)) => {}
                ScaleFuncForm::QuadraticSum(idx1, idx2)
                    if kinematics.iter().any(|&kin| kin == Kinematics::Scale(idx1))
                        && kinematics.iter().any(|&kin| kin == Kinematics::Scale(idx2)) => {}
                _ => return false,
            }
        }

        true
    }
}

/// Error type keeping information if [`Order::from_str`] went wrong.
#[derive(Debug, Error, Eq, PartialEq)]
#[error("{0}")]
pub struct ParseOrderError(String);

// TODO: when possible change the types from `u32` to `u8` to change `try_into` to `into`

/// Coupling powers for each grid.
#[derive(Clone, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
pub struct Order {
    /// Exponent of the strong coupling.
    pub alphas: u32,
    /// Exponent of the electromagnetic coupling.
    pub alpha: u32,
    /// Exponent of the logarithm of the scale factor of the renomalization scale.
    pub logxir: u32,
    /// Exponent of the logarithm of the scale factor of the initial state factorization scale.
    pub logxif: u32,
    /// Exponent of the logarithm of the scale factor of the final state factorization scale
    /// (fragmentation scale).
    pub logxia: u32,
}

impl FromStr for Order {
    type Err = ParseOrderError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut result = Self {
            alphas: 0,
            alpha: 0,
            logxir: 0,
            logxif: 0,
            logxia: 0,
        };

        for tuple in s
            .split(|c: char| c.is_ascii_digit())
            .filter(|s| !s.is_empty())
            .zip(
                s.split(|c: char| !c.is_ascii_digit())
                    .filter(|s| !s.is_empty())
                    .map(str::parse),
            )
        {
            match tuple {
                ("as", Ok(num)) => {
                    result.alphas = num;
                }
                ("a", Ok(num)) => {
                    result.alpha = num;
                }
                ("lr", Ok(num)) => {
                    result.logxir = num;
                }
                ("lf", Ok(num)) => {
                    result.logxif = num;
                }
                ("la", Ok(num)) => {
                    result.logxia = num;
                }
                (label, Err(err)) => {
                    return Err(ParseOrderError(format!(
                        "error while parsing exponent of '{label}': {err}"
                    )));
                }
                (label, Ok(_)) => {
                    return Err(ParseOrderError(format!("unknown coupling: '{label}'")));
                }
            }
        }

        Ok(result)
    }
}

impl Ord for Order {
    fn cmp(&self, other: &Self) -> Ordering {
        // sort leading orders before next-to-leading orders, then the lowest power in alpha, the
        // rest lexicographically
        (self.alphas + self.alpha)
            .cmp(&(other.alphas + other.alpha))
            .then((self.alpha, self.logxir, self.logxif, self.logxia).cmp(&(
                other.alpha,
                other.logxir,
                other.logxif,
                other.logxia,
            )))
    }
}

impl PartialOrd for Order {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Order {
    /// Constructor. This function mainly exists to have a way of constructing `Order` that is less
    /// verbose.
    #[must_use]
    pub const fn new(alphas: u32, alpha: u32, logxir: u32, logxif: u32, logxia: u32) -> Self {
        Self {
            alphas,
            alpha,
            logxir,
            logxif,
            logxia,
        }
    }

    /// Return a mask suitable to pass as the `order_mask` parameter of [`Grid::convolve`],
    /// [`Grid::evolve_with_slice_iter`] or [`Grid::evolve_info`]. The selection of `orders` is
    /// controlled using the `max_as` and `max_al` parameters, for instance setting `max_as = 1`
    /// and `max_al = 0` selects the LO QCD only, `max_as = 2` and `max_al = 0` the NLO QCD;
    /// setting `max_as = 3` and `max_al = 2` would select all NLOs, and the NNLO QCD.
    ///
    /// [`Grid::convolve`]: super::grid::Grid::convolve
    /// [`Grid::evolve_with_slice_iter`]: super::grid::Grid::evolve_with_slice_iter
    /// [`Grid::evolve_info`]: super::grid::Grid::evolve_info
    ///
    /// # Example
    ///
    /// In the case of Drell—Yan, there are the following orders:
    ///
    /// - exactly one leading order (LO),
    /// - two next-to-leading orders (NLO), which are
    ///   - the NLO QCD and
    ///   - the NLO EW, and
    /// - three next-to-next-to-leading orders (NNLO),
    ///   - the NNLO QCD,
    ///   - the NNLO EW, and finally
    ///   - the mixed NNLO QCD—EW.
    ///
    /// ```rust
    /// use pineappl::boc::Order;
    ///
    /// let orders = [
    ///     Order::new(0, 2, 0, 0, 0), //   LO        :          alpha^2
    ///     Order::new(1, 2, 0, 0, 0), //  NLO QCD    : alphas   alpha^2
    ///     Order::new(0, 3, 0, 0, 0), //  NLO  EW    :          alpha^3
    ///     Order::new(2, 2, 0, 0, 0), // NNLO QCD    : alphas^2 alpha^2
    ///     Order::new(1, 3, 0, 0, 0), // NNLO QCD—EW : alphas   alpha^3
    ///     Order::new(0, 4, 0, 0, 0), // NNLO EW     :          alpha^4
    /// ];
    ///
    /// // LO EW
    /// assert_eq!(Order::create_mask(&orders, 0, 1, false), [true, false, false, false, false, false]);
    /// // LO QCD
    /// assert_eq!(Order::create_mask(&orders, 1, 0, false), [true, false, false, false, false, false]);
    /// // LO
    /// assert_eq!(Order::create_mask(&orders, 1, 1, false), [true, false, false, false, false, false]);
    /// // NLO QCD
    /// assert_eq!(Order::create_mask(&orders, 2, 0, false), [true, true, false, false, false, false]);
    /// // NLO EW
    /// assert_eq!(Order::create_mask(&orders, 0, 2, false), [true, false, true, false, false, false]);
    /// // NNLO QCD
    /// assert_eq!(Order::create_mask(&orders, 3, 0, false), [true, true, false, true, false, false]);
    /// // NNLO EW
    /// assert_eq!(Order::create_mask(&orders, 0, 3, false), [true, false, true, false, false, true]);
    /// ```
    ///
    /// Orders containing non-zero powers of logarithms can be selected as well if `logs` is set to
    /// `true`:
    ///
    /// ```rust
    /// use pineappl::boc::Order;
    ///
    /// let orders = [
    ///     Order::new(0, 2, 0, 0, 0), //  LO         :        alpha^2
    ///     Order::new(1, 2, 0, 0, 0), //  NLO QCD    : alphas alpha^2
    ///     Order::new(1, 2, 1, 0, 0), //  NLO QCD    : alphas alpha^2 logxif
    ///     Order::new(0, 3, 0, 0, 0), //  NLO  EW    :        alpha^3
    ///     Order::new(0, 3, 1, 0, 0), //  NLO  EW    :        alpha^3 logxif
    /// ];
    ///
    /// assert_eq!(Order::create_mask(&orders, 0, 2, true), [true, false, false, true, true]);
    /// ```
    ///
    /// For the more complicated example of top-pair production one can see the difference between
    /// the selection for different LOs:
    ///
    /// ```rust
    /// use pineappl::boc::Order;
    ///
    /// let orders = [
    ///     Order::new(2, 0, 0, 0, 0), //   LO QCD    : alphas^2
    ///     Order::new(1, 1, 0, 0, 0), //   LO QCD—EW : alphas   alpha
    ///     Order::new(0, 2, 0, 0, 0), //   LO  EW    :          alpha^2
    ///     Order::new(3, 0, 0, 0, 0), //  NLO QCD    : alphas^3
    ///     Order::new(2, 1, 0, 0, 0), //  NLO QCD—EW : alphas^2 alpha
    ///     Order::new(1, 2, 0, 0, 0), //  NLO QCD—EW : alphas   alpha^2
    ///     Order::new(0, 3, 0, 0, 0), //  NLO EW     :          alpha^3
    /// ];
    ///
    /// // LO EW
    /// assert_eq!(Order::create_mask(&orders, 0, 1, false), [false, false, true, false, false, false, false]);
    /// // LO QCD
    /// assert_eq!(Order::create_mask(&orders, 1, 0, false), [true, false, false, false, false, false, false]);
    /// // LO
    /// assert_eq!(Order::create_mask(&orders, 1, 1, false), [true, true, true, false, false, false, false]);
    /// ```
    #[must_use]
    pub fn create_mask(orders: &[Self], max_as: u32, max_al: u32, logs: bool) -> Vec<bool> {
        // smallest sum of alphas and alpha
        let lo = orders
            .iter()
            .map(|Self { alphas, alpha, .. }| alphas + alpha)
            .min()
            .unwrap_or_default();

        // all leading orders, without logarithms
        let leading_orders: Vec<_> = orders
            .iter()
            .filter(|Self { alphas, alpha, .. }| alphas + alpha == lo)
            .cloned()
            .collect();

        let lo_as = leading_orders
            .iter()
            .map(|Self { alphas, .. }| *alphas)
            .max()
            .unwrap_or_default();
        let lo_al = leading_orders
            .iter()
            .map(|Self { alpha, .. }| *alpha)
            .max()
            .unwrap_or_default();

        let max = max_as.max(max_al);
        let min = max_as.min(max_al);

        orders
            .iter()
            .map(
                |&Self {
                     alphas,
                     alpha,
                     logxir,
                     logxif,
                     logxia,
                 }| {
                    if !logs && (logxir > 0 || logxif > 0 || logxia > 0) {
                        return false;
                    }

                    let pto = alphas + alpha - lo;

                    alphas + alpha < min + lo
                        || (alphas + alpha < max + lo
                            && match max_as.cmp(&max_al) {
                                Ordering::Greater => lo_as + pto == alphas,
                                Ordering::Less => lo_al + pto == alpha,
                                // TODO: when do we hit this condition?
                                Ordering::Equal => false,
                            })
                },
            )
            .collect()
    }
}

/// This structure represents a channel. Each channel consists of a tuple containing in the
/// following order, the particle ID of the first incoming parton, then the particle ID of the
/// second parton, and finally a numerical factor that will multiply the result for this specific
/// combination.
#[derive(Clone, Debug, Deserialize, PartialEq, PartialOrd, Serialize)]
pub struct Channel {
    entry: Vec<(Vec<i32>, f64)>,
}

impl Channel {
    /// Constructor for `Channel`. Note that `entry` must be non-empty, otherwise this function
    /// panics.
    ///
    /// # Examples
    ///
    /// Ordering of the arguments doesn't matter:
    ///
    /// ```rust
    /// use pineappl::boc::Channel;
    ///
    /// let entry1 = Channel::new(vec![(vec![2, 2], 1.0), (vec![4, 4], 1.0)]);
    /// let entry2 = Channel::new(vec![(vec![4, 4], 1.0), (vec![2, 2], 1.0)]);
    ///
    /// // checks that the ordering doesn't matter
    /// assert_eq!(entry1, entry2);
    /// ```
    ///
    /// Same arguments are merged together:
    ///
    /// ```rust
    /// use pineappl::boc::Channel;
    ///
    /// let entry1 = Channel::new(vec![(vec![1, 1], 1.0), (vec![1, 1], 3.0), (vec![3, 3], 1.0), (vec![1, 1], 6.0)]);
    /// let entry2 = Channel::new(vec![(vec![1, 1], 10.0), (vec![3, 3], 1.0)]);
    ///
    /// assert_eq!(entry1, entry2);
    /// ```
    ///
    /// # Panics
    ///
    /// Creating an empty channel panics:
    ///
    /// ```rust,should_panic
    /// use pineappl::boc::Channel;
    ///
    /// let _ = Channel::new(vec![]);
    /// ```
    ///
    /// Creating a channel with entries that have a different number of PIDs panics:
    /// ```rust,should_panic
    /// use pineappl::boc::Channel;
    ///
    /// let _ = Channel::new(vec![(vec![1, 1, 1], 1.0), (vec![1, 1], 1.0)]);
    /// ```
    #[must_use]
    pub fn new(mut entry: Vec<(Vec<i32>, f64)>) -> Self {
        assert!(!entry.is_empty(), "can not create empty channel");
        assert!(
            entry.iter().map(|(pids, _)| pids.len()).all_equal(),
            "can not create channel with a different number of PIDs"
        );

        // sort `entry` because the ordering doesn't matter and because it makes it easier to
        // compare `Channel` objects with each other
        entry.sort_by(|x, y| x.0.cmp(&y.0));

        Self {
            entry: entry
                .into_iter()
                .coalesce(|lhs, rhs| {
                    // sum the factors of repeated elements
                    if lhs.0 == rhs.0 {
                        Ok((lhs.0, lhs.1 + rhs.1))
                    } else {
                        Err((lhs, rhs))
                    }
                })
                // filter zeros
                // TODO: find a better than to hardcode the epsilon limit
                .filter(|&(_, f)| !approx_eq!(f64, f.abs(), 0.0, epsilon = 1e-14))
                .collect(),
        }
    }

    /// Translates `entry` into a different basis using `translator`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use pineappl::boc::Channel;
    /// use pineappl::channel;
    ///
    /// let entry = channel![103, 11, 10.0].translate(&|evol_id| match evol_id {
    ///     103 => vec![(2, 1.0), (-2, -1.0), (1, -1.0), (-1, 1.0)],
    ///     _ => vec![(evol_id, 1.0)],
    /// });
    ///
    /// assert_eq!(entry, channel![2, 11, 10.0; -2, 11, -10.0; 1, 11, -10.0; -1, 11, 10.0]);
    /// ```
    #[must_use]
    pub fn translate(&self, translator: &dyn Fn(i32) -> Vec<(i32, f64)>) -> Self {
        let mut result = Vec::new();

        for (pids, factor) in &self.entry {
            for tuples in pids
                .iter()
                .map(|&pid| translator(pid))
                .multi_cartesian_product()
            {
                result.push((
                    tuples.iter().map(|&(pid, _)| pid).collect(),
                    factor * tuples.iter().map(|(_, f)| f).product::<f64>(),
                ));
            }
        }

        Self::new(result)
    }

    /// Returns a tuple representation of this entry.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use pineappl::channel;
    /// use pineappl::boc::Channel;
    ///
    /// let entry = channel![4, 4, 1.0; 2, 2, 1.0];
    ///
    /// assert_eq!(entry.entry(), [(vec![2, 2], 1.0), (vec![4, 4], 1.0)]);
    /// ```
    #[must_use]
    pub fn entry(&self) -> &[(Vec<i32>, f64)] {
        &self.entry
    }

    /// Create a new object with the PIDs at index `i` and `j` transposed.
    #[must_use]
    pub fn transpose(&self, i: usize, j: usize) -> Self {
        Self::new(
            self.entry
                .iter()
                .map(|(pids, c)| {
                    let mut transposed = pids.clone();
                    transposed.swap(i, j);
                    (transposed, *c)
                })
                .collect(),
        )
    }

    /// If `other` is the same channel when only comparing PIDs and neglecting the factors, return
    /// the number `f1 / f2`, where `f1` is the factor from `self` and `f2` is the factor from
    /// `other`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use pineappl::channel;
    ///
    /// let ch1 = channel![2, 2, 2.0; 4, 4, 2.0];
    /// let ch2 = channel![4, 4, 1.0; 2, 2, 1.0];
    /// let ch3 = channel![3, 4, 1.0; 2, 2, 1.0];
    /// let ch4 = channel![4, 3, 1.0; 2, 3, 2.0];
    /// let ch5 = channel![2, 2, 1.0; 4, 4, 2.0];
    ///
    /// // ch1 is ch2 multiplied by two
    /// assert_eq!(ch1.common_factor(&ch2), Some(2.0));
    /// // ch1 isn't similar to ch3
    /// assert_eq!(ch1.common_factor(&ch3), None);
    /// // ch1 isn't similar to ch4 either
    /// assert_eq!(ch1.common_factor(&ch4), None);
    /// // ch1 is similar to ch5, but they don't share a common factor
    /// assert_eq!(ch1.common_factor(&ch5), None);
    /// ```
    #[must_use]
    pub fn common_factor(&self, other: &Self) -> Option<f64> {
        if self.entry.len() != other.entry.len() {
            return None;
        }

        let result: Option<Vec<_>> = self
            .entry
            .iter()
            .zip(&other.entry)
            .map(|((pids_a, fa), (pids_b, fb))| (pids_a == pids_b).then_some(fa / fb))
            .collect();

        result.and_then(|factors| {
            if factors
                .windows(2)
                .all(|win| approx_eq!(f64, win[0], win[1], ulps = 4))
            {
                factors.first().copied()
            } else {
                None
            }
        })
    }
}

/// Error type keeping information if [`Channel::from_str`] went wrong.
#[derive(Debug, Error)]
#[error("{0}")]
pub struct ParseChannelError(String);

impl FromStr for Channel {
    type Err = ParseChannelError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let result: Vec<_> = s
            .split('+')
            .map(|sub| {
                sub.split_once('*').map_or_else(
                    // TODO: allow a missing numerical factor which then is assumed to be `1`
                    || Err(ParseChannelError(format!("missing '*' in '{sub}'"))),
                    |(factor, pids)| {
                        let vector: Vec<_> = pids
                            .trim()
                            .strip_prefix('(')
                            .ok_or_else(|| ParseChannelError(format!("missing '(' in '{pids}'")))?
                            .strip_suffix(')')
                            .ok_or_else(|| ParseChannelError(format!("missing ')' in '{pids}'")))?
                            .split(',')
                            .map(|pid| {
                                pid.trim().parse::<i32>().map_err(|err| {
                                    ParseChannelError(format!(
                                        "could not parse PID: '{pid}', '{err}'"
                                    ))
                                })
                            })
                            .collect::<Result<_, _>>()?;

                        Ok((
                            vector,
                            str::parse::<f64>(factor.trim())
                                .map_err(|err| ParseChannelError(err.to_string()))?,
                        ))
                    },
                )
            })
            .collect::<Result<_, _>>()?;

        if !result.iter().map(|(pids, _)| pids.len()).all_equal() {
            return Err(ParseChannelError(
                "PID tuples have different lengths".to_owned(),
            ));
        }

        Ok(Self::new(result))
    }
}

/// Helper macro to quickly generate a `Channel` at compile time.
///
/// # Examples
///
/// In the following example `entry1` and `entry2` represent the same values:
///
/// ```rust
/// use pineappl::channel;
///
/// let entry1 = channel![2, 2, 1.0; 4, 4, 1.0];
/// let entry2 = channel![4, 4, 1.0; 2, 2, 1.0];
///
/// assert_eq!(entry1, entry2);
/// ```
#[macro_export]
macro_rules! channel {
    // TODO: generalize this to accept an arbitrary number of PIDs
    ($a:expr, $b:expr, $factor:expr $(; $c:expr, $d:expr, $fac:expr)*) => {
        $crate::boc::Channel::new(vec![(vec![$a, $b], $factor), $((vec![$c, $d], $fac)),*])
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pids;

    #[test]
    fn order_from_str() {
        assert_eq!("as1".parse(), Ok(Order::new(1, 0, 0, 0, 0)));
        assert_eq!("a1".parse(), Ok(Order::new(0, 1, 0, 0, 0)));
        assert_eq!("as1lr1".parse(), Ok(Order::new(1, 0, 1, 0, 0)));
        assert_eq!("as1lf1".parse(), Ok(Order::new(1, 0, 0, 1, 0)));
        assert_eq!("as1la1".parse(), Ok(Order::new(1, 0, 0, 0, 1)));
        assert_eq!(
            "ab12".parse::<Order>().unwrap_err().to_string(),
            "unknown coupling: 'ab'"
        );
        assert_eq!(
            "ab123456789000000"
                .parse::<Order>()
                .unwrap_err()
                .to_string(),
            "error while parsing exponent of 'ab': number too large to fit in target type"
        );
    }

    #[test]
    fn order_cmp() {
        let mut orders = [
            Order::new(1, 2, 1, 0, 0),
            Order::new(1, 2, 0, 1, 0),
            Order::new(1, 2, 0, 0, 0),
            Order::new(0, 3, 1, 0, 0),
            Order::new(0, 3, 0, 1, 0),
            Order::new(0, 3, 0, 0, 0),
            Order::new(0, 2, 0, 0, 0),
        ];

        orders.sort();

        assert_eq!(orders[0], Order::new(0, 2, 0, 0, 0));
        assert_eq!(orders[1], Order::new(1, 2, 0, 0, 0));
        assert_eq!(orders[2], Order::new(1, 2, 0, 1, 0));
        assert_eq!(orders[3], Order::new(1, 2, 1, 0, 0));
        assert_eq!(orders[4], Order::new(0, 3, 0, 0, 0));
        assert_eq!(orders[5], Order::new(0, 3, 0, 1, 0));
        assert_eq!(orders[6], Order::new(0, 3, 1, 0, 0));
    }

    #[test]
    fn order_create_mask() {
        // Drell—Yan orders
        let orders = [
            Order::new(0, 2, 0, 0, 0), //   LO        :          alpha^2
            Order::new(1, 2, 0, 0, 0), //  NLO QCD    : alphas   alpha^2
            Order::new(0, 3, 0, 0, 0), //  NLO  EW    :          alpha^3
            Order::new(2, 2, 0, 0, 0), // NNLO QCD    : alphas^2 alpha^2
            Order::new(1, 3, 0, 0, 0), // NNLO QCD—EW : alphas   alpha^3
            Order::new(0, 4, 0, 0, 0), // NNLO EW     :          alpha^4
        ];

        assert_eq!(
            Order::create_mask(&orders, 0, 0, false),
            [false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 0, 1, false),
            [true, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 0, 2, false),
            [true, false, true, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 0, 3, false),
            [true, false, true, false, false, true]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 0, false),
            [true, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 1, false),
            [true, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 2, false),
            [true, false, true, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 3, false),
            [true, false, true, false, false, true]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 0, false),
            [true, true, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 1, false),
            [true, true, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 2, false),
            [true, true, true, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 3, false),
            [true, true, true, false, false, true]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 0, false),
            [true, true, false, true, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 1, false),
            [true, true, false, true, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 2, false),
            [true, true, true, true, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 3, false),
            [true, true, true, true, true, true]
        );

        // Top-pair production orders
        let orders = [
            Order::new(2, 0, 0, 0, 0), //   LO QCD    : alphas^2
            Order::new(1, 1, 0, 0, 0), //   LO QCD—EW : alphas   alpha
            Order::new(0, 2, 0, 0, 0), //   LO  EW    :          alpha^2
            Order::new(3, 0, 0, 0, 0), //  NLO QCD    : alphas^3
            Order::new(2, 1, 0, 0, 0), //  NLO QCD—EW : alphas^2 alpha
            Order::new(1, 2, 0, 0, 0), //  NLO QCD—EW : alphas   alpha^2
            Order::new(0, 3, 0, 0, 0), //  NLO  EW    :          alpha^3
            Order::new(4, 0, 0, 0, 0), // NNLO QCD    : alphas^4
            Order::new(3, 1, 0, 0, 0), // NNLO QCD—EW : alphas^3 alpha
            Order::new(2, 2, 0, 0, 0), // NNLO QCD—EW : alphas^2 alpha^2
            Order::new(1, 3, 0, 0, 0), // NNLO QCD—EW : alphas   alpha^3
            Order::new(0, 4, 0, 0, 0), // NNLO EW     :          alpha^4
        ];

        assert_eq!(
            Order::create_mask(&orders, 0, 0, false),
            [false, false, false, false, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 0, 1, false),
            [false, false, true, false, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 0, 2, false),
            [false, false, true, false, false, false, true, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 0, 3, false),
            [false, false, true, false, false, false, true, false, false, false, false, true]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 0, false),
            [true, false, false, false, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 1, false),
            [true, true, true, false, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 2, false),
            [true, true, true, false, false, false, true, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 3, false),
            [true, true, true, false, false, false, true, false, false, false, false, true]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 0, false),
            [true, false, false, true, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 1, false),
            [true, true, true, true, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 2, false),
            [true, true, true, true, true, true, true, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 3, false),
            [true, true, true, true, true, true, true, false, false, false, false, true]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 0, false),
            [true, false, false, true, false, false, false, true, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 1, false),
            [true, true, true, true, false, false, false, true, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 2, false),
            [true, true, true, true, true, true, true, true, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 3, false),
            [true, true, true, true, true, true, true, true, true, true, true, true]
        );
    }

    #[test]
    fn channel_translate() {
        let channel = channel![103, 203, 2.0].translate(&pids::evol_to_pdg_mc_ids);

        assert_eq!(
            channel,
            channel![ 2,  2,  2.0;  2, -2, -2.0;  2,  1, -2.0;  2, -1,  2.0;
                     -2,  2,  2.0; -2, -2, -2.0; -2,  1, -2.0; -2, -1,  2.0;
                      1,  2, -2.0;  1, -2,  2.0;  1,  1,  2.0;  1, -1, -2.0;
                     -1,  2, -2.0; -1, -2,  2.0; -1,  1,  2.0; -1, -1, -2.0]
        );
    }

    #[test]
    fn channel_from_str() {
        assert_eq!(
            str::parse::<Channel>(" 1   * (  2 , -2) + 2* (4,-4)").unwrap(),
            channel![2, -2, 1.0; 4, -4, 2.0]
        );

        assert_eq!(
            str::parse::<Channel>("* (  2, -2) + 2* (4,-4)")
                .unwrap_err()
                .to_string(),
            "cannot parse float from empty string"
        );

        assert_eq!(
            str::parse::<Channel>(" 1    (  2 -2) + 2* (4,-4)")
                .unwrap_err()
                .to_string(),
            "missing '*' in ' 1    (  2 -2) '"
        );

        assert_eq!(
            str::parse::<Channel>(" 1   * (  2 -2) + 2* (4,-4)")
                .unwrap_err()
                .to_string(),
            "could not parse PID: '  2 -2', 'invalid digit found in string'"
        );

        assert_eq!(
            str::parse::<Channel>(" 1   *   2, -2) + 2* (4,-4)")
                .unwrap_err()
                .to_string(),
            "missing '(' in '   2, -2) '"
        );

        assert_eq!(
            str::parse::<Channel>(" 1   * (  2, -2 + 2* (4,-4)")
                .unwrap_err()
                .to_string(),
            "missing ')' in ' (  2, -2 '"
        );

        assert_eq!(
            str::parse::<Channel>("1 * (2, 2, 2) + 2 * (4, 4)")
                .unwrap_err()
                .to_string(),
            "PID tuples have different lengths"
        );
    }
}
