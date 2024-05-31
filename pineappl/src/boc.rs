//! Module containing structures for the 3 dimensions of a [`Grid`]: bins, [`Order`] and channels
//! (`boc`).

use float_cmp::approx_eq;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::str::FromStr;
use thiserror::Error;

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
    /// Exponent of the logarithm of the scale factor of the factorization scale.
    pub logxif: u32,
}

impl FromStr for Order {
    type Err = ParseOrderError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut result = Self {
            alphas: 0,
            alpha: 0,
            logxir: 0,
            logxif: 0,
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
            .then((self.alpha, self.logxir, self.logxif).cmp(&(
                other.alpha,
                other.logxir,
                other.logxif,
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
    pub const fn new(alphas: u32, alpha: u32, logxir: u32, logxif: u32) -> Self {
        Self {
            alphas,
            alpha,
            logxir,
            logxif,
        }
    }

    /// Return a mask suitable to pass as the `order_mask` parameter of [`Grid::convolve`],
    /// [`Grid::evolve`] or [`Grid::evolve_info`]. The selection of `orders` is controlled using
    /// the `max_as` and `max_al` parameters, for instance setting `max_as = 1` and `max_al = 0`
    /// selects the LO QCD only, `max_as = 2` and `max_al = 0` the NLO QCD; setting `max_as = 3`
    /// and `max_al = 2` would select all NLOs, and the NNLO QCD.
    ///
    /// [`Grid::convolve`]: super::grid::Grid::convolve
    /// [`Grid::evolve`]: super::grid::Grid::evolve
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
    ///     Order::new(0, 2, 0, 0), //   LO        :          alpha^2
    ///     Order::new(1, 2, 0, 0), //  NLO QCD    : alphas   alpha^2
    ///     Order::new(0, 3, 0, 0), //  NLO  EW    :          alpha^3
    ///     Order::new(2, 2, 0, 0), // NNLO QCD    : alphas^2 alpha^2
    ///     Order::new(1, 3, 0, 0), // NNLO QCD—EW : alphas   alpha^3
    ///     Order::new(0, 4, 0, 0), // NNLO EW     :          alpha^4
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
    ///     Order::new(0, 2, 0, 0), //  LO         :        alpha^2
    ///     Order::new(1, 2, 0, 0), //  NLO QCD    : alphas alpha^2
    ///     Order::new(1, 2, 1, 0), //  NLO QCD    : alphas alpha^2 logxif
    ///     Order::new(0, 3, 0, 0), //  NLO  EW    :        alpha^3
    ///     Order::new(0, 3, 1, 0), //  NLO  EW    :        alpha^3 logxif
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
    ///     Order::new(2, 0, 0, 0), //   LO QCD    : alphas^2
    ///     Order::new(1, 1, 0, 0), //   LO QCD—EW : alphas   alpha
    ///     Order::new(0, 2, 0, 0), //   LO  EW    :          alpha^2
    ///     Order::new(3, 0, 0, 0), //  NLO QCD    : alphas^3
    ///     Order::new(2, 1, 0, 0), //  NLO QCD—EW : alphas^2 alpha
    ///     Order::new(1, 2, 0, 0), //  NLO QCD—EW : alphas   alpha^2
    ///     Order::new(0, 3, 0, 0), //  NLO EW     :          alpha^3
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
                 }| {
                    if !logs && (logxir > 0 || logxif > 0) {
                        return false;
                    }

                    let pto = alphas + alpha - lo;

                    alphas + alpha < min + lo
                        || (alphas + alpha < max + lo
                            && match max_as.cmp(&max_al) {
                                Ordering::Greater => lo_as + pto == alphas,
                                Ordering::Less => lo_al + pto == alpha,
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
    entry: Vec<(i32, i32, f64)>,
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
    /// let entry1 = Channel::new(vec![(2, 2, 1.0), (4, 4, 1.0)]);
    /// let entry2 = Channel::new(vec![(4, 4, 1.0), (2, 2, 1.0)]);
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
    /// let entry1 = Channel::new(vec![(1, 1, 1.0), (1, 1, 3.0), (3, 3, 1.0), (1, 1, 6.0)]);
    /// let entry2 = Channel::new(vec![(1, 1, 10.0), (3, 3, 1.0)]);
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
    #[must_use]
    pub fn new(mut entry: Vec<(i32, i32, f64)>) -> Self {
        assert!(!entry.is_empty());

        // sort `entry` because the ordering doesn't matter and because it makes it easier to
        // compare `Channel` objects with each other
        entry.sort_by(|x, y| (x.0, x.1).cmp(&(y.0, y.1)));

        Self {
            entry: entry
                .into_iter()
                .coalesce(|lhs, rhs| {
                    // sum the factors of repeated elements
                    if (lhs.0, lhs.1) == (rhs.0, rhs.1) {
                        Ok((lhs.0, lhs.1, lhs.2 + rhs.2))
                    } else {
                        Err((lhs, rhs))
                    }
                })
                // filter zeros
                // TODO: find a better than to hardcode the epsilon limit
                .filter(|&(_, _, f)| !approx_eq!(f64, f.abs(), 0.0, epsilon = 1e-14))
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
    /// let entry = Channel::translate(&channel![103, 11, 1.0], &|evol_id| match evol_id {
    ///     103 => vec![(2, 1.0), (-2, -1.0), (1, -1.0), (-1, 1.0)],
    ///     _ => vec![(evol_id, 1.0)],
    /// });
    ///
    /// assert_eq!(entry, channel![2, 11, 1.0; -2, 11, -1.0; 1, 11, -1.0; -1, 11, 1.0]);
    /// ```
    pub fn translate(entry: &Self, translator: &dyn Fn(i32) -> Vec<(i32, f64)>) -> Self {
        let mut tuples = Vec::new();

        for &(a, b, factor) in &entry.entry {
            for (aid, af) in translator(a) {
                for (bid, bf) in translator(b) {
                    tuples.push((aid, bid, factor * af * bf));
                }
            }
        }

        Self::new(tuples)
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
    /// assert_eq!(entry.entry(), [(2, 2, 1.0), (4, 4, 1.0)]);
    /// ```
    #[must_use]
    pub fn entry(&self) -> &[(i32, i32, f64)] {
        &self.entry
    }

    /// Creates a new object with the initial states transposed.
    #[must_use]
    pub fn transpose(&self) -> Self {
        Self::new(self.entry.iter().map(|(a, b, c)| (*b, *a, *c)).collect())
    }

    /// If `other` is the same channel when only comparing PIDs and neglecting the factors, return
    /// the number `f1 / f2`, where `f1` is the factor from `self` and `f2` is the factor from
    /// `other`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use pineappl::boc::Channel;
    ///
    /// let entry1 = Channel::new(vec![(2, 2, 2.0), (4, 4, 2.0)]);
    /// let entry2 = Channel::new(vec![(4, 4, 1.0), (2, 2, 1.0)]);
    /// let entry3 = Channel::new(vec![(3, 4, 1.0), (2, 2, 1.0)]);
    /// let entry4 = Channel::new(vec![(4, 3, 1.0), (2, 3, 2.0)]);
    ///
    /// assert_eq!(entry1.common_factor(&entry2), Some(2.0));
    /// assert_eq!(entry1.common_factor(&entry3), None);
    /// assert_eq!(entry1.common_factor(&entry4), None);
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
            .map(|(a, b)| ((a.0 == b.0) && (a.1 == b.1)).then_some(a.2 / b.2))
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
        Ok(Self::new(
            s.split('+')
                .map(|sub| {
                    sub.split_once('*').map_or_else(
                        || Err(ParseChannelError(format!("missing '*' in '{sub}'"))),
                        |(factor, pids)| {
                            let tuple = pids.split_once(',').map_or_else(
                                || Err(ParseChannelError(format!("missing ',' in '{pids}'"))),
                                |(a, b)| {
                                    Ok((
                                        a.trim()
                                            .strip_prefix('(')
                                            .ok_or_else(|| {
                                                ParseChannelError(format!(
                                                    "missing '(' in '{pids}'"
                                                ))
                                            })?
                                            .trim()
                                            .parse::<i32>()
                                            .map_err(|err| ParseChannelError(err.to_string()))?,
                                        b.trim()
                                            .strip_suffix(')')
                                            .ok_or_else(|| {
                                                ParseChannelError(format!(
                                                    "missing ')' in '{pids}'"
                                                ))
                                            })?
                                            .trim()
                                            .parse::<i32>()
                                            .map_err(|err| ParseChannelError(err.to_string()))?,
                                    ))
                                },
                            )?;

                            Ok((
                                tuple.0,
                                tuple.1,
                                str::parse::<f64>(factor.trim())
                                    .map_err(|err| ParseChannelError(err.to_string()))?,
                            ))
                        },
                    )
                })
                .collect::<Result<_, _>>()?,
        ))
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
    ($a:expr, $b:expr, $factor:expr $(; $c:expr, $d:expr, $fac:expr)*) => {
        $crate::boc::Channel::new(vec![($a, $b, $factor), $(($c, $d, $fac)),*])
    };
}

#[cfg(test)]
mod tests {
    use super::{Channel, Order, ParseOrderError};
    use crate::pids;

    #[test]
    fn order_from_str() {
        assert_eq!("as1".parse(), Ok(Order::new(1, 0, 0, 0)));
        assert_eq!("a1".parse(), Ok(Order::new(0, 1, 0, 0)));
        assert_eq!("as1lr1".parse(), Ok(Order::new(1, 0, 1, 0)));
        assert_eq!("as1lf1".parse(), Ok(Order::new(1, 0, 0, 1)));
        assert_eq!(
            "ab12".parse::<Order>(),
            Err(ParseOrderError("unknown coupling: 'ab'".to_owned()))
        );
        assert_eq!(
            "ab123456789000000".parse::<Order>(),
            Err(ParseOrderError(
                "error while parsing exponent of 'ab': number too large to fit in target type"
                    .to_owned()
            ))
        );
    }

    #[test]
    fn order_cmp() {
        let mut orders = [
            Order::new(1, 2, 1, 0),
            Order::new(1, 2, 0, 1),
            Order::new(1, 2, 0, 0),
            Order::new(0, 3, 1, 0),
            Order::new(0, 3, 0, 1),
            Order::new(0, 3, 0, 0),
            Order::new(0, 2, 0, 0),
        ];

        orders.sort();

        assert_eq!(orders[0], Order::new(0, 2, 0, 0));
        assert_eq!(orders[1], Order::new(1, 2, 0, 0));
        assert_eq!(orders[2], Order::new(1, 2, 0, 1));
        assert_eq!(orders[3], Order::new(1, 2, 1, 0));
        assert_eq!(orders[4], Order::new(0, 3, 0, 0));
        assert_eq!(orders[5], Order::new(0, 3, 0, 1));
        assert_eq!(orders[6], Order::new(0, 3, 1, 0));
    }

    #[test]
    fn order_create_mask() {
        // Drell—Yan orders
        let orders = [
            Order::new(0, 2, 0, 0), //   LO        :          alpha^2
            Order::new(1, 2, 0, 0), //  NLO QCD    : alphas   alpha^2
            Order::new(0, 3, 0, 0), //  NLO  EW    :          alpha^3
            Order::new(2, 2, 0, 0), // NNLO QCD    : alphas^2 alpha^2
            Order::new(1, 3, 0, 0), // NNLO QCD—EW : alphas   alpha^3
            Order::new(0, 4, 0, 0), // NNLO EW     :          alpha^4
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
            Order::new(2, 0, 0, 0), //   LO QCD    : alphas^2
            Order::new(1, 1, 0, 0), //   LO QCD—EW : alphas   alpha
            Order::new(0, 2, 0, 0), //   LO  EW    :          alpha^2
            Order::new(3, 0, 0, 0), //  NLO QCD    : alphas^3
            Order::new(2, 1, 0, 0), //  NLO QCD—EW : alphas^2 alpha
            Order::new(1, 2, 0, 0), //  NLO QCD—EW : alphas   alpha^2
            Order::new(0, 3, 0, 0), //  NLO  EW    :          alpha^3
            Order::new(4, 0, 0, 0), // NNLO QCD    : alphas^4
            Order::new(3, 1, 0, 0), // NNLO QCD—EW : alphas^3 alpha
            Order::new(2, 2, 0, 0), // NNLO QCD—EW : alphas^2 alpha^2
            Order::new(1, 3, 0, 0), // NNLO QCD—EW : alphas   alpha^3
            Order::new(0, 4, 0, 0), // NNLO EW     :          alpha^4
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
        let channel = Channel::translate(&channel![103, 203, 2.0], &pids::evol_to_pdg_mc_ids);

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
            "missing ',' in ' (  2 -2) '"
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
    }
}
