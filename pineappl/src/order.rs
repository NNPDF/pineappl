//! TODO

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

    /// Return a mask suitable to pass as the `order_mask` parameter of [`Grid::convolute`],
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
    /// use pineappl::order::Order;
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
    /// use pineappl::order::Order;
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
    /// use pineappl::order::Order;
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

#[cfg(test)]
mod tests {
    use super::{Order, ParseOrderError};

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
}
