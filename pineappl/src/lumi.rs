//! Module for everything related to luminosity functions.

use serde::{Deserialize, Serialize};

/// This structure represens an entry of a luminosity function. Each entry consists of a tuple,
/// which contains, in the following order, the PDG id of the first incoming parton, then the PDG
/// id of the second parton, and finally a numerical factor that will multiply the result for this
/// specific combination.
#[derive(Clone, Debug, Deserialize, PartialEq, PartialOrd, Serialize)]
pub struct LumiEntry {
    entry: Vec<(i32, i32, f64)>,
}

impl LumiEntry {
    /// Constructor for `LumiEntry`. Note that `entry` must be non-empty, otherwise this function
    /// panics.
    ///
    /// # Examples
    ///
    /// Ordering of the arguments doesn't matter:
    ///
    /// ```rust
    /// use pineappl::lumi::LumiEntry;
    ///
    /// let entry1 = LumiEntry::new(vec![(2, 2, 1.0), (4, 4, 1.0)]);
    /// let entry2 = LumiEntry::new(vec![(4, 4, 1.0), (2, 2, 1.0)]);
    ///
    /// // checks that the ordering doesn't matter
    /// assert_eq!(entry1, entry2);
    /// ```
    ///
    /// # Panics
    ///
    /// Creating an entry with content panics:
    ///
    /// ```rust,should_panic
    /// use pineappl::lumi::LumiEntry;
    ///
    /// let _ = LumiEntry::new(vec![]);
    /// ```
    #[must_use]
    pub fn new(mut entry: Vec<(i32, i32, f64)>) -> Self {
        assert!(!entry.is_empty());

        // sort `entry` because the ordering doesn't matter and because it makes it easier to
        // compare `LumiEntry` objects with each other
        entry.sort_by(|x, y| (x.0, x.1, x.2).partial_cmp(&(y.0, y.1, y.2)).unwrap());

        Self { entry }
    }

    /// Returns a tuple representation of this entry.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use pineappl::lumi_entry;
    /// use pineappl::lumi::LumiEntry;
    ///
    /// let entry = lumi_entry![4, 4, 1.0; 2, 2, 1.0];
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
}

/// Helper macro to quickly generate a LumiEntry at compile time.
///
/// # Examples
///
/// In the following example `entry1` and `entry2` represent the same values:
///
/// ```rust
/// use pineappl::lumi_entry;
///
/// let entry1 = lumi_entry![2, 2, 1.0; 4, 4, 1.0];
/// let entry2 = lumi_entry![4, 4, 1.0; 2, 2, 1.0];
///
/// assert_eq!(entry1, entry2);
/// ```
#[macro_export]
macro_rules! lumi_entry {
    ($a:expr, $b:expr, $factor:expr $(; $c:expr, $d:expr, $fac:expr)*) => {
        $crate::lumi::LumiEntry::new(vec![($a, $b, $factor), $(($c, $d, $fac)),*])
    };
}
