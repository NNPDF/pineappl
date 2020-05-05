//! Module for everything related to luminosity functions.

use serde::{Deserialize, Serialize};

/// This structure represens an entry of a luminosity function. Each entry consists of a tuple,
/// which contains, in the following order, the PDG id of the first incoming parton, then the PDG
/// id of the second parton, and finally a numerical factor that will multiply the result for this
/// specific combination.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct LumiEntry {
    entry: Vec<(i32, i32, f64)>,
}

impl LumiEntry {
    /// Constructor for `LumiEntry`. Note that `entry` must be non-empty, otherwise this function
    /// panics.
    #[must_use]
    pub fn new(entry: Vec<(i32, i32, f64)>) -> Self {
        assert!(!entry.is_empty());

        Self { entry }
    }

    /// Returns a tuple representation of this entry.
    #[must_use]
    pub fn entry(&self) -> &[(i32, i32, f64)] {
        &self.entry
    }
}

/// Helper macro to quickly generate a LumiEntry at compile time. In the following example `entry1`
/// and `entry2` represent the same values:
///
/// ```
/// # use pineappl::lumi_entry;
/// # use pineappl::lumi::LumiEntry;
/// let entry1 = lumi_entry![2, 2, 1.0; 4, 4, 1.0];
/// let entry2 = LumiEntry::new(vec![(2, 2, 1.0), (4, 4, 1.0)]);
/// assert_eq!(entry1, entry2);
/// ```
#[macro_export]
macro_rules! lumi_entry {
    ($a:expr, $b:expr, $factor:expr $(; $c:expr, $d:expr, $fac:expr)*) => {
        $crate::lumi::LumiEntry::new(vec![($a, $b, $factor), $(($c, $d, $fac)),*])
    };
}
