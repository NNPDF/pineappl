//! Module for everything related to luminosity functions.

use noisy_float::types::R64;
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
    pub fn new(mut entry: Vec<(i32, i32, f64)>) -> Self {
        assert!(!entry.is_empty());

        // sort `entry` because the ordering doesn't matter and because it makes it easier to
        // compare `LumiEntry` objects with each other
        entry.sort_by(|x, y| (x.0, x.1, R64::new(x.2)).cmp(&(y.0, y.1, R64::new(y.2))));

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

#[cfg(test)]
mod tests {
    #[test]
    fn lumi_entry() {
        let entry1 = lumi_entry![2, 2, 1.0; 4, 4, 1.0];
        let entry2 = lumi_entry![4, 4, 1.0; 2, 2, 1.0];

        // checks that the ordering doesn't matter
        assert_eq!(entry1, entry2);
    }
}
