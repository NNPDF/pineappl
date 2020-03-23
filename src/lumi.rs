//! Module containing all structures for creating luminosity functions.

/// Structure representing an entry of a luminosity function. Each entry consists of a tuple, which
/// contains, in the following order, the PDG id of the first incoming parton, then the PDG id of
/// the second parton and finally a numerical factor that will multiply the result for this
/// specific combination
#[derive(Debug, PartialEq)]
pub struct LumiEntry {
    entry: Vec<(i32, i32, f64)>,
}

impl LumiEntry {
    /// Constructor for LumiEntry.
    pub fn new(entry: Vec<(i32, i32, f64)>) -> LumiEntry {
        LumiEntry { entry }
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

/// Structure implementing a luminosity function. Each luminosity function is collection of
/// LumiEntries.
pub struct Lumi {
    tuples: Vec<LumiEntry>,
}

impl Lumi {
    /// Constructor for Lumi.
    pub fn new(tuples: Vec<LumiEntry>) -> Lumi {
        Lumi { tuples }
    }

    /// Returns the number of LumiEntries in the Lumi object.
    pub fn len(&self) -> usize {
        self.tuples.len()
    }

    /// Checks if there are no LumiEntries in this Lumi object.
    pub fn is_empty(&self) -> bool {
        self.tuples.is_empty()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test() {
        let lumi = Lumi::new(vec![
            lumi_entry![2, 2, 1.0; 4, 4, 1.0],
            lumi_entry![1, 1, 1.0; 3, 3, 1.0],
            lumi_entry![5, 5, 1.0],
        ]);
        assert_eq!(lumi.len(), 3);
    }
}
