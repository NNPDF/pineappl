//! Module for everything related to luminosity functions.

use super::grid::Grid;
use super::subgrid::Mu2;
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

    /// Translates `entry` into a different basis using `translator`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use pineappl::lumi::LumiEntry;
    /// use pineappl::lumi_entry;
    ///
    /// let entry = LumiEntry::translate(&lumi_entry![103, 11, 1.0], &|evol_id| match evol_id {
    ///     103 => vec![(2, 1.0), (-2, -1.0), (1, -1.0), (-1, 1.0)],
    ///     _ => vec![(evol_id, 1.0)],
    /// });
    ///
    /// assert_eq!(entry, lumi_entry![2, 11, 1.0; -2, 11, -1.0; 1, 11, -1.0; -1, 11, 1.0]);
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

enum Pdfs<'a> {
    Two {
        xfx1: &'a mut dyn FnMut(i32, f64, f64) -> f64,
        xfx2: &'a mut dyn FnMut(i32, f64, f64) -> f64,
    },
    One {
        xfx: &'a mut dyn FnMut(i32, f64, f64) -> f64,
    },
}

/// A cache for evaluating PDFs. Methods like [`Grid::convolute2`] accept instances of this
/// `struct` instead of the PDFs themselves.
pub struct LumiCache<'a> {
    pdfs: Pdfs<'a>,
    alphas: &'a mut dyn FnMut(f64) -> f64,
    mu2_grid: Vec<Mu2>,
    x1_grid: Vec<f64>,
    x2_grid: Vec<f64>,
    cc1: i32,
    cc2: i32,
}

impl<'a> LumiCache<'a> {
    fn new(
        grid: &Grid,
        pdg1: i32,
        pdg2: i32,
        pdfs: Pdfs<'a>,
        alphas: &'a mut dyn FnMut(f64) -> f64,
    ) -> Result<Self, ()> {
        // PDG identifiers of the initial states
        let pdga = grid.key_values().map_or(2212, |kv| {
            kv.get("initial_state_1")
                .map_or(2212, |s| s.parse::<i32>().unwrap())
        });
        let pdgb = grid.key_values().map_or(2212, |kv| {
            kv.get("initial_state_2")
                .map_or(2212, |s| s.parse::<i32>().unwrap())
        });

        // are the initial states hadrons?
        let has_pdfa = grid
            .lumi()
            .iter()
            .all(|entry| entry.entry().iter().all(|&(a, _, _)| a == pdga));
        let has_pdfb = grid
            .lumi()
            .iter()
            .all(|entry| entry.entry().iter().all(|&(_, b, _)| b == pdgb));

        // do we have to charge-conjugate the initial states?
        let cc1 = if pdg1 == pdga {
            1
        } else if pdg1 == -pdga {
            -1
        } else if has_pdfa {
            return Err(());
        } else {
            0
        };
        let cc2 = if pdg2 == pdgb {
            1
        } else if pdg2 == -pdgb {
            -1
        } else if has_pdfb {
            return Err(());
        } else {
            0
        };

        Ok(Self {
            pdfs,
            alphas,
            mu2_grid: vec![],
            x1_grid: vec![],
            x2_grid: vec![],
            cc1,
            cc2,
        })
    }

    /// Construct a luminosity cache with two PDFs, `xfx1` and `xfx2`. The types of hadrons the
    /// PDFs correspond to must be given as `pdg1` and `pdg2`. The function to evaluate the
    /// strong coupling must be given as `alphas`. The grid that the cache will be used with must
    /// be given as `grid`; this parameter determines which of the initial states are hadronic, and
    /// if an initial states is not hadronic the corresponding 'PDF' is set to `xfx = x`. If some
    /// of the PDFs must be charge-conjugated, this is automatically done in this function.
    pub fn with_two(
        grid: &Grid,
        pdg1: i32,
        xfx1: &'a mut dyn FnMut(i32, f64, f64) -> f64,
        pdg2: i32,
        xfx2: &'a mut dyn FnMut(i32, f64, f64) -> f64,
        alphas: &'a mut dyn FnMut(f64) -> f64,
    ) -> Result<Self, ()> {
        Self::new(grid, pdg1, pdg2, Pdfs::Two { xfx1, xfx2 }, alphas)
    }

    /// Construct a luminosity cache with a single PDF `xfx`. The type of hadron the PDF
    /// corresponds to must be given as `pdg`. The function to evaluate the strong coupling must be
    /// given as `alphas`. The grid that the cache should be used with must be given as `grid`;
    /// this parameter determines which of the initial states are hadronic, and if an initial
    /// states is not hadronic the corresponding 'PDF' is set to `xfx = x`. If some of the PDFs
    /// must be charge-conjugated, this is automatically done in this function.
    pub fn with_one(
        grid: &Grid,
        pdg: i32,
        xfx: &'a mut dyn FnMut(i32, f64, f64) -> f64,
        alphas: &'a mut dyn FnMut(f64) -> f64,
    ) -> Result<Self, ()> {
        Self::new(grid, pdg, pdg, Pdfs::One { xfx }, alphas)
    }

    /// Return the PDF (multiplied with `x`) for the first initial state.
    pub fn xfx1(&mut self, pdg_id: i32, ix1: usize, imu2: usize) -> f64 {
        let x = self.x1_grid[ix1];
        if self.cc1 == 0 {
            x
        } else {
            let muf2 = self.mu2_grid[imu2].fac;
            match &mut self.pdfs {
                Pdfs::One { xfx } => xfx(self.cc1 * pdg_id, x, muf2),
                Pdfs::Two { xfx1, .. } => xfx1(self.cc1 * pdg_id, x, muf2),
            }
        }
    }

    /// Return the PDF (multiplied with `x`) for the second initial state.
    pub fn xfx2(&mut self, pdg_id: i32, ix2: usize, imu2: usize) -> f64 {
        let x = self.x2_grid[ix2];
        if self.cc2 == 0 {
            x
        } else {
            let muf2 = self.mu2_grid[imu2].fac;
            match &mut self.pdfs {
                Pdfs::One { xfx } => xfx(self.cc2 * pdg_id, x, muf2),
                Pdfs::Two { xfx2, .. } => xfx2(self.cc2 * pdg_id, x, muf2),
            }
        }
    }

    /// Return the strong coupling for the renormalization scale set with [`LumiCache::set_grids`],
    /// in the grid `mu2_grid` at the index `imu2`.
    pub fn alphas(&mut self, imu2: usize) -> f64 {
        (self.alphas)(self.mu2_grid[imu2].ren)
    }

    /// Clears the cache.
    pub fn clear(&mut self) {
        self.mu2_grid.clear();
        self.x1_grid.clear();
        self.x2_grid.clear();
    }

    /// Set the grids.
    pub fn set_grids(
        &mut self,
        mu2_grid: &[Mu2],
        x1_grid: &[f64],
        x2_grid: &[f64],
        xir: f64,
        xif: f64,
    ) {
        if (x1_grid != self.x1_grid) || (x2_grid != self.x2_grid) {
            self.clear();
            self.x1_grid = x1_grid.to_vec();
            self.x2_grid = x2_grid.to_vec();
        }

        self.mu2_grid = mu2_grid
            .iter()
            .map(|Mu2 { ren, fac }| Mu2 {
                ren: ren * xir * xir,
                fac: fac * xif * xif,
            })
            .collect();
    }
}
