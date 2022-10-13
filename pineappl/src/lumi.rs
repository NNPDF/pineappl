//! Module for everything related to luminosity functions.

use super::grid::Grid;
use super::pids;
use super::subgrid::{Mu2, Subgrid};
use itertools::Itertools;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};

/// This structure represents an entry of a luminosity function. Each entry consists of a tuple,
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
    /// Same arguments are merged together:
    ///
    /// ```rust
    /// use pineappl::lumi::LumiEntry;
    ///
    /// let entry1 = LumiEntry::new(vec![(1, 1, 1.0), (1, 1, 3.0), (3, 3, 1.0), (1, 1, 6.0)]);
    /// let entry2 = LumiEntry::new(vec![(1, 1, 10.0), (3, 3, 1.0)]);
    ///
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
                .collect(),
        }
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

/// Helper macro to quickly generate a `LumiEntry` at compile time.
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
        xfx1_cache: FxHashMap<(i32, usize, usize), f64>,
        xfx2: &'a mut dyn FnMut(i32, f64, f64) -> f64,
        xfx2_cache: FxHashMap<(i32, usize, usize), f64>,
    },
    One {
        xfx: &'a mut dyn FnMut(i32, f64, f64) -> f64,
        xfx_cache: FxHashMap<(i32, usize, usize), f64>,
    },
}

impl<'a> Pdfs<'a> {
    pub fn clear(&mut self) {
        match self {
            Self::One { xfx_cache, .. } => xfx_cache.clear(),
            Self::Two {
                xfx1_cache,
                xfx2_cache,
                ..
            } => {
                xfx1_cache.clear();
                xfx2_cache.clear();
            }
        }
    }
}

/// A cache for evaluating PDFs. Methods like [`Grid::convolute`] accept instances of this `struct`
/// instead of the PDFs themselves.
pub struct LumiCache<'a> {
    pdfs: Pdfs<'a>,
    alphas: &'a mut dyn FnMut(f64) -> f64,
    alphas_cache: Vec<f64>,
    mur2_grid: Vec<f64>,
    muf2_grid: Vec<f64>,
    x_grid: Vec<f64>,
    imur2: Vec<usize>,
    imuf2: Vec<usize>,
    ix1: Vec<usize>,
    ix2: Vec<usize>,
    pdg1: i32,
    pdg2: i32,
    cc1: i32,
    cc2: i32,
}

impl<'a> LumiCache<'a> {
    /// Construct a luminosity cache with two PDFs, `xfx1` and `xfx2`. The types of hadrons the
    /// PDFs correspond to must be given as `pdg1` and `pdg2`. The function to evaluate the
    /// strong coupling must be given as `alphas`. The grid that the cache will be used with must
    /// be given as `grid`; this parameter determines which of the initial states are hadronic, and
    /// if an initial states is not hadronic the corresponding 'PDF' is set to `xfx = x`. If some
    /// of the PDFs must be charge-conjugated, this is automatically done in this function.
    pub fn with_two(
        pdg1: i32,
        xfx1: &'a mut dyn FnMut(i32, f64, f64) -> f64,
        pdg2: i32,
        xfx2: &'a mut dyn FnMut(i32, f64, f64) -> f64,
        alphas: &'a mut dyn FnMut(f64) -> f64,
    ) -> Self {
        Self {
            pdfs: Pdfs::Two {
                xfx1,
                xfx1_cache: FxHashMap::default(),
                xfx2,
                xfx2_cache: FxHashMap::default(),
            },
            alphas,
            alphas_cache: vec![],
            mur2_grid: vec![],
            muf2_grid: vec![],
            x_grid: vec![],
            imur2: Vec::new(),
            imuf2: Vec::new(),
            ix1: Vec::new(),
            ix2: Vec::new(),
            pdg1,
            pdg2,
            cc1: 0,
            cc2: 0,
        }
    }

    /// Construct a luminosity cache with a single PDF `xfx`. The type of hadron the PDF
    /// corresponds to must be given as `pdg`. The function to evaluate the strong coupling must be
    /// given as `alphas`. The grid that the cache should be used with must be given as `grid`;
    /// this parameter determines which of the initial states are hadronic, and if an initial
    /// states is not hadronic the corresponding 'PDF' is set to `xfx = x`. If some of the PDFs
    /// must be charge-conjugated, this is automatically done in this function.
    pub fn with_one(
        pdg: i32,
        xfx: &'a mut dyn FnMut(i32, f64, f64) -> f64,
        alphas: &'a mut dyn FnMut(f64) -> f64,
    ) -> Self {
        Self {
            pdfs: Pdfs::One {
                xfx,
                xfx_cache: FxHashMap::default(),
            },
            alphas,
            alphas_cache: vec![],
            mur2_grid: vec![],
            muf2_grid: vec![],
            x_grid: vec![],
            imur2: Vec::new(),
            imuf2: Vec::new(),
            ix1: Vec::new(),
            ix2: Vec::new(),
            pdg1: pdg,
            pdg2: pdg,
            cc1: 0,
            cc2: 0,
        }
    }

    pub(crate) fn setup(&mut self, grid: &Grid, xi: &[(f64, f64)]) -> Result<(), ()> {
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
        let has_pdfa = !grid
            .lumi()
            .iter()
            .all(|entry| entry.entry().iter().all(|&(a, _, _)| a == pdga));
        let has_pdfb = !grid
            .lumi()
            .iter()
            .all(|entry| entry.entry().iter().all(|&(_, b, _)| b == pdgb));

        // do we have to charge-conjugate the initial states?
        let cc1 = if !has_pdfa {
            0
        } else if self.pdg1 == pdga {
            1
        } else if self.pdg1 == -pdga {
            -1
        } else {
            return Err(());
        };
        let cc2 = if !has_pdfb {
            0
        } else if self.pdg2 == pdgb {
            1
        } else if self.pdg2 == -pdgb {
            -1
        } else {
            return Err(());
        };

        // TODO: try to avoid calling clear
        self.clear();

        let mut x_grid: Vec<_> = grid
            .subgrids()
            .iter()
            .filter_map(|subgrid| {
                if subgrid.is_empty() {
                    None
                } else {
                    let mut vec = subgrid.x1_grid().into_owned();
                    vec.extend_from_slice(&subgrid.x2_grid());
                    Some(vec)
                }
            })
            .flatten()
            .collect();
        x_grid.sort_by(|a, b| a.partial_cmp(b).unwrap_or_else(|| unreachable!()));
        x_grid.dedup();

        let mut mur2_grid: Vec<_> = grid
            .subgrids()
            .iter()
            .filter_map(|subgrid| {
                if subgrid.is_empty() {
                    None
                } else {
                    Some(subgrid.mu2_grid().into_owned())
                }
            })
            .flatten()
            .flat_map(|Mu2 { ren, .. }| {
                xi.iter()
                    .map(|(xir, _)| xir * xir * ren)
                    .collect::<Vec<_>>()
            })
            .collect();
        mur2_grid.sort_by(|a, b| a.partial_cmp(b).unwrap_or_else(|| unreachable!()));
        mur2_grid.dedup();

        let mut muf2_grid: Vec<_> = grid
            .subgrids()
            .iter()
            .filter_map(|subgrid| {
                if subgrid.is_empty() {
                    None
                } else {
                    Some(subgrid.mu2_grid().into_owned())
                }
            })
            .flatten()
            .flat_map(|Mu2 { fac, .. }| {
                xi.iter()
                    .map(|(_, xif)| xif * xif * fac)
                    .collect::<Vec<_>>()
            })
            .collect();
        muf2_grid.sort_by(|a, b| a.partial_cmp(b).unwrap_or_else(|| unreachable!()));
        muf2_grid.dedup();

        self.alphas_cache = mur2_grid.iter().map(|&mur2| (self.alphas)(mur2)).collect();

        self.mur2_grid = mur2_grid;
        self.muf2_grid = muf2_grid;
        self.x_grid = x_grid;
        self.cc1 = cc1;
        self.cc2 = cc2;

        Ok(())
    }

    /// Return the PDF (multiplied with `x`) for the first initial state.
    pub fn xfx1(&mut self, pdg_id: i32, ix1: usize, imu2: usize) -> f64 {
        let ix1 = self.ix1[ix1];
        let x = self.x_grid[ix1];
        if self.cc1 == 0 {
            x
        } else {
            let imuf2 = self.imuf2[imu2];
            let muf2 = self.muf2_grid[imuf2];
            let pid = if self.cc1 == 1 {
                pdg_id
            } else {
                pids::charge_conjugate_pdg_pid(pdg_id)
            };
            let (xfx, xfx_cache) = match &mut self.pdfs {
                Pdfs::One { xfx, xfx_cache, .. } => (xfx, xfx_cache),
                Pdfs::Two {
                    xfx1, xfx1_cache, ..
                } => (xfx1, xfx1_cache),
            };
            *xfx_cache
                .entry((pid, ix1, imuf2))
                .or_insert_with(|| xfx(pid, x, muf2))
        }
    }

    /// Return the PDF (multiplied with `x`) for the second initial state.
    pub fn xfx2(&mut self, pdg_id: i32, ix2: usize, imu2: usize) -> f64 {
        let ix2 = self.ix2[ix2];
        let x = self.x_grid[ix2];
        if self.cc2 == 0 {
            x
        } else {
            let imuf2 = self.imuf2[imu2];
            let muf2 = self.muf2_grid[imuf2];
            let pid = if self.cc2 == 1 {
                pdg_id
            } else {
                pids::charge_conjugate_pdg_pid(pdg_id)
            };
            let (xfx, xfx_cache) = match &mut self.pdfs {
                Pdfs::One { xfx, xfx_cache, .. } => (xfx, xfx_cache),
                Pdfs::Two {
                    xfx2, xfx2_cache, ..
                } => (xfx2, xfx2_cache),
            };
            *xfx_cache
                .entry((pid, ix2, imuf2))
                .or_insert_with(|| xfx(pid, x, muf2))
        }
    }

    /// Return the strong coupling for the renormalization scale set with [`LumiCache::set_grids`],
    /// in the grid `mu2_grid` at the index `imu2`.
    pub fn alphas(&mut self, imu2: usize) -> f64 {
        self.alphas_cache[self.imur2[imu2]]
    }

    /// Clears the cache.
    pub fn clear(&mut self) {
        self.alphas_cache.clear();
        self.pdfs.clear();
        self.mur2_grid.clear();
        self.muf2_grid.clear();
        self.x_grid.clear();
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
        self.imur2 = mu2_grid
            .iter()
            .map(|Mu2 { ren, .. }| {
                self.mur2_grid
                    .iter()
                    .position(|&mur2| mur2 == xir * xir * ren)
                    .unwrap_or_else(|| unreachable!())
            })
            .collect();
        self.imuf2 = mu2_grid
            .iter()
            .map(|Mu2 { fac, .. }| {
                self.muf2_grid
                    .iter()
                    .position(|&muf2| muf2 == xif * xif * fac)
                    .unwrap_or_else(|| unreachable!())
            })
            .collect();
        self.ix1 = x1_grid
            .iter()
            .map(|x1| {
                self.x_grid
                    .iter()
                    .position(|x| x1 == x)
                    .unwrap_or_else(|| unreachable!())
            })
            .collect();

        self.ix2 = x2_grid
            .iter()
            .map(|x2| {
                self.x_grid
                    .iter()
                    .position(|x| x2 == x)
                    .unwrap_or_else(|| unreachable!())
            })
            .collect();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pids;

    #[test]
    fn translate() {
        let lumi = LumiEntry::translate(&lumi_entry![103, 203, 2.0], &pids::evol_to_pdg_mc_ids);

        assert_eq!(
            lumi,
            lumi_entry![ 2,  2,  2.0;  2, -2, -2.0;  2,  1, -2.0;  2, -1,  2.0;
                        -2,  2,  2.0; -2, -2, -2.0; -2,  1, -2.0; -2, -1,  2.0;
                         1,  2, -2.0;  1, -2,  2.0;  1,  1,  2.0;  1, -1, -2.0;
                        -1,  2, -2.0; -1, -2,  2.0; -1,  1,  2.0; -1, -1, -2.0]
        );
    }
}
