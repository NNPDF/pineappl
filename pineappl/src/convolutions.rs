//! Module for everything related to luminosity functions.

use super::grid::Grid;
use super::pids;
use super::subgrid::{Mu2, Subgrid};
use rustc_hash::FxHashMap;

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

/// A cache for evaluating PDFs. Methods like [`Grid::convolve`] accept instances of this `struct`
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
        let convolutions = grid.convolutions();

        // TODO: the following code only works with exactly two convolutions
        assert_eq!(convolutions.len(), 2);

        // do we have to charge-conjugate the initial states?
        let cc1 = if let Some(pid) = convolutions[0].pid() {
            if self.pdg1 == pid {
                1
            } else if self.pdg1 == pids::charge_conjugate_pdg_pid(pid) {
                -1
            } else {
                // TODO: return a proper error
                return Err(());
            }
        } else {
            0
        };
        let cc2 = if let Some(pid) = convolutions[1].pid() {
            if self.pdg2 == pid {
                1
            } else if self.pdg2 == pids::charge_conjugate_pdg_pid(pid) {
                -1
            } else {
                // TODO: return a proper error
                return Err(());
            }
        } else {
            0
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
    #[must_use]
    pub fn alphas(&self, imu2: usize) -> f64 {
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

/// Data type that indentifies different types of convolutions.
#[derive(Debug, Eq, PartialEq)]
pub enum Convolution {
    // TODO: eventually get rid of this value
    /// No convolution.
    None,
    /// Unpolarized parton distribution function. The integer denotes the type of hadron with a PDG
    /// MC ID.
    UnpolPDF(i32),
    /// Polarized parton distribution function. The integer denotes the type of hadron with a PDG
    /// MC ID.
    PolPDF(i32),
    /// Unpolarized fragmentation function. The integer denotes the type of hadron with a PDG MC
    /// ID.
    UnpolFF(i32),
    /// Polarized fragmentation function. The integer denotes the type of hadron with a PDG MC ID.
    PolFF(i32),
}

impl Convolution {
    /// Return the convolution if the PID is charged conjugated.
    #[must_use]
    pub const fn charge_conjugate(&self) -> Self {
        match *self {
            Self::None => Self::None,
            Self::UnpolPDF(pid) => Self::UnpolPDF(pids::charge_conjugate_pdg_pid(pid)),
            Self::PolPDF(pid) => Self::PolPDF(pids::charge_conjugate_pdg_pid(pid)),
            Self::UnpolFF(pid) => Self::UnpolFF(pids::charge_conjugate_pdg_pid(pid)),
            Self::PolFF(pid) => Self::PolFF(pids::charge_conjugate_pdg_pid(pid)),
        }
    }

    /// Return the PID of the convolution if it has any.
    #[must_use]
    pub const fn pid(&self) -> Option<i32> {
        match *self {
            Self::None => None,
            Self::UnpolPDF(pid) | Self::PolPDF(pid) | Self::UnpolFF(pid) | Self::PolFF(pid) => {
                Some(pid)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn convolution_charge_conjugate() {
        assert_eq!(Convolution::None.charge_conjugate(), Convolution::None);
        assert_eq!(
            Convolution::UnpolPDF(2212).charge_conjugate(),
            Convolution::UnpolPDF(-2212)
        );
        assert_eq!(
            Convolution::PolPDF(2212).charge_conjugate(),
            Convolution::PolPDF(-2212)
        );
        assert_eq!(
            Convolution::UnpolFF(2212).charge_conjugate(),
            Convolution::UnpolFF(-2212)
        );
        assert_eq!(
            Convolution::PolFF(2212).charge_conjugate(),
            Convolution::PolFF(-2212)
        );
    }

    #[test]
    fn convolution_pid() {
        assert_eq!(Convolution::None.pid(), None);
        assert_eq!(Convolution::UnpolPDF(2212).pid(), Some(2212));
        assert_eq!(Convolution::PolPDF(2212).pid(), Some(2212));
        assert_eq!(Convolution::UnpolFF(2212).pid(), Some(2212));
        assert_eq!(Convolution::PolFF(2212).pid(), Some(2212));
    }
}
