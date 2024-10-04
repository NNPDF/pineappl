//! Module for everything related to luminosity functions.

use super::boc::Kinematics;
use super::grid::Grid;
use super::pids;
use super::subgrid::{NodeValues, Subgrid};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};

/// A cache for evaluating PDFs. Methods like [`Grid::convolve`] accept instances of this `struct`
/// instead of the PDFs themselves.
pub struct ConvolutionCache<'a> {
    xfx: Vec<&'a mut dyn FnMut(i32, f64, f64) -> f64>,
    xfx_cache: Vec<FxHashMap<(i32, usize, usize), f64>>,
    alphas: &'a mut dyn FnMut(f64) -> f64,
    alphas_cache: Vec<f64>,
    mur2_grid: Vec<f64>,
    muf2_grid: Vec<f64>,
    x_grid: Vec<f64>,
    imur2: Vec<usize>,
    imuf2: Vec<usize>,
    ix: Vec<Vec<usize>>,
    pdg: Vec<i32>,
    perm: Vec<Option<(usize, bool)>>,
}

impl<'a> ConvolutionCache<'a> {
    /// TODO
    pub fn new(
        pdg: Vec<i32>,
        xfx: Vec<&'a mut dyn FnMut(i32, f64, f64) -> f64>,
        alphas: &'a mut dyn FnMut(f64) -> f64,
    ) -> Self {
        Self {
            xfx_cache: vec![FxHashMap::default(); xfx.len()],
            xfx,
            alphas,
            alphas_cache: Vec::new(),
            mur2_grid: Vec::new(),
            muf2_grid: Vec::new(),
            x_grid: Vec::new(),
            imur2: Vec::new(),
            imuf2: Vec::new(),
            ix: Vec::new(),
            pdg,
            perm: Vec::new(),
        }
    }

    pub(crate) fn setup(&mut self, grid: &Grid, xi: &[(f64, f64)]) -> Result<(), ()> {
        let convolutions = grid.convolutions();

        // TODO: the following code only works with exactly two convolutions
        assert_eq!(convolutions.len(), 2);

        self.perm = grid
            .convolutions()
            .iter()
            .enumerate()
            .map(|(max_idx, conv)| {
                conv.pid().map(|pid| {
                    self.pdg
                        .iter()
                        .take(max_idx + 1)
                        .enumerate()
                        .rev()
                        .find_map(|(idx, &pdg)| {
                            if pid == pdg {
                                Some((idx, false))
                            } else if pid == pids::charge_conjugate_pdg_pid(pdg) {
                                Some((idx, true))
                            } else {
                                None
                            }
                        })
                        // TODO: convert `unwrap` to `Err`
                        .unwrap()
                })
            })
            .collect();

        // TODO: try to avoid calling clear
        self.clear();

        let mut x_grid: Vec<_> = grid
            .subgrids()
            .iter()
            .filter_map(|subgrid| {
                if subgrid.is_empty() {
                    None
                } else {
                    let vec: Vec<_> = grid
                        .kinematics()
                        .iter()
                        .zip(subgrid.node_values())
                        .filter(|(kin, _)| matches!(kin, Kinematics::X(_)))
                        .flat_map(|(_, node_values)| node_values.values())
                        .collect();
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
                    Some(
                        grid.kinematics()
                            .iter()
                            .zip(subgrid.node_values())
                            .find_map(|(kin, node_values)| {
                                // TODO: generalize this for arbitrary scales
                                matches!(kin, &Kinematics::Scale(idx) if idx == 0)
                                    .then_some(node_values)
                            })
                            // TODO: convert this into an error
                            .unwrap()
                            .values(),
                    )
                }
            })
            .flatten()
            .flat_map(|ren| {
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
                    Some(
                        grid.kinematics()
                            .iter()
                            .zip(subgrid.node_values())
                            .find_map(|(kin, node_values)| {
                                // TODO: generalize this for arbitrary scales
                                matches!(kin, &Kinematics::Scale(idx) if idx == 0)
                                    .then_some(node_values)
                            })
                            // TODO: convert this into an error
                            .unwrap()
                            .values(),
                    )
                }
            })
            .flatten()
            .flat_map(|fac| {
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

        Ok(())
    }

    /// TODO
    pub fn fx_prod(&mut self, pdg_ids: &[i32], indices: &[usize]) -> f64 {
        // TODO: here we assume that
        // - indices[0] is the (squared) factorization scale,
        // - indices[1] is x1 and
        // - indices[2] is x2.
        // Lift this restriction!
        self.perm
            .iter()
            .zip(pdg_ids)
            .enumerate()
            .filter_map(|(index, (perm, &pdg_id))| {
                perm.map(|(idx, cc)| {
                    let ix = self.ix[index][indices[index + 1]];
                    let imuf2 = self.imuf2[indices[0]];
                    let muf2 = self.muf2_grid[imuf2];
                    let pid = if cc {
                        pids::charge_conjugate_pdg_pid(pdg_id)
                    } else {
                        pdg_id
                    };
                    let xfx = &mut self.xfx[idx];
                    let xfx_cache = &mut self.xfx_cache[idx];
                    *xfx_cache.entry((pid, ix, imuf2)).or_insert_with(|| {
                        let x = self.x_grid[ix];
                        xfx(pid, x, muf2) / x
                    })
                })
            })
            .product()
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
        for xfx_cache in &mut self.xfx_cache {
            xfx_cache.clear();
        }
        self.mur2_grid.clear();
        self.muf2_grid.clear();
        self.x_grid.clear();
    }

    /// Set the grids.
    pub fn set_grids(
        &mut self,
        grid: &Grid,
        node_values: &[NodeValues],
        mu2_grid: &[f64],
        xir: f64,
        xif: f64,
        xia: f64,
    ) {
        // TODO: generalize this for fragmentation functions
        assert_eq!(xia, 1.0);

        self.imur2 = mu2_grid
            .iter()
            .map(|ren| {
                self.mur2_grid
                    .iter()
                    .position(|&mur2| mur2 == xir * xir * ren)
                    .unwrap_or_else(|| unreachable!())
            })
            .collect();
        self.imuf2 = mu2_grid
            .iter()
            .map(|fac| {
                self.muf2_grid
                    .iter()
                    .position(|&muf2| muf2 == xif * xif * fac)
                    .unwrap_or_else(|| unreachable!())
            })
            .collect();

        // TODO: generalize this for arbitrary orderings of x
        self.ix = (0..grid.convolutions().len())
            .map(|idx| {
                grid.kinematics()
                    .iter()
                    .zip(node_values)
                    .find_map(|(kin, node_values)| {
                        matches!(kin, &Kinematics::X(index) if index == idx).then_some(node_values)
                    })
                    // UNWRAP: guaranteed by the grid constructor
                    .unwrap()
                    .values()
                    .iter()
                    .map(|xd| {
                        self.x_grid
                            .iter()
                            .position(|x| xd == x)
                            .unwrap_or_else(|| unreachable!())
                    })
                    .collect()
            })
            .collect();
    }
}

/// Data type that indentifies different types of convolutions.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
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
