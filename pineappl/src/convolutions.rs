//! Module for everything related to convolution functions.

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
    mua2_grid: Vec<f64>,
    x_grid: Vec<f64>,
    imur2: Vec<usize>,
    imuf2: Vec<usize>,
    imua2: Vec<usize>,
    ix: Vec<Vec<usize>>,
    pdg: Vec<Conv>,
    perm: Vec<(usize, bool)>,
}

impl<'a> ConvolutionCache<'a> {
    /// TODO
    pub fn new(
        pdg: Vec<Conv>,
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
            mua2_grid: Vec::new(),
            x_grid: Vec::new(),
            imur2: Vec::new(),
            imuf2: Vec::new(),
            imua2: Vec::new(),
            ix: Vec::new(),
            pdg,
            perm: Vec::new(),
        }
    }

    pub(crate) fn setup(&mut self, grid: &Grid, xi: &[(f64, f64, f64)]) -> Result<(), ()> {
        self.perm = grid
            .convolutions()
            .iter()
            .enumerate()
            .map(|(max_idx, conv)| {
                self.pdg
                    .iter()
                    .take(max_idx + 1)
                    .enumerate()
                    .rev()
                    .find_map(|(idx, pdg)| {
                        if conv == pdg {
                            Some((idx, false))
                        } else if *conv == pdg.cc() {
                            Some((idx, true))
                        } else {
                            None
                        }
                    })
                    // TODO: convert `unwrap` to `Err`
                    .unwrap()
            })
            .collect();

        // TODO: try to avoid calling clear
        self.clear();

        let mut x_grid: Vec<_> = grid
            .subgrids()
            .iter()
            .filter(|subgrid| !subgrid.is_empty())
            .flat_map(|subgrid| {
                grid.kinematics()
                    .iter()
                    .zip(subgrid.node_values())
                    .filter(|(kin, _)| matches!(kin, Kinematics::X(_)))
                    .flat_map(|(_, node_values)| node_values.values())
            })
            .collect();
        x_grid.sort_by(|a, b| a.partial_cmp(b).unwrap_or_else(|| unreachable!()));
        x_grid.dedup();

        let mut mur2_grid: Vec<_> = grid
            .subgrids()
            .iter()
            .filter(|subgrid| !subgrid.is_empty())
            .flat_map(|subgrid| {
                grid.scales()
                    .ren
                    .calc(&subgrid.node_values(), grid.kinematics())
                    .unwrap_or_default()
            })
            .flat_map(|ren| xi.iter().map(move |(xir, _, _)| xir * xir * ren))
            .collect();
        mur2_grid.sort_by(|a, b| a.partial_cmp(b).unwrap_or_else(|| unreachable!()));
        mur2_grid.dedup();

        let mut muf2_grid: Vec<_> = grid
            .subgrids()
            .iter()
            .filter(|subgrid| !subgrid.is_empty())
            .flat_map(|subgrid| {
                grid.scales()
                    .fac
                    .calc(&subgrid.node_values(), grid.kinematics())
                    .unwrap_or_default()
            })
            .flat_map(|fac| xi.iter().map(move |(_, xif, _)| xif * xif * fac))
            .collect();
        muf2_grid.sort_by(|a, b| a.partial_cmp(b).unwrap_or_else(|| unreachable!()));
        muf2_grid.dedup();

        let mut mua2_grid: Vec<_> = grid
            .subgrids()
            .iter()
            .filter(|subgrid| !subgrid.is_empty())
            .flat_map(|subgrid| {
                grid.scales()
                    .frg
                    .calc(&subgrid.node_values(), grid.kinematics())
                    .unwrap_or_default()
            })
            .flat_map(|frg| xi.iter().map(move |(_, _, xia)| xia * xia * frg))
            .collect();
        mua2_grid.sort_by(|a, b| a.partial_cmp(b).unwrap_or_else(|| unreachable!()));
        mua2_grid.dedup();

        self.alphas_cache = mur2_grid.iter().map(|&mur2| (self.alphas)(mur2)).collect();
        self.mur2_grid = mur2_grid;
        self.muf2_grid = muf2_grid;
        self.mua2_grid = mua2_grid;
        self.x_grid = x_grid;

        Ok(())
    }

    /// TODO
    pub fn as_fx_prod(&mut self, pdg_ids: &[i32], as_order: u8, indices: &[usize]) -> f64 {
        // TODO: here we assume that
        // - indices[0] is the (squared) factorization scale,
        // - indices[1] is x1 and
        // - indices[2] is x2.
        // Lift this restriction!
        let fx_prod: f64 = self
            .perm
            .iter()
            .zip(pdg_ids)
            .enumerate()
            .map(|(index, (&(idx, cc), &pdg_id))| {
                let ix = self.ix[index][indices[index + 1]];

                let pid = if cc {
                    pids::charge_conjugate_pdg_pid(pdg_id)
                } else {
                    pdg_id
                };
                let xfx = &mut self.xfx[idx];
                let xfx_cache = &mut self.xfx_cache[idx];
                let (imu2, mu2) = match self.pdg[idx].conv_type() {
                    ConvType::UnpolPDF | ConvType::PolPDF => {
                        let imuf2 = self.imuf2[indices[0]];
                        (imuf2, self.muf2_grid[imuf2])
                    }
                    ConvType::UnpolFF | ConvType::PolFF => {
                        let imua2 = self.imua2[indices[0]];
                        (imua2, self.mua2_grid[imua2])
                    }
                };
                *xfx_cache.entry((pid, ix, imu2)).or_insert_with(|| {
                    let x = self.x_grid[ix];
                    xfx(pid, x, mu2) / x
                })
            })
            .product();
        let alphas_powers = if as_order != 0 {
            self.alphas_cache[self.imur2[indices[0]]].powi(as_order.into())
        } else {
            1.0
        };

        fx_prod * alphas_powers
    }

    /// Clears the cache.
    pub fn clear(&mut self) {
        self.alphas_cache.clear();
        for xfx_cache in &mut self.xfx_cache {
            xfx_cache.clear();
        }
        self.mur2_grid.clear();
        self.muf2_grid.clear();
        self.mua2_grid.clear();
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
        self.imur2 = mu2_grid
            .iter()
            .map(|ren| {
                self.mur2_grid
                    .iter()
                    .position(|&mur2| mur2 == xir * xir * ren)
            })
            .collect::<Option<_>>()
            // if we didn't find a single renormalization scale, we assume we don't need any
            // renormalization scale
            .unwrap_or_default();
        self.imuf2 = mu2_grid
            .iter()
            .map(|fac| {
                self.muf2_grid
                    .iter()
                    .position(|&muf2| muf2 == xif * xif * fac)
            })
            .collect::<Option<_>>()
            // if we didn't find a single factorization scale, we assume we don't need any
            // factorization scale
            .unwrap_or_default();
        self.imua2 = mu2_grid
            .iter()
            .map(|frg| {
                self.mua2_grid
                    .iter()
                    .position(|&mua2| mua2 == xia * xia * frg)
            })
            .collect::<Option<_>>()
            // if we didn't find a single fragmentation scale, we assume we don't need any
            // fragmentation scale
            .unwrap_or_default();

        self.ix = (0..grid.convolutions().len())
            .map(|idx| {
                grid.kinematics()
                    .iter()
                    .zip(node_values)
                    .find_map(|(kin, node_values)| {
                        matches!(kin, &Kinematics::X(index) if index == idx).then_some(node_values)
                    })
                    // UNWRAP: guaranteed by the grid constructor
                    .unwrap_or_else(|| unreachable!())
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

/// TODO
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub enum ConvType {
    /// Unpolarized parton distribution function.
    UnpolPDF,
    /// Polarized parton distribution function.
    PolPDF,
    /// Unpolarized fragmentation function.
    UnpolFF,
    /// Polarized fragmentation function.
    PolFF,
}

impl ConvType {
    /// TODO
    #[must_use]
    pub const fn new(polarized: bool, time_like: bool) -> Self {
        match (polarized, time_like) {
            (false, false) => Self::UnpolPDF,
            (false, true) => Self::UnpolFF,
            (true, false) => Self::PolPDF,
            (true, true) => Self::PolFF,
        }
    }
}

/// Data type that indentifies different types of convolutions.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct Conv {
    conv_type: ConvType,
    pid: i32,
}

impl Conv {
    /// Constructor.
    #[must_use]
    pub const fn new(conv_type: ConvType, pid: i32) -> Self {
        Self { conv_type, pid }
    }

    /// TODO
    #[must_use]
    pub const fn with_pid(&self, pid: i32) -> Self {
        Self {
            conv_type: self.conv_type,
            pid,
        }
    }

    /// Return the convolution if the PID is charged conjugated.
    #[must_use]
    pub const fn cc(&self) -> Self {
        Self {
            conv_type: self.conv_type,
            pid: pids::charge_conjugate_pdg_pid(self.pid),
        }
    }

    /// Return the PID of the convolution.
    #[must_use]
    pub const fn pid(&self) -> i32 {
        self.pid
    }

    /// Return the convolution type of this convolution.
    #[must_use]
    pub const fn conv_type(&self) -> ConvType {
        self.conv_type
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conv_cc() {
        assert_eq!(
            Conv::new(ConvType::UnpolPDF, 2212).cc(),
            Conv::new(ConvType::UnpolPDF, -2212)
        );
        assert_eq!(
            Conv::new(ConvType::PolPDF, 2212).cc(),
            Conv::new(ConvType::PolPDF, -2212)
        );
        assert_eq!(
            Conv::new(ConvType::UnpolFF, 2212).cc(),
            Conv::new(ConvType::UnpolFF, -2212)
        );
        assert_eq!(
            Conv::new(ConvType::PolFF, 2212).cc(),
            Conv::new(ConvType::PolFF, -2212)
        );
    }

    #[test]
    fn conv_pid() {
        assert_eq!(Conv::new(ConvType::UnpolPDF, 2212).pid(), 2212);
        assert_eq!(Conv::new(ConvType::PolPDF, 2212).pid(), 2212);
        assert_eq!(Conv::new(ConvType::UnpolFF, 2212).pid(), 2212);
        assert_eq!(Conv::new(ConvType::PolFF, 2212).pid(), 2212);
    }
}
