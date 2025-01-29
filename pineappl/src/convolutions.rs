//! Module for everything related to convolution functions.

use super::boc::Kinematics;
use super::boc::Scales;
use super::grid::Grid;
use super::pids;
use super::subgrid::{self, Subgrid, SubgridEnum};
use itertools::izip;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};

const REN_IDX: usize = 0;
const FAC_IDX: usize = 1;
const FRG_IDX: usize = 2;
const SCALES_CNT: usize = 3;

struct ConvCache1d<'a> {
    xfx: &'a mut dyn FnMut(i32, f64, f64) -> f64,
    cache: FxHashMap<(i32, usize, usize), f64>,
    conv: Conv,
}

/// A cache for evaluating PDFs. Methods like [`Grid::convolve`] accept instances of this `struct`
/// instead of the PDFs themselves.
pub struct ConvolutionCache<'a> {
    caches: Vec<ConvCache1d<'a>>,
    alphas: &'a mut dyn FnMut(f64) -> f64,
    alphas_cache: Vec<f64>,
    mu2: [Vec<f64>; SCALES_CNT],
    x_grid: Vec<f64>,
}

impl<'a> ConvolutionCache<'a> {
    /// TODO
    pub fn new(
        convolutions: Vec<Conv>,
        xfx: Vec<&'a mut dyn FnMut(i32, f64, f64) -> f64>,
        alphas: &'a mut dyn FnMut(f64) -> f64,
    ) -> Self {
        Self {
            caches: xfx
                .into_iter()
                .zip(convolutions)
                .map(|(xfx, conv)| ConvCache1d {
                    xfx,
                    cache: FxHashMap::default(),
                    conv,
                })
                .collect(),
            alphas,
            alphas_cache: Vec::new(),
            mu2: [const { Vec::new() }; SCALES_CNT],
            x_grid: Vec::new(),
        }
    }

    pub(crate) fn new_grid_conv_cache<'b>(
        &'b mut self,
        grid: &Grid,
        xi: &[(f64, f64, f64)],
    ) -> GridConvCache<'a, 'b> {
        // TODO: try to avoid calling clear
        self.clear();

        let scales: [_; SCALES_CNT] = grid.scales().into();
        let xi: Vec<_> = (0..SCALES_CNT)
            .map(|idx| {
                let mut vars: Vec<_> = xi
                    .iter()
                    .map(|&x| <[_; SCALES_CNT]>::from(x)[idx])
                    .collect();
                vars.sort_by(f64::total_cmp);
                vars.dedup();
                vars
            })
            .collect();

        for (result, scale, xi) in izip!(&mut self.mu2, scales, xi) {
            result.clear();
            result.extend(
                grid.subgrids()
                    .iter()
                    .filter(|subgrid| !subgrid.is_empty())
                    .flat_map(|subgrid| {
                        scale
                            .calc(&subgrid.node_values(), grid.kinematics())
                            .into_owned()
                    })
                    .flat_map(|scale| xi.iter().map(move |&xi| xi * xi * scale)),
            );
            result.sort_by(f64::total_cmp);
            result.dedup();
        }

        let mut x_grid: Vec<_> = grid
            .subgrids()
            .iter()
            .filter(|subgrid| !subgrid.is_empty())
            .flat_map(|subgrid| {
                grid.kinematics()
                    .iter()
                    .zip(subgrid.node_values())
                    .filter(|(kin, _)| matches!(kin, Kinematics::X(_)))
                    .flat_map(|(_, node_values)| node_values)
            })
            .collect();
        x_grid.sort_by(f64::total_cmp);
        x_grid.dedup();

        self.alphas_cache = self.mu2[REN_IDX]
            .iter()
            .map(|&mur2| (self.alphas)(mur2))
            .collect();
        self.x_grid = x_grid;

        let perm = grid
            .convolutions()
            .iter()
            .enumerate()
            .map(|(max_idx, grid_conv)| {
                self.caches
                    .iter()
                    .take(max_idx + 1)
                    .enumerate()
                    .rev()
                    .find_map(|(idx, ConvCache1d { conv, .. })| {
                        if grid_conv == conv {
                            Some((idx, false))
                        } else if *grid_conv == conv.cc() {
                            Some((idx, true))
                        } else {
                            None
                        }
                    })
                    // TODO: convert `unwrap` to `Err`
                    .unwrap()
            })
            .collect();

        GridConvCache {
            cache: self,
            perm,
            imu2: [const { Vec::new() }; SCALES_CNT],
            scales: grid.scales().clone(),
            ix: Vec::new(),
            scale_dims: Vec::new(),
        }
    }

    /// Clears the cache.
    pub fn clear(&mut self) {
        self.alphas_cache.clear();
        for xfx_cache in &mut self.caches {
            xfx_cache.cache.clear();
        }
        for scales in &mut self.mu2 {
            scales.clear();
        }
        self.x_grid.clear();
    }
}

/// TODO
pub struct GridConvCache<'a, 'b> {
    cache: &'b mut ConvolutionCache<'a>,
    perm: Vec<(usize, bool)>,
    imu2: [Vec<usize>; SCALES_CNT],
    scales: Scales,
    ix: Vec<Vec<usize>>,
    scale_dims: Vec<usize>,
}

impl GridConvCache<'_, '_> {
    /// TODO
    pub fn as_fx_prod(&mut self, pdg_ids: &[i32], as_order: u8, indices: &[usize]) -> f64 {
        // TODO: here we assume that
        // - indices[0] is the (squared) factorization scale,
        // - indices[1] is x1 and
        // - indices[2] is x2.
        // Lift this restriction!
        let x_start = indices.len() - pdg_ids.len();
        let indices_scales = &indices[0..x_start];
        let indices_x = &indices[x_start..];

        let ix = self.ix.iter().zip(indices_x).map(|(ix, &index)| ix[index]);
        let idx_pid = self.perm.iter().zip(pdg_ids).map(|(&(idx, cc), &pdg_id)| {
            (
                idx,
                if cc {
                    pids::charge_conjugate_pdg_pid(pdg_id)
                } else {
                    pdg_id
                },
            )
        });

        let fx_prod: f64 = ix
            .zip(idx_pid)
            .map(|(ix, (idx, pid))| {
                let ConvCache1d { xfx, cache, conv } = &mut self.cache.caches[idx];

                let (scale, scale_idx) = match conv.conv_type() {
                    ConvType::UnpolPDF | ConvType::PolPDF => (
                        FAC_IDX,
                        self.scales.fac.idx(indices_scales, &self.scale_dims),
                    ),
                    ConvType::UnpolFF | ConvType::PolFF => (
                        FRG_IDX,
                        self.scales.frg.idx(indices_scales, &self.scale_dims),
                    ),
                };

                let imu2 = self.imu2[scale][scale_idx];
                let mu2 = self.cache.mu2[scale][imu2];

                *cache.entry((pid, ix, imu2)).or_insert_with(|| {
                    let x = self.cache.x_grid[ix];
                    xfx(pid, x, mu2) / x
                })
            })
            .product();
        let alphas_powers = if as_order != 0 {
            let ren_scale_idx = self.scales.ren.idx(indices_scales, &self.scale_dims);
            self.cache.alphas_cache[self.imu2[REN_IDX][ren_scale_idx]].powi(as_order.into())
        } else {
            1.0
        };

        fx_prod * alphas_powers
    }

    /// Set the grids.
    pub fn set_grids(&mut self, grid: &Grid, subgrid: &SubgridEnum, xi: (f64, f64, f64)) {
        let node_values = subgrid.node_values();
        let kinematics = grid.kinematics();
        let scales: [_; SCALES_CNT] = grid.scales().into();
        let xi: [_; SCALES_CNT] = xi.into();

        for (result, values, scale, xi) in izip!(&mut self.imu2, &self.cache.mu2, scales, xi) {
            result.clear();
            result.extend(scale.calc(&node_values, kinematics).iter().map(|s| {
                values
                    .iter()
                    .position(|&value| subgrid::node_value_eq(value, xi * xi * s))
                    // UNWRAP: if this fails, `new_grid_conv_cache` hasn't been called properly
                    .unwrap_or_else(|| unreachable!())
            }));
        }

        self.ix = (0..grid.convolutions().len())
            .map(|idx| {
                kinematics
                    .iter()
                    .zip(&node_values)
                    .find_map(|(kin, node_values)| {
                        matches!(kin, &Kinematics::X(index) if index == idx).then_some(node_values)
                    })
                    // UNWRAP: guaranteed by the grid constructor
                    .unwrap_or_else(|| unreachable!())
                    .iter()
                    .map(|&xd| {
                        self.cache
                            .x_grid
                            .iter()
                            .position(|&x| subgrid::node_value_eq(xd, x))
                            .unwrap_or_else(|| unreachable!())
                    })
                    .collect()
            })
            .collect();

        self.scale_dims = grid
            .kinematics()
            .iter()
            .zip(node_values)
            .filter_map(|(kin, node_values)| {
                matches!(kin, Kinematics::Scale(_)).then_some(node_values.len())
            })
            .collect();
    }
}

/// TODO
#[repr(C)]
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

    /// TODO
    #[must_use]
    pub const fn is_pdf(&self) -> bool {
        matches!(self, Self::UnpolPDF | Self::PolPDF)
    }
}

/// Data type that indentifies different types of convolutions.
#[repr(C)]
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
