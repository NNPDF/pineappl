//! Module containing all traits and supporting structures for grids.

use super::bin::BinLimits;
use super::lumi::Lumi;
use super::ntuple_grid::NtupleSubgrid;
use serde::{Deserialize, Serialize};

/// Coupling powers for each grid.
#[derive(Deserialize, Serialize)]
pub struct Order {
    /// Exponent of the strong coupling.
    pub alphas: u32,
    /// Exponent of the electromagnetic coupling.
    pub alpha: u32,
    /// Exponent of the logarithm of the scale factor of the renomalization scale.
    pub logxir: u32,
    /// Exponent of the logarithm of the scale factor of the factorization scale.
    pub logxif: u32,
}

/// Trait each subgrid must implement.
#[typetag::serde(tag = "type")]
pub trait Subgrid {
    /// Fills the subgrid with `weight` for the parton momentum fractions `x1` and `x2`, and the
    /// scale `q2`.
    fn fill(&mut self, x1: f64, x2: f64, q2: f64, weight: f64);

    /// Scale the subgrid by `factor`.
    fn scale(&mut self, factor: f64);
}

/// Main data structure of `PineAPPL`. This structure contains a `Subgrid` for each `LumiEntry`,
/// bin, and coupling order it was created with.
#[derive(Deserialize, Serialize)]
pub struct Grid {
    // TODO: this should probably be rewritten using something like ndarray
    subgrids: Vec<Vec<Vec<Box<dyn Subgrid>>>>,
    lumi: Lumi,
    bin_limits: BinLimits,
    orders: Vec<Order>,
}

impl Grid {
    /// Constructor.
    #[must_use]
    pub fn new(lumi: Lumi, orders: Vec<Order>, bin_limits: Vec<f64>) -> Self {
        assert!(!bin_limits.is_empty());

        Self {
            // usually we would use vec!, but `Subgrid` does not implement `Clone` (and it probably
            // cannot) so we can't use it in this instance
            subgrids: (0..orders.len())
                .map(|_| {
                    (0..bin_limits.len() - 1)
                        .map(|_| {
                            (0..lumi.len())
                                .map(|_| Box::new(NtupleSubgrid::default()) as Box<dyn Subgrid>)
                                .collect::<Vec<_>>()
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>(),
            orders,
            lumi,
            bin_limits: BinLimits::new(bin_limits),
        }
    }

    /// Fills the grid with events for the parton momentum fractions `x1` and `x2`, the scale `q2`,
    /// and the observable `obs`. The events are stored in `weights` and must be
    /// ordered as the corresponding luminosity function was created. The weights are filled into
    /// the subgrid with the given index.
    pub fn fill(&mut self, x1: f64, x2: f64, q2: f64, obs: f64, weights: &[f64], subgrid: usize) {
        if let Some(bin) = self.bin_limits.index(obs) {
            for (lumi, weight) in weights.iter().enumerate() {
                self.subgrids[subgrid][bin][lumi].fill(x1, x2, q2, *weight);
            }
        }
    }

    /// Returns the luminosity function.
    #[must_use]
    pub const fn lumi(&self) -> &Lumi {
        &self.lumi
    }

    /// Scale all subgrids by `factor`.
    pub fn scale(&mut self, factor: f64) {
        for i in &mut self.subgrids {
            for j in i {
                for k in j {
                    k.scale(factor);
                }
            }
        }
    }

    /// Returns the subgrid parameters.
    #[must_use]
    pub fn orders(&self) -> &[Order] {
        &self.orders
    }
}
