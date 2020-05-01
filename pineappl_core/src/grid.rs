//! Module containing all traits and supporting structures for grids.

use super::bin::BinLimits;
use super::lumi::Lumi;
use super::ntuple_grid::NtupleSubgrid;
use serde::{Deserialize, Serialize};

/// Structure for the metadata of a subgrid.
#[derive(Deserialize, Serialize)]
pub struct SubgridData {
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
    /// Fills the grid with events for the parton momentum fractions `x1` and `x2`, the scale `q2`,
    /// and the observable at the `obs_index`. The events are stored in `weights` and must be
    /// ordered as the corresponding luminosity function was created.
    fn fill(&mut self, x1: f64, x2: f64, q2: f64, obs_index: usize, weights: &[f64]);

    /// Scale the subgrid by `factor`.
    fn scale(&mut self, factor: f64);

    /// Returns the subgrid data.
    fn data(&self) -> &SubgridData;
}

/// A collection of subgrids.
#[derive(Deserialize, Serialize)]
pub struct Grid {
    subgrids: Vec<Box<dyn Subgrid>>,
    lumi: Lumi,
    bin_limits: BinLimits,
}

impl Grid {
    /// Constructor.
    #[must_use]
    pub fn new(lumi: Lumi, subgrid_data: Vec<SubgridData>, bin_limits: Vec<f64>) -> Self {
        Self {
            subgrids: subgrid_data
                .into_iter()
                .map(|d| {
                    let result: Box<dyn Subgrid> =
                        Box::new(NtupleSubgrid::new(lumi.len(), bin_limits.len() - 1, d));
                    result
                })
                .collect::<Vec<Box<dyn Subgrid>>>(),
            lumi,
            bin_limits: BinLimits::new(bin_limits),
        }
    }

    /// Fills the grid with events for the parton momentum fractions `x1` and `x2`, the scale `q2`,
    /// and the observable `obs`. The events are stored in `weights` and must be
    /// ordered as the corresponding luminosity function was created. The weights are filled into
    /// the subgrid with the given index.
    pub fn fill(&mut self, x1: f64, x2: f64, q2: f64, obs: f64, weights: &[f64], subgrid: usize) {
        if let Some(obs_index) = self.bin_limits.index(obs) {
            self.subgrids[subgrid].fill(x1, x2, q2, obs_index, weights);
        }
    }

    /// Returns the luminosity function.
    #[must_use]
    pub const fn lumi(&self) -> &Lumi {
        &self.lumi
    }

    /// Scale all subgrids by `factor`.
    pub fn scale(&mut self, factor: f64) {
        self.subgrids.iter_mut().for_each(|g| g.scale(factor));
    }

    /// Returns the subgrid parameters.
    #[must_use]
    pub fn subgrids(&self) -> &[Box<dyn Subgrid>] {
        &self.subgrids
    }
}
