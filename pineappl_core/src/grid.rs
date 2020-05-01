//! Module containing all traits and supporting structures for grids.

use super::bin::BinLimits;
use super::lumi::Lumi;
use super::ntuple_grid::NtupleSubgrid;

/// Structure for the metadata of a subgrid.
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
pub trait Subgrid {
    /// Fills the grid with events for the parton momentum fractions `x1` and `x2`, the scale `q2`,
    /// and the observable at the `obs_index`. The events are stored in `weights` and must be
    /// ordered as the corresponding luminosity function was created.
    fn fill(&mut self, x1: f64, x2: f64, q2: f64, obs_index: usize, weights: &[f64]);
}

/// A collection of subgrids.
pub struct Grid {
    subgrids: Vec<Box<dyn Subgrid>>,
    lumi: Lumi,
    bin_limits: BinLimits,
}

impl Grid {
    /// Constructor.
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
}
