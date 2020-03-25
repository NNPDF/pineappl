//! Provides an implementation of the `Grid` trait with n-tuples.

use super::bin::BinLimits;
use super::grid::Grid;
use super::grid::Subgrid;
use super::lumi::Lumi;

#[derive(Clone, Debug, PartialEq)]
struct Ntuple {
    x1: f64,
    x2: f64,
    q2: f64,
    weight: f64,
}

impl Ntuple {
    fn new(x1: f64, x2: f64, q2: f64, weight: f64) -> Ntuple {
        Ntuple { x1, x2, q2, weight }
    }
}

/// Structure holding a grid with an n-tuple as the storage method for weights.
pub struct NtupleGrid {
    ntuples: Vec<Vec<Vec<Vec<Ntuple>>>>,
    lumi: Lumi,
    subgrids: Vec<Subgrid>,
    bin_limits: BinLimits,
}

impl NtupleGrid {
    /// Constructor.
    pub fn new(lumi: Lumi, subgrids: Vec<Subgrid>, bin_limits: Vec<f64>) -> NtupleGrid {
        NtupleGrid {
            ntuples: vec![vec![vec![vec![]; lumi.len()]; bin_limits.len() - 1]; subgrids.len()],
            lumi,
            subgrids,
            bin_limits: BinLimits::new(bin_limits),
        }
    }
}

impl Grid for NtupleGrid {
    fn fill(&mut self, x1: f64, x2: f64, q2: f64, obs: f64, weights: &[f64], subgrid: usize) {
        assert!(subgrid < self.subgrids.len());
        assert!(weights.len() == self.lumi.len());

        if let Some(obs_index) = self.bin_limits.index(obs) {
            for (lumi_index, &weight) in weights.iter().enumerate() {
                self.ntuples[subgrid][obs_index][lumi_index].push(Ntuple::new(x1, x2, q2, weight));
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::lumi_entry;

    #[test]
    fn test_fill() {
        // luminosity function with a single entry
        let lumi = Lumi::new(vec![lumi_entry![2, 2, 1.0]]);

        // a single subgrid
        let subgrids = vec![Subgrid {
            alphas: 0,
            alpha: 2,
            logxir: 0,
            logxif: 0,
        }];
        // four bins
        let bin_limits = vec![0.0, 0.25, 0.5, 0.75, 1.0];

        let mut grid = NtupleGrid::new(lumi, subgrids, bin_limits);

        // check first bin
        grid.fill(0.5, 0.75, 1.0, 0.125, &[2.0], 0);
        assert_eq!(grid.ntuples[0][0][0][0], Ntuple::new(0.5, 0.75, 1.0, 2.0));

        // check second bin
        grid.fill(0.75, 0.5, 2.0, 0.375, &[3.0], 0);
        assert_eq!(grid.ntuples[0][0][0][0], Ntuple::new(0.5, 0.75, 1.0, 2.0));
        assert_eq!(grid.ntuples[0][1][0][0], Ntuple::new(0.75, 0.5, 2.0, 3.0));

        // check third bin
        grid.fill(0.125, 0.25, 3.0, 0.625, &[4.0], 0);
        assert_eq!(grid.ntuples[0][0][0][0], Ntuple::new(0.5, 0.75, 1.0, 2.0));
        assert_eq!(grid.ntuples[0][1][0][0], Ntuple::new(0.75, 0.5, 2.0, 3.0));
        assert_eq!(grid.ntuples[0][2][0][0], Ntuple::new(0.125, 0.25, 3.0, 4.0));

        // check fourth bin
        grid.fill(0.5, 0.5, 4.0, 0.875, &[5.0], 0);
        assert_eq!(grid.ntuples[0][0][0][0], Ntuple::new(0.5, 0.75, 1.0, 2.0));
        assert_eq!(grid.ntuples[0][1][0][0], Ntuple::new(0.75, 0.5, 2.0, 3.0));
        assert_eq!(grid.ntuples[0][2][0][0], Ntuple::new(0.125, 0.25, 3.0, 4.0));
        assert_eq!(grid.ntuples[0][3][0][0], Ntuple::new(0.5, 0.5, 4.0, 5.0));

        // check overflow
        grid.fill(0.5, 0.5, 4.0, 2.0, &[99.0], 0);
        assert_eq!(grid.ntuples[0][0][0][0], Ntuple::new(0.5, 0.75, 1.0, 2.0));
        assert_eq!(grid.ntuples[0][1][0][0], Ntuple::new(0.75, 0.5, 2.0, 3.0));
        assert_eq!(grid.ntuples[0][2][0][0], Ntuple::new(0.125, 0.25, 3.0, 4.0));
        assert_eq!(grid.ntuples[0][3][0][0], Ntuple::new(0.5, 0.5, 4.0, 5.0));

        // check underflow
        grid.fill(0.5, 0.5, 4.0, -1.0, &[99.0], 0);
        assert_eq!(grid.ntuples[0][0][0][0], Ntuple::new(0.5, 0.75, 1.0, 2.0));
        assert_eq!(grid.ntuples[0][1][0][0], Ntuple::new(0.75, 0.5, 2.0, 3.0));
        assert_eq!(grid.ntuples[0][2][0][0], Ntuple::new(0.125, 0.25, 3.0, 4.0));
        assert_eq!(grid.ntuples[0][3][0][0], Ntuple::new(0.5, 0.5, 4.0, 5.0));
    }
}
