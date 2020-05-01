//! Provides an implementation of the `Grid` trait with n-tuples.

use super::grid::Subgrid;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
struct Ntuple {
    x1: f64,
    x2: f64,
    q2: f64,
    weight: f64,
}

impl Ntuple {
    const fn new(x1: f64, x2: f64, q2: f64, weight: f64) -> Self {
        Self { x1, x2, q2, weight }
    }
}

/// Structure holding a grid with an n-tuple as the storage method for weights.
#[derive(Default, Deserialize, Serialize)]
pub struct NtupleSubgrid {
    ntuples: Vec<Ntuple>,
}

#[typetag::serde]
impl Subgrid for NtupleSubgrid {
    fn fill(&mut self, x1: f64, x2: f64, q2: f64, weight: f64) {
        self.ntuples.push(Ntuple::new(x1, x2, q2, weight));
    }

    fn scale(&mut self, factor: f64) {
        self.ntuples.iter_mut().for_each(|t| t.weight *= factor);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ntuple_subgrid_fill() {
        let mut grid = NtupleSubgrid::default();

        grid.fill(0.5, 0.75, 1.0, 2.0);
        assert_eq!(grid.ntuples[0], Ntuple::new(0.5, 0.75, 1.0, 2.0));

        grid.fill(0.75, 0.5, 2.0, 3.0);
        assert_eq!(grid.ntuples[0], Ntuple::new(0.5, 0.75, 1.0, 2.0));
        assert_eq!(grid.ntuples[1], Ntuple::new(0.75, 0.5, 2.0, 3.0));

        grid.fill(0.125, 0.25, 3.0, 4.0);
        assert_eq!(grid.ntuples[0], Ntuple::new(0.5, 0.75, 1.0, 2.0));
        assert_eq!(grid.ntuples[1], Ntuple::new(0.75, 0.5, 2.0, 3.0));
        assert_eq!(grid.ntuples[2], Ntuple::new(0.125, 0.25, 3.0, 4.0));

        grid.fill(0.5, 0.5, 4.0, 5.0);
        assert_eq!(grid.ntuples[0], Ntuple::new(0.5, 0.75, 1.0, 2.0));
        assert_eq!(grid.ntuples[1], Ntuple::new(0.75, 0.5, 2.0, 3.0));
        assert_eq!(grid.ntuples[2], Ntuple::new(0.125, 0.25, 3.0, 4.0));
        assert_eq!(grid.ntuples[3], Ntuple::new(0.5, 0.5, 4.0, 5.0));
    }
}
