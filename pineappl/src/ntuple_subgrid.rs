//! Provides an implementation of the `Grid` trait with n-tuples.

use super::grid::{Ntuple, Subgrid};
use serde::{Deserialize, Serialize};
use std::any::Any;

/// Structure holding a grid with an n-tuple as the storage method for weights.
#[derive(Default, Deserialize, Serialize)]
pub struct NtupleSubgridV1 {
    ntuples: Vec<Ntuple<f64>>,
}

impl NtupleSubgridV1 {
    /// Constructor.
    #[must_use]
    pub fn new() -> Self {
        Self { ntuples: vec![] }
    }
}

impl Subgrid for NtupleSubgridV1 {
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }

    fn convolute(&self, lumi: &dyn Fn(f64, f64, f64) -> f64) -> f64 {
        let mut result = 0.0;

        for ntuple in &self.ntuples {
            result += lumi(ntuple.x1, ntuple.x2, ntuple.q2) * ntuple.weight;
        }

        result
    }

    fn fill(&mut self, ntuple: &Ntuple<f64>) {
        self.ntuples.push(ntuple.clone());
    }

    fn is_empty(&self) -> bool {
        self.ntuples.is_empty()
    }

    fn merge(&mut self, other: &mut dyn Subgrid) {
        if let Some(other_grid) = other.as_any_mut().downcast_mut::<Self>() {
            self.ntuples.append(&mut other_grid.ntuples);
        } else {
            todo!();
        }
    }

    fn scale(&mut self, factor: f64) {
        self.ntuples.iter_mut().for_each(|t| t.weight *= factor);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test() {
        let mut subgrid = NtupleSubgridV1::new();
        assert!(subgrid.is_empty());

        subgrid.fill(&Ntuple {
            x1: 0.5,
            x2: 0.5,
            q2: 10.0,
            weight: 1.0,
        });
        assert!(!subgrid.is_empty());

        subgrid.fill(&Ntuple {
            x1: 0.25,
            x2: 0.75,
            q2: 100.0,
            weight: 3.0,
        });
        assert_eq!(subgrid.convolute(&|x1, x2, q2| x1 * x2 * q2), 2.5 + 56.25);

        let mut other_subgrid = NtupleSubgridV1::new();

        other_subgrid.fill(&Ntuple {
            x1: 0.25,
            x2: 0.5,
            q2: 20.0,
            weight: 2.0,
        });
        assert_eq!(other_subgrid.convolute(&|x1, x2, q2| x1 * x2 * q2), 5.0);

        subgrid.merge(&mut other_subgrid);
        assert_eq!(
            subgrid.convolute(&|x1, x2, q2| x1 * x2 * q2),
            2.5 + 56.25 + 5.0
        );

        subgrid.scale(0.5);
        assert_eq!(
            subgrid.convolute(&|x1, x2, q2| x1 * x2 * q2),
            1.25 + 28.125 + 2.5
        );
    }
}
