//! Provides an implementation of the `Grid` trait with n-tuples.

use super::grid::{Ntuple, Subgrid};
use serde::{Deserialize, Serialize};
use std::any::Any;

/// Structure holding a grid with an n-tuple as the storage method for weights.
#[derive(Default, Deserialize, Serialize)]
pub struct NtupleSubgrid {
    ntuples: Vec<Ntuple<f64>>,
}

#[typetag::serde]
impl Subgrid for NtupleSubgrid {
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
