//! Provides an implementation of the `Grid` trait with n-tuples.

use super::grid::{Subgrid, SubgridEntry};
use serde::{Deserialize, Serialize};

/// Structure holding a grid with an n-tuple as the storage method for weights.
#[derive(Default, Deserialize, Serialize)]
pub struct NtupleSubgrid {
    entries: Vec<SubgridEntry<f64>>,
}

#[typetag::serde]
impl Subgrid for NtupleSubgrid {
    fn fill(&mut self, entry: SubgridEntry<f64>) {
        self.entries.push(entry);
    }

    fn scale(&mut self, factor: f64) {
        self.entries.iter_mut().for_each(|t| *t *= factor);
    }
}
