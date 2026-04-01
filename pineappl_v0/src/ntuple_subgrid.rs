//! Provides an implementation of the `Grid` trait with n-tuples.

use super::grid::Ntuple;
use super::subgrid::{Mu2, Subgrid, SubgridIndexedIter};
use serde::Deserialize;
use std::borrow::Cow;

/// Structure holding a grid with an n-tuple as the storage method for weights.
#[derive(Deserialize)]
pub struct NtupleSubgridV1 {
    ntuples: Vec<Ntuple<f64>>,
}

impl Subgrid for NtupleSubgridV1 {
    fn mu2_grid(&self) -> Cow<'_, [Mu2]> {
        Cow::Borrowed(&[])
    }

    fn x1_grid(&self) -> Cow<'_, [f64]> {
        Cow::Borrowed(&[])
    }

    fn x2_grid(&self) -> Cow<'_, [f64]> {
        Cow::Borrowed(&[])
    }

    fn is_empty(&self) -> bool {
        self.ntuples.is_empty()
    }

    fn indexed_iter(&self) -> SubgridIndexedIter<'_> {
        panic!("NtupleSubgridV1 doesn't support the indexed_iter operation");
    }
}
