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

#[cfg(test)]
mod tests {
    use super::*;

    // this subgrid type is now useless, but we have to keep it to keep backwards compatibility

    #[test]
    fn ntuple_subgrid_v1_mu2_grid() {
        assert!(
            NtupleSubgridV1 {
                ntuples: Vec::new()
            }
            .mu2_grid()
            .is_empty()
        );
    }

    #[test]
    fn ntuple_subgrid_v1_x1_grid() {
        assert!(
            NtupleSubgridV1 {
                ntuples: Vec::new()
            }
            .x1_grid()
            .is_empty()
        );
    }

    #[test]
    fn ntuple_subgrid_v1_x2_grid() {
        assert!(
            NtupleSubgridV1 {
                ntuples: Vec::new()
            }
            .x2_grid()
            .is_empty()
        );
    }

    #[test]
    fn ntuple_subgrid_v1_is_empty() {
        assert!(
            NtupleSubgridV1 {
                ntuples: Vec::new()
            }
            .is_empty()
        );
    }

    #[test]
    #[should_panic(expected = "NtupleSubgridV1 doesn't support the indexed_iter operation")]
    fn ntuple_subgrid_v1_indexed_iter() {
        let _ = NtupleSubgridV1 {
            ntuples: Vec::new(),
        }
        .indexed_iter();
    }
}
