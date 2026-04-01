//! Provides an implementation of the `Grid` trait with n-tuples.

use super::grid::Ntuple;
use super::subgrid::{Mu2, Stats, Subgrid, SubgridEnum, SubgridIndexedIter};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::mem;

/// Structure holding a grid with an n-tuple as the storage method for weights.
#[derive(Clone, Default, Deserialize, Serialize)]
pub struct NtupleSubgridV1 {
    ntuples: Vec<Ntuple<f64>>,
}

impl NtupleSubgridV1 {
    /// Constructor.
    #[must_use]
    pub const fn new() -> Self {
        Self { ntuples: vec![] }
    }
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

    fn merge(&mut self, other: &mut SubgridEnum, transpose: bool) {
        assert!(!transpose);

        if let SubgridEnum::NtupleSubgridV1(other_grid) = other {
            self.ntuples.append(&mut other_grid.ntuples);
        } else {
            panic!("NtupleSubgridV1 doesn't support the merge operation with subgrid types other than itself");
        }
    }

    fn scale(&mut self, factor: f64) {
        self.ntuples.iter_mut().for_each(|t| t.weight *= factor);
    }

    fn symmetrize(&mut self) {}

    fn clone_empty(&self) -> SubgridEnum {
        Self::new().into()
    }

    fn indexed_iter(&self) -> SubgridIndexedIter<'_> {
        panic!("NtupleSubgridV1 doesn't support the indexed_iter operation");
    }

    fn stats(&self) -> Stats {
        Stats {
            total: self.ntuples.len(),
            allocated: self.ntuples.len(),
            zeros: 0,
            overhead: 0,
            bytes_per_value: mem::size_of::<Ntuple<f64>>(),
        }
    }

    fn static_scale(&self) -> Option<Mu2> {
        todo!()
    }
}
