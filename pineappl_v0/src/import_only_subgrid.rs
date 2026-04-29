//! TODO

use super::sparse_array3::SparseArray3;
use super::subgrid::{Mu2, Subgrid, SubgridIndexedIter};
use serde::Deserialize;
use std::borrow::Cow;

/// TODO
#[derive(Deserialize)]
pub struct ImportOnlySubgridV1 {
    array: SparseArray3<f64>,
    q2_grid: Vec<f64>,
    x1_grid: Vec<f64>,
    x2_grid: Vec<f64>,
}

impl Subgrid for ImportOnlySubgridV1 {
    fn mu2_grid(&self) -> Cow<'_, [Mu2]> {
        self.q2_grid
            .iter()
            .copied()
            .map(|q2| Mu2 { ren: q2, fac: q2 })
            .collect()
    }

    fn x1_grid(&self) -> Cow<'_, [f64]> {
        Cow::Borrowed(&self.x1_grid)
    }

    fn x2_grid(&self) -> Cow<'_, [f64]> {
        Cow::Borrowed(&self.x2_grid)
    }

    fn is_empty(&self) -> bool {
        self.array.is_empty()
    }

    fn indexed_iter(&self) -> SubgridIndexedIter<'_> {
        Box::new(self.array.indexed_iter())
    }
}

/// TODO
#[derive(Deserialize)]
pub struct ImportOnlySubgridV2 {
    array: SparseArray3<f64>,
    mu2_grid: Vec<Mu2>,
    x1_grid: Vec<f64>,
    x2_grid: Vec<f64>,
}

impl Subgrid for ImportOnlySubgridV2 {
    fn mu2_grid(&self) -> Cow<'_, [Mu2]> {
        Cow::Borrowed(&self.mu2_grid)
    }

    fn x1_grid(&self) -> Cow<'_, [f64]> {
        Cow::Borrowed(&self.x1_grid)
    }

    fn x2_grid(&self) -> Cow<'_, [f64]> {
        Cow::Borrowed(&self.x2_grid)
    }

    fn is_empty(&self) -> bool {
        self.array.is_empty()
    }

    fn indexed_iter(&self) -> SubgridIndexedIter<'_> {
        Box::new(self.array.indexed_iter())
    }
}
