//! TODO

use super::subgrid::{Mu2, Subgrid, SubgridIndexedIter};
use serde::Deserialize;
use std::borrow::Cow;
use std::iter;

/// A subgrid type that is always empty.
#[derive(Deserialize)]
pub struct EmptySubgridV1;

impl Subgrid for EmptySubgridV1 {
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
        true
    }

    fn indexed_iter(&self) -> SubgridIndexedIter<'_> {
        Box::new(iter::empty())
    }
}
