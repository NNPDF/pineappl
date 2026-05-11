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

#[cfg(test)]
mod tests {
    use super::*;

    // the following three methods aren't called if the subgrid is empty

    #[test]
    fn empty_subgrid_v1_x1_grid() {
        assert!(EmptySubgridV1.x1_grid().is_empty());
    }

    #[test]
    fn empty_subgrid_v1_x2_grid() {
        assert!(EmptySubgridV1.x2_grid().is_empty());
    }

    #[test]
    fn empty_subgrid_v1_indexed_iter() {
        assert!(EmptySubgridV1.indexed_iter().collect::<Vec<_>>().is_empty());
    }
}
