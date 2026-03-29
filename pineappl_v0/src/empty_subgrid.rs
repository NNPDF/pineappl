//! TODO

use super::grid::Ntuple;
use super::subgrid::{Mu2, Stats, Subgrid, SubgridEnum, SubgridIndexedIter};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::iter;

/// A subgrid type that is always empty.
#[derive(Clone, Default, Deserialize, Serialize)]
pub struct EmptySubgridV1;

impl Subgrid for EmptySubgridV1 {
    fn convolve(
        &self,
        _: &[f64],
        _: &[f64],
        _: &[Mu2],
        _: &mut dyn FnMut(usize, usize, usize) -> f64,
    ) -> f64 {
        0.0
    }

    fn fill(&mut self, _: &Ntuple<f64>) {
        panic!("EmptySubgridV1 doesn't support the fill operation");
    }

    fn mu2_grid(&self) -> Cow<[Mu2]> {
        Cow::Borrowed(&[])
    }

    fn x1_grid(&self) -> Cow<[f64]> {
        Cow::Borrowed(&[])
    }

    fn x2_grid(&self) -> Cow<[f64]> {
        Cow::Borrowed(&[])
    }

    fn is_empty(&self) -> bool {
        true
    }

    fn merge(&mut self, subgrid: &mut SubgridEnum, _: bool) {
        assert!(
            subgrid.is_empty(),
            "EmptySubgridV1 doesn't support the merge operation for non-empty subgrids"
        );
    }

    fn scale(&mut self, _: f64) {}

    fn symmetrize(&mut self) {}

    fn clone_empty(&self) -> SubgridEnum {
        Self.into()
    }

    fn indexed_iter(&self) -> SubgridIndexedIter {
        Box::new(iter::empty())
    }

    fn stats(&self) -> Stats {
        Stats {
            total: 0,
            allocated: 0,
            zeros: 0,
            overhead: 0,
            bytes_per_value: 0,
        }
    }

    fn static_scale(&self) -> Option<Mu2> {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_empty() {
        let mut subgrid = EmptySubgridV1;
        assert_eq!(subgrid.convolve(&[], &[], &[], &mut |_, _, _| 0.0), 0.0,);
        assert!(subgrid.is_empty());
        subgrid.merge(&mut EmptySubgridV1.into(), false);
        subgrid.scale(2.0);
        subgrid.symmetrize();
        assert!(subgrid.clone_empty().is_empty());
        assert_eq!(
            subgrid.stats(),
            Stats {
                total: 0,
                allocated: 0,
                zeros: 0,
                overhead: 0,
                bytes_per_value: 0,
            }
        );
        assert_eq!(subgrid.static_scale(), None);
    }

    #[test]
    #[should_panic(expected = "EmptySubgridV1 doesn't support the fill operation")]
    fn fill() {
        let mut subgrid = EmptySubgridV1;
        subgrid.fill(&Ntuple {
            x1: 0.0,
            x2: 0.0,
            q2: 0.0,
            weight: 0.0,
        });
    }

    #[test]
    fn q2_grid() {
        assert!(EmptySubgridV1.mu2_grid().is_empty());
    }

    #[test]
    fn x1_grid() {
        assert!(EmptySubgridV1.x1_grid().is_empty());
    }

    #[test]
    fn x2_grid() {
        assert!(EmptySubgridV1.x2_grid().is_empty());
    }
}
