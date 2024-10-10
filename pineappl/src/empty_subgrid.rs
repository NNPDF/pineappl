//! TODO

use super::interpolation::Interp;
use super::subgrid::{Mu2, NodeValues, Stats, Subgrid, SubgridEnum, SubgridIndexedIter};
use serde::{Deserialize, Serialize};
use std::iter;

/// A subgrid type that is always empty.
#[derive(Clone, Default, Deserialize, Serialize)]
pub struct EmptySubgridV1;

impl Subgrid for EmptySubgridV1 {
    fn fill(&mut self, _: &[Interp], _: &[f64], _: f64) {
        panic!("EmptySubgridV1 doesn't support the fill operation");
    }

    fn node_values(&self) -> Vec<NodeValues> {
        Vec::new()
    }

    fn is_empty(&self) -> bool {
        true
    }

    fn merge(&mut self, subgrid: &SubgridEnum, _: Option<(usize, usize)>) {
        assert!(
            subgrid.is_empty(),
            "EmptySubgridV1 doesn't support the merge operation for non-empty subgrids"
        );
    }

    fn scale(&mut self, _: f64) {}

    fn symmetrize(&mut self, _: usize, _: usize) {}

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
    use crate::v0;

    #[test]
    fn create_empty() {
        let mut subgrid = EmptySubgridV1;
        assert!(subgrid.is_empty());
        subgrid.merge(&mut EmptySubgridV1.into(), None);
        subgrid.scale(2.0);
        subgrid.symmetrize(1, 2);
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
        subgrid.fill(&v0::default_interps(2), &[0.0; 3], 0.0);
    }

    #[test]
    fn node_values() {
        assert!(EmptySubgridV1.node_values().is_empty());
    }
}
