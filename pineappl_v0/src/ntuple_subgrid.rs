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
    fn convolve(
        &self,
        _: &[f64],
        _: &[f64],
        _: &[Mu2],
        _: &mut dyn FnMut(usize, usize, usize) -> f64,
    ) -> f64 {
        panic!("NtupleSubgridV1 doesn't support the convolve operation");
    }

    fn fill(&mut self, ntuple: &Ntuple<f64>) {
        if ntuple.weight == 0.0 {
            return;
        }

        self.ntuples.push(ntuple.clone());
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

    fn indexed_iter(&self) -> SubgridIndexedIter {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lagrange_subgrid::LagrangeSubgridV2;
    use crate::subgrid::{ExtraSubgridParams, SubgridParams};

    #[test]
    #[should_panic(expected = "NtupleSubgridV1 doesn't support the convolve operation")]
    fn convolve() {
        NtupleSubgridV1::new().convolve(&[], &[], &[], &mut |_, _, _| 0.0);
    }

    #[test]
    fn fill_zero() {
        let mut subgrid = NtupleSubgridV1::new();

        subgrid.fill(&Ntuple {
            x1: 0.5,
            x2: 0.5,
            q2: 1000.0,
            weight: 0.0,
        });

        assert!(subgrid.is_empty());
    }

    #[test]
    #[should_panic(expected = "NtupleSubgridV1 doesn't support the indexed_iter operation")]
    fn indexed_iter() {
        // `next` isn't called because `indexed_iter` panics, but it suppresses a warning about an
        // unused result
        NtupleSubgridV1::new().indexed_iter().next();
    }

    #[test]
    fn stats() {
        let subgrid = NtupleSubgridV1::new();
        assert_eq!(
            subgrid.stats(),
            Stats {
                total: 0,
                allocated: 0,
                zeros: 0,
                overhead: 0,
                bytes_per_value: 32,
            }
        );
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn static_scale() {
        let subgrid = NtupleSubgridV1::new();
        subgrid.static_scale();
    }

    #[test]
    #[should_panic(
        expected = "NtupleSubgridV1 doesn't support the merge operation with subgrid types other than itself"
    )]
    fn merge_with_lagrange_subgrid() {
        let mut subgrid = NtupleSubgridV1::new();
        let mut other =
            LagrangeSubgridV2::new(&SubgridParams::default(), &ExtraSubgridParams::default())
                .into();
        subgrid.merge(&mut other, false);
    }

    #[test]
    fn test() {
        let mut subgrid1: SubgridEnum = NtupleSubgridV1::new().into();

        assert!(subgrid1.is_empty());

        subgrid1.fill(&Ntuple {
            x1: 0.0,
            x2: 0.0,
            q2: 0.0,
            weight: 1.0,
        });

        assert!(!subgrid1.is_empty());

        assert_eq!(subgrid1.mu2_grid().as_ref(), []);
        assert_eq!(subgrid1.x1_grid().as_ref(), []);
        assert_eq!(subgrid1.x2_grid().as_ref(), []);

        subgrid1.symmetrize();
        subgrid1.scale(2.0);

        let mut subgrid2: SubgridEnum = subgrid1.clone_empty();

        subgrid2.fill(&Ntuple {
            x1: 0.0,
            x2: 0.0,
            q2: 0.0,
            weight: 1.0,
        });

        subgrid2.merge(&mut subgrid1, false);
    }
}
