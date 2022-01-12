//! Provides an implementation of the `Grid` trait with n-tuples.

use super::grid::Ntuple;
use super::subgrid::{Mu2, Stats, Subgrid, SubgridEnum, SubgridIter};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;

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
    fn convolute(
        &self,
        _: &[f64],
        _: &[f64],
        _: &[Mu2],
        _: &mut dyn FnMut(usize, usize, usize) -> f64,
    ) -> f64 {
        todo!();
    }

    fn fill(&mut self, ntuple: &Ntuple<f64>) {
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
            todo!();
        }
    }

    fn scale(&mut self, factor: f64) {
        self.ntuples.iter_mut().for_each(|t| t.weight *= factor);
    }

    fn symmetrize(&mut self) {}

    fn clone_empty(&self) -> SubgridEnum {
        Self::new().into()
    }

    fn iter(&self) -> SubgridIter {
        unimplemented!();
    }

    fn stats(&self) -> Stats {
        unreachable!();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic]
    fn convolute() {
        let _ = NtupleSubgridV1::new().convolute(&[], &[], &[], &mut |_, _, _| 0.0);
    }

    #[test]
    #[should_panic]
    fn iter() {
        let _ = NtupleSubgridV1::new().iter();
    }

    #[test]
    #[should_panic]
    fn stats() {
        let subgrid = NtupleSubgridV1::new();
        subgrid.stats();
    }

    #[test]
    fn test() {
        let mut subgrid1: SubgridEnum = NtupleSubgridV1::new().into();

        assert!(subgrid1.is_empty());

        subgrid1.fill(&Ntuple {
            x1: 0.0,
            x2: 0.0,
            q2: 0.0,
            weight: 0.0,
        });

        assert!(!subgrid1.is_empty());

        assert_eq!(subgrid1.mu2_grid().as_ref(), []);
        assert_eq!(subgrid1.x1_grid().as_ref(), []);
        assert_eq!(subgrid1.x2_grid().as_ref(), []);

        subgrid1.symmetrize();
        subgrid1.scale(2.0);

        let mut subgrid2: SubgridEnum = subgrid1.clone_empty().into();

        subgrid2.fill(&Ntuple {
            x1: 0.0,
            x2: 0.0,
            q2: 0.0,
            weight: 0.0,
        });

        subgrid2.merge(&mut subgrid1, false);
    }
}
