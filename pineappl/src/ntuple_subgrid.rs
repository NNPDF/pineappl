//! Provides an implementation of the `Grid` trait with n-tuples.

use super::grid::Ntuple;
use super::subgrid::{Subgrid, SubgridEnum};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::ops::Range;

/// Structure holding a grid with an n-tuple as the storage method for weights.
#[derive(Default, Deserialize, Serialize)]
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
        _: &[f64],
        _: &dyn Fn(usize, usize, usize) -> f64,
    ) -> f64 {
        todo!();
    }

    fn fill(&mut self, ntuple: &Ntuple<f64>) {
        self.ntuples.push(ntuple.clone());
    }

    fn q2_grid(&self) -> Cow<[f64]> {
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

    fn q2_slice(&self) -> Range<usize> {
        unimplemented!();
    }

    fn fill_q2_slice(&self, _: usize, _: &mut [f64]) {
        unimplemented!();
    }

    fn symmetrize(&mut self) {}

    fn clone_empty(&self) -> SubgridEnum {
        Self::new().into()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = ((usize, usize, usize), &f64)>> {
        unimplemented!();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic]
    fn q2_slice() {
        let subgrid: SubgridEnum = NtupleSubgridV1::new().into();

        subgrid.q2_slice();
    }

    #[test]
    #[should_panic]
    fn fill_q2_slice() {
        let subgrid: SubgridEnum = NtupleSubgridV1::new().into();

        subgrid.fill_q2_slice(0, &mut []);
    }
}
