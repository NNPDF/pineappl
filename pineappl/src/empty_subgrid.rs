//! TODO

use super::grid::Ntuple;
use super::subgrid::{Subgrid, SubgridEnum};
use either::Either;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::iter;
use std::ops::Range;

/// A subgrid type that is always empty.
#[derive(Clone, Default, Deserialize, Serialize)]
pub struct EmptySubgridV1 {}

impl Subgrid for EmptySubgridV1 {
    fn convolute(
        &self,
        _: &[f64],
        _: &[f64],
        _: &[f64],
        _: Either<&dyn Fn(usize, usize, usize) -> f64, &dyn Fn(f64, f64, f64) -> f64>,
    ) -> f64 {
        0.0
    }

    fn fill(&mut self, _: &Ntuple<f64>) {
        unreachable!();
    }

    fn q2_grid(&self) -> Cow<[f64]> {
        unreachable!();
    }

    fn x1_grid(&self) -> Cow<[f64]> {
        unreachable!();
    }

    fn x2_grid(&self) -> Cow<[f64]> {
        unreachable!();
    }

    fn is_empty(&self) -> bool {
        true
    }

    fn merge(&mut self, subgrid: &mut SubgridEnum, _: bool) {
        assert!(subgrid.is_empty());
    }

    fn scale(&mut self, _: f64) {}

    fn q2_slice(&self) -> Range<usize> {
        unreachable!();
    }

    fn fill_q2_slice(&self, _: usize, _: &mut [f64]) {
        unreachable!();
    }

    fn symmetrize(&mut self) {}

    fn clone_empty(&self) -> SubgridEnum {
        Self::default().into()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = ((usize, usize, usize), &f64)> + '_> {
        Box::new(iter::empty())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_empty() {
        let mut subgrid = EmptySubgridV1::default();
        assert_eq!(
            subgrid.convolute(&[], &[], &[], Either::Left(&|_, _, _| 0.0)),
            0.0,
        );
        assert!(subgrid.is_empty());
        subgrid.merge(&mut EmptySubgridV1::default().into(), false);
        subgrid.scale(2.0);
        subgrid.symmetrize();
        assert!(subgrid.clone_empty().is_empty());
    }

    #[test]
    #[should_panic]
    fn fill() {
        let mut subgrid = EmptySubgridV1::default();
        subgrid.fill(&Ntuple {
            x1: 0.0,
            x2: 0.0,
            q2: 0.0,
            weight: 0.0,
        });
    }

    #[test]
    #[should_panic]
    fn q2_grid() {
        let subgrid = EmptySubgridV1::default();
        subgrid.q2_grid();
    }

    #[test]
    #[should_panic]
    fn x1_grid() {
        let subgrid = EmptySubgridV1::default();
        subgrid.x1_grid();
    }

    #[test]
    #[should_panic]
    fn x2_grid() {
        let subgrid = EmptySubgridV1::default();
        subgrid.x2_grid();
    }

    #[test]
    #[should_panic]
    fn q2_slice() {
        let subgrid = EmptySubgridV1::default();
        subgrid.q2_slice();
    }

    #[test]
    #[should_panic]
    fn fill_q2_slice() {
        let subgrid = EmptySubgridV1::default();
        subgrid.fill_q2_slice(0, &mut []);
    }
}
