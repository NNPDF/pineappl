//! TODO

use super::grid::Ntuple;
use super::subgrid::{Subgrid, SubgridEnum};
use either::Either;
use serde::{Deserialize, Serialize};

/// A subgrid type that is always empty.
#[derive(Default, Deserialize, Serialize)]
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

    fn q2_grid(&self) -> Vec<f64> {
        unreachable!();
    }

    fn x1_grid(&self) -> Vec<f64> {
        unreachable!();
    }

    fn x2_grid(&self) -> Vec<f64> {
        unreachable!();
    }

    fn is_empty(&self) -> bool {
        true
    }

    fn merge(&mut self, subgrid: &mut SubgridEnum, _: bool) {
        assert!(subgrid.is_empty());
    }

    fn scale(&mut self, _: f64) {}

    fn q2_slice(&self) -> (usize, usize) {
        unreachable!();
    }

    fn fill_q2_slice(&self, _: usize, _: &mut [f64]) {
        unreachable!();
    }

    fn write_q2_slice(&mut self, _: usize, _: &[f64]) {
        unreachable!();
    }

    fn symmetrize(&mut self) {}

    fn clone_empty(&self) -> SubgridEnum {
        Self::default().into()
    }
}
