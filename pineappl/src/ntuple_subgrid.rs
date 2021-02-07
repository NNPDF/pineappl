//! Provides an implementation of the `Grid` trait with n-tuples.

use super::grid::Ntuple;
use super::subgrid::{Subgrid, SubgridEnum};
use either::Either;
use serde::{Deserialize, Serialize};

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
        lumi: Either<&dyn Fn(usize, usize, usize) -> f64, &dyn Fn(f64, f64, f64) -> f64>,
    ) -> f64 {
        let lumi = lumi.right().unwrap();
        let mut result = 0.0;

        for ntuple in &self.ntuples {
            result += lumi(ntuple.x1, ntuple.x2, ntuple.q2) * ntuple.weight;
        }

        result
    }

    fn fill(&mut self, ntuple: &Ntuple<f64>) {
        self.ntuples.push(ntuple.clone());
    }

    fn q2_grid(&self) -> Vec<f64> {
        vec![]
    }

    fn x1_grid(&self) -> Vec<f64> {
        vec![]
    }

    fn x2_grid(&self) -> Vec<f64> {
        vec![]
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

    fn q2_slice(&self) -> (usize, usize) {
        unimplemented!();
    }

    fn fill_q2_slice(&self, _: usize, _: &mut [f64]) {
        unimplemented!();
    }

    fn write_q2_slice(&mut self, _: usize, _: &[f64]) {
        unimplemented!();
    }

    fn symmetrize(&mut self) {}

    fn clone_empty(&self) -> SubgridEnum {
        Self::new().into()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use either::Either::Right;

    #[test]
    fn test() {
        let mut subgrid: SubgridEnum = NtupleSubgridV1::new().into();
        assert!(subgrid.is_empty());

        subgrid.fill(&Ntuple {
            x1: 0.5,
            x2: 0.5,
            q2: 10.0,
            weight: 1.0,
        });
        assert!(!subgrid.is_empty());

        subgrid.fill(&Ntuple {
            x1: 0.25,
            x2: 0.75,
            q2: 100.0,
            weight: 3.0,
        });
        assert_eq!(
            subgrid.convolute(&[], &[], &[], Right(&|x1, x2, q2| x1 * x2 * q2)),
            2.5 + 56.25
        );

        let mut other_subgrid: SubgridEnum = NtupleSubgridV1::new().into();

        other_subgrid.fill(&Ntuple {
            x1: 0.25,
            x2: 0.5,
            q2: 20.0,
            weight: 2.0,
        });
        assert_eq!(
            other_subgrid.convolute(&[], &[], &[], Right(&|x1, x2, q2| x1 * x2 * q2)),
            5.0
        );

        subgrid.merge(&mut other_subgrid, false);
        assert_eq!(
            subgrid.convolute(&[], &[], &[], Right(&|x1, x2, q2| x1 * x2 * q2)),
            2.5 + 56.25 + 5.0
        );

        subgrid.scale(0.5);
        assert_eq!(
            subgrid.convolute(&[], &[], &[], Right(&|x1, x2, q2| x1 * x2 * q2)),
            1.25 + 28.125 + 2.5
        );
    }

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

    #[test]
    #[should_panic]
    fn write_q2_slice() {
        let mut subgrid: SubgridEnum = NtupleSubgridV1::new().into();

        subgrid.write_q2_slice(0, &[]);
    }
}
