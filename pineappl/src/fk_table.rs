//! Provides the [`FkTable`] type.

use super::grid::{Grid, Order};
use super::subgrid::Subgrid;
use ndarray::Array4;
use std::convert::TryFrom;
use thiserror::Error;

/// Structure implementing FK tables. These are special [`Grid`]s, for which the following
/// additional guarantees are given:
///
/// - all subgrids of the grid evaluate the PDFs at a single scale `q2`. See [`FkTable::q2`].
/// - all subgrids, for both hadronic initial states (if both initial states are hadronic), share
///   the same `x` grid. See [`FkTable::x_grid`].
/// - the luminosity function is *simple*, meaning that every entry consists of a single pair of
///   partons with trivial factor `1.0`, and all tuples are distinct from each other. See
///   [`Grid::lumi`].
/// - the FK table's grid contains only a single [`Order`], whose exponents are all zero.
#[repr(transparent)]
pub struct FkTable {
    grid: Grid,
}

/// The error type returned when a conversion of a [`Grid`] to an [`FkTable`] fails.
#[derive(Debug, Error)]
pub enum TryFromGridError {
    /// Error if the grid contains multiple scales instead of a single one.
    #[error("multiple scales detected")]
    MultipleScales,
    /// Error if the x grids are not the same across all subgrids and hadronic initial states.
    #[error("different x grids detected")]
    NonUniqueGrids,
    /// Error if the luminosity is not simple.
    #[error("complicated luminosity function detected")]
    InvalidLumi,
    /// Error if the order of the grid was not a single one with all zeros in the exponents.
    #[error("multiple orders detected")]
    NonTrivialOrder,
    /// Error if the certain metadata is missing.
    #[error("metadata is missing: expected key `{0}` to have a value")]
    MetadataMissing(String),
}

impl FkTable {
    /// Returns the [`Grid`] object for this `FkTable`.
    pub fn grid(&self) -> &Grid {
        &self.grid
    }

    /// Returns the FK table represented as a four-dimensional array indexed by `bin`, `lumi`,
    /// `x1` and `x2`, in this order.
    pub fn table(&self) -> Array4<f64> {
        let mut subgrid = Array4::zeros((
            self.bins(),
            self.grid.lumi().len(),
            self.grid.subgrid(0, 0, 0).x1_grid().len(),
            self.grid.subgrid(0, 0, 0).x2_grid().len(),
        ));

        for bin in 0..self.bins() {
            for lumi in 0..self.grid.lumi().len() {
                for ((_, ix1, ix2), value) in self.grid().subgrid(0, bin, lumi).iter() {
                    subgrid[[bin, lumi, ix1, ix2]] = *value;
                }
            }
        }

        subgrid
    }

    /// Returns the number of bins for this `FkTable`.
    pub fn bins(&self) -> usize {
        self.grid.bin_info().bins()
    }

    /// Returns the (simplified) luminosity function for this `FkTable`. All factors are `1.0`.
    pub fn lumi(&self) -> Vec<(i32, i32)> {
        self.grid
            .lumi()
            .iter()
            .map(|entry| (entry.entry()[0].0, entry.entry()[0].1))
            .collect()
    }

    ///// Returns the single `q2` scale of this `FkTable`.
    //pub fn q2(&self) -> f64 {
    //    self.q2
    //}

    ///// Returns the x grid that all subgrids for all hadronic initial states share.
    //pub fn x_grid(&self) -> &[f64] {
    //    &self.x_grid
    //}
}

impl TryFrom<Grid> for FkTable {
    type Error = TryFromGridError;

    fn try_from(grid: Grid) -> Result<Self, Self::Error> {
        let mut q2 = -1.0;
        let mut x_grid = Vec::new();

        if grid.orders()
            != [Order {
                alphas: 0,
                alpha: 0,
                logxir: 0,
                logxif: 0,
            }]
        {
            return Err(TryFromGridError::NonTrivialOrder);
        }

        for order in 0..grid.orders().len() {
            for bin in 0..grid.bin_info().bins() {
                for lumi in 0..grid.lumi().len() {
                    let subgrid = grid.subgrid(order, bin, lumi);

                    if subgrid.is_empty() {
                        continue;
                    }

                    let q2_grid = subgrid.q2_grid();
                    let x1_grid = subgrid.x1_grid();
                    let x2_grid = subgrid.x2_grid();

                    if q2_grid.len() > 1 {
                        return Err(TryFromGridError::MultipleScales);
                    }

                    if q2 < 0.0 {
                        q2 = q2_grid[0];

                        if x1_grid.len() == 1 {
                            x_grid = x2_grid.into_owned();
                        } else {
                            x_grid = x1_grid.into_owned();
                        }
                    } else {
                        if q2 != q2_grid[0] {
                            return Err(TryFromGridError::MultipleScales);
                        }

                        if x1_grid.len() != 1 && x1_grid != x_grid {
                            return Err(TryFromGridError::NonUniqueGrids);
                        }

                        if x2_grid.len() != 1 && x2_grid != x_grid {
                            return Err(TryFromGridError::NonUniqueGrids);
                        }
                    }
                }
            }
        }

        for lumi in grid.lumi() {
            let entry = lumi.entry();

            if entry.len() != 1 || entry[0].2 != 1.0 {
                return Err(TryFromGridError::InvalidLumi);
            }
        }

        if (1..grid.lumi().len()).any(|i| grid.lumi()[i..].contains(&grid.lumi()[i - 1])) {
            return Err(TryFromGridError::InvalidLumi);
        }

        if let Some(key_values) = grid.key_values() {
            let keys = vec![
                "initial_state_1".to_string(),
                "initial_state_2".to_string(),
                "lumi_id_types".to_string(),
            ];

            for key in keys.into_iter() {
                if !key_values.contains_key(&key) {
                    return Err(TryFromGridError::MetadataMissing(key));
                }
            }
        } else {
            return Err(TryFromGridError::MetadataMissing(
                "initial_states_1".to_string(),
            ));
        }

        Ok(FkTable {
            grid, /*q2, x_grid*/
        })
    }
}
