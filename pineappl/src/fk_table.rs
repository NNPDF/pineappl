//! Provides the [`FkTable`] type.

use super::grid::{Grid, GridError, Order};
use super::lumi::LumiCache;
use super::subgrid::Subgrid;
use ndarray::Array4;
use std::collections::HashMap;
use std::convert::TryFrom;
use std::fmt::{self, Display, Formatter};
use std::io::Write;
use thiserror::Error;

/// Structure implementing FK tables. These are special [`Grid`]s, for which the following
/// additional guarantees are given:
///
/// - all subgrids of the grid evaluate the PDFs at a single scale `muf2`.
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

/// The optimization assumptions for an [`FkTable`], needed for [`FkTable::optimize`]. Since FK
/// tables are typically stored at very small `Q2 = Q0`, the PDFs `f(x,Q0)` of heavy quarks are
/// typically set to zero at this scale or set to the same value as their anti-quark PDF. This is
/// used to optimize the size of FK tables.
#[derive(Clone, Copy)]
pub enum FkAssumptions {
    /// All quark PDFs are non-zero at the FK table scale and completely independent.
    Nf6Ind,
    /// Like [`Nf6Ind`], but the PDFs of top and anti-top quarks are the same at FK table scale.
    Nf6Sym,
    /// Like [`Nf6Ind`], but the PDFs of top and anti-top quarks are zero at FK table scale.
    Nf5Ind,
    /// Like [`Nf5Ind`], but the PDFs of bottom and anti-bottom quarks are the same at FK table
    /// scale.
    Nf5Sym,
    /// Like [`Nf5Ind`], but the PDFs of bottom and anti-bottom quarks are zero at FK table scale.
    Nf4Ind,
    /// Like [`Nf4Ind`], but the PDFs of charm and anti-charm quarks are the same at FK table
    /// scale. PDF sets that make this assumption are NNPDF4.0 and NNPDF3.1 at fitting scale.
    Nf4Sym,
    /// Like [`Nf4Ind`], but the PDFs of charm and anti-charm quarks are zero at FK table scale.
    /// PDF sets that make this assumption are MSHT20 and NNPDF3.0 at fitting scale.
    Nf3Ind,
    /// Like [`Nf3Ind`], but the PDFs of strange and anti-strange are the same at FK table scale.
    /// A PDF set that makes this assumption is CT18 at fitting scale.
    Nf3Sym,
}

/// Error type when trying to construct [`FkAssumptions`] with a string.
#[derive(Debug, Error)]
#[error("unknown variant for FkAssumptions: {variant}")]
pub struct UnknownFkAssumption {
    variant: String,
}

impl Display for FkAssumptions {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                FkAssumptions::Nf6Ind => "Nf6Ind",
                FkAssumptions::Nf6Sym => "Nf6Sym",
                FkAssumptions::Nf5Ind => "Nf5Ind",
                FkAssumptions::Nf5Sym => "Nf5Sym",
                FkAssumptions::Nf4Ind => "Nf4Ind",
                FkAssumptions::Nf4Sym => "Nf4Sym",
                FkAssumptions::Nf3Ind => "Nf3Ind",
                FkAssumptions::Nf3Sym => "Nf3Sym",
            }
        )
    }
}

impl TryFrom<&str> for FkAssumptions {
    type Error = UnknownFkAssumption;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        Ok(match value {
            "Nf6Ind" => Self::Nf6Ind,
            "Nf6Sym" => Self::Nf6Sym,
            "Nf5Ind" => Self::Nf5Ind,
            "Nf5Sym" => Self::Nf5Sym,
            "Nf4Ind" => Self::Nf4Ind,
            "Nf4Sym" => Self::Nf4Sym,
            "Nf3Ind" => Self::Nf3Ind,
            "Nf3Sym" => Self::Nf3Sym,
            _ => {
                return Err(UnknownFkAssumption {
                    variant: value.to_string(),
                });
            }
        })
    }
}

impl FkTable {
    /// Returns the [`Grid`] object for this `FkTable`.
    #[must_use]
    pub const fn grid(&self) -> &Grid {
        &self.grid
    }

    /// Returns the FK table represented as a four-dimensional array indexed by `bin`, `lumi`,
    /// `x1` and `x2`, in this order.
    #[must_use]
    pub fn table(&self) -> Array4<f64> {
        let first_non_empty_subgrid = self
            .grid
            .subgrids()
            .iter()
            .find(|s| !s.is_empty())
            .unwrap_or_else(|| unreachable!()); // FK tables should never be empty

        let mut subgrid = Array4::zeros((
            self.bins(),
            self.grid.lumi().len(),
            first_non_empty_subgrid.x1_grid().len(),
            first_non_empty_subgrid.x2_grid().len(),
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
    #[must_use]
    pub fn bins(&self) -> usize {
        self.grid.bin_info().bins()
    }

    /// Extract the normalizations for each bin.
    #[must_use]
    pub fn bin_normalizations(&self) -> Vec<f64> {
        self.grid.bin_info().normalizations()
    }

    /// Extract the number of dimensions for bins.
    #[must_use]
    pub fn bin_dimensions(&self) -> usize {
        self.grid.bin_info().dimensions()
    }

    /// Extract the left edges of a specific bin dimension.
    #[must_use]
    pub fn bin_left(&self, dimension: usize) -> Vec<f64> {
        self.grid.bin_info().left(dimension)
    }

    /// Extract the right edges of a specific bin dimension.
    #[must_use]
    pub fn bin_right(&self, dimension: usize) -> Vec<f64> {
        self.grid.bin_info().right(dimension)
    }

    /// Access meta data
    #[must_use]
    pub const fn key_values(&self) -> Option<&HashMap<String, String>> {
        self.grid.key_values()
    }

    /// Returns the (simplified) luminosity function for this `FkTable`. All factors are `1.0`.
    #[must_use]
    pub fn lumi(&self) -> Vec<(i32, i32)> {
        self.grid
            .lumi()
            .iter()
            .map(|entry| (entry.entry()[0].0, entry.entry()[0].1))
            .collect()
    }

    /// Returns the single `muf2` scale of this `FkTable`.
    #[must_use]
    pub fn muf2(&self) -> f64 {
        for bin in 0..self.grid.bin_info().bins() {
            for lumi in 0..self.grid.lumi().len() {
                let subgrid = self.grid.subgrid(0, bin, lumi);

                if subgrid.is_empty() {
                    continue;
                }

                return subgrid.mu2_grid()[0].fac;
            }
        }

        unreachable!();
    }

    /// Returns the x grid that all subgrids for all hadronic initial states share.
    #[must_use]
    pub fn x_grid(&self) -> Vec<f64> {
        for bin in 0..self.grid.bin_info().bins() {
            for lumi in 0..self.grid.lumi().len() {
                let subgrid = self.grid.subgrid(0, bin, lumi);

                if subgrid.is_empty() {
                    continue;
                }

                if subgrid.x1_grid().is_empty() {
                    return subgrid.x2_grid().into_owned();
                }

                return subgrid.x1_grid().into_owned();
            }
        }

        unreachable!();
    }

    /// Propagate write to grid
    ///
    /// # Errors
    ///
    /// TODO
    pub fn write(&self, writer: impl Write) -> Result<(), GridError> {
        self.grid.write(writer)
    }

    /// Propagate `write_lz4` to `Grid`.
    ///
    /// # Errors
    ///
    /// See [`Grid::write_lz4`].
    pub fn write_lz4(&self, writer: impl Write) -> Result<(), GridError> {
        self.grid.write_lz4(writer)
    }

    /// Propagate convolute to grid
    pub fn convolute(
        &self,
        lumi_cache: &mut LumiCache,
        bin_indices: &[usize],
        lumi_mask: &[bool],
    ) -> Vec<f64> {
        self.grid
            .convolute(lumi_cache, &[], bin_indices, lumi_mask, &[(1.0, 1.0)])
    }

    /// Set a metadata key-value pair
    pub fn set_key_value(&mut self, key: &str, value: &str) {
        self.grid.set_key_value(key, value);
    }

    /// Optimizes the storage of FK tables based of assumptions of the PDFs at the FK table's
    /// scale.
    pub fn optimize(&mut self, assumptions: FkAssumptions) {
        let flavor_basis = match self.grid.key_values().unwrap()["lumi_id_types"].as_str() {
            "pdg_mc_ids" => true,
            "evol" => false,
            _ => unimplemented!(),
        };

        let mut add = Vec::new();
        let mut del = Vec::new();

        if flavor_basis {
            match assumptions {
                FkAssumptions::Nf6Ind => {
                    // nothing to do here
                }
                FkAssumptions::Nf6Sym => {
                    add.push((-6, 6));
                }
                FkAssumptions::Nf5Ind => {
                    del.extend_from_slice(&[-6, 6]);
                }
                FkAssumptions::Nf5Sym => {
                    add.push((-5, 5));
                    del.extend_from_slice(&[-6, 6]);
                }
                FkAssumptions::Nf4Ind => {
                    del.extend_from_slice(&[-6, 6, -5, 5]);
                }
                FkAssumptions::Nf4Sym => {
                    add.push((-4, 4));
                    del.extend_from_slice(&[-6, 6, -5, 5]);
                }
                FkAssumptions::Nf3Ind => {
                    del.extend_from_slice(&[-6, 6, -5, 5, -4, 4]);
                }
                FkAssumptions::Nf3Sym => {
                    add.push((-3, 3));
                    del.extend_from_slice(&[-6, 6, -5, 5, -4, 4]);
                }
            }
        } else {
            match assumptions {
                FkAssumptions::Nf6Ind => {
                    // nothing to do here
                }
                FkAssumptions::Nf6Sym => {
                    add.push((235, 200));
                }
                FkAssumptions::Nf5Ind => {
                    add.extend_from_slice(&[(235, 200), (135, 100)]);
                }
                FkAssumptions::Nf5Sym => {
                    add.extend_from_slice(&[(235, 200), (135, 100), (224, 200)]);
                }
                FkAssumptions::Nf4Ind => {
                    add.extend_from_slice(&[(235, 200), (135, 100), (224, 200), (124, 100)]);
                }
                FkAssumptions::Nf4Sym => {
                    add.extend_from_slice(&[
                        (235, 200),
                        (135, 100),
                        (224, 200),
                        (124, 100),
                        (215, 200),
                    ]);
                }
                FkAssumptions::Nf3Ind => {
                    add.extend_from_slice(&[
                        (235, 200),
                        (135, 100),
                        (224, 200),
                        (124, 100),
                        (215, 200),
                        (115, 100),
                    ]);
                }
                FkAssumptions::Nf3Sym => {
                    add.extend_from_slice(&[
                        (235, 200),
                        (135, 100),
                        (224, 200),
                        (124, 100),
                        (215, 200),
                        (115, 100),
                        (208, 200),
                    ]);
                }
            }
        }

        self.grid.rewrite_lumi(&add, &del);

        // store the assumption so that we can check it later on
        self.grid
            .set_key_value("fk_assumptions", &assumptions.to_string());
        self.grid.optimize();
    }
}

impl TryFrom<Grid> for FkTable {
    type Error = TryFromGridError;

    fn try_from(grid: Grid) -> Result<Self, Self::Error> {
        let mut muf2 = -1.0;
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

        for bin in 0..grid.bin_info().bins() {
            for lumi in 0..grid.lumi().len() {
                let subgrid = grid.subgrid(0, bin, lumi);

                if subgrid.is_empty() {
                    continue;
                }

                let mu2_grid = subgrid.mu2_grid();
                let x1_grid = subgrid.x1_grid();
                let x2_grid = subgrid.x2_grid();

                if mu2_grid.len() > 1 {
                    return Err(TryFromGridError::MultipleScales);
                }

                if muf2 < 0.0 {
                    muf2 = mu2_grid[0].fac;

                    if x1_grid.len() == 1 {
                        x_grid = x2_grid.into_owned();
                    } else {
                        x_grid = x1_grid.into_owned();
                    }
                } else {
                    if muf2 != mu2_grid[0].fac {
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

            for key in keys {
                if !key_values.contains_key(&key) {
                    return Err(TryFromGridError::MetadataMissing(key));
                }
            }
        } else {
            return Err(TryFromGridError::MetadataMissing(
                "initial_states_1".to_string(),
            ));
        }

        Ok(Self { grid })
    }
}
