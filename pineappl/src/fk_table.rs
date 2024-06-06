//! Provides the [`FkTable`] type.

use super::boc::Order;
use super::convolutions::{Convolution, LumiCache};
use super::grid::{Grid, GridError};
use super::subgrid::Subgrid;
use float_cmp::approx_eq;
use ndarray::Array4;
use std::collections::HashMap;
use std::fmt::{self, Display, Formatter};
use std::io::Write;
use std::str::FromStr;
use thiserror::Error;

/// Structure implementing FK tables. These are special [`Grid`]s, for which the following
/// additional guarantees are given:
///
/// - all subgrids of the grid evaluate the PDFs at a single factorization scale given by
///   [`FkTable::muf2`].
/// - all subgrids, for both hadronic initial states (if both initial states are hadronic), share
///   the same `x` grid. See [`FkTable::x_grid`].
/// - the channel definitions are *simple*, meaning that every entry consists of a single pair of
///   partons with trivial factor `1.0`, and all tuples are distinct from each other. See
///   [`Grid::channels`].
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
    /// Error if the channels are not simple.
    #[error("complicated channel function detected")]
    InvalidChannel,
    /// Error if the order of the grid was not a single one with all zeros in the exponents.
    #[error("multiple orders detected")]
    NonTrivialOrder,
}

/// The optimization assumptions for an [`FkTable`], needed for [`FkTable::optimize`]. Since FK
/// tables are typically stored at very small `Q2 = Q0`, the PDFs `f(x,Q0)` of heavy quarks are
/// typically set to zero at this scale or set to the same value as their anti-quark PDF. This is
/// used to optimize the size of FK tables.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum FkAssumptions {
    /// All quark PDFs are non-zero at the FK table scale and completely independent.
    Nf6Ind,
    /// Like [`Nf6Ind`](Self::Nf6Ind), but the PDFs of top and anti-top quarks are the same at FK
    /// table scale.
    Nf6Sym,
    /// Like [`Nf6Ind`](Self::Nf6Ind), but the PDFs of top and anti-top quarks are zero at FK table
    /// scale.
    Nf5Ind,
    /// Like [`Nf5Ind`](Self::Nf5Ind), but the PDFs of bottom and anti-bottom quarks are the same
    /// at FK table scale.
    Nf5Sym,
    /// Like [`Nf5Ind`](Self::Nf5Ind), but the PDFs of bottom and anti-bottom quarks are zero at FK
    /// table scale.
    Nf4Ind,
    /// Like [`Nf4Ind`](Self::Nf4Ind), but the PDFs of charm and anti-charm quarks are the same at
    /// FK table scale. PDF sets that make this assumption are NNPDF4.0 and NNPDF3.1 at fitting
    /// scale.
    Nf4Sym,
    /// Like [`Nf4Ind`](Self::Nf4Ind), but the PDFs of charm and anti-charm quarks are zero at FK
    /// table scale. PDF sets that make this assumption are MSHT20 and NNPDF3.0 at fitting scale.
    Nf3Ind,
    /// Like [`Nf3Ind`](Self::Nf3Ind), but the PDFs of strange and anti-strange are the same at FK
    /// table scale. A PDF set that makes this assumption is CT18 at fitting scale.
    Nf3Sym,
}

/// Error type when trying to construct [`FkAssumptions`] with a string.
#[derive(Debug, Eq, Error, PartialEq)]
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
                Self::Nf6Ind => "Nf6Ind",
                Self::Nf6Sym => "Nf6Sym",
                Self::Nf5Ind => "Nf5Ind",
                Self::Nf5Sym => "Nf5Sym",
                Self::Nf4Ind => "Nf4Ind",
                Self::Nf4Sym => "Nf4Sym",
                Self::Nf3Ind => "Nf3Ind",
                Self::Nf3Sym => "Nf3Sym",
            }
        )
    }
}

impl FromStr for FkAssumptions {
    type Err = UnknownFkAssumption;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
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
                    variant: s.to_owned(),
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

    // TODO: when trying to convert the following function to `const` as per clippy's suggestion,
    // the compiler errors out with: 'the destructor for this type cannot be evaluated in constant
    // functions'

    /// Converts the `FkTable` back to a [`Grid`].
    #[must_use]
    pub fn into_grid(self) -> Grid {
        self.grid
    }

    /// Returns the FK table represented as a four-dimensional array indexed by `bin`, `channel`,
    /// `x1` and `x2`, in this order.
    ///
    /// # Panics
    ///
    /// TODO
    #[must_use]
    pub fn table(&self) -> Array4<f64> {
        let has_pdf1 = self.grid.convolutions()[0] != Convolution::None;
        let has_pdf2 = self.grid.convolutions()[1] != Convolution::None;
        let x_grid = self.x_grid();

        let mut result = Array4::zeros((
            self.bins(),
            self.grid.channels().len(),
            if has_pdf1 { x_grid.len() } else { 1 },
            if has_pdf2 { x_grid.len() } else { 1 },
        ));

        for ((_, bin, channel), subgrid) in self.grid().subgrids().indexed_iter() {
            let indices1 = if has_pdf1 {
                subgrid
                    .x1_grid()
                    .iter()
                    .map(|&s| x_grid.iter().position(|&x| approx_eq!(f64, s, x, ulps = 2)))
                    .collect::<Option<_>>()
                    .unwrap()
            } else {
                vec![0]
            };
            let indices2 = if has_pdf2 {
                subgrid
                    .x2_grid()
                    .iter()
                    .map(|&s| x_grid.iter().position(|&x| approx_eq!(f64, s, x, ulps = 2)))
                    .collect::<Option<_>>()
                    .unwrap()
            } else {
                vec![0]
            };

            for ((_, ix1, ix2), value) in subgrid.indexed_iter() {
                result[[bin, channel, indices1[ix1], indices2[ix2]]] = value;
            }
        }

        result
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

    /// Return the channel definition for this `FkTable`. All factors are `1.0`.
    #[must_use]
    pub fn channels(&self) -> Vec<(i32, i32)> {
        self.grid
            .channels()
            .iter()
            .map(|entry| (entry.entry()[0].0, entry.entry()[0].1))
            .collect()
    }

    /// Returns the single `muf2` scale of this `FkTable`.
    #[must_use]
    pub fn muf2(&self) -> f64 {
        if let &[muf2] = &self.grid.evolve_info(&[true]).fac1[..] {
            muf2
        } else {
            // every `FkTable` has only a single factorization scale
            unreachable!()
        }
    }

    /// Returns the x grid that all subgrids for all hadronic initial states share.
    #[must_use]
    pub fn x_grid(&self) -> Vec<f64> {
        self.grid.evolve_info(&[true]).x1
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

    /// Convolve the FK-table. This method has fewer arguments than [`Grid::convolve`], because
    /// FK-tables have all orders merged together and do not support scale variations.
    pub fn convolve(
        &self,
        lumi_cache: &mut LumiCache,
        bin_indices: &[usize],
        channel_mask: &[bool],
    ) -> Vec<f64> {
        self.grid
            .convolve(lumi_cache, &[], bin_indices, channel_mask, &[(1.0, 1.0)])
    }

    /// Set a metadata key-value pair
    pub fn set_key_value(&mut self, key: &str, value: &str) {
        self.grid.set_key_value(key, value);
    }

    /// Optimizes the storage of FK tables based of assumptions of the PDFs at the FK table's
    /// scale.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn optimize(&mut self, assumptions: FkAssumptions) {
        let mut add = Vec::new();

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

        self.grid.rewrite_channels(&add, &[]);

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

        for subgrid in grid.subgrids() {
            if subgrid.is_empty() {
                continue;
            }

            let mu2_grid = subgrid.mu2_grid();

            if mu2_grid.len() > 1 {
                return Err(TryFromGridError::MultipleScales);
            }

            if muf2 < 0.0 {
                muf2 = mu2_grid[0].fac;
            } else if muf2 != mu2_grid[0].fac {
                return Err(TryFromGridError::MultipleScales);
            }
        }

        for channel in grid.channels() {
            let entry = channel.entry();

            if entry.len() != 1 || entry[0].2 != 1.0 {
                return Err(TryFromGridError::InvalidChannel);
            }
        }

        if (1..grid.channels().len())
            .any(|i| grid.channels()[i..].contains(&grid.channels()[i - 1]))
        {
            return Err(TryFromGridError::InvalidChannel);
        }

        Ok(Self { grid })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fk_assumptions_try_from() {
        assert_eq!(FkAssumptions::from_str("Nf6Ind"), Ok(FkAssumptions::Nf6Ind));
        assert_eq!(FkAssumptions::from_str("Nf6Sym"), Ok(FkAssumptions::Nf6Sym));
        assert_eq!(FkAssumptions::from_str("Nf5Ind"), Ok(FkAssumptions::Nf5Ind));
        assert_eq!(FkAssumptions::from_str("Nf5Sym"), Ok(FkAssumptions::Nf5Sym));
        assert_eq!(FkAssumptions::from_str("Nf4Ind"), Ok(FkAssumptions::Nf4Ind));
        assert_eq!(FkAssumptions::from_str("Nf4Sym"), Ok(FkAssumptions::Nf4Sym));
        assert_eq!(FkAssumptions::from_str("Nf3Ind"), Ok(FkAssumptions::Nf3Ind));
        assert_eq!(FkAssumptions::from_str("Nf3Sym"), Ok(FkAssumptions::Nf3Sym));
        assert_eq!(
            FkAssumptions::from_str("XXXXXX"),
            Err(UnknownFkAssumption {
                variant: "XXXXXX".to_owned()
            })
        );
    }

    #[test]
    fn fk_assumptions_display() {
        assert_eq!(format!("{}", FkAssumptions::Nf6Ind), "Nf6Ind");
        assert_eq!(format!("{}", FkAssumptions::Nf6Sym), "Nf6Sym");
        assert_eq!(format!("{}", FkAssumptions::Nf5Ind), "Nf5Ind");
        assert_eq!(format!("{}", FkAssumptions::Nf5Sym), "Nf5Sym");
        assert_eq!(format!("{}", FkAssumptions::Nf4Ind), "Nf4Ind");
        assert_eq!(format!("{}", FkAssumptions::Nf4Sym), "Nf4Sym");
        assert_eq!(format!("{}", FkAssumptions::Nf3Ind), "Nf3Ind");
        assert_eq!(format!("{}", FkAssumptions::Nf3Sym), "Nf3Sym");
    }
}
