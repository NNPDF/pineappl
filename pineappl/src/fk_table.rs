//! Provides the [`FkTable`] type.

use super::boc::{Channel, Kinematics, Order};
use super::convolutions::ConvolutionCache;
use super::error::{Error, Result};
use super::grid::Grid;
use super::pids::OptRules;
use super::subgrid::{self, EmptySubgridV1, Subgrid};
use ndarray::{s, ArrayD};
use std::collections::BTreeMap;
use std::fmt::{self, Display, Formatter};
use std::iter;
use std::str::FromStr;

/// Structure implementing FK tables. These are special [`Grid`]s, for which the following
/// additional guarantees are given:
///
/// - all subgrids of the grid evaluate the convolution functions at a single factorization and
///   fragmentation scale given by [`FkTable::fac0`] and [`FkTable::frg0`], respectively.
/// - all subgrids share the same `x` grid. See [`FkTable::x_grid`].
/// - the channel definitions are *simple*, meaning that every entry consists of a single tuple of
///   partons with trivial factor `1.0`, and all tuples are distinct from each other. See
///   [`Grid::channels`].
/// - the FK table's grid contains only a single [`Order`], whose exponents are all zero.
#[repr(transparent)]
pub struct FkTable {
    grid: Grid,
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
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
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
                return Err(Error::General(format!(
                    "unknown variant for FkAssumptions: {s}",
                )));
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

    /// Return the FK-table represented as an n-dimensional array indexed by `bin`, `channel`
    /// and parton fractions, in this order.
    ///
    /// # Panics
    ///
    /// TODO
    #[must_use]
    pub fn table(&self) -> ArrayD<f64> {
        let x_grid = self.x_grid();

        let mut dim = vec![self.grid.bwfl().len(), self.grid.channels().len()];
        dim.extend(iter::repeat(x_grid.len()).take(self.grid.convolutions().len()));
        let mut idx = vec![0; dim.len()];
        let mut result = ArrayD::zeros(dim);

        for ((_, bin, channel), subgrid) in
            self.grid().subgrids().indexed_iter().filter(|(_, subgrid)|
                // skip empty subgrids, because they have no node values
                !subgrid.is_empty())
        {
            let indices: Vec<Vec<_>> = self
                .grid
                .convolutions()
                .iter()
                .enumerate()
                .map(|(index, _)| {
                    subgrid
                        .node_values()
                        .iter()
                        .zip(self.grid.kinematics())
                        .find_map(|(node_values, kin)| {
                            matches!(kin, Kinematics::X(i) if *i == index).then(|| {
                                node_values
                                    .iter()
                                    .map(|&s| {
                                        x_grid
                                            .iter()
                                            .position(|&x| subgrid::node_value_eq(s, x))
                                            // UNWRAP: must be guaranteed by the grid constructor
                                            .unwrap()
                                    })
                                    .collect()
                            })
                        })
                        // UNWRAP: must be guaranteed by the grid constructor
                        .unwrap()
                })
                .collect();

            for (index, value) in subgrid.indexed_iter() {
                assert_eq!(index[0], 0);
                idx[0] = bin;
                idx[1] = channel;
                for i in 2..result.shape().len() {
                    idx[i] = indices[i - 2][index[i - 1]];
                }
                result[idx.as_slice()] = value;
            }
        }

        result
    }

    /// Return the channel definition for this `FkTable`. All factors are `1.0`.
    #[must_use]
    pub fn channels(&self) -> Vec<Vec<i32>> {
        self.grid
            .channels()
            .iter()
            .map(|entry| entry.entry()[0].0.clone())
            .collect()
    }

    /// Return the squared factorization scale.
    #[must_use]
    pub fn fac0(&self) -> Option<f64> {
        let fac1 = self.grid.evolve_info(&[true]).fac1;

        if let [fac0] = fac1[..] {
            Some(fac0)
        } else {
            // every `FkTable` has either a single factorization scale or none
            assert!(fac1.is_empty());

            None
        }
    }

    /// Return the squared fragmentation scale.
    #[must_use]
    pub fn frg0(&self) -> Option<f64> {
        let frg1 = self.grid.evolve_info(&[true]).frg1;

        if let [frg0] = frg1[..] {
            Some(frg0)
        } else {
            // every `FkTable` has either a single fragmentation scale or none
            assert!(frg1.is_empty());

            None
        }
    }

    /// Return the metadata of this FK-table.
    pub fn metadata_mut(&mut self) -> &mut BTreeMap<String, String> {
        self.grid.metadata_mut()
    }

    /// Returns the x grid that all subgrids for all hadronic initial states share.
    #[must_use]
    pub fn x_grid(&self) -> Vec<f64> {
        self.grid.evolve_info(&[true]).x1
    }

    /// Convolve the FK-table. This method has fewer arguments than [`Grid::convolve`], because
    /// FK-tables have all orders merged together and do not support scale variations.
    pub fn convolve(
        &self,
        convolution_cache: &mut ConvolutionCache,
        bin_indices: &[usize],
        channel_mask: &[bool],
    ) -> Vec<f64> {
        self.grid.convolve(
            convolution_cache,
            &[],
            bin_indices,
            channel_mask,
            &[(1.0, 1.0, 1.0)],
        )
    }

    /// Optimize the size of this FK-table by throwing away heavy quark flavors assumed to be zero
    /// at the FK-table's scales and calling [`Grid::optimize`].
    pub fn optimize(&mut self, assumptions: FkAssumptions) {
        let OptRules(sum, delete) = self.grid.pid_basis().opt_rules(assumptions);

        for idx in 0..self.grid.channels().len() {
            let &[(ref pids, factor)] = self.grid.channels()[idx].entry() else {
                // every FK-table must have a trivial channel definition
                unreachable!()
            };
            let mut pids = pids.clone();

            for pid in &mut pids {
                if delete.contains(pid) {
                    for subgrid in self.grid.subgrids_mut().slice_mut(s![.., .., idx]) {
                        *subgrid = EmptySubgridV1.into();
                    }
                } else if let Some(replace) = sum
                    .iter()
                    .find_map(|&(search, replace)| (*pid == search).then_some(replace))
                {
                    *pid = replace;
                }
            }

            self.grid.channels_mut()[idx] = Channel::new(vec![(pids, factor)]);
        }

        self.grid.optimize();

        // store the assumption so that we can check it later on
        self.grid
            .metadata_mut()
            .insert("fk_assumptions".to_owned(), assumptions.to_string());
    }
}

impl TryFrom<Grid> for FkTable {
    type Error = Error;

    fn try_from(grid: Grid) -> Result<Self> {
        let mut muf2 = -1.0;

        if grid.orders()
            != [Order {
                alphas: 0,
                alpha: 0,
                logxir: 0,
                logxif: 0,
                logxia: 0,
            }]
        {
            return Err(Error::General("multiple orders detected".to_owned()));
        }

        for subgrid in grid.subgrids() {
            if subgrid.is_empty() {
                continue;
            }

            let [fac] = grid
                .scales()
                .fac
                .calc(&subgrid.node_values(), grid.kinematics())[..]
            else {
                return Err(Error::General("multiple scales detected".to_owned()));
            };

            if muf2 < 0.0 {
                muf2 = fac;
            } else if !subgrid::node_value_eq(muf2, fac) {
                return Err(Error::General("multiple scales detected".to_owned()));
            }
        }

        for channel in grid.channels() {
            let entry = channel.entry();

            if entry.len() != 1 || !subgrid::node_value_eq(entry[0].1, 1.0) {
                return Err(Error::General(
                    "complicated channel function detected".to_owned(),
                ));
            }
        }

        if (1..grid.channels().len())
            .any(|i| grid.channels()[i..].contains(&grid.channels()[i - 1]))
        {
            return Err(Error::General(
                "complicated channel function detected".to_owned(),
            ));
        }

        Ok(Self { grid })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fk_assumptions_try_from() {
        assert_eq!(
            FkAssumptions::from_str("Nf6Ind").unwrap(),
            FkAssumptions::Nf6Ind
        );
        assert_eq!(
            FkAssumptions::from_str("Nf6Sym").unwrap(),
            FkAssumptions::Nf6Sym
        );
        assert_eq!(
            FkAssumptions::from_str("Nf5Ind").unwrap(),
            FkAssumptions::Nf5Ind
        );
        assert_eq!(
            FkAssumptions::from_str("Nf5Sym").unwrap(),
            FkAssumptions::Nf5Sym
        );
        assert_eq!(
            FkAssumptions::from_str("Nf4Ind").unwrap(),
            FkAssumptions::Nf4Ind
        );
        assert_eq!(
            FkAssumptions::from_str("Nf4Sym").unwrap(),
            FkAssumptions::Nf4Sym
        );
        assert_eq!(
            FkAssumptions::from_str("Nf3Ind").unwrap(),
            FkAssumptions::Nf3Ind
        );
        assert_eq!(
            FkAssumptions::from_str("Nf3Sym").unwrap(),
            FkAssumptions::Nf3Sym
        );
        assert_eq!(
            FkAssumptions::from_str("XXXXXX").unwrap_err().to_string(),
            "unknown variant for FkAssumptions: XXXXXX"
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
