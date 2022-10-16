//! Module containing all traits and supporting structures for grids.

use super::bin::{BinInfo, BinLimits, BinRemapper};
use super::empty_subgrid::EmptySubgridV1;
use super::evolution::{self, OperatorInfo};
use super::fk_table::FkTable;
use super::import_only_subgrid::ImportOnlySubgridV2;
use super::lagrange_subgrid::{LagrangeSparseSubgridV1, LagrangeSubgridV1, LagrangeSubgridV2};
use super::lumi::{LumiCache, LumiEntry};
use super::lumi_entry;
use super::ntuple_subgrid::NtupleSubgridV1;
use super::pids;
use super::sparse_array3::SparseArray3;
use super::subgrid::{ExtraSubgridParams, Mu2, Subgrid, SubgridEnum, SubgridParams};
use float_cmp::approx_eq;
use git_version::git_version;
use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
use lz4_flex::frame::{FrameDecoder, FrameEncoder};
use ndarray::{s, Array3, Array5, Axis, Dimension};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::convert::{TryFrom, TryInto};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::iter;
use std::mem;
use std::ops::Range;
use thiserror::Error;

// TODO: when possible change the types from `u32` to `u8` to change `try_into` to `into`

/// Coupling powers for each grid.
#[derive(Clone, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
pub struct Order {
    /// Exponent of the strong coupling.
    pub alphas: u32,
    /// Exponent of the electromagnetic coupling.
    pub alpha: u32,
    /// Exponent of the logarithm of the scale factor of the renomalization scale.
    pub logxir: u32,
    /// Exponent of the logarithm of the scale factor of the factorization scale.
    pub logxif: u32,
}

impl Ord for Order {
    fn cmp(&self, other: &Self) -> Ordering {
        // sort leading orders before next-to-leading orders, then the lowest power in alpha, the
        // rest lexicographically
        (self.alphas + self.alpha)
            .cmp(&(other.alphas + other.alpha))
            .then((self.alpha, self.logxir, self.logxif).cmp(&(
                other.alpha,
                other.logxir,
                other.logxif,
            )))
    }
}

impl PartialOrd for Order {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Order {
    /// Constructor. This function mainly exists to have a way of constructing `Order` that is less
    /// verbose.
    #[must_use]
    pub const fn new(alphas: u32, alpha: u32, logxir: u32, logxif: u32) -> Self {
        Self {
            alphas,
            alpha,
            logxir,
            logxif,
        }
    }

    /// Return a mask suitable to pass as the `order_mask` parameter of [`Grid::convolute`]. The
    /// selection of `orders` is controlled using the `max_as` and `max_al` parameters, for
    /// instance setting `max_as = 1` and `max_al = 0` selects the LO QCD only, `max_as = 2` and
    /// `max_al = 0` the NLO QCD; setting `max_as = 3` and `max_al = 2` would select all NLOs, and
    /// the NNLO QCD.
    ///
    /// # Example
    ///
    /// In the case of Drell—Yan, there are the following orders:
    ///
    /// - exactly one leading order (LO),
    /// - two next-to-leading orders (NLO), which are
    ///   - the NLO QCD and
    ///   - the NLO EW, and
    /// - three next-to-next-to-leading orders (NNLO),
    ///   - the NNLO QCD,
    ///   - the NNLO EW, and finally
    ///   - the mixed NNLO QCD—EW.
    ///
    /// ```rust
    /// use pineappl::grid::Order;
    ///
    /// let orders = [
    ///     Order::new(0, 2, 0, 0), //   LO        :          alpha^2
    ///     Order::new(1, 2, 0, 0), //  NLO QCD    : alphas   alpha^2
    ///     Order::new(0, 3, 0, 0), //  NLO  EW    :          alpha^3
    ///     Order::new(2, 2, 0, 0), // NNLO QCD    : alphas^2 alpha^2
    ///     Order::new(1, 3, 0, 0), // NNLO QCD—EW : alphas   alpha^3
    ///     Order::new(0, 4, 0, 0), // NNLO EW     :          alpha^4
    /// ];
    ///
    /// // LO EW
    /// assert_eq!(Order::create_mask(&orders, 0, 1), [true, false, false, false, false, false]);
    /// // LO QCD
    /// assert_eq!(Order::create_mask(&orders, 1, 0), [true, false, false, false, false, false]);
    /// // LO
    /// assert_eq!(Order::create_mask(&orders, 1, 1), [true, false, false, false, false, false]);
    /// // NLO QCD
    /// assert_eq!(Order::create_mask(&orders, 2, 0), [true, true, false, false, false, false]);
    /// // NLO EW
    /// assert_eq!(Order::create_mask(&orders, 0, 2), [true, false, true, false, false, false]);
    /// // NNLO QCD
    /// assert_eq!(Order::create_mask(&orders, 3, 0), [true, true, false, true, false, false]);
    /// // NNLO EW
    /// assert_eq!(Order::create_mask(&orders, 0, 3), [true, false, true, false, false, true]);
    /// ```
    ///
    /// Although not shown in the example above, orders containing non-zero powers of logarithms
    /// are selected as well:
    ///
    /// ```rust
    /// use pineappl::grid::Order;
    ///
    /// let orders = [
    ///     Order::new(0, 2, 0, 0), //  LO         :        alpha^2
    ///     Order::new(1, 2, 0, 0), //  NLO QCD    : alphas alpha^2
    ///     Order::new(1, 2, 1, 0), //  NLO QCD    : alphas alpha^2 logxif
    ///     Order::new(0, 3, 0, 0), //  NLO  EW    :        alpha^3
    ///     Order::new(0, 3, 1, 0), //  NLO  EW    :        alpha^3 logxif
    /// ];
    ///
    /// assert_eq!(Order::create_mask(&orders, 0, 2), [true, false, false, true, true]);
    /// ```
    ///
    /// For the more complicated example of top-pair production one can see the difference between
    /// the selection for different LOs:
    ///
    /// ```rust
    /// use pineappl::grid::Order;
    ///
    /// let orders = [
    ///     Order::new(2, 0, 0, 0), //   LO QCD    : alphas^2
    ///     Order::new(1, 1, 0, 0), //   LO QCD—EW : alphas   alpha
    ///     Order::new(0, 2, 0, 0), //   LO  EW    :          alpha^2
    ///     Order::new(3, 0, 0, 0), //  NLO QCD    : alphas^3
    ///     Order::new(2, 1, 0, 0), //  NLO QCD—EW : alphas^2 alpha
    ///     Order::new(1, 2, 0, 0), //  NLO QCD—EW : alphas   alpha^2
    ///     Order::new(0, 3, 0, 0), //  NLO EW     :          alpha^3
    /// ];
    ///
    /// // LO EW
    /// assert_eq!(Order::create_mask(&orders, 0, 1), [false, false, true, false, false, false, false]);
    /// // LO QCD
    /// assert_eq!(Order::create_mask(&orders, 1, 0), [true, false, false, false, false, false, false]);
    /// // LO
    /// assert_eq!(Order::create_mask(&orders, 1, 1), [true, true, true, false, false, false, false]);
    /// ```
    #[must_use]
    pub fn create_mask(orders: &[Self], max_as: u32, max_al: u32) -> Vec<bool> {
        // smallest sum of alphas and alpha
        let lo = orders
            .iter()
            .map(|Self { alphas, alpha, .. }| alphas + alpha)
            .min()
            .unwrap_or_default();

        // all leading orders, without logarithms
        let leading_orders: Vec<_> = orders
            .iter()
            .filter(|Self { alphas, alpha, .. }| alphas + alpha == lo)
            .cloned()
            .collect();

        let lo_as = leading_orders
            .iter()
            .map(|Self { alphas, .. }| *alphas)
            .max()
            .unwrap_or_default();
        let lo_al = leading_orders
            .iter()
            .map(|Self { alpha, .. }| *alpha)
            .max()
            .unwrap_or_default();

        let max = max_as.max(max_al);
        let min = max_as.min(max_al);

        orders
            .iter()
            .map(|Self { alphas, alpha, .. }| {
                let pto = alphas + alpha - lo;

                alphas + alpha < min + lo
                    || (alphas + alpha < max + lo
                        && match max_as.cmp(&max_al) {
                            Ordering::Greater => lo_as + pto == *alphas,
                            Ordering::Less => lo_al + pto == *alpha,
                            Ordering::Equal => false,
                        })
            })
            .collect()
    }
}

/// This structure represents a position (`x1`, `x2`, `q2`) in a `Subgrid` together with a
/// corresponding `weight`. The type `W` can either be a `f64` or `()`, which is used when multiple
/// weights should be signaled.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct Ntuple<W> {
    /// Momentum fraction of the first parton.
    pub x1: f64,
    /// Momentum fraction of the second parton.
    pub x2: f64,
    /// Squared scale.
    pub q2: f64,
    /// Weight of this entry.
    pub weight: W,
}

/// Error returned when merging two grids fails.
#[derive(Debug, Error)]
pub enum GridError {
    /// Returned when trying to merge two `Grid` objects with incompatible bin limits.
    #[error(transparent)]
    InvalidBinLimits(super::bin::MergeBinError),
    /// Returned if the number of bins in the grid and in the remapper do not agree.
    #[error("the remapper has {remapper_bins} bins, but the grid has {grid_bins}")]
    BinNumberMismatch {
        /// Number of bins in the grid.
        grid_bins: usize,
        /// Number of bins in the remapper.
        remapper_bins: usize,
    },
    /// Returned when it was tried to merge bins that are non-consecutive.
    #[error(transparent)]
    MergeBinError(super::bin::MergeBinError),
    /// Returned when trying to construct a `Grid` using an unknown subgrid type.
    #[error("tried constructing a Grid with unknown Subgrid type `{0}`")]
    UnknownSubgridType(String),
    /// Returned when failed to read a Grid.
    #[error(transparent)]
    ReadFailure(bincode::Error),
    /// Returned when failed to write a Grid.
    #[error(transparent)]
    WriteFailure(bincode::Error),
    /// Returned while performing IO operations.
    #[error(transparent)]
    IoFailure(io::Error),
    /// Returned when trying to read a `PineAPPL` file with file format version that is not
    /// supported.
    #[error("the file version is {file_version}, but supported is only {supported_version}")]
    FileVersionMismatch {
        /// File format version of the file read.
        file_version: u64,
        /// Maximum supported file format version for this library.
        supported_version: u64,
    },
    /// Returned from [`Grid::evolve`] if the evolution failed.
    #[error("failed to evolve grid: {0}")]
    EvolutionFailure(String),
}

#[derive(Clone, Deserialize, Serialize)]
struct Mmv1 {}

#[derive(Clone, Deserialize, Serialize)]
struct Mmv2 {
    remapper: Option<BinRemapper>,
    key_value_db: HashMap<String, String>,
}

#[derive(Clone, Deserialize, Serialize)]
struct Mmv3 {
    remapper: Option<BinRemapper>,
    key_value_db: HashMap<String, String>,
    subgrid_template: SubgridEnum,
}

impl Default for Mmv2 {
    fn default() -> Self {
        Self {
            remapper: None,
            key_value_db: [
                (
                    "pineappl_gitversion".to_owned(),
                    git_version!(
                        args = ["--always", "--dirty", "--long", "--tags"],
                        cargo_prefix = "cargo:",
                        fallback = "unknown"
                    )
                    .to_owned(),
                ),
                // by default we assume there are protons in the initial state
                ("initial_state_1".to_owned(), "2212".to_owned()),
                ("initial_state_2".to_owned(), "2212".to_owned()),
            ]
            .iter()
            .cloned()
            .collect(),
        }
    }
}

impl Mmv3 {
    fn new(subgrid_template: SubgridEnum) -> Self {
        Self {
            remapper: None,
            key_value_db: [
                (
                    "pineappl_gitversion".to_owned(),
                    git_version!(
                        args = ["--always", "--dirty", "--long", "--tags"],
                        cargo_prefix = "cargo:",
                        fallback = "unknown"
                    )
                    .to_owned(),
                ),
                // by default we assume there are protons in the initial state
                ("initial_state_1".to_owned(), "2212".to_owned()),
                ("initial_state_2".to_owned(), "2212".to_owned()),
            ]
            .iter()
            .cloned()
            .collect(),
            subgrid_template,
        }
    }
}

#[derive(Clone, Deserialize, Serialize)]
enum MoreMembers {
    V1(Mmv1),
    V2(Mmv2),
    V3(Mmv3),
}

impl MoreMembers {
    fn upgrade(&mut self) {
        match self {
            Self::V1(_) => {
                *self = Self::V2(Mmv2::default());
            }
            Self::V2(_) | Self::V3(_) => {}
        }
    }
}

/// Information required to calculate the evolution kernel operators (EKO) to perform a conversion
/// of a [`Grid`] using [`Grid::convolute_eko`] to an [`FkTable`].
pub struct GridAxes {
    /// Interpolation grid in x of the `Grid`.
    pub x_grid: Vec<f64>,
    /// Parton IDs used in the grid.
    pub pids: Vec<i32>,
    /// Interpolation grid for the renormalization scale of the `Grid`.
    pub mur2_grid: Vec<f64>,
    /// Interpolation grid for the factorization scale of the `Grid`.
    pub muf2_grid: Vec<f64>,
}

/// Extra information required to perform the conversion of a [`Grid`] to an [`FkTable`] using
/// [`Grid::convolute_eko`].
pub struct EkoInfo {
    /// Scale of the FkTable.
    pub muf2_0: f64,
    /// Strong coupling constants for the renormalization scales in the same ordering as given in
    /// [`GridAxes`].
    pub alphas: Vec<f64>,
    /// Renormalization scale variation.
    pub xir: f64,
    /// Factorization scale variation.
    pub xif: f64,
    /// Interpolation grid in x of the `FkTable`.
    pub target_x_grid: Vec<f64>,
    /// Parton IDs for the `FkTable`.
    pub target_pids: Vec<i32>,
    /// axes shared with the process grid
    pub grid_axes: GridAxes,
    /// TODO: replace this member with the actual data
    pub lumi_id_types: String,
}

/// Main data structure of `PineAPPL`. This structure contains a `Subgrid` for each `LumiEntry`,
/// bin, and coupling order it was created with.
#[derive(Clone, Deserialize, Serialize)]
pub struct Grid {
    subgrids: Array3<SubgridEnum>,
    lumi: Vec<LumiEntry>,
    bin_limits: BinLimits,
    orders: Vec<Order>,
    subgrid_params: SubgridParams,
    more_members: MoreMembers,
}

impl Grid {
    /// Constructor.
    #[must_use]
    pub fn new(
        lumi: Vec<LumiEntry>,
        orders: Vec<Order>,
        bin_limits: Vec<f64>,
        subgrid_params: SubgridParams,
    ) -> Self {
        Self {
            subgrids: Array3::from_shape_simple_fn(
                (orders.len(), bin_limits.len() - 1, lumi.len()),
                || EmptySubgridV1::default().into(),
            ),
            orders,
            lumi,
            bin_limits: BinLimits::new(bin_limits),
            more_members: MoreMembers::V3(Mmv3::new(
                LagrangeSubgridV2::new(&subgrid_params, &ExtraSubgridParams::from(&subgrid_params))
                    .into(),
            )),
            subgrid_params,
        }
    }

    /// Constructor. This function can be used like `new`, but the additional parameter
    /// `subgrid_type` selects the underlying `Subgrid` type. Supported values are:
    /// - `LagrangeSubgrid`
    /// - `LagrangeSparseSubgrid`
    /// - `NtupleSubgrid`
    ///
    /// # Errors
    ///
    /// If `subgrid_type` is none of the values listed above, an error is returned.
    pub fn with_subgrid_type(
        lumi: Vec<LumiEntry>,
        orders: Vec<Order>,
        bin_limits: Vec<f64>,
        subgrid_params: SubgridParams,
        extra: ExtraSubgridParams,
        subgrid_type: &str,
    ) -> Result<Self, GridError> {
        let subgrid_template: SubgridEnum = match subgrid_type {
            "LagrangeSubgrid" | "LagrangeSubgridV2" => {
                LagrangeSubgridV2::new(&subgrid_params, &extra).into()
            }
            "LagrangeSubgridV1" => LagrangeSubgridV1::new(&subgrid_params).into(),
            "NtupleSubgrid" => NtupleSubgridV1::new().into(),
            "LagrangeSparseSubgrid" => LagrangeSparseSubgridV1::new(&subgrid_params).into(),
            _ => return Err(GridError::UnknownSubgridType(subgrid_type.to_string())),
        };

        Ok(Self {
            subgrids: Array3::from_shape_simple_fn(
                (orders.len(), bin_limits.len() - 1, lumi.len()),
                || EmptySubgridV1::default().into(),
            ),
            orders,
            lumi,
            bin_limits: BinLimits::new(bin_limits),
            subgrid_params,
            more_members: MoreMembers::V3(Mmv3::new(subgrid_template)),
        })
    }

    fn pdg_lumi(&self) -> Cow<[LumiEntry]> {
        if let Some(key_values) = self.key_values() {
            if let Some(lumi_id_types) = key_values.get("lumi_id_types") {
                match lumi_id_types.as_str() {
                    "pdg_mc_ids" => {}
                    "evol" => {
                        return self
                            .lumi
                            .iter()
                            .map(|entry| LumiEntry::translate(entry, &pids::evol_to_pdg_mc_ids))
                            .collect();
                    }
                    _ => unimplemented!(),
                }
            }
        }

        Cow::Borrowed(self.lumi())
    }

    /// Perform a convolution using the PDFs and strong coupling in `lumi_cache`, and only
    /// selecting only the orders, bins and luminosities corresponding to `order_mask`,
    /// `bin_indices` and `lumi_mask`. A variation of the scales
    /// is performed using the factors in `xi`; the first factor varies the renormalization scale,
    /// the second the factorization scale. Note that for the variation to be trusted all non-zero
    /// log-grids must be contained.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn convolute(
        &self,
        lumi_cache: &mut LumiCache,
        order_mask: &[bool],
        bin_indices: &[usize],
        lumi_mask: &[bool],
        xi: &[(f64, f64)],
    ) -> Vec<f64> {
        lumi_cache.setup(self, xi).unwrap();

        let bin_indices = if bin_indices.is_empty() {
            (0..self.bin_info().bins()).collect()
        } else {
            bin_indices.to_vec()
        };
        let mut bins: Vec<f64> = vec![0.0; bin_indices.len() * xi.len()];
        let normalizations = self.bin_info().normalizations();
        let self_lumi = self.pdg_lumi();

        for (xi_index, &(xir, xif)) in xi.iter().enumerate() {
            for ((ord, bin, lumi), subgrid) in self.subgrids.indexed_iter() {
                let order = &self.orders[ord];

                if ((order.logxir > 0) && (xir == 1.0)) || ((order.logxif > 0) && (xif == 1.0)) {
                    continue;
                }

                if (!order_mask.is_empty() && !order_mask[ord])
                    || (!lumi_mask.is_empty() && !lumi_mask[lumi])
                {
                    continue;
                }

                let bin_index = match bin_indices.iter().position(|&index| index == bin) {
                    Some(i) => i,
                    None => continue,
                };

                if subgrid.is_empty() {
                    continue;
                }

                let lumi_entry = &self_lumi[lumi];
                let mu2_grid = subgrid.mu2_grid();
                let x1_grid = subgrid.x1_grid();
                let x2_grid = subgrid.x2_grid();

                lumi_cache.set_grids(&mu2_grid, &x1_grid, &x2_grid, xir, xif);

                let mut value =
                    subgrid.convolute(&x1_grid, &x2_grid, &mu2_grid, &mut |ix1, ix2, imu2| {
                        let x1 = x1_grid[ix1];
                        let x2 = x2_grid[ix2];
                        let mut lumi = 0.0;

                        for entry in lumi_entry.entry() {
                            let xfx1 = lumi_cache.xfx1(entry.0, ix1, imu2);
                            let xfx2 = lumi_cache.xfx2(entry.1, ix2, imu2);
                            lumi += xfx1 * xfx2 * entry.2 / (x1 * x2);
                        }

                        let alphas = lumi_cache.alphas(imu2);

                        lumi *= alphas.powi(order.alphas.try_into().unwrap());
                        lumi
                    });

                if order.logxir > 0 {
                    value *= (xir * xir).ln().powi(order.logxir.try_into().unwrap());
                }

                if order.logxif > 0 {
                    value *= (xif * xif).ln().powi(order.logxif.try_into().unwrap());
                }

                bins[xi_index + xi.len() * bin_index] += value / normalizations[bin];
            }
        }

        bins
    }

    /// Convolutes a single subgrid `(order, bin, lumi)` with the PDFs strong coupling given by
    /// `xfx1`, `xfx2` and `alphas`. The convolution result is fully differentially, such that the
    /// axes of the result correspond to the values given by the subgrid `q2`, `x1` and `x2` grid
    /// values.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn convolute_subgrid(
        &self,
        lumi_cache: &mut LumiCache,
        ord: usize,
        bin: usize,
        lumi: usize,
        xir: f64,
        xif: f64,
    ) -> Array3<f64> {
        lumi_cache.setup(self, &[(xir, xif)]).unwrap();

        let normalizations = self.bin_info().normalizations();
        let self_lumi = self.pdg_lumi();

        let subgrid = &self.subgrids[[ord, bin, lumi]];
        let order = &self.orders[ord];

        let lumi_entry = &self_lumi[lumi];
        let mu2_grid = subgrid.mu2_grid();
        let x1_grid = subgrid.x1_grid();
        let x2_grid = subgrid.x2_grid();

        lumi_cache.set_grids(&mu2_grid, &x1_grid, &x2_grid, xir, xif);

        let mut array = Array3::zeros((mu2_grid.len(), x1_grid.len(), x2_grid.len()));

        for ((imu2, ix1, ix2), value) in subgrid.iter() {
            let x1 = x1_grid[ix1];
            let x2 = x2_grid[ix2];
            let mut lumi = 0.0;

            for entry in lumi_entry.entry() {
                let xfx1 = lumi_cache.xfx1(entry.0, ix1, imu2);
                let xfx2 = lumi_cache.xfx2(entry.1, ix2, imu2);
                lumi += xfx1 * xfx2 * entry.2 / (x1 * x2);
            }

            let alphas = lumi_cache.alphas(imu2);

            lumi *= alphas.powi(order.alphas.try_into().unwrap());

            array[[imu2, ix1, ix2]] = lumi * value;
        }

        if order.logxir > 0 {
            array *= (xir * xir).ln().powi(order.logxir.try_into().unwrap());
        }

        if order.logxif > 0 {
            array *= (xif * xif).ln().powi(order.logxif.try_into().unwrap());
        }

        // TODO: what's missing here are possible reweighting factors

        array /= normalizations[bin];
        array
    }

    /// Fills the grid with an ntuple for the given `order`, `observable`, and `lumi`.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn fill(&mut self, order: usize, observable: f64, lumi: usize, ntuple: &Ntuple<f64>) {
        if let Some(bin) = self.bin_limits.index(observable) {
            let subgrid = &mut self.subgrids[[order, bin, lumi]];
            if let SubgridEnum::EmptySubgridV1(_) = subgrid {
                if let MoreMembers::V3(mmv3) = &self.more_members {
                    *subgrid = mmv3.subgrid_template.clone_empty();
                } else {
                    unreachable!();
                }
            }

            subgrid.fill(ntuple);
        }
    }

    /// Construct a `Grid` by deserializing it from `reader`. Reading is buffered.
    ///
    /// # Errors
    ///
    /// If reading from the compressed or uncompressed stream fails an error is returned.
    pub fn read(reader: impl Read) -> Result<Self, GridError> {
        let mut reader = BufReader::new(reader);
        let buffer = reader.fill_buf().map_err(GridError::IoFailure)?;
        let magic_bytes: [u8; 4] = buffer[0..4].try_into().unwrap_or_else(|_| unreachable!());

        if u32::from_le_bytes(magic_bytes) == 0x18_4D_22_04 {
            Self::read_uncompressed(FrameDecoder::new(reader))
        } else {
            Self::read_uncompressed(reader)
        }
    }

    fn read_uncompressed(mut reader: impl BufRead) -> Result<Self, GridError> {
        let magic_bytes: [u8; 16] = reader.fill_buf().map_err(GridError::IoFailure)?[0..16]
            .try_into()
            .unwrap_or_else(|_| unreachable!());

        let file_version = if &magic_bytes[0..8] == b"PineAPPL" {
            reader.consume(16);
            u64::from_le_bytes(
                magic_bytes[8..16]
                    .try_into()
                    .unwrap_or_else(|_| unreachable!()),
            )
        } else {
            0
        };

        if file_version != 0 {
            return Err(GridError::FileVersionMismatch {
                file_version,
                supported_version: 0,
            });
        }

        bincode::deserialize_from(reader).map_err(GridError::ReadFailure)
    }

    /// Serializes `self` into `writer`. Writing is buffered.
    ///
    /// # Errors
    ///
    /// If writing fails an error is returned.
    pub fn write(&self, writer: impl Write) -> Result<(), GridError> {
        let mut writer = BufWriter::new(writer);
        let file_header = b"PineAPPL\0\0\0\0\0\0\0\0";

        // first write PineAPPL file header
        writer.write(file_header).map_err(GridError::IoFailure)?;

        // then serialize
        bincode::serialize_into(writer, self).map_err(GridError::WriteFailure)
    }

    /// Serializes `self` into `writer`, using LZ4 compression. Writing is buffered.
    ///
    /// # Errors
    ///
    /// If writing or compression fails an error is returned.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn write_lz4(&self, writer: impl Write) -> Result<(), GridError> {
        let mut encoder = FrameEncoder::new(writer);
        self.write(&mut encoder)?;
        // TODO: get rid of the unwrap call and return the error
        encoder.try_finish().unwrap();

        Ok(())
    }

    /// Fills the grid with events for the parton momentum fractions `x1` and `x2`, the scale `q2`,
    /// and the `order` and `observable`. The events are stored in `weights` and must be ordered as
    /// the corresponding luminosity function was created.
    pub fn fill_all(
        &mut self,
        order: usize,
        observable: f64,
        ntuple: &Ntuple<()>,
        weights: &[f64],
    ) {
        for (lumi, weight) in weights.iter().enumerate() {
            self.fill(
                order,
                observable,
                lumi,
                &Ntuple {
                    x1: ntuple.x1,
                    x2: ntuple.x2,
                    q2: ntuple.q2,
                    weight: *weight,
                },
            );
        }
    }

    /// Returns the luminosity function.
    #[must_use]
    pub fn lumi(&self) -> &[LumiEntry] {
        &self.lumi
    }

    /// Merges the bins for the corresponding range together in a single one.
    ///
    /// # Errors
    ///
    /// When the given bins are non-consecutive, an error is returned.
    pub fn merge_bins(&mut self, bins: Range<usize>) -> Result<(), GridError> {
        self.bin_limits
            .merge_bins(bins.clone())
            .map_err(GridError::MergeBinError)?;

        if let Some(remapper) = self.remapper_mut() {
            remapper
                .merge_bins(bins.clone())
                .map_err(GridError::MergeBinError)?;
        }

        let bin_count = self.bin_info().bins();
        let mut old_subgrids = mem::replace(
            &mut self.subgrids,
            Array3::from_shape_simple_fn((self.orders.len(), bin_count, self.lumi.len()), || {
                EmptySubgridV1::default().into()
            }),
        );

        for ((order, bin, lumi), subgrid) in old_subgrids.indexed_iter_mut() {
            if subgrid.is_empty() {
                continue;
            }

            if bins.contains(&bin) {
                let new_subgrid = &mut self.subgrids[[order, bins.start, lumi]];

                if new_subgrid.is_empty() {
                    mem::swap(new_subgrid, subgrid);
                } else {
                    new_subgrid.merge(subgrid, false);
                }
            } else {
                let new_bin = if bin > bins.start {
                    bin - (bins.end - bins.start) + 1
                } else {
                    bin
                };

                mem::swap(&mut self.subgrids[[order, new_bin, lumi]], subgrid);
            }
        }

        Ok(())
    }

    /// Merges the non-empty `Subgrid`s contained in `other` into `self`.
    ///
    /// # Errors
    ///
    /// If the bin limits of `self` and `other` are different and if the bin limits of `other` can
    /// not be merged with `self` an error is returned.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn merge(&mut self, mut other: Self) -> Result<(), GridError> {
        let mut new_orders: Vec<Order> = Vec::new();
        let mut new_bins = 0;
        let mut new_entries: Vec<LumiEntry> = Vec::new();

        if self.bin_info() != other.bin_info() {
            let lhs_bins = self.bin_info().bins();
            new_bins = other.bin_info().bins();

            let lhs_remapper = self.remapper_mut();
            let rhs_remapper = other.remapper();

            if let Some(lhs) = lhs_remapper {
                if let Some(rhs) = rhs_remapper {
                    lhs.merge(rhs).map_err(GridError::MergeBinError)?;

                    let a = u32::try_from(lhs_bins).unwrap_or_else(|_| unreachable!());
                    let b = u32::try_from(lhs_bins + new_bins).unwrap_or_else(|_| unreachable!());

                    self.bin_limits = BinLimits::new((0..=b).map(f64::from).collect());
                    other.bin_limits = BinLimits::new((a..=b).map(f64::from).collect());
                } else {
                    // Return an error
                    todo!();
                }
            } else if rhs_remapper.is_none() {
                self.bin_limits
                    .merge(&other.bin_limits)
                    .map_err(GridError::InvalidBinLimits)?;
            } else {
                // Return an error
                todo!();
            }
        }

        for ((i, _, k), _) in other
            .subgrids
            .indexed_iter_mut()
            .filter(|((_, _, _), subgrid)| !subgrid.is_empty())
        {
            let other_order = &other.orders[i];
            let other_entry = &other.lumi[k];

            if !self
                .orders
                .iter()
                .chain(new_orders.iter())
                .any(|x| x == other_order)
            {
                new_orders.push(other_order.clone());
            }

            if !self
                .lumi
                .iter()
                .chain(new_entries.iter())
                .any(|y| y == other_entry)
            {
                new_entries.push(other_entry.clone());
            }
        }

        if !new_orders.is_empty() || !new_entries.is_empty() || (new_bins != 0) {
            self.increase_shape(&(new_orders.len(), new_bins, new_entries.len()));
        }

        self.orders.append(&mut new_orders);
        self.lumi.append(&mut new_entries);

        for ((i, j, k), subgrid) in other
            .subgrids
            .indexed_iter_mut()
            .filter(|((_, _, _), subgrid)| !subgrid.is_empty())
        {
            let other_order = &other.orders[i];
            let other_bin = other.bin_limits.limits()[j];
            let other_entry = &other.lumi[k];

            let self_i = self.orders.iter().position(|x| x == other_order).unwrap();
            let self_j = self
                .bin_limits
                .limits()
                .iter()
                .position(|x| approx_eq!(f64, *x, other_bin, ulps = 4))
                .unwrap();
            let self_k = self.lumi.iter().position(|y| y == other_entry).unwrap();

            if self.subgrids[[self_i, self_j, self_k]].is_empty() {
                mem::swap(&mut self.subgrids[[self_i, self_j, self_k]], subgrid);
            } else {
                self.subgrids[[self_i, self_j, self_k]].merge(&mut *subgrid, false);
            }
        }

        Ok(())
    }

    fn increase_shape(&mut self, new_dim: &(usize, usize, usize)) {
        let old_dim = self.subgrids.raw_dim().into_pattern();
        let mut new_subgrids = Array3::from_shape_simple_fn(
            (
                old_dim.0 + new_dim.0,
                old_dim.1 + new_dim.1,
                old_dim.2 + new_dim.2,
            ),
            || EmptySubgridV1::default().into(),
        );

        for ((i, j, k), subgrid) in self.subgrids.indexed_iter_mut() {
            mem::swap(&mut new_subgrids[[i, j, k]], subgrid);
        }

        mem::swap(&mut self.subgrids, &mut new_subgrids);
    }

    /// Scale all subgrids by `factor`.
    pub fn scale(&mut self, factor: f64) {
        self.subgrids
            .iter_mut()
            .for_each(|subgrid| subgrid.scale(factor));
    }

    /// Scales each subgrid by a factor which is the product of the given values `alphas`, `alpha`,
    /// `logxir`, and `logxif`, each raised to the corresponding powers for each subgrid. In
    /// addition, every subgrid is scaled by a factor `global` independently of its order.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn scale_by_order(
        &mut self,
        alphas: f64,
        alpha: f64,
        logxir: f64,
        logxif: f64,
        global: f64,
    ) {
        for ((i, _, _), subgrid) in self.subgrids.indexed_iter_mut() {
            let order = &self.orders[i];
            let factor = global
                * alphas.powi(order.alphas.try_into().unwrap())
                * alpha.powi(order.alpha.try_into().unwrap())
                * logxir.powi(order.logxir.try_into().unwrap())
                * logxif.powi(order.logxif.try_into().unwrap());

            subgrid.scale(factor);
        }
    }

    /// Scales each subgrid by a bin-dependent factor given in `factors`. If a bin does not have a
    /// corresponding entry in `factors` it is not rescaled. If `factors` has more entries than
    /// there are bins the superfluous entries do not have an effect.
    pub fn scale_by_bin(&mut self, factors: &[f64]) {
        for ((_, bin, _), subgrid) in self.subgrids.indexed_iter_mut() {
            if let Some(&factor) = factors.get(bin) {
                subgrid.scale(factor);
            }
        }
    }

    /// Returns the subgrid parameters.
    #[must_use]
    pub fn orders(&self) -> &[Order] {
        &self.orders
    }

    /// Set the luminosity function for this grid.
    pub fn set_lumis(&mut self, lumis: Vec<LumiEntry>) {
        self.lumi = lumis;
    }

    /// Returns the subgrid with the specified indices `order`, `bin`, and `lumi`.
    #[must_use]
    pub fn subgrid(&self, order: usize, bin: usize, lumi: usize) -> &SubgridEnum {
        &self.subgrids[[order, bin, lumi]]
    }

    /// Returns all subgrids as an `Array3`.
    #[must_use]
    pub const fn subgrids(&self) -> &Array3<SubgridEnum> {
        &self.subgrids
    }

    /// Replaces the subgrid for the specified indices `order`, `bin`, and `lumi` with `subgrid`.
    pub fn set_subgrid(&mut self, order: usize, bin: usize, lumi: usize, subgrid: SubgridEnum) {
        self.subgrids[[order, bin, lumi]] = subgrid;
    }

    /// Sets a remapper. A remapper can change the dimensions and limits of each bin in this grid.
    /// This is useful because many Monte Carlo integrators and also `PineAPPL` do not support
    /// multi-dimensional bins. To work around the problem the multi-dimensional bins can be
    /// projected to one-dimensional bins, and the remapper can be used to restore the multi
    /// dimensionality. Furthermore, it allows to normalize each bin separately, and independently
    /// of the bin widths.
    ///
    /// # Errors
    ///
    /// Returns an error if the number of bins in the grid and in the remapper do not agree.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn set_remapper(&mut self, remapper: BinRemapper) -> Result<(), GridError> {
        if remapper.bins() != self.bin_info().bins() {
            return Err(GridError::BinNumberMismatch {
                grid_bins: self.bin_info().bins(),
                remapper_bins: remapper.bins(),
            });
        }

        self.more_members.upgrade();

        match &mut self.more_members {
            MoreMembers::V1(_) => unreachable!(),
            MoreMembers::V2(mmv2) => mmv2.remapper = Some(remapper),
            MoreMembers::V3(mmv3) => mmv3.remapper = Some(remapper),
        }

        Ok(())
    }

    /// Return the currently set remapper, if there is any.
    #[must_use]
    pub const fn remapper(&self) -> Option<&BinRemapper> {
        match &self.more_members {
            MoreMembers::V1(_) => None,
            MoreMembers::V2(mmv2) => mmv2.remapper.as_ref(),
            MoreMembers::V3(mmv3) => mmv3.remapper.as_ref(),
        }
    }

    fn remapper_mut(&mut self) -> Option<&mut BinRemapper> {
        match &mut self.more_members {
            MoreMembers::V1(_) => None,
            MoreMembers::V2(mmv2) => mmv2.remapper.as_mut(),
            MoreMembers::V3(mmv3) => mmv3.remapper.as_mut(),
        }
    }

    /// Returns all information about the bins in this grid.
    #[must_use]
    pub const fn bin_info(&self) -> BinInfo {
        BinInfo::new(&self.bin_limits, self.remapper())
    }

    /// Optimize the internal datastructures for space efficiency. This changes all subgrids of
    /// type `LagrangeSubgrid` to `LagrangeSparseSubgrid`.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn optimize(&mut self) {
        if self
            .key_values()
            .map_or(true, |map| map["initial_state_1"] == map["initial_state_2"])
        {
            self.symmetrize_lumi();
        }

        self.optimize_orders();
        self.optimize_lumi();

        for subgrid in self.subgrids.iter_mut() {
            if subgrid.is_empty() {
                *subgrid = EmptySubgridV1::default().into();
            } else {
                match subgrid {
                    SubgridEnum::LagrangeSubgridV1(grid) => {
                        let mut new_subgrid = LagrangeSparseSubgridV1::from(&*grid).into();
                        mem::swap(subgrid, &mut new_subgrid);
                    }
                    SubgridEnum::LagrangeSubgridV2(grid) => {
                        let mut new_subgrid = ImportOnlySubgridV2::from(&*grid).into();
                        mem::swap(subgrid, &mut new_subgrid);
                    }
                    SubgridEnum::EmptySubgridV1(_)
                    | SubgridEnum::LagrangeSparseSubgridV1(_)
                    | SubgridEnum::ImportOnlySubgridV1(_)
                    | SubgridEnum::ImportOnlySubgridV2(_) => {
                        // nothing to optimize here
                    }
                    SubgridEnum::NtupleSubgridV1(_) => todo!(),
                }
            }
        }
    }

    fn optimize_lumi(&mut self) {
        let mut indices: Vec<_> = (0..self.lumi.len()).rev().collect();

        // merge luminosities that are the same
        while let Some(index) = indices.pop() {
            if let Some(&other_index) = indices.iter().find(|i| self.lumi[**i] == self.lumi[index])
            {
                let (mut a, mut b) = self
                    .subgrids
                    .multi_slice_mut((s![.., .., other_index], s![.., .., index]));

                // check if in all cases the limits are compatible with merging
                for (lhs, rhs) in a.iter_mut().zip(b.iter_mut()) {
                    if !rhs.is_empty() {
                        if lhs.is_empty() {
                            // we can't merge into an EmptySubgridV1
                            *lhs = rhs.clone_empty();
                            lhs.merge(rhs, false);
                        } else if (lhs.x1_grid() == rhs.x1_grid())
                            && (lhs.x2_grid() == rhs.x2_grid())
                        {
                            lhs.merge(rhs, false);
                        } else {
                            // don't overwrite `rhs`
                            continue;
                        }

                        *rhs = EmptySubgridV1::default().into();
                    }
                }
            }
        }

        let mut keep_lumi_indices = vec![];
        let mut new_lumi_entries = vec![];

        // only keep luminosities that have non-zero factors and for which at least one subgrid is
        // non-empty
        for (lumi, entry) in self.lumi.iter().enumerate() {
            if !entry.entry().iter().all(|&(_, _, factor)| factor == 0.0)
                && !self
                    .subgrids
                    .slice(s![.., .., lumi])
                    .iter()
                    .all(Subgrid::is_empty)
            {
                keep_lumi_indices.push(lumi);
                new_lumi_entries.push(entry.clone());
            }
        }

        // only keep the previously selected subgrids
        let new_subgrids = Array3::from_shape_fn(
            (
                self.orders.len(),
                self.bin_info().bins(),
                keep_lumi_indices.len(),
            ),
            |(order, bin, new_lumi)| {
                mem::replace(
                    &mut self.subgrids[[order, bin, keep_lumi_indices[new_lumi]]],
                    EmptySubgridV1::default().into(),
                )
            },
        );

        self.lumi = new_lumi_entries;
        self.subgrids = new_subgrids;
    }

    fn optimize_orders(&mut self) {
        let mut indices: Vec<_> = (0..self.orders().len()).collect();

        while let Some(index) = indices.pop() {
            if self
                .subgrids
                .slice(s![index, .., ..])
                .iter()
                .all(Subgrid::is_empty)
            {
                self.orders.remove(index);
                self.subgrids.remove_index(Axis(0), index);
            }
        }
    }

    fn symmetrize_lumi(&mut self) {
        let mut indices: Vec<usize> = (0..self.lumi.len()).rev().collect();

        while let Some(index) = indices.pop() {
            let lumi_entry = &self.lumi[index];

            if *lumi_entry == lumi_entry.transpose() {
                // check if in all cases the limits are compatible with merging
                self.subgrids
                    .slice_mut(s![.., .., index])
                    .iter_mut()
                    .for_each(|subgrid| {
                        if !subgrid.is_empty() && (subgrid.x1_grid() == subgrid.x2_grid()) {
                            subgrid.symmetrize();
                        }
                    });
            } else if let Some((j, &other_index)) = indices
                .iter()
                .enumerate()
                .find(|(_, i)| self.lumi[**i] == lumi_entry.transpose())
            {
                indices.remove(j);

                // check if in all cases the limits are compatible with merging
                let (mut a, mut b) = self
                    .subgrids
                    .multi_slice_mut((s![.., .., index], s![.., .., other_index]));

                for (lhs, rhs) in a.iter_mut().zip(b.iter_mut()) {
                    if !rhs.is_empty() {
                        if lhs.is_empty() {
                            // we can't merge into an EmptySubgridV1
                            *lhs = rhs.clone_empty();
                            lhs.merge(rhs, true);
                        } else if (lhs.x1_grid() == rhs.x2_grid())
                            && (lhs.x2_grid() == rhs.x1_grid())
                        {
                            lhs.merge(rhs, true);
                        } else {
                            // don't overwrite `rhs`
                            continue;
                        }

                        *rhs = EmptySubgridV1::default().into();
                    }
                }
            }
        }
    }

    /// Upgrades the internal data structures to their latest versions.
    pub fn upgrade(&mut self) {
        self.more_members.upgrade();
    }

    /// Returns a map with key-value pairs, if there are any stored in this grid.
    #[must_use]
    pub const fn key_values(&self) -> Option<&HashMap<String, String>> {
        match &self.more_members {
            MoreMembers::V3(mmv3) => Some(&mmv3.key_value_db),
            MoreMembers::V2(mmv2) => Some(&mmv2.key_value_db),
            MoreMembers::V1(_) => None,
        }
    }

    /// Returns a map with key-value pairs, if there are any stored in this grid.
    ///
    /// # Panics
    ///
    /// TODO
    #[must_use]
    pub fn key_values_mut(&mut self) -> &mut HashMap<String, String> {
        self.more_members.upgrade();

        match &mut self.more_members {
            MoreMembers::V1(_) => unreachable!(),
            MoreMembers::V2(mmv2) => &mut mmv2.key_value_db,
            MoreMembers::V3(mmv3) => &mut mmv3.key_value_db,
        }
    }

    /// Sets a specific key-value pair in this grid.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn set_key_value(&mut self, key: &str, value: &str) {
        self.key_values_mut()
            .insert(key.to_owned(), value.to_owned());
    }

    /// Provide information used to compute a suitable EKO for the current grid.
    /// More specific, the `x_grid` and `muf2_grid` are extracted and checked.
    ///
    /// # Panics
    ///
    /// TODO
    #[must_use]
    pub fn axes(&self) -> Option<GridAxes> {
        // are the initial states hadrons?
        let has_pdf1 = self.has_pdf1();
        let has_pdf2 = self.has_pdf1();

        let mut mur2_grid = Vec::new();
        let mut muf2_grid = Vec::new();
        let mut x_grid = Vec::new();
        let pids = Vec::new();

        // Within each lane, that is for a specific combination of (order, bin) ...
        for lane in self.subgrids().lanes(Axis(2)) {
            // for all luminosities ...

            // the renormalization and factorization grid must be the same, ...
            if !lane
                .iter()
                .filter_map(|subgrid| (!subgrid.is_empty()).then(|| subgrid.mu2_grid()))
                .all_equal()
            {
                return None;
            }

            // the x1 grid must be the same and finally ...
            if has_pdf1
                && !lane
                    .iter()
                    .filter_map(|subgrid| (!subgrid.is_empty()).then(|| subgrid.x1_grid()))
                    .all_equal()
            {
                return None;
            }

            // the x2 grid must be the same
            if has_pdf2
                && !lane
                    .iter()
                    .filter_map(|subgrid| (!subgrid.is_empty()).then(|| subgrid.x2_grid()))
                    .all_equal()
            {
                return None;
            }

            // not all luminosities are equal (some appear only at higher orders)
            for subgrid in lane.iter() {
                mur2_grid.append(&mut subgrid.mu2_grid().iter().map(|mu2| mu2.ren).collect());
                muf2_grid.append(&mut subgrid.mu2_grid().iter().map(|mu2| mu2.fac).collect());
                if has_pdf1 {
                    x_grid.extend_from_slice(&subgrid.x1_grid());
                }
                if has_pdf2 {
                    x_grid.extend_from_slice(&subgrid.x2_grid());
                }
            }
        }

        // make grids unique
        x_grid.sort_by(|a, b| a.partial_cmp(b).unwrap());
        x_grid.dedup_by(|a, b| approx_eq!(f64, *a, *b, ulps = 64));
        mur2_grid.sort_by(|a, b| a.partial_cmp(b).unwrap());
        mur2_grid.dedup_by(|a, b| approx_eq!(f64, *a, *b, ulps = 64));
        muf2_grid.sort_by(|a, b| a.partial_cmp(b).unwrap());
        muf2_grid.dedup_by(|a, b| approx_eq!(f64, *a, *b, ulps = 64));

        Some(GridAxes {
            x_grid,
            pids, // TODO: for the time being they are just empty, but we might use them for slicing the eko
            mur2_grid,
            muf2_grid,
        })
    }

    /// Applies an evolution kernel operator (EKO) to the grids to evolve them from different
    /// values of the factorization scale to a single one given by the parameter `q2`.
    /// Using `xir` and `xif` you can trigger renormalization and factorization scale
    /// variations respectively in the grid.
    ///
    /// # Panics
    ///
    /// Panics if the parameters do not match with the given grid.
    #[must_use]
    pub fn convolute_eko(
        &self,
        operator: Array5<f64>,
        eko_info: EkoInfo,
        order_mask: &[bool],
    ) -> Option<FkTable> {
        // Check operator layout
        let dim = operator.shape();

        assert_eq!(dim[0], eko_info.grid_axes.mur2_grid.len());
        assert_eq!(dim[0], eko_info.grid_axes.muf2_grid.len());
        assert_eq!(dim[1], eko_info.target_pids.len());
        assert_eq!(dim[3], eko_info.grid_axes.pids.len());

        // swap axes around to optimize convolution
        let operator = operator.permuted_axes([3, 1, 4, 0, 2]);
        let operator = operator.as_standard_layout();

        // determine what and how many hadrons are in the initial state
        let initial_state_1 = self.initial_state_1();
        let initial_state_2 = self.initial_state_2();

        // are the initial states hadrons?
        let has_pdf1 = self.has_pdf1();
        let has_pdf2 = self.has_pdf2();

        let pids1 = if has_pdf1 {
            eko_info.grid_axes.pids.clone()
        } else {
            vec![initial_state_1]
        };
        let pids2 = if has_pdf2 {
            eko_info.grid_axes.pids.clone()
        } else {
            vec![initial_state_2]
        };
        // create target luminosities
        let tgt_pids1 = if has_pdf1 {
            eko_info.target_pids.clone()
        } else {
            vec![initial_state_1]
        };
        let tgt_pids2 = if has_pdf2 {
            eko_info.target_pids.clone()
        } else {
            vec![initial_state_2]
        };
        let lumi: Vec<_> = tgt_pids1
            .iter()
            .cartesian_product(tgt_pids2.iter())
            .map(|(a, b)| lumi_entry![*a, *b, 1.0])
            .collect();

        // create target subgrid dimensions
        let tgt_q2_grid = vec![Mu2 {
            ren: eko_info.muf2_0,
            fac: eko_info.muf2_0,
        }];
        let tgt_x1_grid = if has_pdf1 {
            eko_info.target_x_grid.clone()
        } else {
            vec![1.0]
        };
        let tgt_x2_grid = if has_pdf2 {
            eko_info.target_x_grid.clone()
        } else {
            vec![1.0]
        };

        // create target grid
        let mut result = Self {
            subgrids: Array3::from_shape_simple_fn((1, self.bin_info().bins(), lumi.len()), || {
                EmptySubgridV1::default().into()
            }),
            lumi,
            bin_limits: self.bin_limits.clone(),
            orders: vec![Order {
                alphas: 0,
                alpha: 0,
                logxir: 0,
                logxif: 0,
            }],
            subgrid_params: SubgridParams::default(),
            more_members: self.more_members.clone(),
        };
        // write additional metadata
        result.set_key_value("lumi_id_types", &eko_info.lumi_id_types);

        // collect source grid informations
        let grid_axes = self.axes()?;

        // Setup progress bar
        let bar = ProgressBar::new(
            (self.bin_info().bins() * self.lumi.len() * pids1.len() * pids2.len()) as u64,
        );
        bar.set_style(ProgressStyle::default_bar().template(
            "[{elapsed_precise}] {bar:50.cyan/blue} {pos:>7}/{len:7} - ETA: {eta_precise} {msg}",
        ));

        // which (tgt_pid, src_pid) tuples are non-zero in general?
        let non_zero_pid_indices: Vec<_> = (0..operator.dim().0)
            .cartesian_product(0..operator.dim().1)
            .filter(|&(tgt_pid_idx, src_pid_idx)| {
                operator
                    .slice(s![tgt_pid_idx, src_pid_idx, .., .., ..])
                    .iter()
                    .any(|&value| value != 0.0)
            })
            .collect();

        // iterate over all bins, which are mapped one-to-one from the target to the source grid
        for bin in 0..self.bin_info().bins() {
            // iterate over the source grid luminosities
            for (src_lumi, src_entries) in self.lumi.iter().enumerate() {
                // create a sorted and unique vector with the `q2` for all orders
                let mut src_array_q2_grid: Vec<f64> = (0..self.orders.len())
                    .flat_map(|order| {
                        self.subgrids[[order, bin, src_lumi]]
                            .mu2_grid()
                            .iter()
                            .map(|mu2| mu2.fac)
                            .collect::<Vec<f64>>()
                    })
                    .collect();
                src_array_q2_grid.sort_by(|a, b| a.partial_cmp(b).unwrap());
                src_array_q2_grid.dedup();
                let src_array_q2_grid = src_array_q2_grid;

                let mut src_array = SparseArray3::<f64>::new(
                    src_array_q2_grid.len(),
                    if has_pdf1 { grid_axes.x_grid.len() } else { 1 },
                    if has_pdf2 { grid_axes.x_grid.len() } else { 1 },
                );

                // iterate over the source grid orders and add all of them together into
                // `src_array`, using the right powers of alphas
                for (order, powers) in self.orders.iter().enumerate() {
                    if !order_mask.is_empty() && !order_mask[order] {
                        continue;
                    }

                    let logs = if (eko_info.xir, eko_info.xif) == (1.0, 1.0) {
                        if (powers.logxir > 0) || (powers.logxif > 0) {
                            continue;
                        }

                        1.0
                    } else {
                        (eko_info.xir * eko_info.xir)
                            .ln()
                            .powi(powers.logxir.try_into().unwrap())
                            * (eko_info.xif * eko_info.xif)
                                .ln()
                                .powi(powers.logxif.try_into().unwrap())
                    };

                    let src_subgrid = &self.subgrids[[order, bin, src_lumi]];

                    // source x1/x2 grid might differ and be differently sorted than the operator
                    let x1_grid = if has_pdf1 {
                        src_subgrid
                            .x1_grid()
                            .iter()
                            .map(|x| {
                                eko_info
                                    .grid_axes
                                    .x_grid
                                    .iter()
                                    .position(|xi| approx_eq!(f64, *xi, *x, ulps = 64))
                                    .unwrap_or_else(|| unreachable!())
                            })
                            .collect()
                    } else {
                        Vec::new()
                    };
                    let x2_grid = if has_pdf2 {
                        src_subgrid
                            .x2_grid()
                            .iter()
                            .map(|x| {
                                eko_info
                                    .grid_axes
                                    .x_grid
                                    .iter()
                                    .position(|xi| approx_eq!(f64, *xi, *x, ulps = 64))
                                    .unwrap_or_else(|| unreachable!())
                            })
                            .collect()
                    } else {
                        Vec::new()
                    };

                    for ((iq2, ix1, ix2), &value) in src_subgrid.iter() {
                        let scale = src_subgrid.mu2_grid()[iq2].fac;
                        let src_iq2 = src_array_q2_grid
                            .iter()
                            .position(|&q2| q2 == scale)
                            .unwrap();
                        let als_iq2 = eko_info
                            .grid_axes
                            .mur2_grid
                            .iter()
                            .position(|&q2| {
                                approx_eq!(f64, q2, eko_info.xir * eko_info.xir * scale, ulps = 64)
                            })
                            .unwrap_or_else(|| {
                                panic!(
                                    "Couldn't find mur2: {:?} with xir: {:?} and mur2_grid: {:?}",
                                    scale, eko_info.xir, eko_info.grid_axes.mur2_grid
                                )
                            });

                        let ix1 = if has_pdf1 { x1_grid[ix1] } else { ix1 };
                        let ix2 = if has_pdf2 { x2_grid[ix2] } else { ix2 };

                        src_array[[src_iq2, ix1, ix2]] += eko_info.alphas[als_iq2]
                            .powi(powers.alphas.try_into().unwrap())
                            * logs
                            * value;
                    }
                }

                // Now we have our final grid
                let src_array = src_array;

                if src_array.is_empty() {
                    bar.inc((pids1.len() * pids2.len()) as u64);
                    continue;
                }

                // Next we need to apply the tensor
                let eko_src_q2_indices: Vec<_> = src_array_q2_grid
                    .iter()
                    .map(|&src_q2| {
                        eko_info
                            .grid_axes
                            .muf2_grid
                            .iter()
                            .position(|&q2| {
                                approx_eq!(f64, q2, eko_info.xif * eko_info.xif * src_q2, ulps = 64)
                            })
                            .unwrap_or_else(|| {
                                panic!(
                                    "Couldn't find muf2: {:?} with xif: {:?} and muf2_grid: {:?}",
                                    src_q2, eko_info.xif, eko_info.grid_axes.muf2_grid
                                )
                            })
                    })
                    .collect();
                // Iterate target lumis
                for (tgt_lumi, (tgt_pid1_idx, tgt_pid2_idx)) in (0..pids1.len())
                    .cartesian_product(0..pids2.len())
                    .enumerate()
                {
                    for (src_pid1, src_pid2, factor) in src_entries.entry().iter() {
                        // find source lumi position
                        let src_pid1_idx = if has_pdf1 {
                            eko_info
                                .grid_axes
                                .pids
                                .iter()
                                .position(|x| {
                                    // if `pid == 0` the gluon is meant
                                    if *src_pid1 == 0 {
                                        *x == 21
                                    } else {
                                        x == src_pid1
                                    }
                                })
                                .unwrap()
                        } else {
                            0
                        };
                        let src_pid2_idx = if has_pdf2 {
                            eko_info
                                .grid_axes
                                .pids
                                .iter()
                                .position(|x| {
                                    // `pid == 0` is the gluon exception, which might be 0 or 21
                                    if *src_pid2 == 0 {
                                        *x == 21
                                    } else {
                                        x == src_pid2
                                    }
                                })
                                .unwrap()
                        } else {
                            0
                        };

                        // if `op1` and `op2` below are zero there's no work to do
                        // TODO: ideally we change the for loops instead of vetoing here
                        if (has_pdf1
                            && !non_zero_pid_indices
                                .iter()
                                .any(|&tuple| tuple == (tgt_pid1_idx, src_pid1_idx)))
                            || (has_pdf2
                                && !non_zero_pid_indices
                                    .iter()
                                    .any(|&tuple| tuple == (tgt_pid2_idx, src_pid2_idx)))
                        {
                            continue;
                        }

                        // create target subgrid
                        let mut tgt_array =
                            SparseArray3::new(1, tgt_x1_grid.len(), tgt_x2_grid.len());

                        // slice the operater (which has already been reshuffled in the beginning)
                        let op1 = operator.slice(s![tgt_pid1_idx, src_pid1_idx, .., .., ..]);
                        let op2 = operator.slice(s![tgt_pid2_idx, src_pid2_idx, .., .., ..]);

                        // -- this is by far the slowest section, and has to be optimized

                        // iterate the target x position
                        for (tgt_x1_idx, tgt_x2_idx) in
                            (0..tgt_x1_grid.len()).cartesian_product(0..tgt_x2_grid.len())
                        {
                            for ((src_q2_idx, src_x1_idx, src_x2_idx), value) in
                                src_array.indexed_iter()
                            {
                                // do the linear algebra
                                let mut value = factor * value;
                                let eko_src_q2_idx = eko_src_q2_indices[src_q2_idx];

                                if has_pdf1 {
                                    value *= op1[[tgt_x1_idx, eko_src_q2_idx, src_x1_idx]];
                                }

                                // it's possible that at least one of the operators is zero - so skip, if possible
                                if value == 0.0 {
                                    continue;
                                }

                                if has_pdf2 {
                                    value *= op2[[tgt_x2_idx, eko_src_q2_idx, src_x2_idx]];
                                }

                                // it's possible that at least one of the operators is zero - so skip, if possible
                                if value == 0.0 {
                                    continue;
                                }

                                tgt_array[[0, tgt_x1_idx, tgt_x2_idx]] += value;
                            }
                        }

                        // --

                        // Now transfer the computed subgrid into the target grid
                        if !tgt_array.is_empty() {
                            let mut tgt_subgrid = mem::replace(
                                &mut result.subgrids[[0, bin, tgt_lumi]],
                                EmptySubgridV1::default().into(),
                            );

                            let mut subgrid = match tgt_subgrid {
                                SubgridEnum::EmptySubgridV1(_) => ImportOnlySubgridV2::new(
                                    tgt_array,
                                    tgt_q2_grid.clone(),
                                    tgt_x1_grid.clone(),
                                    tgt_x2_grid.clone(),
                                )
                                .into(),
                                SubgridEnum::ImportOnlySubgridV2(ref mut array) => {
                                    let array = array.array_mut();

                                    for ((_, tgt_x1_idx, tgt_x2_idx), &value) in
                                        tgt_array.indexed_iter()
                                    {
                                        array[[0, tgt_x1_idx, tgt_x2_idx]] += value;
                                    }

                                    tgt_subgrid
                                }
                                _ => unreachable!(),
                            };

                            mem::swap(&mut subgrid, &mut result.subgrids[[0, bin, tgt_lumi]]);
                        }
                    }

                    bar.inc(1);
                }
            }
        }

        bar.finish();

        result.optimize();
        FkTable::try_from(result).ok()
    }

    /// Converts this `Grid` into an [`FkTable`] using an evolution kernel operator (EKO) given as
    /// `operator`. The dimensions and properties of this operator must be described using `info`.
    /// The parameter `order_mask` can be used to include or exclude orders from this operation,
    /// and must correspond to the ordering given by [`Grid::orders`]. Orders that are not given
    /// are enabled, and in particular if `order_mask` is empty all orders are activated.
    ///
    /// # Errors
    ///
    /// Returns a [`GridError::EvolutionFailure`] if either the `operator` or its `info` is
    /// incompatible with this `Grid`.
    pub fn evolve(
        &self,
        operator: &Array5<f64>,
        info: &OperatorInfo,
        order_mask: &[bool],
    ) -> Result<FkTable, GridError> {
        let op_info_dim = (
            info.fac1.len(),
            info.pids1.len(),
            info.x1.len(),
            info.pids0.len(),
            info.x0.len(),
        );

        if operator.dim() != op_info_dim {
            return Err(GridError::EvolutionFailure(format!(
                "operator information {:?} does not match the operator's dimensions: {:?}",
                op_info_dim,
                operator.dim(),
            )));
        }

        let (subgrids, lumi) = if self.has_pdf1() && self.has_pdf2() {
            evolution::evolve_with_two(self, operator, info, order_mask)
        } else {
            evolution::evolve_with_one(self, operator, info, order_mask)
        }?;

        let mut grid = Self {
            subgrids,
            lumi,
            bin_limits: self.bin_limits.clone(),
            orders: vec![Order::new(0, 0, 0, 0)],
            subgrid_params: SubgridParams::default(),
            more_members: self.more_members.clone(),
        };

        // write additional metadata
        grid.set_key_value("lumi_id_types", &info.lumi_id_types);

        Ok(FkTable::try_from(grid).unwrap())
    }

    /// Deletes bins with the corresponding `bin_indices`. Repeated indices and indices larger or
    /// equal the bin length are ignored.
    pub fn delete_bins(&mut self, bin_indices: &[usize]) {
        let mut bin_indices: Vec<_> = bin_indices
            .iter()
            .copied()
            // ignore indices corresponding to bin that don't exist
            .filter(|&index| index < self.bin_info().bins())
            .collect();

        // sort and remove repeated indices
        bin_indices.sort_unstable();
        bin_indices.dedup();
        let bin_indices = bin_indices;

        let mut bin_ranges: Vec<Range<_>> = Vec::new();

        // convert indices into consecutive ranges
        for &bin_index in &bin_indices {
            match bin_ranges.last_mut() {
                Some(range) if range.end == bin_index => range.end += 1,
                _ => bin_ranges.push(bin_index..(bin_index + 1)),
            }
        }

        let bin_ranges = bin_ranges;
        let mut ranges = bin_ranges.as_slice();
        let old_limits = self.bin_limits.limits();

        // remove the bins from the right first, so as not to invalidate any indices
        if let Some((range, remainder)) = ranges.split_last() {
            if range.end == self.bin_info().bins() {
                self.bin_limits.delete_bins_right(range.end - range.start);
                ranges = remainder;
            }
        }

        // indices on the left aren't affected by removal of bins to their right
        if let Some((range, remainder)) = ranges.split_first() {
            if range.start == 0 {
                self.bin_limits.delete_bins_left(range.end);
                ranges = remainder;
            }
        }

        if !ranges.is_empty() {
            // if there's no remapper we need to store the bin limits in a new remapper
            if self.remapper_mut().is_none() {
                self.set_remapper(
                    BinRemapper::new(
                        old_limits.windows(2).map(|win| win[1] - win[0]).collect(),
                        old_limits.windows(2).map(|win| (win[0], win[1])).collect(),
                    )
                    .unwrap_or_else(|_| unreachable!()),
                )
                .unwrap_or_else(|_| unreachable!());
            }

            // the following should not be needed, but let's set these limits to integer values
            self.bin_limits = BinLimits::new(
                iter::successors(Some(0.0), |x| Some(x + 1.0))
                    .take(old_limits.len() - bin_indices.len())
                    .collect(),
            );
        }

        if let Some(remapper) = self.remapper_mut() {
            remapper.delete_bins(&bin_ranges);
        }

        for &bin_index in bin_indices.iter().rev() {
            self.subgrids.remove_index(Axis(1), bin_index);
        }
    }

    pub(crate) fn rewrite_lumi(&mut self, add: &[(i32, i32)], del: &[i32]) {
        self.lumi = self
            .lumi
            .iter()
            .map(|entry| {
                LumiEntry::new(
                    entry
                        .entry()
                        .iter()
                        .map(|(a, b, f)| {
                            (
                                // if `a` is to be added to another pid replace it with this pid
                                add.iter().fold(
                                    *a,
                                    |id, &(source, target)| if id == source { target } else { id },
                                ),
                                // if `b` is to be added to another pid replace it with this pid
                                add.iter().fold(
                                    *b,
                                    |id, &(source, target)| if id == source { target } else { id },
                                ),
                                // if any of the pids `a` or `b` are to b deleted set the factor to
                                // zero
                                if del.iter().any(|id| id == a || id == b) {
                                    0.0
                                } else {
                                    *f
                                },
                            )
                        })
                        .collect(),
                )
            })
            .collect();
    }

    /// Returns `true` if the first initial state needs a convolution, `false` otherwise.
    pub fn has_pdf1(&self) -> bool {
        let initial_state_1 = self.initial_state_1();

        !self
            .lumi()
            .iter()
            .all(|entry| entry.entry().iter().all(|&(a, _, _)| a == initial_state_1))
    }

    /// Returns `true` if the second initial state needs a convolution, `false` otherwise.
    pub fn has_pdf2(&self) -> bool {
        let initial_state_2 = self.initial_state_2();

        !self
            .lumi()
            .iter()
            .all(|entry| entry.entry().iter().all(|&(_, b, _)| b == initial_state_2))
    }

    /// Returns the particle identifier of the first initial state. This is usually but not always
    /// a proton, which is represented by the PDG ID `2212`.
    pub fn initial_state_1(&self) -> i32 {
        self.key_values()
            .map_or(Some("2212"), |kv| {
                kv.get("initial_state_1").map(String::as_str)
            })
            .map(str::parse)
            .unwrap()
            .unwrap()
    }

    /// Returns the particle identifier of the second initial state. This is usually but not always
    /// a proton, which is represented by the PDG ID `2212`.
    pub fn initial_state_2(&self) -> i32 {
        self.key_values()
            .map_or(Some("2212"), |kv| {
                kv.get("initial_state_2").map(String::as_str)
            })
            .map(str::parse)
            .unwrap()
            .unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::import_only_subgrid::ImportOnlySubgridV1;
    use crate::lumi_entry;

    #[test]
    fn order_cmp() {
        let mut orders = vec![
            Order::new(1, 2, 1, 0),
            Order::new(1, 2, 0, 1),
            Order::new(1, 2, 0, 0),
            Order::new(0, 3, 1, 0),
            Order::new(0, 3, 0, 1),
            Order::new(0, 3, 0, 0),
            Order::new(0, 2, 0, 0),
        ];

        orders.sort();

        assert_eq!(orders[0], Order::new(0, 2, 0, 0));
        assert_eq!(orders[1], Order::new(1, 2, 0, 0));
        assert_eq!(orders[2], Order::new(1, 2, 0, 1));
        assert_eq!(orders[3], Order::new(1, 2, 1, 0));
        assert_eq!(orders[4], Order::new(0, 3, 0, 0));
        assert_eq!(orders[5], Order::new(0, 3, 0, 1));
        assert_eq!(orders[6], Order::new(0, 3, 1, 0));
    }

    #[test]
    fn order_create_mask() {
        // Drell—Yan orders
        let orders = [
            Order::new(0, 2, 0, 0), //   LO        :          alpha^2
            Order::new(1, 2, 0, 0), //  NLO QCD    : alphas   alpha^2
            Order::new(0, 3, 0, 0), //  NLO  EW    :          alpha^3
            Order::new(2, 2, 0, 0), // NNLO QCD    : alphas^2 alpha^2
            Order::new(1, 3, 0, 0), // NNLO QCD—EW : alphas   alpha^3
            Order::new(0, 4, 0, 0), // NNLO EW     :          alpha^4
        ];

        assert_eq!(
            Order::create_mask(&orders, 0, 0),
            [false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 0, 1),
            [true, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 0, 2),
            [true, false, true, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 0, 3),
            [true, false, true, false, false, true]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 0),
            [true, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 1),
            [true, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 2),
            [true, false, true, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 3),
            [true, false, true, false, false, true]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 0),
            [true, true, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 1),
            [true, true, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 2),
            [true, true, true, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 3),
            [true, true, true, false, false, true]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 0),
            [true, true, false, true, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 1),
            [true, true, false, true, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 2),
            [true, true, true, true, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 3),
            [true, true, true, true, true, true]
        );

        // Top-pair production orders
        let orders = [
            Order::new(2, 0, 0, 0), //   LO QCD    : alphas^2
            Order::new(1, 1, 0, 0), //   LO QCD—EW : alphas   alpha
            Order::new(0, 2, 0, 0), //   LO  EW    :          alpha^2
            Order::new(3, 0, 0, 0), //  NLO QCD    : alphas^3
            Order::new(2, 1, 0, 0), //  NLO QCD—EW : alphas^2 alpha
            Order::new(1, 2, 0, 0), //  NLO QCD—EW : alphas   alpha^2
            Order::new(0, 3, 0, 0), //  NLO  EW    :          alpha^3
            Order::new(4, 0, 0, 0), // NNLO QCD    : alphas^4
            Order::new(3, 1, 0, 0), // NNLO QCD—EW : alphas^3 alpha
            Order::new(2, 2, 0, 0), // NNLO QCD—EW : alphas^2 alpha^2
            Order::new(1, 3, 0, 0), // NNLO QCD—EW : alphas   alpha^3
            Order::new(0, 4, 0, 0), // NNLO EW     :          alpha^4
        ];

        assert_eq!(
            Order::create_mask(&orders, 0, 0),
            [false, false, false, false, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 0, 1),
            [false, false, true, false, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 0, 2),
            [false, false, true, false, false, false, true, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 0, 3),
            [false, false, true, false, false, false, true, false, false, false, false, true]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 0),
            [true, false, false, false, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 1),
            [true, true, true, false, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 2),
            [true, true, true, false, false, false, true, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 3),
            [true, true, true, false, false, false, true, false, false, false, false, true]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 0),
            [true, false, false, true, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 1),
            [true, true, true, true, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 2),
            [true, true, true, true, true, true, true, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 3),
            [true, true, true, true, true, true, true, false, false, false, false, true]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 0),
            [true, false, false, true, false, false, false, true, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 1),
            [true, true, true, true, false, false, false, true, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 2),
            [true, true, true, true, true, true, true, true, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 3),
            [true, true, true, true, true, true, true, true, true, true, true, true]
        );
    }

    #[test]
    fn grid_with_subgrid_type() {
        let subgrid_type = String::from("Idontexist");
        let result = Grid::with_subgrid_type(
            vec![],
            vec![],
            vec![],
            SubgridParams::default(),
            ExtraSubgridParams::default(),
            &subgrid_type,
        );

        matches!(result, Err(GridError::UnknownSubgridType(x)) if x == subgrid_type);
    }

    #[test]
    fn grid_merge_empty_subgrids() {
        let mut grid = Grid::new(
            vec![
                lumi_entry![2, 2, 1.0; 4, 4, 1.0],
                lumi_entry![1, 1, 1.0; 3, 3, 1.0],
            ],
            vec![Order::new(0, 2, 0, 0)],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            SubgridParams::default(),
        );

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.lumi().len(), 2);
        assert_eq!(grid.orders().len(), 1);

        let other = Grid::new(
            vec![
                // differently ordered than `grid`
                lumi_entry![1, 1, 1.0; 3, 3, 1.0],
                lumi_entry![2, 2, 1.0; 4, 4, 1.0],
            ],
            vec![Order::new(1, 2, 0, 0), Order::new(1, 2, 0, 1)],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            SubgridParams::default(),
        );

        // merging with empty subgrids should not change the grid
        assert!(grid.merge(other).is_ok());

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.lumi().len(), 2);
        assert_eq!(grid.orders().len(), 1);
    }

    #[test]
    fn grid_merge_orders() {
        let mut grid = Grid::new(
            vec![
                lumi_entry![2, 2, 1.0; 4, 4, 1.0],
                lumi_entry![1, 1, 1.0; 3, 3, 1.0],
            ],
            vec![Order::new(0, 2, 0, 0)],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            SubgridParams::default(),
        );

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.lumi().len(), 2);
        assert_eq!(grid.orders().len(), 1);

        let mut other = Grid::new(
            vec![
                lumi_entry![2, 2, 1.0; 4, 4, 1.0],
                lumi_entry![1, 1, 1.0; 3, 3, 1.0],
            ],
            vec![
                Order::new(1, 2, 0, 0),
                Order::new(1, 2, 0, 1),
                Order::new(0, 2, 0, 0),
            ],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            SubgridParams::default(),
        );

        other.fill_all(
            0,
            0.1,
            &Ntuple {
                x1: 0.1,
                x2: 0.2,
                q2: 90.0_f64.powi(2),
                weight: (),
            },
            &[1.0, 2.0],
        );
        other.fill_all(
            1,
            0.1,
            &Ntuple {
                x1: 0.1,
                x2: 0.2,
                q2: 90.0_f64.powi(2),
                weight: (),
            },
            &[1.0, 2.0],
        );

        // merge with four non-empty subgrids
        assert!(grid.merge(other).is_ok());

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.lumi().len(), 2);
        assert_eq!(grid.orders().len(), 3);
    }

    #[test]
    fn grid_merge_lumi_entries() {
        let mut grid = Grid::new(
            vec![
                lumi_entry![2, 2, 1.0; 4, 4, 1.0],
                lumi_entry![1, 1, 1.0; 3, 3, 1.0],
            ],
            vec![Order::new(0, 2, 0, 0)],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            SubgridParams::default(),
        );

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.lumi().len(), 2);
        assert_eq!(grid.orders().len(), 1);

        let mut other = Grid::new(
            vec![lumi_entry![22, 22, 1.0], lumi_entry![2, 2, 1.0; 4, 4, 1.0]],
            vec![Order::new(0, 2, 0, 0)],
            vec![0.0, 0.25, 0.5, 0.75, 1.0],
            SubgridParams::default(),
        );

        // fill the photon-photon entry
        other.fill(
            0,
            0.1,
            0,
            &Ntuple {
                x1: 0.1,
                x2: 0.2,
                q2: 90.0_f64.powi(2),
                weight: 3.0,
            },
        );

        assert!(grid.merge(other).is_ok());

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.lumi().len(), 3);
        assert_eq!(grid.orders().len(), 1);
    }

    #[test]
    fn grid_merge_bins() {
        let mut grid = Grid::new(
            vec![
                lumi_entry![2, 2, 1.0; 4, 4, 1.0],
                lumi_entry![1, 1, 1.0; 3, 3, 1.0],
            ],
            vec![Order::new(0, 2, 0, 0)],
            vec![0.0, 0.25, 0.5],
            SubgridParams::default(),
        );

        assert_eq!(grid.bin_info().bins(), 2);
        assert_eq!(grid.lumi().len(), 2);
        assert_eq!(grid.orders().len(), 1);

        let mut other = Grid::new(
            vec![
                // luminosity function is differently sorted
                lumi_entry![1, 1, 1.0; 3, 3, 1.0],
                lumi_entry![2, 2, 1.0; 4, 4, 1.0],
            ],
            vec![Order::new(0, 2, 0, 0)],
            vec![0.5, 0.75, 1.0],
            SubgridParams::default(),
        );

        other.fill_all(
            0,
            0.1,
            &Ntuple {
                x1: 0.1,
                x2: 0.2,
                q2: 90.0_f64.powi(2),
                weight: (),
            },
            &[2.0, 3.0],
        );

        assert!(grid.merge(other).is_ok());

        assert_eq!(grid.bin_info().bins(), 4);
        assert_eq!(grid.lumi().len(), 2);
        assert_eq!(grid.orders().len(), 1);
    }

    // TODO: convolute_subgrid, merge_bins, subgrid, set_subgrid

    #[test]
    fn grid_key_value() {
        let mut grid = Grid::new(
            vec![lumi_entry![21, 21, 1.0]],
            vec![Order {
                alphas: 0,
                alpha: 0,
                logxir: 0,
                logxif: 0,
            }],
            vec![0.0, 1.0],
            SubgridParams::default(),
        );

        assert_eq!(
            grid.key_values().unwrap().get("initial_state_1").unwrap(),
            "2212"
        );

        grid.key_values_mut()
            .insert("initial_state_1".into(), "-2212".into());
        grid.set_key_value("initial_state_2", "-2212");

        assert_eq!(
            grid.key_values()
                .unwrap()
                .get("initial_state_1".into())
                .unwrap(),
            "-2212"
        );
        assert_eq!(
            grid.key_values()
                .unwrap()
                .get("initial_state_2".into())
                .unwrap(),
            "-2212"
        );
    }

    // TODO: properly test axes returned

    fn simple_grid() -> (Grid, GridAxes) {
        let mur2_grid = vec![20.];
        let muf2_grid = vec![20.];
        let x_grid = vec![0.1, 0.5, 1.];

        let mut subgrid_params = SubgridParams::default();
        subgrid_params.set_x_order(1);
        subgrid_params.set_x_bins(1);
        subgrid_params.set_q2_order(1);
        subgrid_params.set_q2_bins(1);

        let mut array = Array3::zeros((1, 3, 3));
        array[[0, 0, 0]] = 1.;
        let sparse_array = SparseArray3::from_ndarray(&array, 0, 3);
        let subgrid = ImportOnlySubgridV1::new(
            sparse_array,
            muf2_grid.clone(),
            x_grid.clone(),
            x_grid.clone(),
        )
        .into();

        let pids = vec![21, 1, 2];
        let mut grid = Grid::new(
            vec![lumi_entry![21, 21, 1.0], lumi_entry![1, 2, 1.0]],
            vec![Order {
                alphas: 0,
                alpha: 0,
                logxir: 0,
                logxif: 0,
            }],
            vec![0.0, 1.0],
            subgrid_params,
        );

        grid.set_subgrid(0, 0, 0, subgrid);

        (
            grid,
            GridAxes {
                x_grid,
                pids,
                mur2_grid,
                muf2_grid,
            },
        )
    }

    #[test]
    fn grid_axes() {
        let (grid, axes) = simple_grid();

        let ret_axes = grid.axes().unwrap();
        assert_eq!(ret_axes.x_grid, axes.x_grid);
        assert_eq!(ret_axes.mur2_grid, axes.mur2_grid);
        assert_eq!(ret_axes.muf2_grid, axes.muf2_grid);
        assert_eq!(ret_axes.pids, vec![]);
    }

    #[test]
    fn grid_convolute_eko() {
        let (grid, axes) = simple_grid();
        let target_x_grid = vec![1e-7, 1e-2, 1.];
        let target_pids = vec![21, 1, 2];

        let lumi_id_types = "pdg_mc_ids".to_owned();

        let eko_info = EkoInfo {
            muf2_0: 1.,
            alphas: vec![1.],
            xir: 1.,
            xif: 1.,
            target_x_grid: target_x_grid.clone(),
            target_pids: target_pids.clone(),
            grid_axes: GridAxes {
                x_grid: axes.x_grid.clone(),
                pids: axes.pids.clone(),
                mur2_grid: axes.mur2_grid.clone(),
                muf2_grid: axes.muf2_grid.clone(),
            },
            lumi_id_types,
        };
        let operator = ndarray::Array::from_shape_vec(
            (1, 3, 3, 3, 3),
            (0..4)
                .map(|_| (0..3))
                .multi_cartesian_product()
                .map(|v| if v[0] == v[2] && v[1] == v[3] { 1. } else { 0. })
                .collect(),
        )
        .unwrap();
        let fk = grid.convolute_eko(operator, eko_info, &[]).unwrap();

        assert_eq!(fk.bins(), 1);
    }

    fn execute_evolve_test(
        setname: &str,
        grid: &str,
        metadata: &str,
        alphas: &str,
        operator: &str,
        ref_results: &[f64],
        ref_evolved_results: &[f64],
    ) {
        use float_cmp::assert_approx_eq;
        use lhapdf::Pdf;
        use ndarray::Array1;
        use std::fs::File;

        #[derive(Deserialize)]
        struct Metadata {
            #[serde(rename = "Q2grid")]
            q2_grid: Vec<f64>,
            #[allow(dead_code)]
            eko_version: String,
            inputgrid: Vec<f64>,
            inputpids: Vec<i32>,
            #[allow(dead_code)]
            interpolation_is_log: bool,
            #[allow(dead_code)]
            interpolation_polynomial_degree: usize,
            #[allow(dead_code)]
            interpolation_xgrid: Vec<f64>,
            q2_ref: f64,
            targetgrid: Vec<f64>,
            targetpids: Vec<i32>,
        }

        let grid = Grid::read(File::open(grid).unwrap()).unwrap();

        let lhapdf = Pdf::with_setname_and_nmem(setname).unwrap();
        let mut pdf = |id, x, q2| lhapdf.xfx_q2(id, x, q2);
        let mut als = |q2| lhapdf.alphas_q2(q2);
        let mut cache = LumiCache::with_one(2212, &mut pdf, &mut als);
        let results = grid.convolute(&mut cache, &[], &[], &[], &[(1.0, 1.0)]);

        assert_eq!(results, ref_results);

        let metadata: Metadata = serde_yaml::from_reader(File::open(metadata).unwrap()).unwrap();
        let alphas: Array1<f64> = ndarray_npy::read_npy(alphas).unwrap();
        let alphas = alphas.to_vec();
        let operator: Array5<f64> = ndarray_npy::read_npy(operator).unwrap();

        let info = OperatorInfo {
            fac1: metadata.q2_grid.clone(),
            pids0: metadata.inputpids,
            x0: metadata.inputgrid,
            pids1: metadata.targetpids,
            x1: metadata.targetgrid,
            fac0: metadata.q2_ref,
            ren1: metadata.q2_grid, // TODO: check whether this is true in the general case
            alphas,
            xir: 1.0,
            xif: 1.0,
            lumi_id_types: "pdg_mc_ids".to_string(),
        };

        //let mut grid_axes = grid.axes().unwrap();
        //grid_axes.pids = metadata.targetpids.clone(); // TODO: which one them?
        //let eko_info = EkoInfo {
        //    muf2_0: metadata.q2_ref,
        //    alphas,
        //    xir: 1.0,
        //    xif: 1.0,
        //    target_x_grid: metadata.inputgrid,
        //    target_pids: metadata.targetpids,
        //    grid_axes,
        //    lumi_id_types: "pdg_mc_ids".to_string(),
        //};

        //let fk_table = grid.convolute_eko(operator, eko_info, &[]).unwrap();
        let fk_table = grid.evolve(&operator, &info, &[]).unwrap();
        let evolved_results = fk_table.convolute(&mut cache, &[], &[]);

        assert_eq!(evolved_results.len(), ref_evolved_results.len());

        for (bin, diff) in results
            .iter()
            .zip(evolved_results.iter())
            .map(|(a, b)| (b / a - 1.0) * 1000.0)
            .enumerate()
            .filter(|&(_, diff)| diff > 1.0)
        {
            println!("WARNING: bin {} diff {}", bin, diff);
        }

        for (&result, &ref_result) in evolved_results.iter().zip(ref_evolved_results.iter()) {
            assert_approx_eq!(f64, result, ref_result, ulps = 16);
        }
    }

    #[test]
    #[ignore]
    fn evolve_hera_cc_318gev_ep_sigmared() {
        execute_evolve_test(
            "NNPDF40_nlo_as_01180",
            "../HERA_CC_318GEV_EP_SIGMARED.pineappl.lz4",
            "../HERA_CC_318GEV_EP_SIGMARED/metadata.yaml",
            "../HERA_CC_318GEV_EP_SIGMARED/alphas.npy",
            "../HERA_CC_318GEV_EP_SIGMARED/operators.npy",
            &[
                1.158872296528366,
                1.1204108571716151,
                0.8875902564141837,
                0.6064124901549366,
                0.46730446763293626,
                0.9175036665313503,
                0.9625683080629993,
                0.832097616454224,
                0.5840101875241862,
                0.45055068848064506,
                0.6688041946734723,
                0.6889352106823404,
                0.5297495101909169,
                0.41609405897210827,
                0.23143317308882982,
                0.5649694171604525,
                0.4792399077439775,
                0.3861828122239765,
                0.21793994357039714,
                0.08432403516120096,
                0.467724373428524,
                0.43253481860290705,
                0.3588386557913573,
                0.20658438492949507,
                0.36145545861402606,
                0.350280931434649,
                0.3095655884296923,
                0.18715074960942407,
                0.0737622486172105,
                0.2292088799939892,
                0.2277411334118036,
                0.15507169325765838,
                0.06400693681274111,
                0.1403126590329724,
                0.11639162212711233,
                0.052756148828243,
                0.05464232942385792,
                0.033578480958376296,
                0.01095350422009362,
            ],
            &[
                1.1586791793006517,
                1.1202273184348848,
                0.8874599897872428,
                0.6063712625073608,
                0.4672874245885059,
                0.9173403727229149,
                0.962399631849818,
                0.8319677648956049,
                0.5839675146069822,
                0.4505325132086283,
                0.6686768094263615,
                0.6888210100292264,
                0.5297074425880604,
                0.41607535867486134,
                0.23142390114495495,
                0.5648727887895276,
                0.47919973616456246,
                0.38616444708721975,
                0.21793107752932925,
                0.08431941687497238,
                0.46764290121503427,
                0.43249701889959064,
                0.35882091245690373,
                0.20657579928724892,
                0.36139047088472026,
                0.35024730127068276,
                0.3095490328079981,
                0.1871428094003754,
                0.07375775279995289,
                0.22918160031673335,
                0.22772688522902013,
                0.15506498963894486,
                0.0640029182481115,
                0.14030150245012324,
                0.11638652740576161,
                0.052752620249130584,
                0.05463985560251649,
                0.03357575618674528,
                0.010951799924697785,
            ],
        );
    }

    #[test]
    #[ignore]
    fn evolve_atlaswzrap36pb() {
        execute_evolve_test(
            "NNPDF40_nlo_as_01180",
            "../ATLASWZRAP36PB-ATLAS-arXiv:1109.5141-Z0_eta34.pineappl.lz4",
            "../ATLASWZRAP36PB-ATLAS-arXiv:1109.5141-Z0_eta34/metadata.yaml",
            "../ATLASWZRAP36PB-ATLAS-arXiv:1109.5141-Z0_eta34/alphas.npy",
            "../ATLASWZRAP36PB-ATLAS-arXiv:1109.5141-Z0_eta34/operators.npy",
            &[
                134632.4480167966,
                133472.53856965082,
                130776.9197998554,
                126991.17650299768,
                121168.15135839234,
                112761.75563082621,
                98484.12455342771,
                57843.82009866679,
            ],
            &[
                134500.65897999398,
                133344.24845948347,
                130663.70226447149,
                126890.6156547337,
                121090.7595443343,
                112698.92115188896,
                98421.19065191328,
                57813.500329070215,
            ],
        );
    }

    #[test]
    #[ignore]
    fn evolve_lhcb_wp_8tev() {
        execute_evolve_test(
            "NNPDF40_nlo_as_01180",
            "../LHCB_WP_8TEV.pineappl.lz4",
            "../LHCB_WP_8TEV/metadata.yaml",
            "../LHCB_WP_8TEV/alphas.npy",
            "../LHCB_WP_8TEV/operators.npy",
            &[
                877.9498750860583,
                823.5123400865052,
                735.605093586326,
                616.3465722662226,
                478.63703336207277,
                341.06729874517384,
                174.3688634669724,
                48.27440593682665,
            ],
            &[
                877.710496991437,
                823.2923843876903,
                735.4119248950592,
                616.18687672059,
                478.5155078388758,
                340.9838587318582,
                174.33023801119046,
                48.2671526729747,
            ],
        );
    }

    #[test]
    #[ignore]
    fn evolve_atlas_wm_jet_8tev_pt() {
        execute_evolve_test(
            "NNPDF40_nlo_as_01180",
            "../ATLAS_WM_JET_8TEV_PT-atlas-atlas-wjets-arxiv-1711.03296-xsec003.pineappl.lz4",
            "../ATLAS_WM_JET_8TEV_PT-atlas-atlas-wjets-arxiv-1711.03296-xsec003/metadata.yaml",
            "../ATLAS_WM_JET_8TEV_PT-atlas-atlas-wjets-arxiv-1711.03296-xsec003/alphas.npy",
            "../ATLAS_WM_JET_8TEV_PT-atlas-atlas-wjets-arxiv-1711.03296-xsec003/operators.npy",
            &[
                1469.0819829706113,
                4364.509970805006,
                2072.647269805865,
                742.6160478667135,
                312.195449854513,
                145.782031788469,
                72.64484245290579,
                38.88392569051013,
                17.093513927247354,
                6.157887576517187,
                2.4621639421455903,
                1.0817834362838417,
                0.5084171526510098,
                0.2520459057372801,
                0.10211930488819178,
                0.021141855492915994,
            ],
            &[
                1467.6009449768221,
                4359.87917180896,
                2070.4278159799974,
                741.9298754488171,
                311.9073865166957,
                145.65671641953438,
                72.58340237308579,
                38.85208071336316,
                17.080194246318936,
                6.153298093496777,
                2.4604627649049604,
                1.0810950528772425,
                0.5081137620517796,
                0.25190465608989626,
                0.10206534388970377,
                0.02113211970400249,
            ],
        );
    }

    #[test]
    #[ignore]
    fn evolve_nnpdf_pos_anti_up_40() {
        execute_evolve_test(
            "NNPDF40_nnlo_as_01180",
            "../NNPDF_POS_ANTI_UP_40.pineappl.lz4",
            "../NNPDF_POS_ANTI_UP_40/metadata.yaml",
            "../HERA_CC_318GEV_EP_SIGMARED/alphas.npy",
            "../NNPDF_POS_ANTI_UP_40/operators.npy",
            &[
                1.24053558491532,
                1.0556043292537103,
                0.8982894541272358,
                0.764653254688176,
                0.6497645161935843,
                0.5455683653190757,
                0.44285618029811136,
                0.34236858231394923,
                0.23851111837499242,
                0.10826740305077837,
                0.053312425688123236,
                0.022300403914252444,
                0.007689357135047722,
                0.0030482582072498296,
                0.001778467203699489,
                0.0016533253586389347,
                0.001889830982843056,
                0.0014804650577204803,
                0.0006649270383861285,
                0.00014006028310118777,
            ],
            &[
                1.2207842700252995,
                1.045048066777534,
                0.8929558020715721,
                0.7616886315970761,
                0.6473998210787656,
                0.5429673354151222,
                0.4397140650493807,
                0.3388398929365999,
                0.2360361559509584,
                0.10965325070809523,
                0.0565065025445373,
                0.025950818992108653,
                0.010936975130524709,
                0.005500955939196705,
                0.003402371551444896,
                0.0025976596425020862,
                0.002333911364590294,
                0.001628487038508355,
                0.0006948686757361142,
                0.00014172972322202258,
            ],
        );
    }

    #[test]
    #[ignore]
    fn evolve_chorus_cc_nb_pb_sigmared() {
        execute_evolve_test(
            "NNPDF40_nlo_as_01180",
            "../CHORUS_CC_NB_PB_SIGMARED.pineappl.lz4",
            "../CHORUS_CC_NB_PB_SIGMARED/metadata.yaml",
            "../CHORUS_CC_NB_PB_SIGMARED/alphas.npy",
            "../CHORUS_CC_NB_PB_SIGMARED/operators.npy",
            &[
                0.17508418054652689,
                0.2695273897938339,
                0.23765585105922482,
                0.25169985216985463,
                0.3441821600653684,
                0.3036860496902356,
                0.23681995980798465,
                0.42195630412804347,
                0.34270475962218155,
                0.40787940477815376,
                0.3132748894380318,
                0.2577807128549423,
                0.4300117620481792,
                0.43819261599661424,
                0.40225026401919767,
                0.5845553253061538,
                0.3593543733351632,
                0.27124972491279786,
                0.5536384436876284,
                0.42654101852129855,
                0.6640122148055845,
                0.5768259324369139,
                0.586191593826013,
                0.39093894579856103,
                0.7031461953731108,
                0.543154270152632,
                0.48743175575755404,
                0.48144934089502267,
                0.5929176631119987,
                0.41005526792552716,
                0.8318257880948878,
                0.7230182190621122,
                0.8515005407535482,
                0.7831785444888832,
                0.5261111093752395,
                0.614809628004991,
                0.7163592480665573,
                0.7596144675276578,
                0.5265857916523599,
                0.9945270628537345,
                0.6916650300420164,
                0.7773664402348859,
                0.547006181338868,
                0.9333286158343123,
                0.6936364669791876,
                0.9719407903384734,
                0.656861275827149,
                0.5439580179066053,
                1.1348808424845245,
                0.7975487381974588,
                0.9350502432439056,
                0.5558717571090206,
                1.139608410113371,
                0.9749761310145566,
                0.9473294807096161,
                1.2924478860309754,
                0.6763055089242359,
                1.1487852454645615,
                0.7380675797619043,
                0.5418437412727434,
                0.7643540773865097,
                1.1836666045484023,
                0.8825207229758716,
                0.8378868368156698,
                0.9458420189262449,
                1.0717847698282192,
                1.1317647499274746,
                0.679484167316982,
                1.3505945482579,
                1.1560798218905395,
                1.015468845847213,
                0.9373101587601604,
                1.1614862224122455,
                0.738431687759785,
                0.8452181806489453,
                1.3857792588478974,
                1.036177724313531,
                0.7583172253276966,
                0.9049773173534823,
                1.4478265180633512,
                1.2966739403264655,
                1.5135613517314832,
                1.0337045110443905,
                1.1252182889654128,
                1.1888702235492254,
                0.7092876883785939,
                0.8290400203545237,
                1.350013991301613,
                0.967448732398068,
                1.4878146949120024,
                1.2530927139381476,
                1.530228415350425,
                0.7011085242557515,
                1.3445830308171336,
                0.8626417692613555,
                1.1915801548973077,
                1.0086814810585862,
                0.8920774478992273,
                0.9847734708434532,
                1.1776858163648567,
                1.2729869484971787,
                1.0972889655410152,
                0.988412982615934,
                1.1325311352467975,
                1.2177953834502937,
                1.460937502136669,
                1.6804277138399746,
                1.576247446853074,
                1.3279307482373712,
                0.8737981370721432,
                0.6168649653436493,
                0.7870608360027186,
                0.9257553976220181,
                1.4511510723212608,
                1.0433665740952258,
                1.4109408679024265,
                0.9143982683572018,
                1.6300600868911526,
                1.054494230406241,
                1.1887441685274223,
                0.7482989409408192,
                1.0685501763604146,
                1.097104591267553,
                1.4698679072389773,
                1.218960950634183,
                1.2410987909994937,
                0.7129852945507955,
                1.5884922274685551,
                1.2525099961522117,
                0.749992569612053,
                1.1168068534203335,
                0.8761441221105105,
                0.8218718329015258,
                1.3415655577236016,
                0.9216633349011075,
                0.7978099753009874,
                1.449914680736598,
                1.4768574326730393,
                1.0643666911240752,
                0.8986610109026638,
                0.9931605057645142,
                0.611909644587469,
                0.9587418333694464,
                1.4655537300020285,
                1.2153816186491477,
                0.9511722729956501,
                0.6470100275601826,
                1.2832897617135286,
                0.7469847754293129,
                1.1175530181004367,
                1.322063983692295,
                1.0814404270774267,
                1.1238631252990405,
                0.7464032635717058,
                1.4975034508853429,
                1.6001798413799317,
                1.2806206728943808,
                0.7726802678146308,
                0.7593136740492077,
                0.5018542523697478,
                0.6285650525184223,
                0.9176468626653805,
                0.8636077728845849,
                1.0888724586268557,
                0.6758126769417545,
                1.3509038027384432,
                0.9297974499651537,
                1.4616942865839717,
                1.5099712863567014,
                1.3172747295603051,
                0.9897254333927065,
                0.8602656242958056,
                1.0036866991746578,
                1.210801321026721,
                0.7950043225996786,
                0.7306772036442314,
                1.1553078694407823,
                1.3082308332888373,
                1.0700709437093272,
                1.1508593937386442,
                0.9307773207273198,
                0.7871843081182327,
                0.6686964066774572,
                1.1313991064024929,
                0.6226175413667078,
                0.48911763992073976,
                0.5626797887982953,
                1.301874007142872,
                0.7902143928142813,
                1.528717300988022,
                1.3573117482938861,
                0.9356426747990623,
                0.5588958418980426,
                0.9395033314179385,
                0.7666937986475366,
                0.5110051440927045,
                0.9057735106618026,
                1.1130045579045142,
                0.5887318238011767,
                0.9866413295914536,
                0.725836209651192,
                1.1855835835560118,
                1.2970432635941802,
                1.0381924269928706,
                1.4546096807856124,
                1.1319609699293744,
                0.9155145612483051,
                0.8495706756372534,
                1.3692865321955805,
                0.5127815290378634,
                0.627399176593509,
                0.3789663321519362,
                0.6852968908549385,
                1.3185532293864661,
                1.056993701095617,
                0.665324593722263,
                0.46777024528292577,
                1.205346037778319,
                0.7924999144790976,
                0.8873298512070742,
                0.709459280784969,
                0.8276159355969919,
                1.3618625223069039,
                0.8101756783666472,
                1.1379919536339014,
                0.6297152151085044,
                0.37781204993042583,
                0.35560079573069786,
                0.7719036007569298,
                0.407314923597201,
                0.937990215357849,
                0.9417295550225834,
                0.40841083934057804,
                1.0657294306252212,
                1.1308969315823405,
                0.6048407945325175,
                1.1168346609759434,
                0.9606822425660947,
                0.9824349108378805,
                0.48960123994054394,
                0.5418283022915388,
                0.8412707228195404,
                0.7117366128115412,
                1.2792981029926656,
                1.200431062908973,
                0.7901172745397984,
                0.8981280632102917,
                0.5452045610826048,
                0.8663590225461809,
                0.6929932191811116,
                0.2588265140046349,
                0.7039942227014601,
                0.2742217803379954,
                0.504459105427192,
                0.5199522596034711,
                0.62597953911183,
                1.0435944866226887,
                0.6577885116272207,
                0.5650099068699643,
                0.4490846780324496,
                1.3430785271464998,
                1.142262454236171,
                0.6346986016397684,
                0.2865953219020086,
                0.19530973756051623,
                1.1447112538933928,
                1.367662564206789,
                0.38102312219626017,
                0.8292656866455711,
                0.777267398358757,
                0.28600508105457895,
                0.9761630463697777,
                0.622849462267499,
                0.2522951919540529,
                0.9468827960067999,
                0.43794101157049925,
                0.7010075046093125,
                0.1189023400415132,
                1.1959293031659233,
                0.849477208474895,
                1.0934360288790417,
                0.9773565293132962,
                0.48927275636415846,
                0.8313096442494344,
                0.3744253232382924,
                0.5268564446869717,
                0.33802865605582194,
                0.5249430384660696,
                0.6243972273466691,
                0.39888029626494703,
                0.19384371833083308,
                1.1450998582660805,
                0.08919399238257084,
                0.7866933793052904,
                0.6739693625567023,
                0.39447603928009467,
                1.032870458279164,
                0.6517181805823248,
                0.5462278142839979,
                0.43468982370668574,
                0.8807666292883145,
                0.535583589662179,
                0.30725075320625195,
                0.24400838454393775,
                0.7217194179540221,
                0.49943753457265216,
                0.9879848350716192,
                0.6397539331394319,
                0.19004703184019897,
                0.8430736944274474,
                1.1646912559844074,
                0.9500374153431089,
                0.4431011097650916,
                0.1358517225098266,
                0.3816983503496567,
                0.2734027662570491,
                0.7816697296690134,
                0.17264102023330116,
                0.6397732061925886,
                0.9727134547847681,
                0.4884817541453219,
                0.5153534913286036,
                0.6884565112874529,
                0.2821889370605011,
                0.18171114275940467,
                1.0238646640222304,
                0.5312036757045824,
                0.13708018059931698,
                0.10908203203508102,
                1.187844330026233,
                0.8234678520136302,
                0.867250037405063,
                0.5280545412466593,
                0.29298765566069757,
                0.2331119684994854,
                0.8207829720762455,
                0.3714881071395831,
                0.32529299461879857,
                0.5299622465209836,
                0.6218875130630758,
                0.27761677944132257,
                0.7344117376273341,
                0.4952243643300778,
                0.3926592728605549,
                0.8535099583622792,
                0.6442397533753098,
                0.4180782333604763,
                0.2409433056320422,
                0.7824022300069466,
                0.3027927077223924,
                0.6550368059503859,
                0.38444900185065894,
                1.1482329600700933,
                0.9519789512388735,
                0.060580882253723545,
                0.3816923021118775,
                0.19257650866758338,
                0.6438276289212486,
                0.11148017924313727,
                0.9684765768161862,
                0.08028721882900387,
                0.7842496487219783,
                0.09193695829466855,
                0.4482665888923368,
                0.2721624183103696,
                0.8561621559270862,
                0.28192108332583904,
                0.12812530881831327,
                0.1022029969166139,
                0.4869392994935116,
                0.502027655826873,
                0.1859802642666868,
                0.8121893548777561,
                0.3688192118317807,
                0.3154642616461757,
                0.6756161354062145,
                0.2786542931017233,
                0.17239019126750388,
                0.38764494496793006,
                1.0091693564839619,
                0.5082684416248578,
                0.7784557875594005,
                0.3034748623167903,
                0.5190743285546534,
                0.640507464448423,
                0.37666056704261075,
                0.22098632205470153,
                0.12755794902012385,
                0.48984474758340746,
                0.16601973443100057,
                0.5339786133263489,
                0.6185863797916844,
                0.27139680627083085,
                0.6361388553735544,
                0.20029926006564597,
                0.40163959663315296,
                0.232581415954885,
                0.7857430698310202,
                0.04064886014154843,
                0.6461618251759489,
                0.12143335100382315,
                0.9537647916055934,
                0.3810457264475656,
                0.11662315870272344,
                0.2707718275731995,
                0.19113966967629645,
                0.8048973447935833,
                0.07456956855621856,
                0.9609457914314185,
                0.4524011968793949,
                0.6654620677686393,
                0.27564461841532323,
                0.1654134978047173,
                0.09489139248107531,
                0.05538962624230017,
                0.7748124999748075,
                0.628763751660673,
                0.8385779629812827,
                0.26552252771486734,
                0.4846269413293837,
                0.3652365816241996,
                0.48863177584806516,
                0.18269227749188724,
                0.3041157980083019,
                0.10846041314090618,
                0.5362714057914194,
                0.6154381782036658,
                0.2664502831484185,
                0.38151419595140723,
                0.08726360691339036,
                0.6474703306111875,
                0.3036088813982311,
                0.36751657132302845,
                0.6294436265357662,
                0.1999183403172808,
                0.3891576520469917,
                0.22606778328384097,
                0.5097187036589822,
                0.12855778273226898,
                0.20936169478422986,
                0.12168225162688405,
                0.4838776323831572,
                0.16095027922066846,
                0.1897001511733311,
                0.07049009502330647,
                0.6570565217492003,
                0.3796506374040504,
                0.26867643894376625,
                0.11555836381974315,
                0.11193452442863465,
                0.6124842654624286,
                0.03778605693724704,
                0.7929431374080002,
                0.482283545810485,
                0.4782938175685489,
                0.17995502424142656,
                0.27181577867628687,
                0.15758542000860468,
                0.0570704314503577,
                0.6237387262108742,
                0.37924948437600586,
                0.08822802545228833,
                0.05196026936795638,
                0.7682964482449042,
                0.6105628167239305,
                0.502247236341904,
                0.1282941652548714,
                0.20080098532881438,
                0.11723141671996572,
                0.3611032638350682,
                0.2928802761150611,
                0.10587639125688059,
                0.2605471585193985,
                0.47886783088384827,
                0.15692300068049003,
                0.37499375353688147,
                0.08378971160423875,
                0.19897434129308406,
                0.218483801202603,
                0.3029957359976254,
                0.3582195132912819,
                0.3883811178827541,
                0.11628629425078416,
                0.19238242893327612,
                0.49003417113436026,
                0.4709207945548197,
                0.06675558654324969,
                0.025843406464136125,
                0.6474140875944007,
                0.2729289452234853,
                0.06531514680112217,
                0.4999824947905057,
                0.19450942266987922,
                0.08409306192210333,
                0.0501535763684259,
                0.48090644016409956,
                0.03652336440852907,
                0.17813366010665932,
                0.36400284358607593,
                0.28423070003244144,
                0.10402113897321216,
                0.27180502351511204,
                0.15049751807593523,
                0.05606619217111426,
                0.37370508841136557,
                0.08174997490641407,
                0.6201553441862678,
                0.3638735923280538,
                0.13189451721353904,
                0.35224849290304794,
                0.07012011109529989,
                0.11272110606715918,
                0.15236037922378579,
                0.25583632878103524,
                0.20370906009562176,
                0.21102859757421724,
                0.08042145551388168,
                0.11450808703206192,
                0.26940817822191465,
                0.06446058846330588,
                0.1890298072014,
                0.03516964099665347,
                0.3599724405816431,
                0.27750749150097737,
                0.06254585379774286,
                0.02463075395239731,
                0.45778222082838066,
                0.03490060843046899,
                0.3687153235755149,
                0.2676057440791352,
                0.1448204818805799,
                0.054476341129201536,
                0.0472651730338768,
                0.34624971790901565,
                0.4891481722249711,
                0.18450967020066855,
                0.17425792234569298,
                0.10133976826252562,
                0.25085886777306265,
                0.07848085739266623,
                0.06964059808513534,
                0.20529122630381336,
                0.12990903317526245,
                0.10766326009632754,
                0.14754395256487632,
                0.1862446006894729,
                0.03444040744713824,
                0.05948520227285855,
                0.023683713683422248,
                0.26406670855428793,
                0.14031929958532208,
                0.11243166358818434,
                0.06323600590529867,
                0.07507376362968197,
                0.24669876664761947,
                0.03307057280674311,
                0.2670564213867311,
                0.1710786276382489,
                0.05254919080349757,
                0.20059387111543683,
                0.044513590783645134,
                0.336679961089545,
                0.1038077051293141,
                0.0984517782779149,
                0.14371966777063166,
                0.0751842143887107,
                0.01559259099045098,
                0.03674273528485569,
                0.06885967441174741,
                0.057104823205109886,
                0.11063730228417895,
                0.03346893151645903,
                0.16837710274585108,
                0.02257444055398897,
                0.13343891534204733,
                0.06178162633117294,
                0.1007073467373077,
                0.042471442225472444,
                0.14055238322835578,
                0.031280948168247946,
                0.09609437285090261,
                0.050536303503963736,
                0.0726264120189156,
                0.1931878072042266,
                0.014980558154686598,
                0.036290411970244844,
                0.05355538573472793,
                0.04086359640260389,
                0.060520864856410284,
                0.03238209739944476,
                0.09410569269679359,
                0.021457164768928578,
                0.0299296085548296,
                0.07054524627613962,
                0.04893630055839461,
                0.09591344047647017,
                0.019902630830190192,
                0.014331642222965658,
                0.028851213767815274,
                0.03147715582928427,
                0.020593742773379987,
                0.0476119046649274,
                0.038431402299134304,
                0.01381143270758339,
                0.019893155135724604,
            ],
            &[
                0.9846578039527537,
                1.1428160730074823,
                1.0157400591515289,
                0.9373047460961609,
                1.2103085838284144,
                1.07380338424002,
                0.7684689940081988,
                1.279756177030059,
                1.0152905000655528,
                1.1642827825647268,
                0.9054185397460988,
                0.7151559537863046,
                1.2343265971792428,
                1.0889723034278065,
                1.081155071172186,
                1.4074796700349008,
                0.8576943116723308,
                0.6600110308254856,
                1.272255932823374,
                0.9872479868239301,
                1.48432129739738,
                1.292121409884949,
                1.1839710135536405,
                0.7993933215340576,
                1.3624497500745392,
                1.0573263376094482,
                0.9340313158568363,
                0.9789908470120022,
                1.1496188075755744,
                0.7371970725332893,
                1.5406978615602802,
                1.3409488760478188,
                1.4355946161127415,
                1.2822835235784886,
                0.8676653081121799,
                0.9958318374399315,
                1.1410369005167196,
                1.3182029010165506,
                0.8742080239576873,
                1.5889255647986937,
                1.1285775441215042,
                1.2050474209514361,
                0.8006725383236973,
                1.4023158074558837,
                1.0441470853256825,
                1.3624965351669247,
                0.9256393810106471,
                0.7758755383370003,
                1.5551757094613372,
                1.0746134793831859,
                1.2332230932981874,
                0.7397398556906186,
                1.5862474289981363,
                1.3582616089812056,
                1.2521507744659113,
                1.6692556717199287,
                0.8570755114661468,
                1.4337975388185895,
                0.9366317385884309,
                0.6892153354401825,
                0.9823830541760954,
                1.4694780271336392,
                1.0972377744467945,
                1.0024251445418058,
                1.1724516252065078,
                1.3092148177448175,
                1.3130180505985327,
                0.7944817393203242,
                1.6197648163527454,
                1.3871576038465498,
                1.1648716107615382,
                1.0933977811302296,
                1.3105855965735682,
                0.8357125958674643,
                0.9308538130211269,
                1.5240726616291256,
                1.1412385927939876,
                0.8430627416646357,
                0.9865171167043952,
                1.5542168162471084,
                1.381798359764537,
                1.6458954233797203,
                1.0882356398263828,
                1.2024386157072753,
                1.2404177637653306,
                0.7465751526915569,
                0.865106603000001,
                1.4196734539347613,
                1.0180342736330041,
                1.5706646267023336,
                1.3231267552483184,
                1.5678119587857564,
                0.7201741664521962,
                1.3633041525327652,
                0.877309242102593,
                1.2240152798288562,
                1.020741963826803,
                0.9095039353275298,
                0.993604922759415,
                1.1869976325063978,
                1.275979707723701,
                1.1069931861796032,
                0.9980747423260864,
                1.1347014742230916,
                1.224588890765672,
                1.4591906825445102,
                1.6842968910850602,
                1.580170115049403,
                1.331250106457184,
                0.8743784539888663,
                0.6173438080570497,
                0.7867406929390685,
                0.9256099944955801,
                1.4511603892005631,
                1.0433756041894395,
                1.4110659168537616,
                0.9144772304267657,
                1.63025852876531,
                1.0546189730754731,
                1.1888734107091454,
                0.7483261051844592,
                1.0686505480953365,
                1.0970869647500705,
                1.4698861489035782,
                1.2189723306112468,
                1.2411587048698787,
                0.7130206006744954,
                1.5885268972221829,
                1.2525099247497116,
                0.7499926314230843,
                1.1168067383099014,
                0.876141993434363,
                0.8218642824510636,
                1.3415625634379356,
                0.9216610372261911,
                0.7978053893083532,
                1.449899126372724,
                1.476847705332237,
                1.0643599107324762,
                0.8986551842761088,
                0.9931546498056035,
                0.6119077542415317,
                0.958730452956305,
                1.4655467087600806,
                1.2153748480769544,
                0.9511690438337452,
                0.6470062788930065,
                1.2832716066971843,
                0.7469753236357717,
                1.1175355594302643,
                1.3220588863624538,
                1.0814358527904078,
                1.1238587216600255,
                0.7464007136663318,
                1.4974883808022084,
                1.6001740492661205,
                1.280607710647327,
                0.7726731654727896,
                0.7593100064751936,
                0.5018511539303819,
                0.628556346863019,
                0.9176382441374444,
                0.8635928996889607,
                1.088862050999555,
                0.6758071866581457,
                1.3508994977747555,
                0.9297936838034961,
                1.4616832908047759,
                1.5099417817784118,
                1.317248522833804,
                0.9897145568742565,
                0.8602587005242931,
                1.0036652631960274,
                1.2107885670764857,
                0.7949929843297018,
                0.7306809216448841,
                1.155281584464726,
                1.3082212623693805,
                1.0700622574801901,
                1.1508563395440718,
                0.9307748891420934,
                0.7871674899326594,
                0.6686885422849159,
                1.1313931149408478,
                0.62261286017232,
                0.48910638298216746,
                0.5626842094670441,
                1.3018566581004074,
                0.7902048244877213,
                1.5286966610628414,
                1.357305719686446,
                0.9356373989108703,
                0.5588945325201345,
                0.9394911004201509,
                0.7666886089854141,
                0.5109959784084821,
                0.9057518975700104,
                1.1129879282777555,
                0.588724232779656,
                0.9866267010492198,
                0.7258346814276995,
                1.1855512905528838,
                1.2970296895337297,
                1.0381642497510992,
                1.4545924111511648,
                1.1319581558817173,
                0.9155125549602462,
                0.8495586859659635,
                1.36924665875997,
                0.5127759686453649,
                0.6273850877628259,
                0.37895345058066043,
                0.6852864849654823,
                1.3185309443011708,
                1.0569801115861845,
                0.6653130207412206,
                0.46775802407243866,
                1.2053292541014717,
                0.7924848019354067,
                0.8873330176345925,
                0.709463728750861,
                0.8275941505981728,
                1.3618545740162549,
                0.8101625400764274,
                1.137984443388146,
                0.6297091701401701,
                0.37780105841922934,
                0.3555889579323264,
                0.7718967693611626,
                0.4073147972841343,
                0.9379627144439735,
                0.9417215253305974,
                0.40841887927808795,
                1.0656964423715485,
                1.1308745373618423,
                0.6048295878380807,
                1.116831432639867,
                0.9606624861471456,
                0.9824167577119173,
                0.4895856564395406,
                0.5418360230853796,
                0.8412559814492886,
                0.7117348654791876,
                1.2792802912893553,
                1.2004116653663646,
                0.790099763177358,
                0.8981249071897095,
                0.5452031847098657,
                0.8663640824453414,
                0.6929997931282208,
                0.2588219479216775,
                0.7039791205320579,
                0.2742098731754484,
                0.5044456800268781,
                0.5199457102788779,
                0.6259626435336597,
                1.0435773559241601,
                0.6577737778175918,
                0.5649906498610809,
                0.4490682528462317,
                1.3430483894075582,
                1.1422534120289813,
                0.6346914252671366,
                0.2865947310563475,
                0.1953062739024531,
                1.1446845911181278,
                1.3676514520743444,
                0.3810067177635386,
                0.8292443239788192,
                0.7772584922116912,
                0.28601703850852,
                0.9761373392785603,
                0.6228324356458697,
                0.25228338634547637,
                0.9468728734844488,
                0.43793381623448996,
                0.7010041274302387,
                0.11890725155914843,
                1.1959078247835222,
                0.8494843519713849,
                1.093429820886188,
                0.9773360735547205,
                0.4892548842151652,
                0.8312922239345566,
                0.37441121542620937,
                0.5268667040132265,
                0.3380138730026699,
                0.5249354739703403,
                0.624378470809524,
                0.39887809500471666,
                0.1938428421936126,
                1.1450895383784212,
                0.0891971965278184,
                0.7866736607280926,
                0.6739781420218406,
                0.39448788345605695,
                1.0328511371871365,
                0.6517015758940885,
                0.546206443863707,
                0.434671546446116,
                0.8807609660496004,
                0.5355802222926985,
                0.3072459351055423,
                0.244004851993911,
                0.7216983581765535,
                0.49942189332714976,
                0.9879551615628381,
                0.6397449531622061,
                0.19006204680893152,
                0.8430472055480123,
                1.1646546766435881,
                0.9500261049935138,
                0.44309295746946314,
                0.1358508368580552,
                0.3816798883703529,
                0.27338858345420075,
                0.7816594872622985,
                0.17262952798760126,
                0.6397509940282039,
                0.9726910710714768,
                0.4884622629659889,
                0.5153657418426606,
                0.6884509239552447,
                0.28218636625825777,
                0.18170891100341025,
                1.0238439287023562,
                0.5311811833193738,
                0.13708520263317062,
                0.10908610193612234,
                1.1878206615564184,
                0.8234762987257431,
                0.8672424312807067,
                0.5280495815684301,
                0.29298241979657735,
                0.2331079650356981,
                0.820763509785017,
                0.37147228895234624,
                0.3252756847422742,
                0.5299536379931126,
                0.6218668446139202,
                0.27763204938804253,
                0.7343860672871082,
                0.4952074058083206,
                0.3926552058993296,
                0.8534805047331688,
                0.644221492210066,
                0.4180578029425472,
                0.240928330943711,
                0.7823803899835161,
                0.3027738899198869,
                0.655047593734567,
                0.3844632038265053,
                1.1482202916378605,
                0.9519668772892353,
                0.06058343369547032,
                0.38167253670303064,
                0.19257382120816555,
                0.6438172091646727,
                0.11146891090483704,
                0.9684526984060637,
                0.08029077342204308,
                0.784238199792356,
                0.09193586513484762,
                0.4482575170750042,
                0.2721469210969564,
                0.8561533830245281,
                0.2819164066773412,
                0.12813046873071404,
                0.10220710512214042,
                0.48691833789062233,
                0.5020412820874053,
                0.18599812753581396,
                0.8121687517652941,
                0.36880251167109485,
                0.3154454270385834,
                0.6756084201143353,
                0.2786501203117945,
                0.1723867699420482,
                0.387639637609794,
                1.0091468735648987,
                0.508242653901033,
                0.7784327379188821,
                0.30345497136360333,
                0.5190677038536843,
                0.640519131922536,
                0.376675995333094,
                0.22098157991684783,
                0.1275551497757849,
                0.4898267085245536,
                0.16600542275125038,
                0.5339688360762072,
                0.6185641988968787,
                0.27141389525603604,
                0.6361194803550576,
                0.2002839895666459,
                0.40161785054677474,
                0.23256506030896762,
                0.7857309561801054,
                0.04065046670367107,
                0.6461505408858991,
                0.12143921594383017,
                0.9537508522799212,
                0.3810247878512162,
                0.11664293801399515,
                0.2707557030633519,
                0.19113601535066962,
                0.8048760900903663,
                0.07457312461046703,
                0.9609205429919944,
                0.4523912217352152,
                0.665453663970266,
                0.2756399732493783,
                0.16540985711143172,
                0.09489567150271352,
                0.055392149520247014,
                0.7747885225815303,
                0.6287763081282781,
                0.8385676335778723,
                0.26551657566497,
                0.4846048820850437,
                0.3652190692671813,
                0.4886462194400875,
                0.1827114511884209,
                0.304095413528203,
                0.10844705733961149,
                0.5362608855801766,
                0.6154151421132362,
                0.2664682045031514,
                0.3815082477808231,
                0.08726084978323194,
                0.6474589415643739,
                0.30358803667695256,
                0.3675329660581807,
                0.6294234133875171,
                0.19990241405851114,
                0.38913394620793723,
                0.22604998240670895,
                0.5097114350346291,
                0.12855446997582357,
                0.20935662686035392,
                0.12167883056579297,
                0.4838588368279654,
                0.1609349028517144,
                0.1896934654087842,
                0.0704866079064006,
                0.6570342740496697,
                0.3796227170901338,
                0.2686507535162928,
                0.11557361441637819,
                0.11191723380548069,
                0.6124495482858422,
                0.037777293458026176,
                0.7928897019411747,
                0.48224943390257385,
                0.47826025307398207,
                0.17995845261781546,
                0.2717920010683322,
                0.15754196870719447,
                0.057054577712252835,
                0.6236890174560398,
                0.3791564610675984,
                0.08820066551006676,
                0.0519446513225905,
                0.7682570978380291,
                0.6105065905277058,
                0.502204521901495,
                0.12828104113130173,
                0.20074218906532487,
                0.11719705552286333,
                0.3610660680354089,
                0.29280769673410617,
                0.10584569115912862,
                0.2605410204316948,
                0.47882891946942724,
                0.15688498146540839,
                0.37496548143410047,
                0.08376904903991991,
                0.19895111534561313,
                0.21844432385379636,
                0.3029691060288074,
                0.35821569226438876,
                0.3883537708422474,
                0.11630035287251823,
                0.19236766318623888,
                0.4900043001447056,
                0.47089459397751593,
                0.06673990712378386,
                0.025838069420445763,
                0.6473659837086239,
                0.27289944730778637,
                0.0652972997883429,
                0.4999458645441742,
                0.19446263958430418,
                0.08407243442038238,
                0.05014174910893161,
                0.48087086803902285,
                0.03651659529935935,
                0.1781423051559216,
                0.3639734316244291,
                0.2841810936875253,
                0.1039990891343954,
                0.2717878667136994,
                0.150470552840648,
                0.05605613940226065,
                0.37368633362850345,
                0.08173837509096951,
                0.620123808420381,
                0.3638271009783325,
                0.1318876928617264,
                0.35225359836990383,
                0.07014080377359379,
                0.1127095422717575,
                0.15234170453324902,
                0.25585073477479553,
                0.20369121020089143,
                0.21101041848590366,
                0.08042702010861291,
                0.11452904880661925,
                0.26938952867501825,
                0.06444846126207517,
                0.18902512048581,
                0.035168215493436586,
                0.3599526410857824,
                0.27748698962877094,
                0.06255033639861617,
                0.024632615880804905,
                0.4577961549745699,
                0.0349031735372127,
                0.3687082462467413,
                0.26760007827270466,
                0.14481846478948573,
                0.05447485419554104,
                0.0472685671673231,
                0.3462658925090089,
                0.4891394297692007,
                0.18450728997176272,
                0.17427770797339237,
                0.10132569893693144,
                0.2508769475567074,
                0.07847918438097932,
                0.06966274513561732,
                0.20527245136019176,
                0.1299055869072704,
                0.1076613674931596,
                0.1475275362546755,
                0.1862405092051298,
                0.03443920650065171,
                0.0594897375351457,
                0.023685523090850816,
                0.2640614762666253,
                0.14031786384774558,
                0.11245293194014107,
                0.06322288701673176,
                0.07507968155152947,
                0.246716946916963,
                0.033073183194136144,
                0.2670346358960148,
                0.17109854489275034,
                0.05254813042832415,
                0.2005747243295308,
                0.0445170914555632,
                0.336695887032949,
                0.10380654121052656,
                0.09843676344528281,
                0.1437025649636455,
                0.07518308738816534,
                0.015593708397882409,
                0.03673027551266197,
                0.06888206764359746,
                0.05710916700940627,
                0.11065863757864756,
                0.033468000516455834,
                0.16839696705440485,
                0.022576171498287372,
                0.13343838317306927,
                0.06176774676398228,
                0.10070683579684733,
                0.04247506594430676,
                0.1405350156684872,
                0.031283580385173926,
                0.09607887943495196,
                0.050535646007401284,
                0.07262585482828961,
                0.1931678793222537,
                0.014981629285630805,
                0.036277308422216936,
                0.05355986280268147,
                0.04086694668560617,
                0.06050650034922825,
                0.0323814211886023,
                0.09408973918545771,
                0.021458813392771746,
                0.02993210736141958,
                0.07054484943455624,
                0.04893586040241371,
                0.09591350510261923,
                0.019902048716216418,
                0.014332669765393097,
                0.028853689269610825,
                0.031476696135032706,
                0.020595366637073804,
                0.0476115689696842,
                0.03843480451999081,
                0.013812447488514188,
                0.019894872558384832,
            ],
        );
    }

    #[test]
    #[ignore]
    fn evolve_atlaspht15() {
        execute_evolve_test(
            "NNPDF40_nlo_as_01180",
            "../ATLASPHT15-ATLASPHT15_Et_3bin.pineappl.lz4",
            "../ATLASPHT15-ATLASPHT15_Et_3bin/metadata.yaml",
            "../ATLASPHT15-ATLASPHT15_Et_3bin/alphas.npy",
            "../ATLASPHT15-ATLASPHT15_Et_3bin/operators.npy",
            &[
                847.6443451935788,
                389.2518461733433,
                197.52262245783854,
                85.5341661018376,
                30.87251252799057,
                12.857919383339262,
                5.881543988434407,
                2.5226447467912503,
                0.9852473741631275,
                0.33567120929573163,
                0.1103227781235781,
                0.03105085344647281,
                0.005857170880761543,
            ],
            &[
                846.9388599419221,
                388.9452263757685,
                197.35723798133498,
                85.45945643409014,
                30.84568815110067,
                12.847951715560423,
                5.876789744699159,
                2.5205570958215677,
                0.9845342100753295,
                0.33542931711292495,
                0.11024603718181779,
                0.031033712901365054,
                0.005854961951848535,
            ],
        );
    }
}
