//! Module containing all traits and supporting structures for grids.

use super::bin::{BinInfo, BinLimits, BinRemapper};
use super::empty_subgrid::EmptySubgridV1;
use super::evolution::{self, EvolveInfo, OperatorInfo};
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
use ndarray::{s, Array3, Array5, ArrayView5, Axis, Dimension};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::convert::{TryFrom, TryInto};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::iter;
use std::mem;
use std::ops::Range;
use std::slice;
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

    /// Return a mask suitable to pass as the `order_mask` parameter of [`Grid::convolute`],
    /// [`Grid::evolve`] or [`Grid::evolve_info`]. The selection of `orders` is controlled using
    /// the `max_as` and `max_al` parameters, for instance setting `max_as = 1` and `max_al = 0`
    /// selects the LO QCD only, `max_as = 2` and `max_al = 0` the NLO QCD; setting `max_as = 3`
    /// and `max_al = 2` would select all NLOs, and the NNLO QCD.
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
    /// assert_eq!(Order::create_mask(&orders, 0, 1, false), [true, false, false, false, false, false]);
    /// // LO QCD
    /// assert_eq!(Order::create_mask(&orders, 1, 0, false), [true, false, false, false, false, false]);
    /// // LO
    /// assert_eq!(Order::create_mask(&orders, 1, 1, false), [true, false, false, false, false, false]);
    /// // NLO QCD
    /// assert_eq!(Order::create_mask(&orders, 2, 0, false), [true, true, false, false, false, false]);
    /// // NLO EW
    /// assert_eq!(Order::create_mask(&orders, 0, 2, false), [true, false, true, false, false, false]);
    /// // NNLO QCD
    /// assert_eq!(Order::create_mask(&orders, 3, 0, false), [true, true, false, true, false, false]);
    /// // NNLO EW
    /// assert_eq!(Order::create_mask(&orders, 0, 3, false), [true, false, true, false, false, true]);
    /// ```
    ///
    /// Orders containing non-zero powers of logarithms can be selected as well if `logs` is set to
    /// `true`:
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
    /// assert_eq!(Order::create_mask(&orders, 0, 2, true), [true, false, false, true, true]);
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
    /// assert_eq!(Order::create_mask(&orders, 0, 1, false), [false, false, true, false, false, false, false]);
    /// // LO QCD
    /// assert_eq!(Order::create_mask(&orders, 1, 0, false), [true, false, false, false, false, false, false]);
    /// // LO
    /// assert_eq!(Order::create_mask(&orders, 1, 1, false), [true, true, true, false, false, false, false]);
    /// ```
    #[must_use]
    pub fn create_mask(orders: &[Self], max_as: u32, max_al: u32, logs: bool) -> Vec<bool> {
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
            .map(
                |&Self {
                     alphas,
                     alpha,
                     logxir,
                     logxif,
                 }| {
                    if !logs && (logxir > 0 || logxif > 0) {
                        return false;
                    }

                    let pto = alphas + alpha - lo;

                    alphas + alpha < min + lo
                        || (alphas + alpha < max + lo
                            && match max_as.cmp(&max_al) {
                                Ordering::Greater => lo_as + pto == alphas,
                                Ordering::Less => lo_al + pto == alpha,
                                Ordering::Equal => false,
                            })
                },
            )
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
#[deprecated(since = "0.6.0", note = "use EvolveInfo instead")]
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
#[deprecated(since = "0.6.0", note = "use OperatorInfo instead")]
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
    #[allow(deprecated)]
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
        let mut bins = vec![0.0; bin_indices.len() * xi.len()];
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

                // TODO: use let-else statement when MSRV is 1.65
                let bin_index =
                    if let Some(bin) = bin_indices.iter().position(|&index| index == bin) {
                        bin
                    } else {
                        continue;
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

        for ((imu2, ix1, ix2), value) in subgrid.indexed_iter() {
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

        let bin_indices: Vec<_> = (0..other.bin_info().bins())
            .map(|bin| {
                self.bin_info()
                    .find_bin(&other.bin_info().bin_limits(bin))
                    .unwrap_or_else(|| panic!("failed for {bin}"))
            })
            .collect();

        for ((i, j, k), subgrid) in other
            .subgrids
            .indexed_iter_mut()
            .filter(|((_, _, _), subgrid)| !subgrid.is_empty())
        {
            let other_order = &other.orders[i];
            let other_entry = &other.lumi[k];

            let self_i = self.orders.iter().position(|x| x == other_order).unwrap();
            let self_j = bin_indices[j];
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
        // first convert everything into `ImportOnlySubgridV2`

        for subgrid in self.subgrids.iter_mut() {
            if subgrid.is_empty() {
                *subgrid = EmptySubgridV1::default().into();
            } else {
                match subgrid {
                    // can't be reach because we already caught empty grids above
                    SubgridEnum::EmptySubgridV1(_) => unreachable!(),
                    // can't be optimized without losing information
                    SubgridEnum::NtupleSubgridV1(_) => continue,
                    _ => {
                        let mut new_subgrid = ImportOnlySubgridV2::from(&*subgrid).into();
                        mem::swap(subgrid, &mut new_subgrid);
                    }
                }
            }
        }

        if self
            .key_values()
            .map_or(true, |map| map["initial_state_1"] == map["initial_state_2"])
        {
            self.symmetrize_lumi();
        }

        self.optimize_orders();
        self.optimize_lumi();
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
                        }
                        lhs.merge(rhs, false);

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
                        }

                        lhs.merge(rhs, true);
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
    #[deprecated(since = "0.6.0", note = "use evolve_info instead")]
    #[allow(deprecated)]
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
    #[deprecated(since = "0.6.0", note = "use evolve instead")]
    #[allow(deprecated)]
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
            &eko_info.grid_axes.pids
        } else {
            slice::from_ref(&initial_state_1)
        };
        let pids2 = if has_pdf2 {
            &eko_info.grid_axes.pids
        } else {
            slice::from_ref(&initial_state_2)
        };
        // create target luminosities
        let tgt_pids1 = if has_pdf1 {
            &eko_info.target_pids
        } else {
            slice::from_ref(&initial_state_1)
        };
        let tgt_pids2 = if has_pdf2 {
            &eko_info.target_pids
        } else {
            slice::from_ref(&initial_state_2)
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
                let mut src_array_q2_grid: Vec<_> = (0..self.orders.len())
                    .flat_map(|order| {
                        self.subgrids[[order, bin, src_lumi]]
                            .mu2_grid()
                            .iter()
                            .map(|mu2| mu2.fac)
                            .collect::<Vec<_>>()
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

                    for ((iq2, ix1, ix2), value) in src_subgrid.indexed_iter() {
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

                                    for ((_, tgt_x1_idx, tgt_x2_idx), value) in
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

    /// Returns information for the generation of evolution operators that are being used in
    /// [`Grid::evolve`] with the parameter `order_mask`.
    #[must_use]
    pub fn evolve_info(&self, order_mask: &[bool]) -> EvolveInfo {
        use super::evolution::EVOLVE_INFO_TOL_ULPS;

        let has_pdf1 = self.has_pdf1();
        let has_pdf2 = self.has_pdf2();

        let mut ren1 = Vec::new();
        let mut fac1 = Vec::new();
        let mut x1 = Vec::new();
        let mut pids1 = Vec::new();

        for (lumi, subgrid) in self
            .subgrids()
            .indexed_iter()
            .filter_map(|(tuple, subgrid)| {
                (!subgrid.is_empty() && (order_mask.is_empty() || order_mask[tuple.0]))
                    .then_some((tuple.2, subgrid))
            })
        {
            ren1.extend(subgrid.mu2_grid().iter().map(|Mu2 { ren, .. }| *ren));
            ren1.sort_by(f64::total_cmp);
            ren1.dedup_by(|a, b| approx_eq!(f64, *a, *b, ulps = EVOLVE_INFO_TOL_ULPS));

            fac1.extend(subgrid.mu2_grid().iter().map(|Mu2 { fac, .. }| *fac));
            fac1.sort_by(f64::total_cmp);
            fac1.dedup_by(|a, b| approx_eq!(f64, *a, *b, ulps = EVOLVE_INFO_TOL_ULPS));

            if has_pdf1 {
                x1.extend(subgrid.x1_grid().iter().copied());
            }
            if has_pdf2 {
                x1.extend(subgrid.x2_grid().iter());
            }

            x1.sort_by(f64::total_cmp);
            x1.dedup_by(|a, b| approx_eq!(f64, *a, *b, ulps = EVOLVE_INFO_TOL_ULPS));

            if has_pdf1 {
                pids1.extend(self.lumi()[lumi].entry().iter().map(|(a, _, _)| a));
            }
            if has_pdf2 {
                pids1.extend(self.lumi()[lumi].entry().iter().map(|(_, b, _)| b));
            }

            pids1.sort_unstable();
            pids1.dedup();
        }

        EvolveInfo {
            fac1,
            pids1,
            x1,
            ren1,
        }
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
        operator: ArrayView5<f64>,
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
            evolution::evolve_with_two(self, &operator, info, order_mask)
        } else {
            evolution::evolve_with_one(self, &operator, info, order_mask)
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

        Ok(FkTable::try_from(grid).unwrap_or_else(|_| unreachable!()))
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
    #[must_use]
    pub fn has_pdf1(&self) -> bool {
        let initial_state_1 = self.initial_state_1();

        !self
            .lumi()
            .iter()
            .all(|entry| entry.entry().iter().all(|&(a, _, _)| a == initial_state_1))
    }

    /// Returns `true` if the second initial state needs a convolution, `false` otherwise.
    #[must_use]
    pub fn has_pdf2(&self) -> bool {
        let initial_state_2 = self.initial_state_2();

        !self
            .lumi()
            .iter()
            .all(|entry| entry.entry().iter().all(|&(_, b, _)| b == initial_state_2))
    }

    /// Returns the particle identifier of the first initial state. This is usually but not always
    /// a proton, which is represented by the PDG ID `2212`.
    ///
    /// # Panics
    ///
    /// TODO
    #[must_use]
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
    ///
    /// # Panics
    ///
    /// TODO
    #[must_use]
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
    use float_cmp::assert_approx_eq;
    use std::fs::File;

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
            Order::create_mask(&orders, 0, 0, false),
            [false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 0, 1, false),
            [true, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 0, 2, false),
            [true, false, true, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 0, 3, false),
            [true, false, true, false, false, true]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 0, false),
            [true, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 1, false),
            [true, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 2, false),
            [true, false, true, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 3, false),
            [true, false, true, false, false, true]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 0, false),
            [true, true, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 1, false),
            [true, true, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 2, false),
            [true, true, true, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 3, false),
            [true, true, true, false, false, true]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 0, false),
            [true, true, false, true, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 1, false),
            [true, true, false, true, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 2, false),
            [true, true, true, true, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 3, false),
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
            Order::create_mask(&orders, 0, 0, false),
            [false, false, false, false, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 0, 1, false),
            [false, false, true, false, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 0, 2, false),
            [false, false, true, false, false, false, true, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 0, 3, false),
            [false, false, true, false, false, false, true, false, false, false, false, true]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 0, false),
            [true, false, false, false, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 1, false),
            [true, true, true, false, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 2, false),
            [true, true, true, false, false, false, true, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 1, 3, false),
            [true, true, true, false, false, false, true, false, false, false, false, true]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 0, false),
            [true, false, false, true, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 1, false),
            [true, true, true, true, false, false, false, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 2, false),
            [true, true, true, true, true, true, true, false, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 2, 3, false),
            [true, true, true, true, true, true, true, false, false, false, false, true]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 0, false),
            [true, false, false, true, false, false, false, true, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 1, false),
            [true, true, true, true, false, false, false, true, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 2, false),
            [true, true, true, true, true, true, true, true, false, false, false, false]
        );
        assert_eq!(
            Order::create_mask(&orders, 3, 3, false),
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
            grid.key_values().unwrap().get("initial_state_1").unwrap(),
            "-2212"
        );
        assert_eq!(
            grid.key_values().unwrap().get("initial_state_2").unwrap(),
            "-2212"
        );
    }

    // TODO: properly test axes returned

    #[allow(deprecated)]
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
        let sparse_array = SparseArray3::from_ndarray(array.view(), 0, 3);
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
    #[allow(deprecated)]
    fn grid_axes() {
        let (grid, axes) = simple_grid();

        let ret_axes = grid.axes().unwrap();
        assert_eq!(ret_axes.x_grid, axes.x_grid);
        assert_eq!(ret_axes.mur2_grid, axes.mur2_grid);
        assert_eq!(ret_axes.muf2_grid, axes.muf2_grid);
        assert_eq!(ret_axes.pids, vec![]);
    }

    #[test]
    #[allow(deprecated)]
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
            target_x_grid,
            target_pids,
            grid_axes: GridAxes {
                x_grid: axes.x_grid,
                pids: axes.pids,
                mur2_grid: axes.mur2_grid.clone(),
                muf2_grid: axes.muf2_grid,
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

    #[test]
    fn evolve_info() {
        let grid =
            Grid::read(File::open("../pineappl_cli/data/LHCB_WP_7TEV.pineappl.lz4").unwrap())
                .unwrap();
        let info = grid.evolve_info(&[]);

        assert_eq!(info.fac1.len(), 1);
        assert_approx_eq!(f64, info.fac1[0], 6456.443904000001, ulps = 64);

        assert_eq!(info.pids1, [-3, -1, 0, 2, 4, 22]);

        assert_eq!(info.x1.len(), 50);
        assert_approx_eq!(f64, info.x1[0], 1.9999999999999954e-7, ulps = 64);

        assert_approx_eq!(f64, info.x1[1], 3.034304765867952e-7, ulps = 64);
        assert_approx_eq!(f64, info.x1[2], 4.6035014748963906e-7, ulps = 64);
        assert_approx_eq!(f64, info.x1[3], 6.984208530700364e-7, ulps = 64);
        assert_approx_eq!(f64, info.x1[4], 1.0596094959101024e-6, ulps = 64);
        assert_approx_eq!(f64, info.x1[5], 1.607585498470808e-6, ulps = 64);
        assert_approx_eq!(f64, info.x1[6], 2.438943292891682e-6, ulps = 64);
        assert_approx_eq!(f64, info.x1[7], 3.7002272069854957e-6, ulps = 64);
        assert_approx_eq!(f64, info.x1[8], 5.613757716930151e-6, ulps = 64);
        assert_approx_eq!(f64, info.x1[9], 8.516806677573355e-6, ulps = 64);
        assert_approx_eq!(f64, info.x1[10], 1.292101569074731e-5, ulps = 64);
        assert_approx_eq!(f64, info.x1[11], 1.9602505002391748e-5, ulps = 64);
        assert_approx_eq!(f64, info.x1[12], 2.97384953722449e-5, ulps = 64);
        assert_approx_eq!(f64, info.x1[13], 4.511438394964044e-5, ulps = 64);
        assert_approx_eq!(f64, info.x1[14], 6.843744918967896e-5, ulps = 64);
        assert_approx_eq!(f64, info.x1[15], 0.00010381172986576898, ulps = 64);
        assert_approx_eq!(f64, info.x1[16], 0.00015745605600841445, ulps = 64);
        assert_approx_eq!(f64, info.x1[17], 0.00023878782918561914, ulps = 64);
        assert_approx_eq!(f64, info.x1[18], 0.00036205449638139736, ulps = 64);
        assert_approx_eq!(f64, info.x1[19], 0.0005487795323670796, ulps = 64);
        assert_approx_eq!(f64, info.x1[20], 0.0008314068836488144, ulps = 64);
        assert_approx_eq!(f64, info.x1[21], 0.0012586797144272762, ulps = 64);
        assert_approx_eq!(f64, info.x1[22], 0.0019034634022867384, ulps = 64);
        assert_approx_eq!(f64, info.x1[23], 0.0028738675812817515, ulps = 64);
        assert_approx_eq!(f64, info.x1[24], 0.004328500638820811, ulps = 64);
        assert_approx_eq!(f64, info.x1[25], 0.006496206194633799, ulps = 64);
        assert_approx_eq!(f64, info.x1[26], 0.009699159574043398, ulps = 64);
        assert_approx_eq!(f64, info.x1[27], 0.014375068581090129, ulps = 64);
        assert_approx_eq!(f64, info.x1[28], 0.02108918668378717, ulps = 64);
        assert_approx_eq!(f64, info.x1[29], 0.030521584007828916, ulps = 64);
        assert_approx_eq!(f64, info.x1[30], 0.04341491741702269, ulps = 64);
        assert_approx_eq!(f64, info.x1[31], 0.060480028754447364, ulps = 64);
        assert_approx_eq!(f64, info.x1[32], 0.08228122126204893, ulps = 64);
        assert_approx_eq!(f64, info.x1[33], 0.10914375746330703, ulps = 64);
        assert_approx_eq!(f64, info.x1[34], 0.14112080644440345, ulps = 64);
        assert_approx_eq!(f64, info.x1[35], 0.17802566042569432, ulps = 64);
        assert_approx_eq!(f64, info.x1[36], 0.2195041265003886, ulps = 64);
        assert_approx_eq!(f64, info.x1[37], 0.2651137041582823, ulps = 64);
        assert_approx_eq!(f64, info.x1[38], 0.31438740076927585, ulps = 64);
        assert_approx_eq!(f64, info.x1[39], 0.3668753186482242, ulps = 64);
        assert_approx_eq!(f64, info.x1[40], 0.4221667753589648, ulps = 64);
        assert_approx_eq!(f64, info.x1[41], 0.4798989029610255, ulps = 64);
        assert_approx_eq!(f64, info.x1[42], 0.5397572337880445, ulps = 64);
        assert_approx_eq!(f64, info.x1[43], 0.601472197967335, ulps = 64);
        assert_approx_eq!(f64, info.x1[44], 0.6648139482473823, ulps = 64);
        assert_approx_eq!(f64, info.x1[45], 0.7295868442414312, ulps = 64);
        assert_approx_eq!(f64, info.x1[46], 0.7956242522922756, ulps = 64);
        assert_approx_eq!(f64, info.x1[47], 0.8627839323906108, ulps = 64);
        assert_approx_eq!(f64, info.x1[48], 0.9309440808717544, ulps = 64);
        assert_approx_eq!(f64, info.x1[49], 1.0, ulps = 64);

        assert_eq!(info.ren1.len(), 1);
        assert_approx_eq!(f64, info.ren1[0], 6456.443904000001, ulps = 64);
    }
}
