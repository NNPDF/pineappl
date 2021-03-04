//! Module containing all traits and supporting structures for grids.

use super::bin::{BinInfo, BinLimits, BinRemapper};
use super::lagrange_subgrid::{LagrangeSparseSubgridV1, LagrangeSubgridV1, LagrangeSubgridV2};
use super::lumi::LumiEntry;
use super::lumi_entry;
use super::ntuple_subgrid::NtupleSubgridV1;
use super::read_only_sparse_subgrid::ReadOnlySparseSubgridV1;
use super::sparse_array3::SparseArray3;
use super::subgrid::{ExtraSubgridParams, Subgrid, SubgridEnum, SubgridParams};
use either::Either::{Left, Right};
use float_cmp::approx_eq;
use git_version::git_version;
use itertools::Itertools;
use lz_fear::{framed::DecompressionError::WrongMagic, LZ4FrameReader};
use ndarray::{Array3, Dimension};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::convert::TryInto;
use std::io::{Read, Seek, SeekFrom, Write};
use std::mem;
use std::ptr;
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

    /// Compares two vectors of `Order` for equality after sorting them.
    #[must_use]
    pub fn equal_after_sort(lhs: &[Self], rhs: &[Self]) -> bool {
        let mut lhs = lhs.to_vec();
        let mut rhs = rhs.to_vec();

        lhs.sort();
        rhs.sort();

        lhs == rhs
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
pub enum GridMergeError {
    /// Returned when trying to merge two `Grid` objects with different bin limits and different
    /// orders.
    #[error("the merged grid has different orders")]
    DifferentOrders,
    /// Returned when trying to merge two `Grid` Objects with different bin limits and different
    /// luminosity functions.
    #[error("the merged grid has a different luminosity function")]
    DifferentLumi,
    /// Returned when trying to merge two `Grid` objects with incompatible bin limits.
    #[error(transparent)]
    DifferentBins(super::bin::MergeBinError),
}

/// Error returned when trying to construct a `Grid` using an unknown subgrid type.
#[derive(Debug, Error)]
#[error("tried constructing a Grid with unknown Subgrid type `{0}`")]
pub struct UnknownSubgrid(String);

#[derive(Clone, Deserialize, Serialize)]
struct Mmv1 {}

#[derive(Clone, Deserialize, Serialize)]
struct Mmv2 {
    remapper: Option<BinRemapper>,
    key_value_db: HashMap<String, String>,
}

impl Default for Mmv2 {
    fn default() -> Self {
        Self {
            remapper: None,
            key_value_db: [
                (
                    "pineappl_gitversion".to_owned(),
                    git_version!(args = ["--always", "--dirty", "--long", "--tags"]).to_owned(),
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

#[derive(Clone, Deserialize, Serialize)]
enum MoreMembers {
    V1(Mmv1),
    V2(Mmv2),
}

impl MoreMembers {
    fn upgrade(&mut self) {
        match self {
            Self::V1(_) => {
                *self = Self::V2(Mmv2::default());
            }
            Self::V2(_) => {}
        }
    }
}

/// Error type returned by `Grid::set_remapper`.
#[derive(Debug, Error)]
pub enum GridSetBinRemapperError {
    /// Returned if the number of bins in the grid and in the remapper do not agree.
    #[error("the number of bins in the remapper, {remapper_bins}, does not agree with the number of bins in the grid: {grid_bins}")]
    BinMismatch {
        /// Number of bins in the grid.
        grid_bins: usize,
        /// Number of bins in the remapper.
        remapper_bins: usize,
    },
}

pub struct EkoInfo {
    pub x_grid: Vec<f64>,
    pub q2_grid: Vec<f64>,
}

/// Main data structure of `PineAPPL`. This structure contains a `Subgrid` for each `LumiEntry`,
/// bin, and coupling order it was created with.
#[derive(Deserialize, Serialize)]
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
                || {
                    LagrangeSubgridV2::new(
                        &subgrid_params,
                        &ExtraSubgridParams::from(&subgrid_params),
                    )
                    .into()
                },
            ),
            orders,
            lumi,
            bin_limits: BinLimits::new(bin_limits),
            subgrid_params,
            more_members: MoreMembers::V2(Mmv2::default()),
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
    ) -> Result<Self, UnknownSubgrid> {
        let subgrid_maker: Box<dyn Fn() -> SubgridEnum> = match subgrid_type {
            "LagrangeSubgrid" | "LagrangeSubgridV2" => {
                Box::new(|| LagrangeSubgridV2::new(&subgrid_params, &extra).into())
            }
            "LagrangeSubgridV1" => Box::new(|| LagrangeSubgridV1::new(&subgrid_params).into()),
            "NtupleSubgrid" => Box::new(|| NtupleSubgridV1::new().into()),
            "LagrangeSparseSubgrid" => {
                Box::new(|| LagrangeSparseSubgridV1::new(&subgrid_params).into())
            }
            _ => return Err(UnknownSubgrid(subgrid_type.to_string())),
        };

        Ok(Self {
            subgrids: Array3::from_shape_simple_fn(
                (orders.len(), bin_limits.len() - 1, lumi.len()),
                subgrid_maker,
            ),
            orders,
            lumi,
            bin_limits: BinLimits::new(bin_limits),
            subgrid_params,
            more_members: MoreMembers::V2(Mmv2::default()),
        })
    }

    /// Performs a convolution of the contained subgrids with the given PDFs, `xfx1` for the first
    /// parton and `xfx2` for the second parton, `alphas` for the evaluation of the strong
    /// coupling. The parameters `order_mask` and `lumi_mask` can be used to selectively enable
    /// perturbative orders and luminosities; they must either be empty (everything enabled) or as
    /// large as the orders and luminosity function, respectively. If the corresponding entry is
    /// `true` the order/luminosity is enable, `false` disables the entry. The tuple `xi` can be
    /// used to independently vary the renormalization (first element) and factorization scale
    /// (second element) from their central value `(1.0, 1.0)`.
    pub fn convolute(
        &self,
        xfx1: &dyn Fn(i32, f64, f64) -> f64,
        xfx2: &dyn Fn(i32, f64, f64) -> f64,
        alphas: &dyn Fn(f64) -> f64,
        order_mask: &[bool],
        bin_indices: &[usize],
        lumi_mask: &[bool],
        xi: &[(f64, f64)],
    ) -> Vec<f64> {
        let bin_indices = if bin_indices.is_empty() {
            (0..self.bin_limits.bins()).collect()
        } else {
            bin_indices.to_vec()
        };
        let mut bins: Vec<f64> = vec![0.0; bin_indices.len() * xi.len()];
        let bin_sizes = self.bin_info().normalizations();

        let pdf_cache1 = RefCell::new(FxHashMap::default());
        let pdf_cache2 = RefCell::new(FxHashMap::default());
        let alphas_cache = RefCell::new(FxHashMap::default());
        let mut last_xif = 0.0;

        let mut q2_grid = self.subgrids[[0, 0, 0]].q2_grid();
        let mut x1_grid = self.subgrids[[0, 0, 0]].x1_grid();
        let mut x2_grid = self.subgrids[[0, 0, 0]].x2_grid();
        let use_cache = !q2_grid.is_empty() && !x1_grid.is_empty() && !x2_grid.is_empty();
        let two_caches = !ptr::eq(&xfx1, &xfx2);

        let mut xir_values: Vec<_> = xi.iter().map(|xi| xi.0).collect();
        xir_values.sort_by(|lhs, rhs| lhs.partial_cmp(rhs).unwrap());
        xir_values.dedup();
        let xir_indices: Vec<_> = xi
            .iter()
            .map(|xi| xir_values.iter().position(|xir| xi.0 == *xir).unwrap())
            .collect();

        // iterate over the elements of `xi` and a corresponding index, but sorted using the
        // factorisation value of `xi`
        for (l, (&(xir, xif), &xir_index)) in xi
            .iter()
            .zip(xir_indices.iter())
            .enumerate()
            .sorted_by(|lhs, rhs| (lhs.1).1.partial_cmp((rhs.1).1).unwrap())
        {
            // whenever the value `xif` changes we can clear the PDF cache
            if xif != last_xif {
                pdf_cache1.borrow_mut().clear();
                pdf_cache2.borrow_mut().clear();
                last_xif = xif;
            }

            for ((i, j, k), subgrid) in self.subgrids.indexed_iter() {
                let order = &self.orders[i];

                if ((order.logxir > 0) && (xir == 1.0)) || ((order.logxif > 0) && (xif == 1.0)) {
                    continue;
                }

                if (!order_mask.is_empty() && !order_mask[i])
                    || (!lumi_mask.is_empty() && !lumi_mask[k])
                {
                    continue;
                }

                let bin_index = match bin_indices.iter().position(|&index| index == j) {
                    Some(i) => i,
                    None => continue,
                };

                let lumi_entry = &self.lumi[k];

                let mut value = if subgrid.is_empty() {
                    0.0
                } else if use_cache {
                    let new_q2_grid = subgrid.q2_grid();
                    let new_x1_grid = subgrid.x1_grid();
                    let new_x2_grid = subgrid.x2_grid();
                    let q2_grid_changed = new_q2_grid != q2_grid;

                    if q2_grid_changed {
                        q2_grid = new_q2_grid;
                        alphas_cache.borrow_mut().clear();
                    }

                    if q2_grid_changed || (new_x1_grid != x1_grid) || (new_x2_grid != x2_grid) {
                        x1_grid = new_x1_grid;
                        x2_grid = new_x2_grid;
                        pdf_cache1.borrow_mut().clear();
                        pdf_cache2.borrow_mut().clear();
                    }

                    subgrid.convolute(
                        &x1_grid,
                        &x2_grid,
                        &q2_grid,
                        Left(&|ix1, ix2, iq2| {
                            let mut pdf_cache1 = pdf_cache1.borrow_mut();
                            let mut pdf_cache2 = pdf_cache2.borrow_mut();
                            let x1 = x1_grid[ix1];
                            let x2 = x2_grid[ix2];
                            let q2 = q2_grid[iq2];
                            let q2f = xif * xif * q2;

                            let mut lumi = 0.0;

                            for entry in lumi_entry.entry() {
                                let xfx1 = *pdf_cache1
                                    .entry((entry.0, ix1, iq2))
                                    .or_insert_with(|| xfx1(entry.0, x1, q2f));
                                let xfx2 = if two_caches {
                                    *pdf_cache2
                                        .entry((entry.1, ix2, iq2))
                                        .or_insert_with(|| xfx2(entry.1, x2, q2f))
                                } else {
                                    *pdf_cache1
                                        .entry((entry.1, ix2, iq2))
                                        .or_insert_with(|| xfx2(entry.1, x2, q2f))
                                };
                                lumi += xfx1 * xfx2 * entry.2 / (x1 * x2);
                            }

                            let mut alphas_cache = alphas_cache.borrow_mut();
                            let alphas = alphas_cache
                                .entry(xir_values.len() * iq2 + xir_index)
                                .or_insert_with(|| alphas(xir * xir * q2));

                            lumi *= alphas.powi(order.alphas.try_into().unwrap());
                            lumi
                        }),
                    )
                } else {
                    subgrid.convolute(
                        &x1_grid,
                        &x2_grid,
                        &q2_grid,
                        Right(&|x1, x2, q2| {
                            let mut lumi = 0.0;
                            let q2f = xif * xif * q2;

                            for entry in lumi_entry.entry() {
                                let xfx1 = xfx1(entry.0, x1, q2f);
                                let xfx2 = xfx2(entry.1, x2, q2f);
                                lumi += xfx1 * xfx2 * entry.2 / (x1 * x2);
                            }

                            lumi *= alphas(xir * xir * q2).powi(order.alphas.try_into().unwrap());
                            lumi
                        }),
                    )
                };

                if order.logxir > 0 {
                    value *= (xir * xir).ln().powi(order.logxir.try_into().unwrap());
                }

                if order.logxif > 0 {
                    value *= (xif * xif).ln().powi(order.logxif.try_into().unwrap());
                }

                bins[l + xi.len() * bin_index] += value / bin_sizes[j];
            }
        }

        bins
    }

    /// Fills the grid with an ntuple for the given `order`, `observable`, and `lumi`.
    pub fn fill(&mut self, order: usize, observable: f64, lumi: usize, ntuple: &Ntuple<f64>) {
        if let Some(bin) = self.bin_limits.index(observable) {
            self.subgrids[[order, bin, lumi]].fill(ntuple);
        }
    }

    /// Constructs a `Grid` by deserializing it from `reader`. Reading is not buffered.
    ///
    /// # Errors
    ///
    /// If reading from the compressed or uncompressed stream fails an error is returned.
    pub fn read(mut reader: impl Read + Seek) -> anyhow::Result<Self> {
        match LZ4FrameReader::new(&mut reader) {
            Ok(reader) => Ok(bincode::deserialize_from(reader.into_read())?),
            Err(WrongMagic(_)) => {
                reader.seek(SeekFrom::Start(0))?;
                Ok(bincode::deserialize_from(reader)?)
            }
            Err(e) => Err(anyhow::Error::new(e)),
        }
    }

    /// Serializes `self` into `writer`. Writing is not buffered.
    ///
    /// # Errors
    ///
    /// If writing fails an error is returned.
    pub fn write(&self, writer: impl Write) -> anyhow::Result<()> {
        Ok(bincode::serialize_into(writer, self)?)
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

    /// Merges the non-empty `Subgrid`s contained in `other` into `self`. This performs one of two
    /// possible operations:
    /// 1. If the bin limits of `self` and `other` are different and can be concatenated with each
    ///    other the bins are merged. In this case both grids are assumed to have the same orders
    ///    and the same luminosity functions. If this is not the case, an error is returned.
    /// 2. If the bin limits of `self` and `other` are the same, the luminosity functions and
    ///    perturbative orders of `self` and `other` may be different.
    ///
    /// # Errors
    ///
    /// If in the first case describe above the perturbative orders or the luminosity function is
    /// different an error is returned.
    pub fn merge(&mut self, mut other: Self) -> Result<(), GridMergeError> {
        if self.bin_limits == other.bin_limits {
            let mut new_orders: Vec<Order> = Vec::new();
            let mut new_entries: Vec<LumiEntry> = Vec::new();

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

            if !new_orders.is_empty() || !new_entries.is_empty() {
                self.increase_shape(&(new_orders.len(), 0, new_entries.len()));
            }

            self.orders.append(&mut new_orders);
            self.lumi.append(&mut new_entries);
        } else {
            if !Order::equal_after_sort(&self.orders, &other.orders) {
                return Err(GridMergeError::DifferentOrders);
            }

            if !LumiEntry::equal_after_sort(&self.lumi, &other.lumi) {
                return Err(GridMergeError::DifferentLumi);
            }

            if let Err(e) = self.bin_limits.merge(&other.bin_limits) {
                return Err(GridMergeError::DifferentBins(e));
            }

            let new_bins = other.bin_limits.bins();

            self.increase_shape(&(0, new_bins, 0));

            // TODO: figure out a better strategy than removing the remapper
            match &mut self.more_members {
                MoreMembers::V1(_) => {}
                MoreMembers::V2(mmv2) => {
                    mmv2.remapper = None;
                }
            }
        }

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
            || self.subgrids[[0, 0, 0]].clone_empty(),
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

    /// Returns the subgrid parameters.
    #[must_use]
    pub fn orders(&self) -> &[Order] {
        &self.orders
    }

    /// Returns the subgrid with the specified indices `order`, `bin`, and `lumi`.
    #[must_use]
    pub fn subgrid(&self, order: usize, bin: usize, lumi: usize) -> &SubgridEnum {
        &self.subgrids[[order, bin, lumi]]
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
    pub fn set_remapper(&mut self, remapper: BinRemapper) -> Result<(), GridSetBinRemapperError> {
        if remapper.bins() != self.bin_limits.bins() {
            return Err(GridSetBinRemapperError::BinMismatch {
                grid_bins: self.bin_limits.bins(),
                remapper_bins: remapper.bins(),
            });
        }

        self.more_members.upgrade();

        match &mut self.more_members {
            MoreMembers::V1(_) => unreachable!(),
            MoreMembers::V2(mmv2) => mmv2.remapper = Some(remapper),
        }

        Ok(())
    }

    /// Returns all information about the bins in this grid.
    #[must_use]
    pub const fn bin_info(&self) -> BinInfo {
        BinInfo::new(
            &self.bin_limits,
            match &self.more_members {
                MoreMembers::V1(_) => None,
                MoreMembers::V2(mmv2) => mmv2.remapper.as_ref(),
            },
        )
    }

    /// Optimize the internal datastructures for space efficiency. This changes all subgrids of
    /// type `LagrangeSubgrid` to `LagrangeSparseSubgrid`.
    pub fn optimize(&mut self) {
        if self
            .key_values()
            .map_or(true, |map| map["initial_state_1"] == map["initial_state_2"])
        {
            self.symmetrize();
        }

        for subgrid in self.subgrids.iter_mut() {
            match subgrid {
                SubgridEnum::LagrangeSubgridV1(grid) => {
                    let mut new_subgrid = LagrangeSparseSubgridV1::from(&*grid).into();
                    mem::swap(subgrid, &mut new_subgrid);
                }
                SubgridEnum::LagrangeSubgridV2(grid) => {
                    let mut new_subgrid = ReadOnlySparseSubgridV1::from(&*grid).into();
                    mem::swap(subgrid, &mut new_subgrid);
                }
                SubgridEnum::LagrangeSparseSubgridV1(_)
                | SubgridEnum::ReadOnlySparseSubgridV1(_) => {
                    // nothing to optimize here
                }
                SubgridEnum::NtupleSubgridV1(_) => todo!(),
            }
        }
    }

    fn symmetrize(&mut self) {
        let mut indices: Vec<usize> = (0..self.lumi.len()).rev().collect();
        let mut pairs: Vec<(usize, usize)> = Vec::new();

        while let Some(index) = indices.pop() {
            let lumi_entry = &self.lumi[index];

            if *lumi_entry == lumi_entry.transpose() {
                pairs.push((index, index));
            } else {
                let (j, &other_index) = indices
                    .iter()
                    .enumerate()
                    .find(|(_, i)| self.lumi[**i] == lumi_entry.transpose())
                    .unwrap();

                pairs.push((index, other_index));
                indices.remove(j);
            }
        }

        let i_size = self.orders.len();
        let j_size = self.bin_limits.bins();
        let k_size = self.lumi.len();

        let subgrids = mem::replace(
            &mut self.subgrids,
            Array3::from_shape_vec((0, 0, 0), vec![]).unwrap(),
        );
        let mut subgrids = Array3::from_shape_vec(
            (i_size, j_size, k_size),
            subgrids.into_raw_vec().into_iter().map(Some).collect(),
        )
        .unwrap();

        for i in 0..i_size {
            for j in 0..j_size {
                for &(k1, k2) in &pairs {
                    if k1 == k2 {
                        subgrids[[i, j, k1]].as_mut().unwrap().symmetrize();
                    } else {
                        let subgrid = mem::replace(&mut subgrids[[i, j, k2]], None);
                        subgrids[[i, j, k1]]
                            .as_mut()
                            .unwrap()
                            .merge(&mut subgrid.unwrap(), true);
                    }
                }
            }
        }

        self.subgrids = Array3::from_shape_vec(
            (i_size, j_size, pairs.len()),
            subgrids
                .into_raw_vec()
                .into_iter()
                .filter_map(|s| s)
                .collect(),
        )
        .unwrap();

        self.lumi = pairs
            .iter()
            .map(|(index, _)| self.lumi[*index].clone())
            .collect();
    }

    /// Upgrades the internal data structures to their latest versions.
    pub fn upgrade(&mut self) {
        self.more_members.upgrade();
    }

    /// Returns a map with key-value pairs, if there are any stored in this grid.
    #[must_use]
    pub const fn key_values(&self) -> Option<&HashMap<String, String>> {
        match &self.more_members {
            MoreMembers::V2(mmv2) => Some(&mmv2.key_value_db),
            MoreMembers::V1(_) => None,
        }
    }

    /// Returns a map with key-value pairs, if there are any stored in this grid.
    #[must_use]
    pub fn key_values_mut(&mut self) -> &mut HashMap<String, String> {
        self.more_members.upgrade();

        match &mut self.more_members {
            MoreMembers::V1(_) => unreachable!(),
            MoreMembers::V2(mmv2) => &mut mmv2.key_value_db,
        }
    }

    /// Sets a specific key-value pair in this grid.
    pub fn set_key_value(&mut self, key: &str, value: &str) {
        self.more_members.upgrade();

        let mmv2 = match &mut self.more_members {
            MoreMembers::V1(_) => unreachable!(),
            MoreMembers::V2(mmv2) => mmv2,
        };

        mmv2.key_value_db.insert(key.to_owned(), value.to_owned());
    }

    pub fn eko_info(&self) -> Option<EkoInfo> {
        let mut q2_grid = Vec::<f64>::new();
        let mut x_grid = Vec::<f64>::new();

        for subgrid in &self.subgrids {
            if !subgrid.is_empty() {
                if q2_grid.is_empty() {
                    q2_grid = subgrid.q2_grid();
                    let x = subgrid.x1_grid();

                    // if subgrid.x2_grid() != x {
                    // return None;
                    // }

                    x_grid = x;
                }
                // } else if (subgrid.x1_grid() != x_grid) || (subgrid.x2_grid() != x_grid) {
                // return None;
                // }

                q2_grid.append(&mut subgrid.q2_grid());
                q2_grid.sort_by(|a, b| a.partial_cmp(b).unwrap());
                q2_grid.dedup();
            }
        }

        Some(EkoInfo { x_grid, q2_grid })
    }

    /// Applies an evolution kernel operator (EKO) to the grids to evolve them from different
    /// values of the factorization scale to a single one given by the parameter `q2`.
    #[must_use]
    pub fn convolute_eko(
        &self,
        q2: f64,
        alphas: &[f64],
        (xir, xif): (f64, f64),
        pids: &[i32],
        x_grid: Vec<f64>,
        q2_grid: Vec<f64>,
        operator: Vec<Vec<Vec<Vec<Vec<f64>>>>>,
    ) -> Option<Self> {
        // let EkoInfo { x_grid, q2_grid } = if let Some(eko_info) = self.eko_info() {
        // eko_info
        // } else {
        // return None;
        // };
        println!("{:?}\n{:?}", x_grid, q2_grid);

        let initial_state_1 = self.key_values().map_or(2212, |map| {
            map.get("initial_state_1").unwrap().parse::<i32>().unwrap()
        });
        let initial_state_2 = self.key_values().map_or(2212, |map| {
            map.get("initial_state_2").unwrap().parse::<i32>().unwrap()
        });
        let has_pdf1 = match initial_state_1 {
            2212 | -2212 => true,
            11 | 13 | -11 | -13 => false,
            _ => unimplemented!(),
        };
        let has_pdf2 = match initial_state_2 {
            2212 | -2212 => true,
            11 | 13 | -11 | -13 => false,
            _ => unimplemented!(),
        };

        let pids1 = if has_pdf1 {
            pids.to_vec()
        } else {
            vec![initial_state_1]
        };
        let pids2 = if has_pdf2 {
            pids.to_vec()
        } else {
            vec![initial_state_2]
        };

        let lumi: Vec<_> = pids1
            .iter()
            .cartesian_product(pids2.iter())
            .map(|(a, b)| lumi_entry![*a, *b, 1.0])
            .collect();

        let q2low_grid = vec![q2];
        let x1low_grid = if has_pdf1 { x_grid.clone() } else { vec![1.] };
        let x2low_grid = if has_pdf2 { x_grid.clone() } else { vec![1.] };
        let x1_len = x1low_grid.len();
        let x2_len = x2low_grid.len();

        let mut result = Self {
            subgrids: Array3::from_shape_simple_fn(self.subgrids.dim(), || {
                ReadOnlySparseSubgridV1::new(
                    SparseArray3::new(1, x1_len, x2_len),
                    q2low_grid.clone(),
                    x1low_grid.clone(),
                    x2low_grid.clone(),
                )
                .into()
            }),
            lumi: lumi.clone(),
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

        // TODO: put original perturbative orders and order of the EKO inside new metadata

        // println!("{:?}, {:?}", pids1, pids2);

        // Iterate over RESULT = LOW
        for ((_, bin, low_lumi), subgrid) in result.subgrids.indexed_iter_mut() {
            let mut array = SparseArray3::<f64>::new(1, x1_len, x2_len);

            let pid_low1_idx = pids1
                .iter()
                .position(|pid| *pid == lumi[low_lumi].entry()[0].0)
                .unwrap();
            let pid_low2_idx = pids2
                .iter()
                .position(|pid| *pid == lumi[low_lumi].entry()[0].1)
                .unwrap();

            // Iterate over SELF = HIGH
            for (order_idx, order) in self.orders.iter().enumerate() {
                // if order_idx > 0 {
                // continue;
                // }
                // skip log grids if we don't want the scale varied from the central choice
                if ((order.logxir > 0) && (xir == 1.0)) || ((order.logxif > 0) && (xif == 1.0)) {
                    continue;
                }

                // Iterate over SELF = HIGH
                for (high_lumi_idx, high_lumi) in self.lumi.iter().enumerate() {
                    // println!("ll: {}, o: {}, hl: {}", low_lumi, order_idx, high_lumi_idx);

                    let subgrid_high = &self.subgrids[[order_idx, bin, high_lumi_idx]];
                    if subgrid_high.is_empty() {
                        continue;
                    }
                    let x1high_grid = self.subgrids[[order_idx, bin, high_lumi_idx]].x1_grid();
                    let x2high_grid = self.subgrids[[order_idx, bin, high_lumi_idx]].x2_grid();
                    let q2high_grid = self.subgrids[[order_idx, bin, high_lumi_idx]].q2_grid();

                    // Iterate over SELF = HIGH
                    for (pid_high1, pid_high2, factor) in high_lumi.entry() {
                        // TODO: call the gluon 21 instead of 0
                        if pid_high1 == &0 || pid_high2 == &0 {
                            continue;
                        }
                        let pid_high1_idx = pids1.iter().position(|pid| pid == pid_high1).unwrap();
                        let pid_high2_idx = pids2.iter().position(|pid| pid == pid_high2).unwrap();

                        // Iterate over RESULT = LOW
                        for (x1_low, x2_low) in (0..x1_len).cartesian_product(0..x2_len) {
                            let eko_x1_low_idx = x_grid
                                .iter()
                                .position(|x| x == &x1low_grid[x1_low])
                                .unwrap();
                            let eko_x2_low_idx = x_grid
                                .iter()
                                .position(|x| x == &x2low_grid[x2_low])
                                .unwrap();
                            // Iterate over SELF = HIGH
                            let convoluted = subgrid_high.convolute(
                                &x1high_grid,
                                &x2high_grid,
                                &[],
                                Left(&|x1_high_idx, x2_high_idx, q2_index| {
                                    // index for the EK operator
                                    let eko_x1_high_idx = x_grid
                                        .iter()
                                        .position(|x| x == &x1high_grid[x1_high_idx])
                                        .unwrap();
                                    let eko_x2_high_idx = x_grid
                                        .iter()
                                        .position(|x| x == &x2high_grid[x2_high_idx])
                                        .unwrap();
                                    let eko_q2_index = q2_grid
                                        .iter()
                                        .position(|q2| q2 == &q2high_grid[q2_index])
                                        .unwrap_or_else(|| {
                                            eprintln!("{}: {}", q2_index, q2high_grid[q2_index]);
                                            panic!();
                                        });
                                    // assert_eq!(eko_x1_high_idx + x1_high_idx, x_grid.len() - 1);
                                    // assert_eq!(eko_x1_low_idx + x1_low, x_grid.len() - 1);
                                    // assert_eq!(eko_x2_high_idx + x2_high_idx, x_grid.len() - 1);
                                    // assert_eq!(eko_x2_low_idx + x2_low_idx, x_grid.len() - 1);
                                    let op1 = if has_pdf1 {
                                        operator[eko_q2_index][pid_high1_idx]
                                            [x_grid.len() - 1 - x1_high_idx][pid_low1_idx]
                                            [x_grid.len() - 1 - x1_low]
                                    } else {
                                        1.
                                    };
                                    let op2 = if has_pdf2 {
                                        operator[eko_q2_index][pid_high2_idx][eko_x2_high_idx]
                                            [pid_low2_idx][eko_x2_low_idx]
                                    } else {
                                        1.
                                    };
                                    assert_eq!(op2, 1.);
                                    // let mut value = alphas[q2_index]
                                    let value = alphas[eko_q2_index]
                                        .powi(order.alphas.try_into().unwrap())
                                        * op1
                                        * op2;

                                    // TODO: implement scale variations

                                    //if order.logxir > 0 {
                                    //    value *= (xir * xir)
                                    //        .ln()
                                    //        .powi(order.logxir.try_into().unwrap());
                                    //}

                                    //if order.logxif > 0 {
                                    //    value *= (xif * xif)
                                    //        .ln()
                                    //        .powi(order.logxif.try_into().unwrap());
                                    //}

                                    value
                                }),
                            );

                            let saved = factor * convoluted;
                            if saved != 0. {
                                array[[0, x1_low, x2_low]] += saved;
                            }
                        }
                    }
                }
            }

            *subgrid = ReadOnlySparseSubgridV1::new(
                array,
                q2low_grid.clone(),
                x1low_grid.clone(),
                x2low_grid.clone(),
            )
            .into();
        }

        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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

        matches!(result, Err(UnknownSubgrid(x)) if x == subgrid_type);
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
}
