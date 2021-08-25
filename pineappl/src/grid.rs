//! Module containing all traits and supporting structures for grids.

use super::bin::{BinInfo, BinLimits, BinRemapper};
use super::empty_subgrid::EmptySubgridV1;
use super::import_only_subgrid::ImportOnlySubgridV1;
use super::lagrange_subgrid::{LagrangeSparseSubgridV1, LagrangeSubgridV1, LagrangeSubgridV2};
use super::lumi::LumiEntry;
use super::ntuple_subgrid::NtupleSubgridV1;
use super::subgrid::{ExtraSubgridParams, Subgrid, SubgridEnum, SubgridParams};
use float_cmp::approx_eq;
use git_version::git_version;
use itertools::Itertools;
use lz_fear::{framed::DecompressionError::WrongMagic, LZ4FrameReader};
use ndarray::{s, Array3, Dimension};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::convert::TryInto;
use std::io::{Read, Seek, SeekFrom, Write};
use std::mem;
use std::ops::Range;
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
    /// Returned when trying to merge two `Grid` objects with incompatible bin limits.
    #[error(transparent)]
    DifferentBins(super::bin::MergeBinError),
}

/// Error returned when trying to construct a `Grid` using an unknown subgrid type.
#[derive(Debug, Error)]
#[error("tried constructing a Grid with unknown Subgrid type `{0}`")]
pub struct UnknownSubgrid(String);

#[derive(Deserialize, Serialize)]
struct Mmv1 {}

#[derive(Deserialize, Serialize)]
struct Mmv2 {
    remapper: Option<BinRemapper>,
    key_value_db: HashMap<String, String>,
}

#[derive(Deserialize, Serialize)]
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

#[derive(Deserialize, Serialize)]
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
    ) -> Result<Self, UnknownSubgrid> {
        let subgrid_template: SubgridEnum = match subgrid_type {
            "LagrangeSubgrid" | "LagrangeSubgridV2" => {
                LagrangeSubgridV2::new(&subgrid_params, &extra).into()
            }
            "LagrangeSubgridV1" => LagrangeSubgridV1::new(&subgrid_params).into(),
            "NtupleSubgrid" => NtupleSubgridV1::new().into(),
            "LagrangeSparseSubgrid" => LagrangeSparseSubgridV1::new(&subgrid_params).into(),
            _ => return Err(UnknownSubgrid(subgrid_type.to_string())),
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

    /// Performs a convolution of the contained subgrids with the given PDFs, `xfx1` for the first
    /// parton and `xfx2` for the second parton, `alphas` for the evaluation of the strong
    /// coupling. The parameters `order_mask` and `lumi_mask` can be used to selectively enable
    /// perturbative orders and luminosities; they must either be empty (everything enabled) or as
    /// large as the orders and luminosity function, respectively. If the corresponding entry is
    /// `true` the order/luminosity is enable, `false` disables the entry. The tuple `xi` can be
    /// used to independently vary the renormalization (first element) and factorization scale
    /// (second element) from their central value `(1.0, 1.0)`.
    ///
    /// # Panics
    ///
    /// TODO
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

        let (mut q2_grid, mut x1_grid, mut x2_grid) = self
            .subgrids
            .iter()
            .find(|subgrid| !subgrid.is_empty())
            .map_or_else(
                || (Cow::default(), Cow::default(), Cow::default()),
                |grid| (grid.q2_grid(), grid.x1_grid(), grid.x2_grid()),
            );
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

                    subgrid.convolute(&x1_grid, &x2_grid, &q2_grid, &|ix1, ix2, iq2| {
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
                    })
                } else {
                    todo!();
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
        xfx1: &dyn Fn(i32, f64, f64) -> f64,
        xfx2: &dyn Fn(i32, f64, f64) -> f64,
        alphas: &dyn Fn(f64) -> f64,
        order: usize,
        bin: usize,
        lumi: usize,
        xir: f64,
        xif: f64,
    ) -> Array3<f64> {
        let normalization = self.bin_info().normalizations()[bin];

        let pdf_cache1 = RefCell::new(FxHashMap::default());
        let pdf_cache2 = RefCell::new(FxHashMap::default());
        let alphas_cache = RefCell::new(FxHashMap::default());

        let subgrid = &self.subgrids[[order, bin, lumi]];
        let order = &self.orders[order];

        let mut array = if subgrid.is_empty() {
            Array3::zeros((0, 0, 0))
        } else {
            let q2_grid = subgrid.q2_grid();
            let x1_grid = subgrid.x1_grid();
            let x2_grid = subgrid.x2_grid();

            let use_cache = !q2_grid.is_empty() && !x1_grid.is_empty() && !x2_grid.is_empty();
            let two_caches = !ptr::eq(&xfx1, &xfx2);

            let lumi_entry = &self.lumi[lumi];

            if use_cache {
                let mut array = Array3::zeros((q2_grid.len(), x1_grid.len(), x2_grid.len()));

                for ((iq2, ix1, ix2), value) in subgrid.iter() {
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
                        .entry(iq2)
                        .or_insert_with(|| alphas(xir * xir * q2));

                    lumi *= alphas.powi(order.alphas.try_into().unwrap());
                    array[[iq2, ix1, ix2]] = lumi * value;
                }

                array
            } else {
                todo!();
            }
        };

        if order.logxir > 0 {
            array *= (xir * xir).ln().powi(order.logxir.try_into().unwrap());
        }

        if order.logxif > 0 {
            array *= (xif * xif).ln().powi(order.logxif.try_into().unwrap());
        }

        array /= normalization;
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

    /// Merges the bins for the corresponding range together in a single one.
    ///
    /// # Errors
    ///
    /// TODO
    pub fn merge_bins(&mut self, bins: Range<usize>) -> Result<(), ()> {
        if (bins.start >= self.bin_limits.bins()) || (bins.end > self.bin_limits.bins()) {
            return Err(());
        }

        self.bin_limits.merge_bins(bins.clone());

        match &mut self.more_members {
            MoreMembers::V1(_) => {}
            MoreMembers::V2(mmv2) => {
                if let Some(remapper) = &mut mmv2.remapper {
                    remapper.merge_bins(bins.clone())?;
                }
            }
            MoreMembers::V3(mmv3) => {
                if let Some(remapper) = &mut mmv3.remapper {
                    remapper.merge_bins(bins.clone())?;
                }
            }
        }

        let mut old_subgrids = mem::replace(
            &mut self.subgrids,
            Array3::from_shape_simple_fn(
                (self.orders.len(), self.bin_limits.bins(), self.lumi.len()),
                || EmptySubgridV1::default().into(),
            ),
        );

        for ((order, bin, lumi), subgrid) in old_subgrids.indexed_iter_mut() {
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
    pub fn merge(&mut self, mut other: Self) -> Result<(), GridMergeError> {
        let mut new_orders: Vec<Order> = Vec::new();
        let mut new_bins = 0;
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

        if self.bin_limits != other.bin_limits {
            if let Err(e) = self.bin_limits.merge(&other.bin_limits) {
                return Err(GridMergeError::DifferentBins(e));
            }

            new_bins = other.bin_limits.bins();

            // TODO: figure out a better strategy than removing the remapper
            match &mut self.more_members {
                MoreMembers::V1(_) => {}
                MoreMembers::V2(mmv2) => {
                    mmv2.remapper = None;
                }
                MoreMembers::V3(mmv3) => {
                    mmv3.remapper = None;
                }
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
    ///
    /// # Panics
    ///
    /// TODO
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
            MoreMembers::V3(mmv3) => mmv3.remapper = Some(remapper),
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
                MoreMembers::V3(mmv3) => mmv3.remapper.as_ref(),
            },
        )
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
            self.optimize_lumi();
        }

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
                        let mut new_subgrid = ImportOnlySubgridV1::from(&*grid).into();
                        mem::swap(subgrid, &mut new_subgrid);
                    }
                    SubgridEnum::EmptySubgridV1(_)
                    | SubgridEnum::LagrangeSparseSubgridV1(_)
                    | SubgridEnum::ImportOnlySubgridV1(_) => {
                        // nothing to optimize here
                    }
                    SubgridEnum::NtupleSubgridV1(_) => todo!(),
                }
            }
        }
    }

    fn optimize_lumi(&mut self) {
        let mut keep_lumi_indices = vec![];
        let mut new_lumi_entries = vec![];

        for (lumi, entry) in self.lumi.iter().enumerate() {
            if !self
                .subgrids
                .slice(s![.., .., lumi])
                .iter()
                .all(|subgrid| subgrid.is_empty())
            {
                keep_lumi_indices.push(lumi);
                new_lumi_entries.push(entry.clone());
            }
        }

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

    // TODO: simplify the method, because `optimize_lumi` already removes empty entries
    fn symmetrize_lumi(&mut self) {
        let mut indices: Vec<usize> = (0..self.lumi.len()).rev().collect();
        let mut pairs: Vec<(usize, usize)> = Vec::new();
        let mut not_symmetrized: Vec<usize> = Vec::new();

        'looop: while let Some(index) = indices.pop() {
            let lumi_entry = &self.lumi[index];

            if *lumi_entry == lumi_entry.transpose() {
                for order in 0..self.orders.len() {
                    for bin in 0..self.bin_limits.bins() {
                        let subgrid = &self.subgrids[[order, bin, index]];

                        // check if in all cases the limits are compatible with merging
                        if !subgrid.is_empty() && (subgrid.x1_grid() != subgrid.x2_grid()) {
                            not_symmetrized.push(index);

                            continue 'looop;
                        }
                    }
                }

                pairs.push((index, index));
            } else if let Some((j, &other_index)) = indices
                .iter()
                .enumerate()
                .find(|(_, i)| self.lumi[**i] == lumi_entry.transpose())
            {
                indices.remove(j);

                for order in 0..self.orders.len() {
                    for bin in 0..self.bin_limits.bins() {
                        let lhs = &self.subgrids[[order, bin, index]];
                        let rhs = &self.subgrids[[order, bin, other_index]];

                        // check if in all cases the limits are compatible with merging
                        if !lhs.is_empty()
                            && !rhs.is_empty()
                            && ((lhs.x1_grid() != rhs.x2_grid())
                                || (lhs.x2_grid() != rhs.x1_grid()))
                        {
                            not_symmetrized.push(index);
                            not_symmetrized.push(other_index);

                            continue 'looop;
                        }
                    }
                }

                pairs.push((index, other_index));
            } else {
                not_symmetrized.push(index);
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
                        let mut lhs = mem::replace(&mut subgrids[[i, j, k1]], None).unwrap();
                        let mut rhs = mem::replace(&mut subgrids[[i, j, k2]], None).unwrap();

                        subgrids[[i, j, k1]] = Some(if rhs.is_empty() {
                            lhs
                        } else if lhs.is_empty() {
                            let mut new_lhs = rhs.clone_empty();
                            new_lhs.merge(&mut rhs, true);
                            new_lhs
                        } else {
                            lhs.merge(&mut rhs, true);
                            lhs
                        });
                    }
                }
            }
        }

        self.subgrids = Array3::from_shape_vec(
            (i_size, j_size, pairs.len() + not_symmetrized.len()),
            subgrids.into_raw_vec().into_iter().flatten().collect(),
        )
        .unwrap();

        let mut new_lumi_indices: Vec<_> = pairs
            .iter()
            .map(|(index, _)| index)
            .chain(not_symmetrized.iter())
            .copied()
            .collect();
        new_lumi_indices.sort_unstable();

        self.lumi = new_lumi_indices
            .iter()
            .map(|i| self.lumi[*i].clone())
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
}
