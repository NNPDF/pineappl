//! Module containing structures for the 3 dimensions of a [`Grid`]: [`Bin`], [`Order`] and
//! channels (`boc`).
//!
//! [`Grid`]: super::grid::Grid

// use super::grid::GridError;
use super::convert;
use super::error::{Error, Result};
use float_cmp::approx_eq;
use itertools::{izip, Itertools};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::cmp::Ordering;
use std::ops::Range;
use std::str::FromStr;

/// TODO
#[repr(C)]
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub enum Kinematics {
    /// TODO
    Scale(usize),
    /// TODO
    X(usize),
}

/// TODO
#[repr(C)]
#[derive(Clone, Deserialize, Eq, PartialEq, Serialize)]
pub enum ScaleFuncForm {
    /// TODO
    NoScale,
    /// TODO
    Scale(usize),
    /// TODO
    QuadraticSum(usize, usize),
    /// TODO
    QuadraticMean(usize, usize),
    /// TODO
    QuadraticSumOver4(usize, usize),
    /// TODO
    LinearMean(usize, usize),
    /// TODO
    LinearSum(usize, usize),
    /// TODO
    ScaleMax(usize, usize),
    /// TODO
    ScaleMin(usize, usize),
    /// TODO
    Prod(usize, usize),
    /// TODO
    S2plusS1half(usize, usize),
    /// TODO
    Pow4Sum(usize, usize),
    /// TODO
    WgtAvg(usize, usize),
    /// TODO
    S2plusS1fourth(usize, usize),
    /// TODO
    ExpProd2(usize, usize),
}

impl ScaleFuncForm {
    /// TODO
    #[must_use]
    pub fn calc<'a>(
        &self,
        node_values: &'a [Vec<f64>],
        kinematics: &[Kinematics],
    ) -> Cow<'a, [f64]> {
        match self.clone() {
            Self::NoScale => Cow::Borrowed(&[]),
            Self::Scale(index) => {
                if node_values.is_empty() {
                    // TODO: empty subgrid should have as many node values as dimensions
                    Cow::Borrowed(&[])
                } else {
                    Cow::Borrowed(
                        &node_values[kinematics
                            .iter()
                            .position(|&kin| kin == Kinematics::Scale(index))
                            // UNWRAP: this should be guaranteed by `Grid::new`
                            .unwrap_or_else(|| unreachable!())],
                    )
                }
            }
            Self::QuadraticSum(idx1, idx2)
            | Self::QuadraticMean(idx1, idx2)
            | Self::QuadraticSumOver4(idx1, idx2)
            | Self::LinearMean(idx1, idx2)
            | Self::LinearSum(idx1, idx2)
            | Self::ScaleMax(idx1, idx2)
            | Self::ScaleMin(idx1, idx2)
            | Self::Prod(idx1, idx2)
            | Self::S2plusS1half(idx1, idx2)
            | Self::Pow4Sum(idx1, idx2)
            | Self::WgtAvg(idx1, idx2)
            | Self::S2plusS1fourth(idx1, idx2)
            | Self::ExpProd2(idx1, idx2) => {
                let calc_scale: fn((f64, f64)) -> f64 = match self.clone() {
                    Self::QuadraticSum(_, _) => |(s1, s2)| s1 + s2,
                    Self::QuadraticMean(_, _) => |(s1, s2)| 0.5 * (s1 + s2),
                    Self::QuadraticSumOver4(_, _) => |(s1, s2)| 0.25 * (s1 + s2),
                    Self::LinearMean(_, _) => |(s1, s2)| 0.25 * (s1.sqrt() + s2.sqrt()).powi(2),
                    Self::LinearSum(_, _) => |(s1, s2)| (s1.sqrt() + s2.sqrt()).powi(2),
                    Self::ScaleMax(_, _) => |(s1, s2)| s1.max(s2),
                    Self::ScaleMin(_, _) => |(s1, s2)| s1.min(s2),
                    Self::Prod(_, _) => |(s1, s2)| s1 * s2,
                    Self::S2plusS1half(_, _) => |(s1, s2)| 0.5 * s2.mul_add(2.0, s1),
                    Self::Pow4Sum(_, _) => |(s1, s2)| s1.hypot(s2),
                    Self::WgtAvg(_, _) => |(s1, s2)| s1.mul_add(s1, s2 * s2) / (s1 + s2),
                    Self::S2plusS1fourth(_, _) => |(s1, s2)| s1.mul_add(0.25, s2),
                    Self::ExpProd2(_, _) => {
                        |(s1, s2)| (s1.sqrt() * (0.3 * s2.sqrt()).exp()).powi(2)
                    }
                    _ => unreachable!(),
                };

                let scales1 = &node_values[kinematics
                    .iter()
                    .position(|&kin| kin == Kinematics::Scale(idx1))
                    // UNWRAP: this should be guaranteed by `Grid::new`
                    .unwrap_or_else(|| unreachable!())];
                let scales2 = &node_values[kinematics
                    .iter()
                    .position(|&kin| kin == Kinematics::Scale(idx2))
                    // UNWRAP: this should be guaranteed by `Grid::new`
                    .unwrap_or_else(|| unreachable!())];

                Cow::Owned(
                    scales1
                        .iter()
                        .copied()
                        .cartesian_product(scales2.iter().copied())
                        .map(calc_scale)
                        .collect(),
                )
            }
        }
    }

    /// TODO
    #[must_use]
    pub fn idx(&self, indices: &[usize], scale_dims: &[usize]) -> usize {
        match self.clone() {
            Self::NoScale => unreachable!(),
            Self::Scale(index) => indices[index],
            Self::QuadraticSum(idx1, idx2)
            | Self::QuadraticMean(idx1, idx2)
            | Self::QuadraticSumOver4(idx1, idx2)
            | Self::LinearMean(idx1, idx2)
            | Self::LinearSum(idx1, idx2)
            | Self::ScaleMax(idx1, idx2)
            | Self::ScaleMin(idx1, idx2)
            | Self::Prod(idx1, idx2)
            | Self::S2plusS1half(idx1, idx2)
            | Self::Pow4Sum(idx1, idx2)
            | Self::WgtAvg(idx1, idx2)
            | Self::S2plusS1fourth(idx1, idx2)
            | Self::ExpProd2(idx1, idx2) => indices[idx1] * scale_dims[1] + indices[idx2],
        }
    }
}

/// TODO
#[derive(Clone, Deserialize, Eq, PartialEq, Serialize)]
pub struct Scales {
    /// TODO
    pub ren: ScaleFuncForm,
    /// TODO
    pub fac: ScaleFuncForm,
    /// TODO
    pub frg: ScaleFuncForm,
}

impl<'a> From<&'a Scales> for [&'a ScaleFuncForm; 3] {
    fn from(scales: &'a Scales) -> [&'a ScaleFuncForm; 3] {
        [&scales.ren, &scales.fac, &scales.frg]
    }
}

impl Scales {
    /// TODO
    pub fn compatible_with(&self, kinematics: &[Kinematics]) -> bool {
        for scale in [&self.ren, &self.fac, &self.frg].map(Clone::clone) {
            match scale {
                ScaleFuncForm::NoScale => {}
                ScaleFuncForm::Scale(index)
                    if kinematics
                        .iter()
                        .any(|&kin| kin == Kinematics::Scale(index)) => {}
                ScaleFuncForm::QuadraticSum(idx1, idx2)
                | ScaleFuncForm::QuadraticMean(idx1, idx2)
                | ScaleFuncForm::QuadraticSumOver4(idx1, idx2)
                | ScaleFuncForm::LinearMean(idx1, idx2)
                | ScaleFuncForm::LinearSum(idx1, idx2)
                | ScaleFuncForm::ScaleMax(idx1, idx2)
                | ScaleFuncForm::ScaleMin(idx1, idx2)
                | ScaleFuncForm::Prod(idx1, idx2)
                | ScaleFuncForm::S2plusS1half(idx1, idx2)
                | ScaleFuncForm::Pow4Sum(idx1, idx2)
                | ScaleFuncForm::WgtAvg(idx1, idx2)
                | ScaleFuncForm::S2plusS1fourth(idx1, idx2)
                | ScaleFuncForm::ExpProd2(idx1, idx2)
                    if kinematics.iter().any(|&kin| kin == Kinematics::Scale(idx1))
                        && kinematics.iter().any(|&kin| kin == Kinematics::Scale(idx2)) => {}
                _ => return false,
            }
        }

        true
    }
}

/// TODO
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Bin {
    limits: Vec<(f64, f64)>,
    normalization: f64,
}

impl Bin {
    /// TODO
    ///
    /// # Panics
    ///
    /// TODO
    #[must_use]
    pub fn new(limits: Vec<(f64, f64)>, normalization: f64) -> Self {
        for limits in &limits {
            assert!(limits.1 >= limits.0);
        }

        Self {
            limits,
            normalization,
        }
    }

    /// TODO
    #[must_use]
    pub fn dimensions(&self) -> usize {
        self.limits.len()
    }

    /// TODO
    #[must_use]
    pub const fn normalization(&self) -> f64 {
        self.normalization
    }

    /// TODO
    #[must_use]
    pub fn limits(&self) -> &[(f64, f64)] {
        &self.limits
    }

    /// TODO
    #[must_use]
    pub fn partial_eq_with_ulps(&self, other: &Self, ulps: i64) -> bool {
        self.limits.iter().zip(other.limits()).all(|(&lhs, &rhs)| {
            approx_eq!(f64, lhs.0, rhs.0, ulps = ulps) && approx_eq!(f64, lhs.1, rhs.1, ulps = ulps)
        }) && approx_eq!(f64, self.normalization, other.normalization, ulps = ulps)
    }
}

/// TODO
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct BinsWithFillLimits {
    bins: Vec<Bin>,
    fill_limits: Vec<f64>,
}

impl BinsWithFillLimits {
    /// TODO
    ///
    /// # Errors
    ///
    /// TODO
    pub fn new(bins: Vec<Bin>, fill_limits: Vec<f64>) -> Result<Self> {
        // TODO: validate the bins

        // - there must be at least one bin
        // - all fill limits must be ascending
        // - all dimensions must be the same
        // - limits must not overlap

        let fill_limits_len = fill_limits.len();
        let bins_len_p1 = bins.len() + 1;

        if fill_limits_len != bins_len_p1 {
            return Err(Error::General(
                "number of bins must agree with the number of fill limits plus 1".to_owned(),
            ));
        }

        Ok(Self { bins, fill_limits })
    }

    /// TODO
    ///
    /// # Errors
    ///
    /// TODO
    pub fn from_fill_limits(fill_limits: Vec<f64>) -> Result<Self> {
        let bins = fill_limits
            .windows(2)
            .map(|win| Bin::new(vec![(win[0], win[1])], win[1] - win[0]))
            .collect();

        Self::new(bins, fill_limits)
    }

    /// TODO
    ///
    /// # Errors
    ///
    /// TODO
    ///
    /// # Panics
    ///
    /// TODO
    pub fn from_limits_and_normalizations(
        limits: Vec<Vec<(f64, f64)>>,
        normalizations: Vec<f64>,
    ) -> Result<Self> {
        let limits_len = limits.len();
        let normalizations_len = normalizations.len();

        if limits_len != normalizations_len {
            return Err(Error::General(
                "number of limits be the same as the number of normalizations".to_owned(),
            ));
        }

        let fill_limits = (0..=limits.len()).map(convert::f64_from_usize).collect();
        let bins = limits
            .into_iter()
            .zip(normalizations)
            .map(|(limits, normalization)| Bin::new(limits, normalization))
            .collect();

        Self::new(bins, fill_limits)
    }

    /// TODO
    pub fn slices(&self) -> Vec<Range<usize>> {
        if self.dimensions() == 1 {
            // TODO: check that bins are contiguous
            vec![0..self.len()]
        } else {
            self.bins()
                .iter()
                .flat_map(Bin::limits)
                .enumerate()
                .filter_map(|(index, x)| {
                    ((index % self.dimensions()) != (self.dimensions() - 1)).then_some(x)
                })
                .collect::<Vec<_>>()
                .chunks_exact(self.dimensions() - 1)
                .enumerate()
                .dedup_by_with_count(|(_, x), (_, y)| x == y)
                .map(|(count, (index, _))| index..index + count)
                .collect()
        }
    }

    /// TODO
    #[must_use]
    pub fn bins(&self) -> &[Bin] {
        &self.bins
    }

    /// TODO
    #[must_use]
    pub fn len(&self) -> usize {
        self.bins.len()
    }

    /// TODO
    #[must_use]
    pub fn dimensions(&self) -> usize {
        self.bins
            .first()
            // UNWRAP: `Bin::new` should guarantee that there's at least one bin
            .unwrap_or_else(|| unreachable!())
            .dimensions()
    }

    /// TODO
    #[must_use]
    pub fn fill_index(&self, value: f64) -> Option<usize> {
        match self
            .fill_limits
            .binary_search_by(|left| left.total_cmp(&value))
        {
            Err(0) => None,
            Err(index) if index == self.fill_limits.len() => None,
            Ok(index) if index == (self.fill_limits.len() - 1) => None,
            Ok(index) => Some(index),
            Err(index) => Some(index - 1),
        }
    }

    /// TODO
    #[must_use]
    pub fn fill_limits(&self) -> &[f64] {
        &self.fill_limits
    }

    /// TODO
    pub fn normalizations(&self) -> Vec<f64> {
        self.bins.iter().map(Bin::normalization).collect()
    }

    /// TODO
    ///
    /// # Errors
    ///
    /// TODO
    // TODO: change range to `RangeBounds<usize>`
    pub fn merge(&self, range: Range<usize>) -> Result<Self> {
        // TODO: allow more flexible merging
        if !self
            .slices()
            .iter()
            .any(|&Range { start, end }| (start <= range.start) && (range.end <= end))
        {
            // TODO: implement proper error handling
            return Err(Error::General("bins are not simply connected".to_string()));
        }

        let mut limits: Vec<_> = self.bins.iter().map(|bin| bin.limits().to_vec()).collect();
        let mut normalizations = self.normalizations();
        let mut fill_limits = self.fill_limits.clone();
        let dim = self.dimensions();

        limits[range.start][dim - 1].1 = limits[range.end - 1][dim - 1].1;
        normalizations[range.start] += normalizations[range.start + 1..range.end]
            .iter()
            .sum::<f64>();

        limits.drain(range.start + 1..range.end);
        normalizations.drain(range.start + 1..range.end);
        fill_limits.drain(range.start + 1..range.end);

        Self::new(
            limits
                .into_iter()
                .zip(normalizations)
                .map(|(limits, normalization)| Bin::new(limits, normalization))
                .collect(),
            fill_limits,
        )
    }

    /// TODO
    ///
    /// # Panics
    ///
    /// TODO
    pub fn remove(&mut self, index: usize) -> Bin {
        assert!(self.len() > 1);

        self.fill_limits.pop().unwrap();
        self.bins.remove(index)
    }

    /// TODO
    #[must_use]
    pub fn bins_partial_eq_with_ulps(&self, other: &Self, ulps: i64) -> bool {
        self.bins
            .iter()
            .zip(other.bins())
            .all(|(lhs, rhs)| lhs.partial_eq_with_ulps(rhs, ulps))
    }
}

impl FromStr for BinsWithFillLimits {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        let remaps: Result<Vec<Vec<Vec<f64>>>> = s
            .split(';')
            .map(|string| {
                string
                    .split('|')
                    .map(|string| {
                        string
                            .split_once(':')
                            .map_or(Ok(string), |(lhs, rhs)| {
                                match (lhs.trim().parse::<usize>(), rhs.trim().parse::<usize>()) {
                                    (Err(lhs), Err(rhs)) => Err(Error::General(format!(
                                        "unable to parse 'N:M' syntax from: '{string}' (N: '{lhs}', M: '{rhs}')"
                                    ))),
                                    // skip :N specification
                                    (Err(_), Ok(_)) => Ok(lhs),
                                    // skip N: specification
                                    (Ok(_), Err(_)) => Ok(rhs),
                                    // skip N:M specification
                                    (Ok(_), Ok(_)) => Ok(""),
                                }
                            })?
                            .split(',')
                            .filter_map(|string| {
                                let string = string.trim();
                                if string.is_empty() {
                                    None
                                } else {
                                    Some(string.parse::<f64>().map_err(|err| {
                                        Error::General(format!(
                                            "unable to parse limit '{string}': '{err}')"
                                        ))
                                    }))
                                }
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect();
        let mut remaps = remaps?;

        if let Some(first) = remaps.first() {
            if first.len() != 1 {
                return Err(Error::General(
                    "'|' syntax not meaningful for first dimension".to_owned(),
                ));
            }
        }

        // go over `remaps` again, and repeat previous entries as requested with the `|` syntax
        for vec in &mut remaps {
            for i in 1..vec.len() {
                if vec[i].is_empty() {
                    if vec[i - 1].is_empty() {
                        return Err(Error::General("empty repetition with '|'".to_owned()));
                    }

                    vec[i] = vec[i - 1].clone();
                }
            }
        }

        // go over `remaps` again, this time remove bin as requested with the `:N` or `N:` syntax
        for (vec, string) in remaps.iter_mut().zip(s.split(';')) {
            for (vec, string) in vec.iter_mut().zip(string.split('|')) {
                let (lhs, rhs) = {
                    if let Some((lhs, rhs)) = string.split_once(':') {
                        (lhs.parse::<usize>(), rhs.parse::<usize>())
                    } else {
                        // there's no colon
                        continue;
                    }
                };

                if let Ok(num) = rhs {
                    vec.truncate(vec.len() - num);
                }

                if let Ok(num) = lhs {
                    vec.drain(0..num);
                }

                if vec.len() <= 1 {
                    return Err(Error::General("no limits due to ':' syntax".to_owned()));
                }
            }
        }

        let dimensions = remaps.len();
        let mut normalizations = Vec::new();
        let mut limits = Vec::new();
        let mut buffer = Vec::with_capacity(dimensions);
        let mut pipe_indices = vec![0; dimensions];
        let mut last_indices = vec![0; dimensions];

        'looop: for indices in remaps
            .iter()
            .map(|vec| 0..vec.iter().map(|vec| vec.len() - 1).max().unwrap())
            .multi_cartesian_product()
        {
            // calculate `pipe_indices`, which stores the indices for the second dimension of `remaps`
            for d in 0..dimensions - 1 {
                if indices[d] > last_indices[d] {
                    for dp in d + 1..dimensions {
                        if remaps[dp].len() != 1 {
                            pipe_indices[dp] += 1;
                        }
                    }
                }
            }

            last_indices.clone_from(&indices);

            let mut normalization = 1.0;

            for (remap, &pipe_index, &i) in izip!(&remaps, &pipe_indices, &indices) {
                if let Some(r) = remap.get(pipe_index) {
                    if r.len() <= (i + 1) {
                        buffer.clear();

                        // this index doesn't exist
                        continue 'looop;
                    }

                    let left = r[i];
                    let right = r[i + 1];

                    buffer.push((left, right));
                    normalization *= right - left;
                } else {
                    return Err(Error::General(
                        "missing '|' specification: number of variants too small".to_owned(),
                    ));
                }
            }

            // TODO: rewrite the code to avoid `buffer.clear()`
            limits.push(buffer.clone());
            buffer.clear();
            normalizations.push(normalization);
        }

        Ok(Self::from_limits_and_normalizations(limits, normalizations)
            // TODO: implement proper error handling
            .unwrap())
    }
}

/// Coupling powers for each grid.
#[derive(Clone, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
pub struct Order {
    /// Exponent of the strong coupling.
    pub alphas: u8,
    /// Exponent of the electromagnetic coupling.
    pub alpha: u8,
    /// Exponent of the logarithm of the scale factor of the renomalization scale.
    pub logxir: u8,
    /// Exponent of the logarithm of the scale factor of the initial state factorization scale.
    pub logxif: u8,
    /// Exponent of the logarithm of the scale factor of the final state factorization scale
    /// (fragmentation scale).
    pub logxia: u8,
    // /// Reserved for future usage.
    // pub other: [u8; 3],
}

impl FromStr for Order {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        let mut result = Self {
            alphas: 0,
            alpha: 0,
            logxir: 0,
            logxif: 0,
            logxia: 0,
        };

        for tuple in s
            .split(|c: char| c.is_ascii_digit())
            .filter(|s| !s.is_empty())
            .zip(
                s.split(|c: char| !c.is_ascii_digit())
                    .filter(|s| !s.is_empty())
                    .map(str::parse),
            )
        {
            match tuple {
                ("as", Ok(num)) => {
                    result.alphas = num;
                }
                ("a", Ok(num)) => {
                    result.alpha = num;
                }
                ("lr", Ok(num)) => {
                    result.logxir = num;
                }
                ("lf", Ok(num)) => {
                    result.logxif = num;
                }
                ("la", Ok(num)) => {
                    result.logxia = num;
                }
                (label, Err(err)) => {
                    return Err(Error::General(format!(
                        "error while parsing exponent of '{label}': {err}"
                    )));
                }
                (label, Ok(_)) => {
                    return Err(Error::General(format!("unknown coupling: '{label}'")));
                }
            }
        }

        Ok(result)
    }
}

impl Ord for Order {
    fn cmp(&self, other: &Self) -> Ordering {
        // sort leading orders before next-to-leading orders, then the lowest power in alpha, the
        // rest lexicographically
        (self.alphas + self.alpha)
            .cmp(&(other.alphas + other.alpha))
            .then((self.alpha, self.logxir, self.logxif, self.logxia).cmp(&(
                other.alpha,
                other.logxir,
                other.logxif,
                other.logxia,
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
    pub const fn new(alphas: u8, alpha: u8, logxir: u8, logxif: u8, logxia: u8) -> Self {
        Self {
            alphas,
            alpha,
            logxir,
            logxif,
            logxia,
        }
    }

    /// Return a mask suitable to pass as the `order_mask` parameter of [`Grid::convolve`],
    /// [`Grid::evolve`] or [`Grid::evolve_info`]. The selection of `orders` is controlled using
    /// the `max_as` and `max_al` parameters, for instance setting `max_as = 1` and `max_al = 0`
    /// selects the LO QCD only, `max_as = 2` and `max_al = 0` the NLO QCD; setting `max_as = 3`
    /// and `max_al = 2` would select all NLOs, and the NNLO QCD.
    ///
    /// [`Grid::convolve`]: super::grid::Grid::convolve
    /// [`Grid::evolve`]: super::grid::Grid::evolve
    /// [`Grid::evolve_info`]: super::grid::Grid::evolve_info
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
    /// use pineappl::boc::Order;
    ///
    /// let orders = [
    ///     Order::new(0, 2, 0, 0, 0), //   LO        :          alpha^2
    ///     Order::new(1, 2, 0, 0, 0), //  NLO QCD    : alphas   alpha^2
    ///     Order::new(0, 3, 0, 0, 0), //  NLO  EW    :          alpha^3
    ///     Order::new(2, 2, 0, 0, 0), // NNLO QCD    : alphas^2 alpha^2
    ///     Order::new(1, 3, 0, 0, 0), // NNLO QCD—EW : alphas   alpha^3
    ///     Order::new(0, 4, 0, 0, 0), // NNLO EW     :          alpha^4
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
    /// use pineappl::boc::Order;
    ///
    /// let orders = [
    ///     Order::new(0, 2, 0, 0, 0), //  LO         :        alpha^2
    ///     Order::new(1, 2, 0, 0, 0), //  NLO QCD    : alphas alpha^2
    ///     Order::new(1, 2, 1, 0, 0), //  NLO QCD    : alphas alpha^2 logxif
    ///     Order::new(0, 3, 0, 0, 0), //  NLO  EW    :        alpha^3
    ///     Order::new(0, 3, 1, 0, 0), //  NLO  EW    :        alpha^3 logxif
    /// ];
    ///
    /// assert_eq!(Order::create_mask(&orders, 0, 2, true), [true, false, false, true, true]);
    /// ```
    ///
    /// For the more complicated example of top-pair production one can see the difference between
    /// the selection for different LOs:
    ///
    /// ```rust
    /// use pineappl::boc::Order;
    ///
    /// let orders = [
    ///     Order::new(2, 0, 0, 0, 0), //   LO QCD    : alphas^2
    ///     Order::new(1, 1, 0, 0, 0), //   LO QCD—EW : alphas   alpha
    ///     Order::new(0, 2, 0, 0, 0), //   LO  EW    :          alpha^2
    ///     Order::new(3, 0, 0, 0, 0), //  NLO QCD    : alphas^3
    ///     Order::new(2, 1, 0, 0, 0), //  NLO QCD—EW : alphas^2 alpha
    ///     Order::new(1, 2, 0, 0, 0), //  NLO QCD—EW : alphas   alpha^2
    ///     Order::new(0, 3, 0, 0, 0), //  NLO EW     :          alpha^3
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
    pub fn create_mask(orders: &[Self], max_as: u8, max_al: u8, logs: bool) -> Vec<bool> {
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
                     logxia,
                 }| {
                    if !logs && (logxir > 0 || logxif > 0 || logxia > 0) {
                        return false;
                    }

                    let pto = alphas + alpha - lo;

                    alphas + alpha < min + lo
                        || (alphas + alpha < max + lo
                            && match max_as.cmp(&max_al) {
                                Ordering::Greater => lo_as + pto == alphas,
                                Ordering::Less => lo_al + pto == alpha,
                                // TODO: when do we hit this condition?
                                Ordering::Equal => false,
                            })
                },
            )
            .collect()
    }
}

/// This structure represents a channel. Each channel consists of a tuple containing in the
/// following order, the particle ID of the first incoming parton, then the particle ID of the
/// second parton, and finally a numerical factor that will multiply the result for this specific
/// combination.
#[derive(Clone, Debug, Deserialize, PartialEq, PartialOrd, Serialize)]
pub struct Channel {
    entry: Vec<(Vec<i32>, f64)>,
}

impl Channel {
    /// Constructor for `Channel`. Note that `entry` must be non-empty, otherwise this function
    /// panics.
    ///
    /// # Examples
    ///
    /// Ordering of the arguments doesn't matter:
    ///
    /// ```rust
    /// use pineappl::boc::Channel;
    ///
    /// let entry1 = Channel::new(vec![(vec![2, 2], 1.0), (vec![4, 4], 1.0)]);
    /// let entry2 = Channel::new(vec![(vec![4, 4], 1.0), (vec![2, 2], 1.0)]);
    ///
    /// // checks that the ordering doesn't matter
    /// assert_eq!(entry1, entry2);
    /// ```
    ///
    /// Same arguments are merged together:
    ///
    /// ```rust
    /// use pineappl::boc::Channel;
    ///
    /// let entry1 = Channel::new(vec![(vec![1, 1], 1.0), (vec![1, 1], 3.0), (vec![3, 3], 1.0), (vec![1, 1], 6.0)]);
    /// let entry2 = Channel::new(vec![(vec![1, 1], 10.0), (vec![3, 3], 1.0)]);
    ///
    /// assert_eq!(entry1, entry2);
    /// ```
    ///
    /// # Panics
    ///
    /// Creating an empty channel panics:
    ///
    /// ```rust,should_panic
    /// use pineappl::boc::Channel;
    ///
    /// let _ = Channel::new(vec![]);
    /// ```
    ///
    /// Creating a channel with entries that have a different number of PIDs panics:
    /// ```rust,should_panic
    /// use pineappl::boc::Channel;
    ///
    /// let _ = Channel::new(vec![(vec![1, 1, 1], 1.0), (vec![1, 1], 1.0)]);
    /// ```
    #[must_use]
    pub fn new(mut entry: Vec<(Vec<i32>, f64)>) -> Self {
        assert!(!entry.is_empty(), "can not create empty channel");
        assert!(
            entry.iter().map(|(pids, _)| pids.len()).all_equal(),
            "can not create channel with a different number of PIDs"
        );

        // sort `entry` because the ordering doesn't matter and because it makes it easier to
        // compare `Channel` objects with each other
        entry.sort_by(|x, y| x.0.cmp(&y.0));

        Self {
            entry: entry
                .into_iter()
                .coalesce(|lhs, rhs| {
                    // sum the factors of repeated elements
                    if lhs.0 == rhs.0 {
                        Ok((lhs.0, lhs.1 + rhs.1))
                    } else {
                        Err((lhs, rhs))
                    }
                })
                // filter zeros
                // TODO: find a better than to hardcode the epsilon limit
                .filter(|&(_, f)| !approx_eq!(f64, f.abs(), 0.0, epsilon = 1e-14))
                .collect(),
        }
    }

    /// Translates `entry` into a different basis using `translator`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use pineappl::boc::Channel;
    /// use pineappl::channel;
    ///
    /// let entry = channel![10.0 * (103, 11)].translate(&|evol_id| match evol_id {
    ///     103 => vec![(2, 1.0), (-2, -1.0), (1, -1.0), (-1, 1.0)],
    ///     _ => vec![(evol_id, 1.0)],
    /// });
    ///
    /// assert_eq!(entry, channel![10.0 * (2, 11) + -10.0 * (-2, 11) + -10.0 * (1, 11) + 10.0 * (-1, 11)]);
    /// ```
    #[must_use]
    pub fn translate(&self, translator: &dyn Fn(i32) -> Vec<(i32, f64)>) -> Self {
        let mut result = Vec::new();

        for (pids, factor) in &self.entry {
            for tuples in pids
                .iter()
                .map(|&pid| translator(pid))
                .multi_cartesian_product()
            {
                result.push((
                    tuples.iter().map(|&(pid, _)| pid).collect(),
                    factor * tuples.iter().map(|(_, f)| f).product::<f64>(),
                ));
            }
        }

        Self::new(result)
    }

    /// Returns a tuple representation of this entry.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use pineappl::channel;
    /// use pineappl::boc::Channel;
    ///
    /// let entry = channel![1.0 * (4, 4) + 1.0 * (2, 2)];
    ///
    /// assert_eq!(entry.entry(), [(vec![2, 2], 1.0), (vec![4, 4], 1.0)]);
    /// ```
    #[must_use]
    pub fn entry(&self) -> &[(Vec<i32>, f64)] {
        &self.entry
    }

    /// Create a new object with the PIDs at index `i` and `j` transposed.
    #[must_use]
    pub fn transpose(&self, i: usize, j: usize) -> Self {
        Self::new(
            self.entry
                .iter()
                .map(|(pids, c)| {
                    let mut transposed = pids.clone();
                    transposed.swap(i, j);
                    (transposed, *c)
                })
                .collect(),
        )
    }

    /// If `other` is the same channel when only comparing PIDs and neglecting the factors, return
    /// the number `f1 / f2`, where `f1` is the factor from `self` and `f2` is the factor from
    /// `other`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use pineappl::channel;
    ///
    /// let ch1 = channel![2.0 * (2, 2) + 2.0 * (4, 4)];
    /// let ch2 = channel![1.0 * (4, 4) + 1.0 * (2, 2)];
    /// let ch3 = channel![1.0 * (3, 4) + 1.0 * (2, 2)];
    /// let ch4 = channel![1.0 * (4, 3) + 2.0 * (2, 3)];
    /// let ch5 = channel![1.0 * (2, 2) + 2.0 * (4, 4)];
    ///
    /// // ch1 is ch2 multiplied by two
    /// assert_eq!(ch1.common_factor(&ch2), Some(2.0));
    /// // ch1 isn't similar to ch3
    /// assert_eq!(ch1.common_factor(&ch3), None);
    /// // ch1 isn't similar to ch4 either
    /// assert_eq!(ch1.common_factor(&ch4), None);
    /// // ch1 is similar to ch5, but they don't share a common factor
    /// assert_eq!(ch1.common_factor(&ch5), None);
    /// ```
    #[must_use]
    pub fn common_factor(&self, other: &Self) -> Option<f64> {
        if self.entry.len() != other.entry.len() {
            return None;
        }

        let result: Option<Vec<_>> = self
            .entry
            .iter()
            .zip(&other.entry)
            .map(|((pids_a, fa), (pids_b, fb))| (pids_a == pids_b).then_some(fa / fb))
            .collect();

        result.and_then(|factors| {
            if factors
                .windows(2)
                .all(|win| approx_eq!(f64, win[0], win[1], ulps = 4))
            {
                factors.first().copied()
            } else {
                None
            }
        })
    }
}

impl FromStr for Channel {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        let result: Vec<_> = s
            .split('+')
            .map(|sub| {
                sub.split_once('*').map_or_else(
                    // TODO: allow a missing numerical factor which then is assumed to be `1`
                    || Err(Error::General(format!("missing '*' in '{sub}'"))),
                    |(factor, pids)| {
                        let vector: Vec<_> = pids
                            .trim()
                            .strip_prefix('(')
                            .ok_or_else(|| Error::General(format!("missing '(' in '{pids}'")))?
                            .strip_suffix(')')
                            .ok_or_else(|| Error::General(format!("missing ')' in '{pids}'")))?
                            .split(',')
                            .map(|pid| {
                                pid.trim().parse::<i32>().map_err(|err| {
                                    Error::General(format!("could not parse PID: '{pid}', '{err}'"))
                                })
                            })
                            .collect::<Result<_>>()?;

                        Ok((
                            vector,
                            str::parse::<f64>(factor.trim())
                                .map_err(|err| Error::General(err.to_string()))?,
                        ))
                    },
                )
            })
            .collect::<Result<_>>()?;

        if !result.iter().map(|(pids, _)| pids.len()).all_equal() {
            return Err(Error::General(
                "PID tuples have different lengths".to_owned(),
            ));
        }

        Ok(Self::new(result))
    }
}

/// Helper macro to quickly generate a `Channel` at compile time.
///
/// # Examples
///
/// In the following example `entry1` and `entry2` represent the same values:
///
/// ```rust
/// use pineappl::channel;
///
/// let entry1 = channel![1.0 * (2, 2) + 1.0 * (4, 4)];
/// let entry2 = channel![1.0 * (4, 4) + 1.0 * (2, 2)];
///
/// assert_eq!(entry1, entry2);
/// ```
#[macro_export]
macro_rules! channel {
    ($factor:literal * ($($pids:expr),+) $(+ $more_factors:literal * ($($more_pids:expr),+))*) => {
        $crate::boc::Channel::new(
            vec![
                (vec![$($pids),+], $factor) $(,
                (vec![$($more_pids),+], $more_factors)
                )*
            ]
        )
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::assert_approx_eq;

    #[test]
    fn order_from_str() {
        assert_eq!("as1".parse::<Order>().unwrap(), Order::new(1, 0, 0, 0, 0));
        assert_eq!("a1".parse::<Order>().unwrap(), Order::new(0, 1, 0, 0, 0));
        assert_eq!(
            "as1lr1".parse::<Order>().unwrap(),
            Order::new(1, 0, 1, 0, 0)
        );
        assert_eq!(
            "as1lf1".parse::<Order>().unwrap(),
            Order::new(1, 0, 0, 1, 0)
        );
        assert_eq!(
            "as1la1".parse::<Order>().unwrap(),
            Order::new(1, 0, 0, 0, 1)
        );
        assert_eq!(
            "ab12".parse::<Order>().unwrap_err().to_string(),
            "unknown coupling: 'ab'"
        );
        assert_eq!(
            "ab123456789000000"
                .parse::<Order>()
                .unwrap_err()
                .to_string(),
            "error while parsing exponent of 'ab': number too large to fit in target type"
        );
    }

    #[test]
    fn order_cmp() {
        let mut orders = [
            Order::new(1, 2, 1, 0, 0),
            Order::new(1, 2, 0, 1, 0),
            Order::new(1, 2, 0, 0, 0),
            Order::new(0, 3, 1, 0, 0),
            Order::new(0, 3, 0, 1, 0),
            Order::new(0, 3, 0, 0, 0),
            Order::new(0, 2, 0, 0, 0),
        ];

        orders.sort();

        assert_eq!(orders[0], Order::new(0, 2, 0, 0, 0));
        assert_eq!(orders[1], Order::new(1, 2, 0, 0, 0));
        assert_eq!(orders[2], Order::new(1, 2, 0, 1, 0));
        assert_eq!(orders[3], Order::new(1, 2, 1, 0, 0));
        assert_eq!(orders[4], Order::new(0, 3, 0, 0, 0));
        assert_eq!(orders[5], Order::new(0, 3, 0, 1, 0));
        assert_eq!(orders[6], Order::new(0, 3, 1, 0, 0));
    }

    #[test]
    fn order_create_mask() {
        // Drell—Yan orders
        let orders = [
            Order::new(0, 2, 0, 0, 0), //   LO        :          alpha^2
            Order::new(1, 2, 0, 0, 0), //  NLO QCD    : alphas   alpha^2
            Order::new(0, 3, 0, 0, 0), //  NLO  EW    :          alpha^3
            Order::new(2, 2, 0, 0, 0), // NNLO QCD    : alphas^2 alpha^2
            Order::new(1, 3, 0, 0, 0), // NNLO QCD—EW : alphas   alpha^3
            Order::new(0, 4, 0, 0, 0), // NNLO EW     :          alpha^4
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
            Order::new(2, 0, 0, 0, 0), //   LO QCD    : alphas^2
            Order::new(1, 1, 0, 0, 0), //   LO QCD—EW : alphas   alpha
            Order::new(0, 2, 0, 0, 0), //   LO  EW    :          alpha^2
            Order::new(3, 0, 0, 0, 0), //  NLO QCD    : alphas^3
            Order::new(2, 1, 0, 0, 0), //  NLO QCD—EW : alphas^2 alpha
            Order::new(1, 2, 0, 0, 0), //  NLO QCD—EW : alphas   alpha^2
            Order::new(0, 3, 0, 0, 0), //  NLO  EW    :          alpha^3
            Order::new(4, 0, 0, 0, 0), // NNLO QCD    : alphas^4
            Order::new(3, 1, 0, 0, 0), // NNLO QCD—EW : alphas^3 alpha
            Order::new(2, 2, 0, 0, 0), // NNLO QCD—EW : alphas^2 alpha^2
            Order::new(1, 3, 0, 0, 0), // NNLO QCD—EW : alphas   alpha^3
            Order::new(0, 4, 0, 0, 0), // NNLO EW     :          alpha^4
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
    fn channel_from_str() {
        assert_eq!(
            str::parse::<Channel>(" 1   * (  2 , -2) + 2* (4,-4)").unwrap(),
            channel![1.0 * (2, -2) + 2.0 * (4, -4)]
        );

        assert_eq!(
            str::parse::<Channel>("* (  2, -2) + 2* (4,-4)")
                .unwrap_err()
                .to_string(),
            "cannot parse float from empty string"
        );

        assert_eq!(
            str::parse::<Channel>(" 1    (  2 -2) + 2* (4,-4)")
                .unwrap_err()
                .to_string(),
            "missing '*' in ' 1    (  2 -2) '"
        );

        assert_eq!(
            str::parse::<Channel>(" 1   * (  2 -2) + 2* (4,-4)")
                .unwrap_err()
                .to_string(),
            "could not parse PID: '  2 -2', 'invalid digit found in string'"
        );

        assert_eq!(
            str::parse::<Channel>(" 1   *   2, -2) + 2* (4,-4)")
                .unwrap_err()
                .to_string(),
            "missing '(' in '   2, -2) '"
        );

        assert_eq!(
            str::parse::<Channel>(" 1   * (  2, -2 + 2* (4,-4)")
                .unwrap_err()
                .to_string(),
            "missing ')' in ' (  2, -2 '"
        );

        assert_eq!(
            str::parse::<Channel>("1 * (2, 2, 2) + 2 * (4, 4)")
                .unwrap_err()
                .to_string(),
            "PID tuples have different lengths"
        );
    }

    #[test]
    fn scale_func_form() {
        let node_values = [vec![1.0, 2.0, 3.0], vec![4.0, 5.0]];
        let kinematics = [Kinematics::Scale(0), Kinematics::Scale(1)];
        let sff = ScaleFuncForm::QuadraticSum(0, 1);

        let ref_calc = [5.0, 6.0, 6.0, 7.0, 7.0, 8.0];
        let calc = sff.calc(&node_values, &kinematics).into_owned();

        assert_eq!(calc.len(), ref_calc.len());

        for (&calc, ref_calc) in calc.iter().zip(ref_calc) {
            assert_approx_eq!(f64, calc, ref_calc, ulps = 2);
        }

        let scale_dims = [3, 2];

        assert_eq!(sff.idx(&[0, 0, 1], &scale_dims), 0);
        assert_eq!(sff.idx(&[0, 1, 1], &scale_dims), 1);
        assert_eq!(sff.idx(&[1, 0, 1], &scale_dims), 2);
        assert_eq!(sff.idx(&[1, 1, 1], &scale_dims), 3);
        assert_eq!(sff.idx(&[2, 0, 1], &scale_dims), 4);
        assert_eq!(sff.idx(&[2, 1, 1], &scale_dims), 5);
    }

    #[test]
    fn bwfl_limit_parsing_failure() {
        assert_eq!(
            BinsWithFillLimits::from_str("0,1,2,x")
                .unwrap_err()
                .to_string(),
            "unable to parse limit 'x': 'invalid float literal')"
        );
    }

    #[test]
    fn bwfl_pipe_syntax_first_dimension() {
        assert_eq!(
            BinsWithFillLimits::from_str("|0,1,2")
                .unwrap_err()
                .to_string(),
            "'|' syntax not meaningful for first dimension"
        );
    }

    #[test]
    fn bwfl_pipe_syntax_first_empty() {
        assert_eq!(
            BinsWithFillLimits::from_str("0,1,2;0,2,4;||")
                .unwrap_err()
                .to_string(),
            "empty repetition with '|'"
        );
    }

    #[test]
    fn bwfl_colon_syntax_bad_string() {
        assert_eq!(
            BinsWithFillLimits::from_str("0,1,2;0,2,4;1,2,3,4,5|::")
                .unwrap_err()
                .to_string(),
            "unable to parse 'N:M' syntax from: '::' (N: 'cannot parse integer from empty string', M: 'invalid digit found in string')"
        );
    }

    #[test]
    fn bwfl_colon_syntax_bad_lhs() {
        assert_eq!(
            BinsWithFillLimits::from_str("0,1,2;0,2,4;1,2,3,4,5|2.5:|:3|:3")
                .unwrap_err()
                .to_string(),
            "unable to parse 'N:M' syntax from: '2.5:' (N: 'invalid digit found in string', M: 'cannot parse integer from empty string')"
        );
    }

    #[test]
    fn bwfl_colon_syntax_bad_rhs() {
        assert_eq!(
            BinsWithFillLimits::from_str("0,1,2;0,2,4;1,2,3,4,5|:2.5|:3|:3")
                .unwrap_err()
                .to_string(),
            "unable to parse 'N:M' syntax from: ':2.5' (N: 'cannot parse integer from empty string', M: 'invalid digit found in string')"
        );
    }

    #[test]
    fn bwfl_colon_syntax_no_limits() {
        assert_eq!(
            BinsWithFillLimits::from_str("0,1,2;0,2,4;1,2,3,4,5|:4|:3|:3")
                .unwrap_err()
                .to_string(),
            "no limits due to ':' syntax"
        );
    }

    #[test]
    fn bwfl_pipe_syntax_too_few_pipes() {
        assert_eq!(
            BinsWithFillLimits::from_str("0,1,2;0,2,4;1,2,3|4,5,6|7,8,9")
                .unwrap_err()
                .to_string(),
            "missing '|' specification: number of variants too small"
        );
    }
}
