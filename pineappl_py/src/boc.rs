//! Interface for bins, orders and channels.

use numpy::{IntoPyArray, PyArray1};
use pineappl::boc::{Bin, BinsWithFillLimits, Channel, Kinematics, Order, ScaleFuncForm, Scales};
use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::boc::Bin <boc/struct.Bin.html>`.
#[pyclass(name = "Bin")]
#[repr(transparent)]
#[derive(Clone)]
pub struct PyBin {
    pub(crate) bin: Bin,
}

#[pymethods]
impl PyBin {
    /// Constructor.
    ///
    /// Parameters
    /// ----------
    /// bin_limits: list(tuple(float, float))
    ///     edges of the bins
    /// normalization: float
    ///     normalization factor
    #[new]
    #[must_use]
    pub fn new(bin_limits: Vec<(f64, f64)>, normalization: f64) -> Self {
        Self {
            bin: Bin::new(bin_limits, normalization),
        }
    }

    /// Get the dimension for a given bin.
    ///
    /// Returns
    /// -------
    /// int:
    ///     dimension on which the observable is defined
    #[getter]
    #[must_use]
    pub fn dimensions(&self) -> usize {
        self.bin.dimensions()
    }

    /// Get the normalization factor for this given bin.
    ///
    /// Returns
    /// -------
    /// float:
    ///     normalization factor
    #[getter]
    #[must_use]
    pub const fn normalization(&self) -> f64 {
        self.bin.normalization()
    }

    /// Get the edges of this given bin.
    ///
    /// Returns
    /// -------
    /// list(tuple(float, float)):
    ///     edges of the current bin
    #[getter]
    #[must_use]
    pub fn bin_limits(&self) -> Vec<(f64, f64)> {
        self.bin.limits().to_vec()
    }
}

/// PyO3 wrapper to :rustdoc:`pineappl::boc::Bin <boc/struct.Bin.html>`.
#[pyclass(name = "BinsWithFillLimits")]
#[repr(transparent)]
#[derive(Clone)]
pub struct PyBinsWithFillLimits {
    pub(crate) bins_fill_limits: BinsWithFillLimits,
}

#[pymethods]
impl PyBinsWithFillLimits {
    /// Constructor for `BinsWithFillLimits`.
    ///
    /// # Panics
    /// TODO
    ///
    /// # Errors
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// bins: Bin
    ///     a list containing the bin specifications
    /// fill_limits: list(float)
    ///     edges of the bins
    #[new]
    #[must_use]
    pub fn new(bins: Vec<PyRef<PyBin>>, fill_limits: Vec<f64>) -> Self {
        Self {
            bins_fill_limits: BinsWithFillLimits::new(
                bins.into_iter().map(|b| b.bin.clone()).collect(),
                fill_limits,
            )
            .unwrap(),
        }
    }

    /// Construct the bin specifications using the bin edges.
    ///
    /// # Panics
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// fill_limits: list(float)
    ///     edges of the bins
    #[must_use]
    #[staticmethod]
    pub fn from_fill_limits(fill_limits: Vec<f64>) -> Self {
        Self {
            bins_fill_limits: BinsWithFillLimits::from_fill_limits(fill_limits).unwrap(),
        }
    }

    /// Construct the bin specifications using the edges and the normalizations.
    ///
    /// # Panics
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// limits: list(list(float, float))
    ///     edges of the bins
    /// normalizations: list(float)
    ///     a list of normalization factors
    #[must_use]
    #[staticmethod]
    pub fn from_limits_and_normalizations(
        limits: Vec<Vec<(f64, f64)>>,
        normalizations: Vec<f64>,
    ) -> Self {
        Self {
            bins_fill_limits: BinsWithFillLimits::from_limits_and_normalizations(
                limits,
                normalizations,
            )
            .unwrap(),
        }
    }

    /// Get the list of `Bin` specifications.
    ///
    /// Returns
    /// -------
    /// list(Bin):
    ///     a list of `Bin` objects
    #[must_use]
    pub fn bins(&self) -> Vec<PyBin> {
        self.bins_fill_limits
            .bins()
            .iter()
            .map(|b| PyBin { bin: b.clone() })
            .collect()
    }

    /// Get the size/length of the bins.
    ///
    /// Returns
    /// -------
    /// int:
    ///     the size of the bins
    #[must_use]
    pub fn len(&self) -> usize {
        self.bins_fill_limits.len()
    }

    /// Get the dimensions of the bins.
    ///
    /// Returns
    /// -------
    /// int:
    ///     dimension of the defined observable
    #[must_use]
    pub fn dimensions(&self) -> usize {
        self.bins_fill_limits.dimensions()
    }

    /// Get the removed bin using the index.
    ///
    /// Parameters
    /// ----------
    /// index: int
    ///     index of the bin to be removed
    /// Returns
    /// -------
    /// Bin:
    ///     a `Bin` object from the removed index
    #[must_use]
    pub fn removed_index(&self, index: usize) -> PyBin {
        PyBin {
            bin: self.bins_fill_limits.clone().remove(index),
        }
    }

    /// Get the bin slices.
    ///
    /// Returns
    /// -------
    /// list(list(int)):
    ///     a list of indices representing the slices
    #[must_use]
    pub fn slices(&self) -> Vec<Vec<usize>> {
        self.bins_fill_limits
            .slices()
            .iter()
            .map(|slice| slice.clone().collect())
            .collect()
    }

    /// Get the bin limits/edges.
    ///
    /// Returns
    /// -------
    /// list(list(tuple(float, float))):
    ///     the bin edges with shape (n_bins, n_dimension, 2)
    #[must_use]
    pub fn bin_limits(&self) -> Vec<Vec<(f64, f64)>> {
        self.bins_fill_limits
            .bins()
            .iter()
            .map(|b| b.limits().to_vec())
            .collect()
    }

    /// Get the normalizations of all bins.
    ///
    /// Returns
    /// -------
    /// list(float):
    ///     a list containing the normalizations
    #[must_use]
    pub fn bin_normalizations(&self) -> Vec<f64> {
        self.bins_fill_limits
            .bins()
            .iter()
            .map(Bin::normalization)
            .collect()
    }
}

/// PyO3 wrapper to :rustdoc:`pineappl::boc::Channel <boc/struct.Channel.html>`.
///
/// Each entry consists of a tuple, which contains, in the following order:
///
/// 1. a list containing the PDG value of the 1st, 2nd, and etc. of the incoming parton
/// 2. a numerical factor that will multiply the result for this specific combination.
#[pyclass(name = "Channel")]
#[repr(transparent)]
pub struct PyChannel {
    pub(crate) entry: Channel,
}

#[pymethods]
impl PyChannel {
    /// Constructor.
    ///
    /// Parameters
    /// ----------
    /// entry: list(tuple(list(int),float))
    ///     channel configuration
    #[new]
    #[must_use]
    pub fn new(entry: Vec<(Vec<i32>, f64)>) -> Self {
        Self {
            entry: Channel::new(entry),
        }
    }

    /// Get list representation.
    ///
    /// Returns
    /// -------
    /// list(tuple(list(int),float)) :
    ///     list representation
    #[must_use]
    pub fn into_array(&self) -> Vec<(Vec<i32>, f64)> {
        self.entry.entry().to_vec()
    }
}

/// PyO3 wrapper to :rustdoc:`pineappl::boc::Kinematics <boc/enum.Kinematics.html>`.
#[pyclass(eq, name = "Kinematics")]
#[derive(Clone, PartialEq, Eq)]
pub enum PyKinematics {
    /// map to Kinematics::Scale
    Scale(usize),
    /// map to Kinematics::X
    X(usize),
}

impl From<PyKinematics> for Kinematics {
    fn from(item: PyKinematics) -> Self {
        match item {
            PyKinematics::X(v) => Self::X(v),
            PyKinematics::Scale(v) => Self::Scale(v),
        }
    }
}

/// PyO3 wrapper to :rustdoc:`pineappl::boc::ScaleFuncForm <boc/enum.ScaleFuncForm.html>`.
#[pyclass(eq, name = "ScaleFuncForm")]
#[derive(Clone, PartialEq, Eq)]
pub enum PyScaleFuncForm {
    /// map to ScaleFuncForm::NoScale
    /// NOTE No variant is not supported in complex enums
    NoScale(usize),
    /// map to ScaleFuncForm::Scale
    Scale(usize),
    /// map to ScaleFuncForm::QuadraticSum
    QuadraticSum(usize, usize),
    /// map to ScaleFuncForm::QuadraticMean
    QuadraticMean(usize, usize),
    /// map to ScaleFuncForm::QuadraticSumOver4
    QuadraticSumOver4(usize, usize),
    /// map to ScaleFuncForm::LinearMean
    LinearMean(usize, usize),
    /// map to ScaleFuncForm::LinearSum
    LinearSum(usize, usize),
    /// map to ScaleFuncForm::ScaleMax
    ScaleMax(usize, usize),
    /// map to ScaleFuncForm::ScaleMin
    ScaleMin(usize, usize),
    /// map to ScaleFuncForm::Prod
    Prod(usize, usize),
    /// map to ScaleFuncForm::S2plusS1half
    S2plusS1half(usize, usize),
    /// map to ScaleFuncForm::Pow4Sum
    Pow4Sum(usize, usize),
    /// map to ScaleFuncForm::WgtAvg
    WgtAvg(usize, usize),
    /// map to ScaleFuncForm::S2plusS1fourth
    S2plusS1fourth(usize, usize),
    /// map to ScaleFuncForm::ExpProd2
    ExpProd2(usize, usize),
}

impl From<PyScaleFuncForm> for ScaleFuncForm {
    fn from(item: PyScaleFuncForm) -> Self {
        match item {
            PyScaleFuncForm::NoScale(_) => Self::NoScale,
            PyScaleFuncForm::Scale(v) => Self::Scale(v),
            PyScaleFuncForm::QuadraticSum(v1, v2) => Self::QuadraticSum(v1, v2),
            PyScaleFuncForm::QuadraticMean(v1, v2) => Self::QuadraticMean(v1, v2),
            PyScaleFuncForm::QuadraticSumOver4(v1, v2) => Self::QuadraticSumOver4(v1, v2),
            PyScaleFuncForm::LinearMean(v1, v2) => Self::LinearMean(v1, v2),
            PyScaleFuncForm::LinearSum(v1, v2) => Self::LinearSum(v1, v2),
            PyScaleFuncForm::ScaleMax(v1, v2) => Self::ScaleMax(v1, v2),
            PyScaleFuncForm::ScaleMin(v1, v2) => Self::ScaleMin(v1, v2),
            PyScaleFuncForm::Prod(v1, v2) => Self::Prod(v1, v2),
            PyScaleFuncForm::S2plusS1half(v1, v2) => Self::S2plusS1half(v1, v2),
            PyScaleFuncForm::Pow4Sum(v1, v2) => Self::Pow4Sum(v1, v2),
            PyScaleFuncForm::WgtAvg(v1, v2) => Self::WgtAvg(v1, v2),
            PyScaleFuncForm::S2plusS1fourth(v1, v2) => Self::S2plusS1fourth(v1, v2),
            PyScaleFuncForm::ExpProd2(v1, v2) => Self::ExpProd2(v1, v2),
        }
    }
}

/// PyO3 wrapper to :rustdoc:`pineappl::boc::Scales <boc/struct.Scales.html>`.
#[pyclass(name = "Scales")]
pub struct PyScales {
    pub(crate) scales: Scales,
}

impl PyScales {
    pub(crate) const fn new(scales: Scales) -> Self {
        Self { scales }
    }
}

#[pymethods]
impl PyScales {
    /// Constructor for `Scales`
    #[new]
    #[must_use]
    pub fn news_scales(
        ren: PyRef<PyScaleFuncForm>,
        fac: PyRef<PyScaleFuncForm>,
        frg: PyRef<PyScaleFuncForm>,
    ) -> Self {
        let ren = ren.clone().into();
        let fac = fac.clone().into();
        let frg = frg.clone().into();
        Self::new(Scales { ren, fac, frg })
    }
}

/// PyO3 wrapper to :rustdoc:`pineappl::boc::Order <boc/struct.Order.html>`.
#[pyclass(name = "Order")]
#[repr(transparent)]
pub struct PyOrder {
    pub(crate) order: Order,
}

impl PyOrder {
    pub(crate) const fn new(order: Order) -> Self {
        Self { order }
    }
}

#[pymethods]
impl PyOrder {
    /// Constructor.
    ///
    /// Parameters
    /// ----------
    /// alphas : int
    ///     power of :math:`\alpha_s`
    /// alpha : int
    ///     power of :math:`\alpha`
    /// logxir : int
    ///     power of :math:`\ln(\xi_r)`
    /// logxif : int
    ///     power of :math:`\ln(\xi_f)`
    /// logxia : int
    ///     power of :math:`\ln(\xi_a)`
    #[new]
    #[must_use]
    pub const fn new_order(alphas: u8, alpha: u8, logxir: u8, logxif: u8, logxia: u8) -> Self {
        Self::new(Order::new(alphas, alpha, logxir, logxif, logxia))
    }

    /// Tuple representation.
    ///
    /// Returns
    /// -------
    /// alphas : int
    ///     power of :math:`\alpha_s`
    /// alpha : int
    ///     power of :math:`\alpha`
    /// logxir : int
    ///     power of :math:`\ln(\xi_r)`
    /// logxif : int
    ///     power of :math:`\ln(\xi_f)`
    /// logxia : int
    ///     power of :math:`\ln(\xi_a)`
    #[must_use]
    pub const fn as_tuple(&self) -> (u8, u8, u8, u8, u8) {
        (
            self.order.alphas,
            self.order.alpha,
            self.order.logxir,
            self.order.logxif,
            self.order.logxia,
        )
    }

    /// Return a mask suitable to pass as the `order_mask` parameter of [`Grid::convolve`].
    ///
    /// The selection of `orders` is controlled using the `max_as` and `max_al` parameters, for
    /// instance setting `max_as = 1` and `max_al = 0` selects the LO QCD only, `max_as = 2` and
    /// `max_al = 0` the NLO QCD; setting `max_as = 3` and `max_al = 2` would select all NLOs, and
    /// the NNLO QCD.
    ///
    /// See `pineappl` crate docs for more examples.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray(bool)
    ///     boolean array, to be used as orders' mask
    #[staticmethod]
    #[must_use]
    pub fn create_mask<'py>(
        orders: Vec<PyRef<Self>>,
        max_as: u8,
        max_al: u8,
        logs: bool,
        py: Python<'py>,
    ) -> Bound<'py, PyArray1<bool>> {
        Order::create_mask(
            &orders.iter().map(|o| o.order.clone()).collect::<Vec<_>>(),
            max_as,
            max_al,
            logs,
        )
        .into_pyarray(py)
    }
}

/// Register submodule in parent.
///
/// # Errors
///
/// Raises an error if (sub)module is not found.
pub fn register(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "boc")?;
    m.setattr(
        pyo3::intern!(m.py(), "__doc__"),
        "Interface for bins, orders and channels.",
    )?;
    pyo3::py_run!(
        parent_module.py(),
        m,
        "import sys; sys.modules['pineappl.boc'] = m"
    );
    m.add_class::<PyBin>()?;
    m.add_class::<PyBinsWithFillLimits>()?;
    m.add_class::<PyChannel>()?;
    m.add_class::<PyOrder>()?;
    m.add_class::<PyKinematics>()?;
    m.add_class::<PyScaleFuncForm>()?;
    m.add_class::<PyScales>()?;
    parent_module.add_submodule(&m)
}
