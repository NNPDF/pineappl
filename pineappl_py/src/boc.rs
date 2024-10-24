//! Interface for bins, orders and channels.

use numpy::{IntoPyArray, PyArray1};
use pineappl::boc::{Channel, Kinematics, Order, ScaleFuncForm, Scales};
use pyo3::prelude::*;

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
#[pyclass(name = "Kinematics")]
#[repr(transparent)]
pub struct PyKinematics {
    pub(crate) kinematics: Kinematics,
}

impl PyKinematics {
    pub(crate) const fn new(kinematics: Kinematics) -> Self {
        Self { kinematics }
    }
}

#[pymethods]
impl PyKinematics {
    /// Constructor.
    ///
    /// Parameters
    /// ----------
    /// kintype: str
    ///     the type of the kinematics, can be either `Scale` or `X`
    /// value: int
    ///     the index mapping to the value, eg. for `X` the integer `0` => `X1`
    #[new]
    #[must_use]
    pub fn new_kin(kintype: &str, value: usize) -> Self {
        let kins = match kintype {
            "X" => Kinematics::X(value),
            "Scale" => Kinematics::Scale(value),
            _ => todo!(),
        };
        Self::new(kins)
    }
}

/// PyO3 wrapper to :rustdoc:`pineappl::boc::ScaleFuncForm <boc/enum.ScaleFuncForm.html>`.
#[pyclass(name = "ScaleFuncForm")]
#[repr(transparent)]
pub struct PyScaleFuncForm {
    pub(crate) scale_func: ScaleFuncForm,
}

impl PyScaleFuncForm {
    pub(crate) const fn new(scale_func: ScaleFuncForm) -> Self {
        Self { scale_func }
    }
}

#[pymethods]
impl PyScaleFuncForm {
    /// Constructor.
    ///
    /// Parameters
    /// ----------
    /// scaletype: str
    ///     the type of the scales, can be `NoScale`, `Scale(int)`, etc.
    /// value: list(int)
    ///     list containing the index mapping to the corresponding value
    #[new]
    #[must_use]
    pub fn new_scale(scaletype: &str, value: Vec<usize>) -> Self {
        let scale = match scaletype {
            "NoScale" => ScaleFuncForm::NoScale,
            "Scale" => ScaleFuncForm::Scale(value[0]),
            "QuadraticSum" => ScaleFuncForm::QuadraticSum(value[0], value[1]),
            "QuadraticMean" => ScaleFuncForm::QuadraticMean(value[0], value[1]),
            "QuadraticSumOver4" => ScaleFuncForm::QuadraticSumOver4(value[0], value[1]),
            "LinearMean" => ScaleFuncForm::LinearMean(value[0], value[1]),
            "LinearSum" => ScaleFuncForm::LinearSum(value[0], value[1]),
            "ScaleMax" => ScaleFuncForm::ScaleMax(value[0], value[1]),
            "ScaleMin" => ScaleFuncForm::ScaleMin(value[0], value[1]),
            "Prod" => ScaleFuncForm::Prod(value[0], value[1]),
            "S2plusS1half" => ScaleFuncForm::S2plusS1half(value[0], value[1]),
            "Pow4Sum" => ScaleFuncForm::Pow4Sum(value[0], value[1]),
            "WgtAvg" => ScaleFuncForm::WgtAvg(value[0], value[1]),
            "S2plusS1fourth" => ScaleFuncForm::S2plusS1fourth(value[0], value[1]),
            "ExpProd2" => ScaleFuncForm::ExpProd2(value[0], value[1]),
            _ => todo!(),
        };
        Self::new(scale)
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

impl Default for PyScales {
    fn default() -> Self {
        Self::new(Scales {
            ren: ScaleFuncForm::Scale(0),
            fac: ScaleFuncForm::Scale(0),
            frg: ScaleFuncForm::NoScale,
        })
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
        let ren = ren.scale_func.clone();
        let fac = fac.scale_func.clone();
        let frg = frg.scale_func.clone();
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
        .into_pyarray_bound(py)
    }
}

/// Register submodule in parent.
/// # Errors
///
/// Raises an error if (sub)module is not found.
pub fn register(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new_bound(parent_module.py(), "boc")?;
    m.setattr(
        pyo3::intern!(m.py(), "__doc__"),
        "Interface for bins, orders and channels.",
    )?;
    pyo3::py_run!(
        parent_module.py(),
        m,
        "import sys; sys.modules['pineappl.boc'] = m"
    );
    m.add_class::<PyChannel>()?;
    m.add_class::<PyOrder>()?;
    m.add_class::<PyKinematics>()?;
    m.add_class::<PyScaleFuncForm>()?;
    m.add_class::<PyScales>()?;
    parent_module.add_submodule(&m)
}
