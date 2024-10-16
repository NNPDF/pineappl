//! Interface for bins, orders and channels.

use numpy::{IntoPyArray, PyArray1};
use pineappl::boc::{Channel, Order};
use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::boc::Channel <boc/struct.Channel.html>`.
///
/// Each entry consists of a tuple, which contains, in the following order:
///
/// 1. the PDG id of the first incoming parton
/// 2. the PDG id of the second parton
/// 3. a numerical factor that will multiply the result for this specific combination.
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
    pub const fn new_order(alphas: u32, alpha: u32, logxir: u32, logxif: u32, logxia: u32) -> Self {
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
    #[must_use]
    pub const fn as_tuple(&self) -> (u32, u32, u32, u32) {
        (
            self.order.alphas,
            self.order.alpha,
            self.order.logxir,
            self.order.logxif,
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
    #[allow(clippy::needless_pass_by_value)]
    pub fn create_mask<'py>(
        orders: Vec<PyRef<Self>>,
        max_as: u32,
        max_al: u32,
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
        "import sys; sys.modules['pineappl.channel'] = m"
    );
    m.add_class::<PyChannel>()?;
    m.add_class::<PyOrder>()?;
    parent_module.add_submodule(&m)
}
