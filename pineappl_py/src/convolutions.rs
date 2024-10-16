//! Convolution interface.

use pineappl::convolutions::{Conv, ConvType};
use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::`.
#[pyclass(name = "ConvType")]
#[repr(transparent)]
pub struct PyConvType {
    pub(crate) convtype: ConvType,
}

impl PyConvType {
    pub(crate) const fn new(convtype: ConvType) -> Self {
        Self { convtype }
    }
}

#[pymethods]
impl PyConvType {
    /// Constructor.
    #[new]
    #[must_use]
    pub const fn new_convtype(polarized: bool, time_like: bool) -> Self {
        Self::new(ConvType::new(polarized, time_like))
    }
}

/// PyO3 wrapper to :rustdoc:`pineappl::`.
#[pyclass(name = "Conv")]
#[repr(transparent)]
pub struct PyConv {
    pub(crate) conv: Conv,
}

impl PyConv {
    pub(crate) const fn new(conv: Conv) -> Self {
        Self { conv }
    }
}

#[pymethods]
impl PyConv {
    /// Constructor.
    #[new]
    #[must_use]
    #[allow(clippy::needless_pass_by_value)]
    pub fn new_conv(conv_type: PyRef<PyConvType>, pid: i32) -> Self {
        Self::new(Conv::new(conv_type.convtype, pid))
    }
}

/// Register submodule in parent.
/// # Errors
///
/// Raises an error if (sub)module is not found.
pub fn register(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new_bound(parent_module.py(), "convolutions")?;
    m.setattr(
        pyo3::intern!(m.py(), "__doc__"),
        "Define the type of convolutions.",
    )?;
    pyo3::py_run!(
        parent_module.py(),
        m,
        "import sys; sys.modules['pineappl.convolutions'] = m"
    );
    m.add_class::<PyConvType>()?;
    m.add_class::<PyConv>()?;
    parent_module.add_submodule(&m)
}
