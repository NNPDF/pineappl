//! Convolution interface.

use pineappl::convolutions::{Conv, ConvType};
use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::convolutions::ConvType <convolutions/enum.ConvType.html>`.
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

    /// Returns a boolean on whether or not the convolution type is polarized
    #[getter]
    #[must_use]
    pub const fn polarized(&self) -> bool {
        matches!(self.convtype, ConvType::PolPDF | ConvType::PolFF)
    }

    /// Returns a boolean on whether or not the convolution type is timelike
    #[getter]
    #[must_use]
    pub const fn time_like(&self) -> bool {
        matches!(self.convtype, ConvType::UnpolFF | ConvType::PolFF)
    }
}

/// PyO3 wrapper to :rustdoc:`pineappl::convolutions::Conv <convolutions/struct.Conv.html>`.
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
    pub fn new_conv(conv_type: PyRef<PyConvType>, pid: i32) -> Self {
        Self::new(Conv::new(conv_type.convtype, pid))
    }

    /// Return the convolution type of this convolution.
    #[getter]
    #[must_use]
    pub const fn conv_type(&self) -> PyConvType {
        PyConvType {
            convtype: self.conv.conv_type(),
        }
    }
}

/// Register submodule in parent.
///
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
