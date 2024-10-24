//! Interpolation interface.

use pineappl::interpolation::{Interp, InterpMeth, Map, ReweightMeth};
use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::interpolation::Interp <interpolation/struct.Interp.html>`.
#[pyclass(name = "Interp")]
#[repr(transparent)]
pub struct PyInterp {
    pub(crate) interp: Interp,
}

impl PyInterp {
    pub(crate) const fn new(interp: Interp) -> Self {
        Self { interp }
    }
}

/// TODO
#[pyclass(eq, eq_int, name = "ReweightingMethod")]
#[derive(Clone, PartialEq, Eq)]
pub enum PyReweightingMethod {
    /// map to ReweightMeth::NoReweight
    NoReweight,
    /// map to ReweightMeth::ApplGridX
    ApplGridX,
}

impl From<PyReweightingMethod> for ReweightMeth {
    fn from(value: PyReweightingMethod) -> Self {
        match value {
            PyReweightingMethod::NoReweight => Self::NoReweight,
            PyReweightingMethod::ApplGridX => Self::ApplGridX,
        }
    }
}

/// TODO
#[pyclass(eq, eq_int, name = "MappingMethod")]
#[derive(Clone, PartialEq, Eq)]
pub enum PyMappingMethod {
    /// map to Map::ApplGridF2
    ApplGridF2,
    /// map to Map::ApplGridH0
    ApplGridH0,
}

impl From<PyMappingMethod> for Map {
    fn from(value: PyMappingMethod) -> Self {
        match value {
            PyMappingMethod::ApplGridF2 => Self::ApplGridF2,
            PyMappingMethod::ApplGridH0 => Self::ApplGridH0,
        }
    }
}

/// TODO
#[pyclass(eq, eq_int, name = "InterpolationMethod")]
#[derive(Clone, PartialEq, Eq)]
pub enum PyInterpolationMethod {
    /// map to InterpMeth::Lagrange
    Lagrange,
}

impl From<PyInterpolationMethod> for InterpMeth {
    fn from(value: PyInterpolationMethod) -> Self {
        match value {
            PyInterpolationMethod::Lagrange => Self::Lagrange,
        }
    }
}

#[pymethods]
impl PyInterp {
    /// Constructor.
    ///
    /// TODO: Exposes ReweightMeth,reweight Map, and InterpMeth
    ///
    /// Parameteters
    /// ------------
    /// min : float
    ///     minimum value of the node
    /// max : float
    ///     maximum value of the node
    /// nodes : int
    ///     number of nodes
    /// order : int
    ///     order of the interpolation
    /// reweght_meth : Optional[PyReweightingMethod]
    ///     re-weighting method to be used
    /// map : Optional[PyMappingMethod]
    ///     the type of mapping to be used
    /// interpolation_meth : Optional[PyInterpolationMethod]
    ///     the type of interpolation to be used
    #[new]
    #[must_use]
    #[pyo3(signature = (min, max, nodes, order, reweight_meth = None, map = None, interpolation_meth = None))]
    pub fn new_interp(
        min: f64,
        max: f64,
        nodes: usize,
        order: usize,
        reweight_meth: Option<PyReweightingMethod>,
        map: Option<PyMappingMethod>,
        interpolation_meth: Option<PyInterpolationMethod>,
    ) -> Self {
        let reweight = reweight_meth.unwrap_or(PyReweightingMethod::NoReweight);
        let mapping = map.unwrap_or(PyMappingMethod::ApplGridF2);
        let interp_method = interpolation_meth.unwrap_or(PyInterpolationMethod::Lagrange);

        Self::new(Interp::new(
            min,
            max,
            nodes,
            order,
            reweight.into(),
            mapping.into(),
            interp_method.into(),
        ))
    }
}

/// Register submodule in parent.
/// # Errors
///
/// Raises an error if (sub)module is not found.
pub fn register(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new_bound(parent_module.py(), "interpolation")?;
    m.setattr(pyo3::intern!(m.py(), "__doc__"), "Interpolation submodule.")?;
    pyo3::py_run!(
        parent_module.py(),
        m,
        "import sys; sys.modules['pineappl.interpolation'] = m"
    );
    m.add_class::<PyInterp>()?;
    m.add_class::<PyReweightingMethod>()?;
    m.add_class::<PyMappingMethod>()?;
    m.add_class::<PyInterpolationMethod>()?;
    parent_module.add_submodule(&m)
}
