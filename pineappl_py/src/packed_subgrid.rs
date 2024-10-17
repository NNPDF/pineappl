//! PyPackedSubgrid* interface.

use super::subgrid::PySubgridEnum;
use numpy::PyReadonlyArray3;
use pineappl::packed_array::PackedArray;
use pineappl::packed_subgrid::PackedQ1X2SubgridV1;
use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl`.
#[pyclass(name = "PackedSubgrid")]
#[derive(Clone)]
#[repr(transparent)]
pub struct PyPackedSubgrid {
    pub(crate) packed_subgrid: PackedQ1X2SubgridV1,
}

#[pymethods]
impl PyPackedSubgrid {
    /// Constructor.
    /// Constructor.
    ///
    /// Parameters
    /// ----------
    /// array : numpy.ndarray(float)
    ///     3D array with all weights
    /// scales : list(float)
    ///     scales grid
    /// x1_grid : list(float)
    ///     first momentum fraction grid
    /// x2_grid : list(float)
    ///     second momentum fraction grid
    #[new]
    pub fn new(
        array: PyReadonlyArray3<f64>,
        scales: Vec<f64>,
        x1_grid: Vec<f64>,
        x2_grid: Vec<f64>,
    ) -> Self {
        let node_values: Vec<Vec<f64>> = vec![scales, x1_grid, x2_grid];
        let mut sparse_array: PackedArray<f64> =
            PackedArray::new(node_values.iter().map(Vec::len).collect());

        for ((iscale, ix1, ix2), value) in array
            .as_array()
            .indexed_iter()
            .filter(|((_, _, _), value)| **value != 0.0)
        {
            sparse_array[[iscale, ix1, ix2]] = *value;
        }

        Self {
            packed_subgrid: PackedQ1X2SubgridV1::new(sparse_array, node_values),
        }
    }

    /// TODO
    #[must_use]
    pub fn into(&self) -> PySubgridEnum {
        PySubgridEnum {
            subgrid_enum: self.packed_subgrid.clone().into(),
        }
    }
}

/// Register submodule in parent.
/// # Errors
///
/// Raises an error if (sub)module is not found.
pub fn register(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new_bound(parent_module.py(), "packed_subgrid")?;
    m.setattr(
        pyo3::intern!(m.py(), "__doc__"),
        "Interface for packed subgrid specs.",
    )?;
    pyo3::py_run!(
        parent_module.py(),
        m,
        "import sys; sys.modules['pineappl.packed_subgrid'] = m"
    );
    m.add_class::<PyPackedSubgrid>()?;
    parent_module.add_submodule(&m)
}
