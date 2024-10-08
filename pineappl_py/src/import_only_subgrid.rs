//! PyImportOnlySubgrid* interface.

use super::subgrid::PySubgridEnum;
use numpy::PyReadonlyArray3;
use pineappl::import_only_subgrid::ImportOnlySubgridV1;
use pineappl::import_only_subgrid::ImportOnlySubgridV2;
use pineappl::sparse_array3::SparseArray3;
use pineappl::subgrid::Mu2;
use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::import_only_subgrid::ImportOnlySubgridV2 <import_only_subgrid/struct.ImportOnlySubgridV2.html>`.
#[pyclass(name = "ImportOnlySubgridV2")]
#[derive(Clone)]
#[repr(transparent)]
pub struct PyImportOnlySubgridV2 {
    pub(crate) import_only_subgrid: ImportOnlySubgridV2,
}

#[pymethods]
impl PyImportOnlySubgridV2 {
    /// Constructor.
    ///
    /// Parameters
    /// ----------
    /// array : numpy.ndarray(float)
    ///     3D array with all weights
    /// mu2_grid : list(tuple(float, float))
    ///     scales grid
    /// x1_grid : list(float)
    ///     first momentum fraction grid
    /// x2_grid : list(float)
    ///     second momentum fraction grid
    #[new]
    pub fn new(
        array: PyReadonlyArray3<f64>,
        mu2_grid: Vec<(f64, f64)>,
        x1_grid: Vec<f64>,
        x2_grid: Vec<f64>,
    ) -> Self {
        let mut sparse_array = SparseArray3::new(mu2_grid.len(), x1_grid.len(), x2_grid.len());

        for ((imu2, ix1, ix2), value) in array
            .as_array()
            .indexed_iter()
            .filter(|((_, _, _), value)| **value != 0.0)
        {
            sparse_array[[imu2, ix1, ix2]] = *value;
        }
        Self {
            import_only_subgrid: ImportOnlySubgridV2::new(
                sparse_array,
                mu2_grid
                    .iter()
                    .map(|(ren, fac)| Mu2 {
                        ren: *ren,
                        fac: *fac,
                    })
                    .collect(),
                x1_grid,
                x2_grid,
            ),
        }
    }

    /// Wrapper to match :meth:`pineappl.pineappl.PyGrid.set_subgrid()`.
    ///
    /// Returns
    /// -------
    /// PySubgridEnum :
    ///     casted object
    pub fn into(&self) -> PySubgridEnum {
        PySubgridEnum {
            subgrid_enum: self.import_only_subgrid.clone().into(),
        }
    }
}

/// PyO3 wrapper to :rustdoc:`pineappl::import_only_subgrid::ImportOnlySubgridV1 <import_only_subgrid/struct.ImportOnlySubgridV1.html>`.
#[pyclass(name = "ImportOnlySubgridV1")]
#[derive(Clone)]
#[repr(transparent)]
pub struct PyImportOnlySubgridV1 {
    pub(crate) import_only_subgrid: ImportOnlySubgridV1,
}

impl PyImportOnlySubgridV1 {
    pub(crate) fn new(import_only_subgrid: ImportOnlySubgridV1) -> Self {
        Self {
            import_only_subgrid,
        }
    }
}

#[pymethods]
impl PyImportOnlySubgridV1 {
    /// Constructor.
    ///
    /// Parameters
    /// ----------
    /// array : numpy.ndarray(float)
    ///     3D array with all weights
    /// mu2_grid : list(tuple(float, float))
    ///     scales grid
    /// x1_grid : list(float)
    ///     first momentum fraction grid
    /// x2_grid : list(float)
    ///     second momentum fraction grid
    #[new]
    pub fn new_import_only_subgrid(
        array: PyReadonlyArray3<f64>,
        q2_grid: Vec<f64>,
        x1_grid: Vec<f64>,
        x2_grid: Vec<f64>,
    ) -> Self {
        let mut sparse_array = SparseArray3::new(q2_grid.len(), x1_grid.len(), x2_grid.len());

        for ((iq2, ix1, ix2), value) in array
            .as_array()
            .indexed_iter()
            .filter(|((_, _, _), value)| **value != 0.0)
        {
            sparse_array[[iq2, ix1, ix2]] = *value;
        }

        Self::new(ImportOnlySubgridV1::new(
            sparse_array,
            q2_grid,
            x1_grid,
            x2_grid,
        ))
    }

    /// Wrapper to match :meth:`pineappl.pineappl.PyGrid.set_subgrid()`.
    ///
    /// Returns
    /// -------
    /// PySubgridEnum :
    ///     casted object
    pub fn into(&self) -> PySubgridEnum {
        PySubgridEnum {
            subgrid_enum: self.import_only_subgrid.clone().into(),
        }
    }
}

/// Register submodule in parent.
pub fn register(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new_bound(parent_module.py(), "import_only_subgrid")?;
    m.setattr(
        pyo3::intern!(m.py(), "__doc__"),
        "ImportOnlySubgrid* interface.",
    )?;
    pyo3::py_run!(
        parent_module.py(),
        m,
        "import sys; sys.modules['pineappl.import_only_subgrid'] = m"
    );
    m.add_class::<PyImportOnlySubgridV1>()?;
    m.add_class::<PyImportOnlySubgridV2>()?;
    parent_module.add_submodule(&m)
}
