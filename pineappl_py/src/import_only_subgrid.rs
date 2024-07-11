use super::subgrid::PySubgridEnum;

use numpy::{PyArrayMethods, PyReadonlyArray1, PyReadonlyArray3};
use pineappl::import_only_subgrid::ImportOnlySubgridV1;
use pineappl::import_only_subgrid::ImportOnlySubgridV2;
use pineappl::sparse_array3::SparseArray3;
use pineappl::subgrid::Mu2;
use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::import_only_subgrid::ImportOnlySubgridV2 <import_only_subgrid/struct.ImportOnlySubgridV2.html>`
///
/// **Usage**: `pineko`
#[pyclass]
#[derive(Clone)]
#[repr(transparent)]
pub struct PyImportOnlySubgridV2 {
    pub(crate) import_only_subgrid: ImportOnlySubgridV2,
}

#[pymethods]
impl PyImportOnlySubgridV2 {
    #[new]
    pub fn new(
        array: PyReadonlyArray3<f64>,
        mu2_grid: Vec<(f64, f64)>,
        x1_grid: PyReadonlyArray1<f64>,
        x2_grid: PyReadonlyArray1<f64>,
    ) -> Self {
        let mut sparse_array = SparseArray3::new(
            mu2_grid.len(),
            x1_grid.as_slice().unwrap().len(),
            x2_grid.as_slice().unwrap().len(),
        );

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
                x1_grid.to_vec().unwrap(),
                x2_grid.to_vec().unwrap(),
            ),
        }
    }

    /// Wrapper to match :meth:`pineappl.pineappl.PyGrid.set_subgrid()`
    ///
    /// Returns
    /// -------
    ///     PySubgridEnum :
    ///         casted object
    pub fn into(&self) -> PySubgridEnum {
        PySubgridEnum {
            subgrid_enum: self.import_only_subgrid.clone().into(),
        }
    }
}

/// PyO3 wrapper to :rustdoc:`pineappl::import_only_subgrid::ImportOnlySubgridV1 <import_only_subgrid/struct.ImportOnlySubgridV1.html>`
///
/// **Usage**: `yadism`
#[pyclass]
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
    #[new]
    pub fn new_import_only_subgrid(
        array: PyReadonlyArray3<f64>,
        q2_grid: PyReadonlyArray1<f64>,
        x1_grid: PyReadonlyArray1<f64>,
        x2_grid: PyReadonlyArray1<f64>,
    ) -> Self {
        let mut sparse_array = SparseArray3::new(
            q2_grid.as_slice().unwrap().len(),
            x1_grid.as_slice().unwrap().len(),
            x2_grid.as_slice().unwrap().len(),
        );

        for ((iq2, ix1, ix2), value) in array
            .as_array()
            .indexed_iter()
            .filter(|((_, _, _), value)| **value != 0.0)
        {
            sparse_array[[iq2, ix1, ix2]] = *value;
        }

        Self::new(ImportOnlySubgridV1::new(
            sparse_array,
            q2_grid.to_vec().unwrap(),
            x1_grid.to_vec().unwrap(),
            x2_grid.to_vec().unwrap(),
        ))
    }

    /// Wrapper to match :meth:`pineappl.pineappl.PyGrid.set_subgrid()`
    ///
    /// Returns
    /// -------
    ///     PySubgridEnum :
    ///         casted object
    pub fn into(&self) -> PySubgridEnum {
        PySubgridEnum {
            subgrid_enum: self.import_only_subgrid.clone().into(),
        }
    }
}
