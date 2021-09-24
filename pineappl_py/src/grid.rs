use pineappl::grid::{EkoInfo, Grid, Ntuple, Order};

use super::bin::PyBinRemapper;
use super::fk_table::PyFkTable;
use super::lumi::PyLumiEntry;
use super::subgrid::{PySubgridEnum, PySubgridParams};

use ndarray::{Array, Ix5};

use std::fs::File;
use std::io::BufReader;

use pyo3::prelude::*;

/// PyO3 wrapper to :rustdoc:`pineappl::grid::Order <grid/struct.Order.html>`
///
/// **Usage**: `yadism`
#[pyclass]
#[repr(transparent)]
pub struct PyOrder {
    pub(crate) order: Order,
}

impl PyOrder {
    pub(crate) fn new(order: Order) -> Self {
        Self { order }
    }
}

#[pymethods]
impl PyOrder {
    #[new]
    pub fn new_order(alphas: u32, alpha: u32, logxir: u32, logxif: u32) -> Self {
        Self::new(Order::new(alphas, alpha, logxir, logxif))
    }
}

/// PyO3 wrapper to :rustdoc:`pineappl::grid::Grid <grid/struct.Grid.html>`
///
/// **Usage**: `yadism`, `pineko`, FKTable interface
#[pyclass]
#[repr(transparent)]
pub struct PyGrid {
    pub(crate) grid: Grid,
}

impl PyGrid {
    pub(crate) fn new(grid: Grid) -> Self {
        Self { grid }
    }
}

#[pymethods]
impl PyGrid {
    #[new]
    pub fn new_grid(
        lumi: Vec<PyRef<PyLumiEntry>>,
        orders: Vec<PyRef<PyOrder>>,
        bin_limits: Vec<f64>,
        subgrid_params: PySubgridParams,
    ) -> Self {
        Self::new(Grid::new(
            lumi.iter().map(|pyl| pyl.lumi_entry.clone()).collect(),
            orders.iter().map(|pyo| pyo.order.clone()).collect(),
            bin_limits,
            subgrid_params.subgrid_params,
        ))
    }

    /// Add a point to the grid.
    ///
    /// Parameters
    /// ----------
    ///     x1 : float
    ///         first momentum fraction
    ///     x2 : float
    ///         second momentum fraction
    ///     q2 : float
    ///         process scale
    ///     order : int
    ///         order index
    ///     observable : float
    ///         reference point (to be binned)
    ///     lumi : int
    ///         luminosity index
    ///     weight : float
    ///         cross section weight
    pub fn fill(
        &mut self,
        x1: f64,
        x2: f64,
        q2: f64,
        order: usize,
        observable: f64,
        lumi: usize,
        weight: f64,
    ) {
        self.grid.fill(
            order,
            observable,
            lumi,
            &Ntuple::<f64> { x1, x2, q2, weight },
        );
    }

    /// Set a metadata key-value pair in the grid.
    ///
    /// **Usage:** `yadism`
    ///
    /// Parameters
    /// ----------
    ///     key : str
    ///         key
    ///     value : str
    ///         value
    pub fn set_key_value(&mut self, key: &str, value: &str) {
        self.grid.set_key_value(key, value);
    }

    /// Set a subgrid.
    ///
    /// **Usage:** `yadism`
    pub fn set_subgrid(&mut self, order: usize, bin: usize, lumi: usize, subgrid: PySubgridEnum) {
        self.grid
            .set_subgrid(order, bin, lumi, subgrid.subgrid_enum);
    }

    /// Set the normalizations.
    ///
    /// **Usage:** `yadism`
    ///
    /// Parameters
    /// ----------
    ///     remapper: BinRemapper
    ///         Remapper object
    pub fn set_remapper(&mut self, remapper: PyBinRemapper) {
        self.grid.set_remapper(remapper.bin_remapper).unwrap();
    }

    /// Extract the necessary informations for EKO.
    ///
    /// **Usage:** `pineko`
    ///
    /// Returns
    /// -------
    ///     x_grid: list(float)
    ///         interpolation grid
    ///     muf2_grid : list(float)
    ///         factorization scale list
    pub fn eko_info(&self) -> (Vec<f64>, Vec<f64>) {
        let EkoInfo { x_grid, muf2_grid } = self.grid.eko_info().unwrap();
        (x_grid, muf2_grid)
    }

    /// Convolute grid with pdf.
    ///
    /// **Usage:** `pineko`
    ///
    /// Parameters
    /// ----------
    ///     xfx1 : callable
    ///         lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x_1`-grid
    ///     xfx2 : callable
    ///         lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x_2`-grid
    ///     alphas : callable
    ///         lhapdf like callable with arguments `Q2` returning :math:`\alpha_s`
    ///
    /// Returns
    /// -------
    ///     list(float) :
    ///         cross sections for all bins
    pub fn convolute(&self, xfx1: &PyAny, xfx2: &PyAny, alphas: &PyAny) -> Vec<f64> {
        self.grid.convolute(
            &|id, x, q2| f64::extract(xfx1.call1((id, x, q2)).unwrap()).unwrap(),
            &|id, x, q2| f64::extract(xfx2.call1((id, x, q2)).unwrap()).unwrap(),
            &|q2| f64::extract(alphas.call1((q2,)).unwrap()).unwrap(),
            &[],
            &[],
            &[],
            &[(1.0, 1.0)],
        )
    }

    /// Convolute with with an evolution operator.
    ///
    /// **Usage:** `pineko`
    ///
    /// Parameters
    /// ----------
    ///     q2 : float
    ///         reference scale
    ///     alphas : list(float)
    ///         list with :math:`\alpha_s(Q2)` for the process scales
    ///     pids : list(int)
    ///         sorting of the particles in the tensor
    ///     x_grid : list(float)
    ///         interpolation grid
    ///     q2_grid : list(float)
    ///         list of process scales
    ///     operator_flattened : list(float)
    ///         evolution tensor as a flat list
    ///     operator_shape : list(int)
    ///         shape of the evolution tensor
    ///
    /// Returns
    /// -------
    ///     PyFkTable :
    ///         produced FK table
    pub fn convolute_eko(
        &self,
        q2: f64,
        alphas: Vec<f64>,
        pids: Vec<i32>,
        x_grid: Vec<f64>,
        q2_grid: Vec<f64>,
        operator_flattened: Vec<f64>,
        operator_shape: Vec<usize>,
    ) -> PyFkTable {
        let operator = Array::from_shape_vec(operator_shape, operator_flattened).unwrap();
        let evolved_grid = self
            .grid
            .convolute_eko(
                q2,
                &alphas,
                (1., 1.),
                &pids,
                x_grid,
                q2_grid,
                operator.into_dimensionality::<Ix5>().unwrap(),
            )
            .unwrap();
        PyFkTable {
            fk_table: evolved_grid,
        }
    }

    /// Load grid from file.
    ///
    /// **Usage:** `pineko`, FKTable generation
    ///
    /// Parameters
    /// ----------
    ///     path : str
    ///         file path
    ///
    /// Returns
    /// -------
    ///     PyGrid :
    ///         grid
    #[staticmethod]
    pub fn read(path: &str) -> Self {
        Self::new(Grid::read(BufReader::new(File::open(path).unwrap())).unwrap())
    }

    /// Write grid to file.
    ///
    /// **Usage:** `yadism`
    ///
    /// Parameters
    /// ----------
    ///     path : str
    ///         file path
    pub fn write(&self, path: &str) {
        self.grid.write(File::create(path).unwrap()).unwrap();
    }

    /// Optimize grid content.
    ///
    /// **Usage:** `yadism`
    pub fn optimize(&mut self) {
        self.grid.optimize();
    }

    /// Extract the number of dimensions for bins.
    ///
    /// **Usage:** `pineko`
    ///
    /// E.g.: two differential cross-sections will return 2.
    ///
    /// Returns
    /// -------
    ///     int :
    ///         bin dimension
    pub fn bin_dimensions(&self) -> usize {
        self.grid.bin_info().dimensions()
    }

    /// Extract the left edges of a specific bin dimension.
    ///
    /// **Usage:** `pineko`
    ///
    /// Parameters
    /// ----------
    ///     dimension : int
    ///         bin dimension
    ///
    /// Returns
    /// -------
    ///     list(float) :
    ///         left edges of bins
    pub fn bin_left(&self, dimension: usize) -> Vec<f64> {
        self.grid.bin_info().left(dimension)
    }

    /// Extract the right edges of a specific bin dimension.
    ///
    /// **Usage:** `pineko`
    ///
    /// Parameters
    /// ----------
    ///     dimension : int
    ///         bin dimension
    ///
    /// Returns
    /// -------
    ///     list(float) :
    ///         right edges of bins
    pub fn bin_right(&self, dimension: usize) -> Vec<f64> {
        self.grid.bin_info().right(dimension)
    }
}
