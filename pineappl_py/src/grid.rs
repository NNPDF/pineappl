use pineappl::grid::{EkoInfo, Grid, GridAxes, Ntuple, Order};
use pineappl::lumi::LumiCache;

use super::bin::PyBinRemapper;
use super::fk_table::PyFkTable;
use super::lumi::PyLumiEntry;
use super::subgrid::{PySubgridEnum, PySubgridParams};

use ndarray::{Array, Ix5};
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray2};

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;

use pyo3::exceptions::PyValueError;
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

    /// Tuple representation.
    ///
    /// Returns
    /// -------
    ///     alphas : int
    ///         power of :math:`\alpha_s`
    ///     alpha : int
    ///         power of :math:`\alpha`
    ///     logxir : int
    ///         power of :math:` \ln(\xi_r)`
    ///     logxif : int
    ///         power of :math:` \ln(\xi_f)`
    pub fn as_tuple(&self) -> (u32, u32, u32, u32) {
        (
            self.order.alphas,
            self.order.alpha,
            self.order.logxir,
            self.order.logxif,
        )
    }
    /// Return a mask suitable to pass as the `order_mask` parameter of [`Grid::convolute`]. The
    /// selection of `orders` is controlled using the `max_as` and `max_al` parameters, for
    /// instance setting `max_as = 1` and `max_al = 0` selects the LO QCD only, `max_as = 2` and
    /// `max_al = 0` the NLO QCD; setting `max_as = 3` and `max_al = 2` would select all NLOs, and
    /// the NNLO QCD.
    ///
    /// See `pineappl` crate docs for relevant examples
    #[staticmethod]
    pub fn create_mask(orders: Vec<PyRef<Self>>, max_as: u32, max_al: u32) -> Vec<bool> {
        Order::create_mask(
            &orders.iter().map(|o| o.order.clone()).collect::<Vec<_>>(),
            max_as,
            max_al,
        )
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

    /// Add an array to the grid.
    ///
    /// Useful to avoid multiple python calls, leading to performance improvement.
    // TODO: signature
    pub fn fill_with_array(
        &mut self,
        ntuples: PyReadonlyArray2<f64>,
        order: usize,
        observable: f64,
        lumi: usize,
    ) {
        for nt in ntuples.as_array().outer_iter() {
            let (x1, x2, q2, weight) = (nt[0], nt[1], nt[2], nt[3]);
            self.grid.fill(
                order,
                observable,
                lumi,
                &Ntuple::<f64> { x1, x2, q2, weight },
            );
        }
    }

    /// Add a point to the grid for all lumis.
    // TODO: signature
    pub fn fill_all(
        &mut self,
        x1: f64,
        x2: f64,
        q2: f64,
        order: usize,
        observable: f64,
        weights: Vec<f64>,
    ) {
        self.grid.fill_all(
            order,
            observable,
            &Ntuple::<()> {
                x1,
                x2,
                q2,
                weight: (),
            },
            &weights,
        );
    }

    /// Get metadata values stored in the grid.
    ///
    ///
    /// Returns
    /// -------
    ///     dict :
    ///         key, value map
    pub fn key_values(&self) -> HashMap<String, String> {
        self.grid.key_values().unwrap().clone()
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

    /// Retrieve a subgrid.
    ///
    /// **Usage:** `yadism`
    pub fn subgrid(&self, order: usize, bin: usize, lumi: usize) -> PySubgridEnum {
        PySubgridEnum {
            subgrid_enum: self.grid.subgrid(order, bin, lumi).clone(),
        }
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
    ///     pids: list(int)
    ///         particle ids
    ///     muf2_grid : list(float)
    ///         factorization scale list
    pub fn axes(&self) -> (Vec<f64>, Vec<i32>, Vec<f64>) {
        let GridAxes {
            x_grid,
            pids,
            muf2_grid,
        } = self.grid.axes().unwrap();
        (x_grid, pids, muf2_grid)
    }

    /// Convolute grid with pdf.
    ///
    /// **Usage:** `pineko`
    ///
    /// Parameters
    /// ----------
    ///     pdg_id : int
    ///         PDG Monte Carlo ID of the hadronic particle `xfx` is the PDF for
    ///     xfx : callable
    ///         lhapdf like callable with arguments `pid, x, Q2` returning x*pdf for :math:`x`-grid
    ///     alphas : callable
    ///         lhapdf like callable with arguments `Q2` returning :math:`\alpha_s`
    ///     order_mask : list(bool)
    ///         Mask for selecting specific orders. The value `True` means the corresponding order
    ///         is included. An empty list corresponds to all orders being enabled.
    ///     bin_indices : list(int)
    ///         A list with the indices of the corresponding bins that should be calculated. An
    ///         empty list means that all orders should be calculated.
    ///     lumi_mask : list(bool)
    ///         Mask for selecting specific luminosity channels. The value `True` means the
    ///         corresponding channel is included. An empty list corresponds to all channels being
    ///         enabled.
    ///     xi : list((float, float))
    ///         A list with the scale variation factors that should be used to calculate
    ///         scale-varied results. The first entry of a tuple corresponds to the variation of
    ///         the renormalization scale, the second entry to the variation of the factorization
    ///         scale. If only results for the central scale are need the list should contain
    ///         `(1.0, 1.0)`.
    ///
    /// Returns
    /// -------
    ///     list(float) :
    ///         cross sections for all bins, for each scale-variation tuple (first all bins, then
    ///         the scale variation)
    pub fn convolute_with_one(
        &self,
        pdg_id: i32,
        xfx: &PyAny,
        alphas: &PyAny,
        order_mask: Vec<bool>,
        bin_indices: Vec<usize>,
        lumi_mask: Vec<bool>,
        xi: Vec<(f64, f64)>,
    ) -> Vec<f64> {
        let mut xfx = |id, x, q2| f64::extract(xfx.call1((id, x, q2)).unwrap()).unwrap();
        let mut alphas = |q2| f64::extract(alphas.call1((q2,)).unwrap()).unwrap();
        let mut lumi_cache = LumiCache::with_one(pdg_id, &mut xfx, &mut alphas);
        self.grid
            .convolute(&mut lumi_cache, &order_mask, &bin_indices, &lumi_mask, &xi)
    }

    /// Convolute with with an evolution operator.
    ///
    /// **Usage:** `pineko`
    ///
    /// Parameters
    /// ----------
    ///     muf2_0 : float
    ///         reference scale
    ///     alphas : list(float)
    ///         list with :math:`\alpha_s(Q2)` for the process scales
    ///     pids : list(int)
    ///         sorting of the particles in the tensor
    ///     x_grid : list(float)
    ///         interpolation grid
    ///     target_pids : list(int)
    ///         sorting of the particles in the tensor for final FkTable
    ///     target_x_grid : list(float)
    ///         final FKTable interpolation grid
    ///     muf2_grid : list(float)
    ///         list of process scales
    ///     operator_flattened : list(float)
    ///         evolution tensor as a flat list
    ///     operator_shape : list(int)
    ///         shape of the evolution tensor
    ///     additional_metadata : dict(str: str)
    ///         further metadata
    ///
    /// Returns
    /// -------
    ///     PyFkTable :
    ///         produced FK table
    pub fn convolute_eko(
        &self,
        muf2_0: f64,
        alphas: Vec<f64>,
        pids: Vec<i32>,
        x_grid: Vec<f64>,
        target_pids: Vec<i32>,
        target_x_grid: Vec<f64>,
        muf2_grid: Vec<f64>,
        operator_flattened: Vec<f64>,
        operator_shape: Vec<usize>,
        lumi_id_types: String,
        order_mask: Vec<bool>,
    ) -> PyFkTable {
        let operator = Array::from_shape_vec(operator_shape, operator_flattened).unwrap();
        let eko_info = EkoInfo {
            muf2_0,
            alphas,
            xir: 1.,
            xif: 1.,
            target_pids,
            target_x_grid,
            grid_axes: GridAxes {
                x_grid,
                pids,
                muf2_grid,
            },
            lumi_id_types,
        };

        let evolved_grid = self
            .grid
            .convolute_eko(
                operator.into_dimensionality::<Ix5>().unwrap(),
                eko_info,
                &order_mask,
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

    /// Write grid to compressed file.
    ///
    /// Parameters
    /// ----------
    ///     path : str
    ///         file path
    pub fn write_lz4(&self, path: &str) {
        self.grid.write_lz4(File::create(path).unwrap()).unwrap();
    }

    /// Optimize grid content.
    ///
    /// **Usage:** `yadism`
    pub fn optimize(&mut self) {
        self.grid.optimize();
    }

    /// Merge grid with another one, loaded from file
    ///
    /// Note
    /// ----
    /// For a current limitation with the implementation of the bound object `Grid` is not possible
    /// to operate with two `Grid`s in memory, since is not possible to pass a `Grid` by argument
    pub fn merge_from_file(&mut self, path: &str) -> PyResult<()> {
        match self.grid.merge(Self::read(path).grid) {
            Ok(()) => Ok(()),
            Err(x) => Err(PyValueError::new_err(format!("{:?}", x))),
        }
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

    /// Extract the normalizations for each bin.
    ///
    /// **Usage:** `runcards`
    ///
    /// Returns
    /// -------
    ///     np.array
    ///         bin normalizations
    pub fn bin_normalizations<'py>(&self, py: Python<'py>) -> &'py PyArray1<f64> {
        self.grid.bin_info().normalizations().into_pyarray(py)
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

    /// Extract the available perturbative orders and scale variations.
    ///
    /// Returns
    /// -------
    ///     list(PyOrder) :
    ///         list with perturbative orders and scale variations
    pub fn orders(&self) -> Vec<PyOrder> {
        self.grid
            .orders()
            .iter()
            .map(|order| PyOrder {
                order: order.clone(),
            })
            .collect()
    }

    /// Scale all subgrids.
    ///
    /// Parameters
    /// ----------
    /// factor : float
    ///     scalar factor by which scaling
    pub fn scale(&mut self, factor: f64) {
        self.grid.scale(factor);
    }

    /// Delete bins.
    ///
    /// Repeated bins and those exceeding length are ignored.
    ///
    /// Parameters
    /// ----------
    /// bin_indices : list(int)
    ///     list of indices of bins to removed
    pub fn delete_bins(&mut self, bin_indices: Vec<usize>) {
        self.grid.delete_bins(&bin_indices)
    }
}
