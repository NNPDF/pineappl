//! Grid interface.

use super::bin::PyBinRemapper;
use super::boc::{PyChannel, PyKinematics, PyOrder, PyScales};
use super::convolutions::PyConv;
use super::evolution::{PyEvolveInfo, PyOperatorSliceInfo};
use super::fk_table::PyFkTable;
use super::interpolation::PyInterp;
use super::pids::PyPidBasis;
use super::subgrid::PySubgridEnum;
use itertools::izip;
use ndarray::CowArray;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray4};
use pineappl::convolutions::ConvolutionCache;
use pineappl::evolution::AlphasTable;
use pineappl::grid::Grid;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyIterator;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

/// PyO3 wrapper to :rustdoc:`pineappl::grid::Grid <grid/struct.Grid.html>`.
#[pyclass(name = "Grid")]
#[repr(transparent)]
#[derive(Clone)]
pub struct PyGrid {
    pub(crate) grid: Grid,
}

#[pymethods]
impl PyGrid {
    /// Constructor to instantiate a new Grid.
    ///
    /// TODO: Exposes `Scales`
    ///
    /// # Panics
    ///
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// pid_basis : PidBasis
    ///     choice of basis which can be `Evol` or `Pdg`
    /// channels : list(PyChannel)
    ///     channels
    /// orders : list(PyOrder)
    ///     orders
    /// bin_limits : list(float)
    ///     bin configurations
    /// convolutions : list(PyConv)
    ///     contains the types of convolution
    /// interpolations : list(PyInterp)
    ///     types of interpolations required by each kinematic
    /// kinematics : list(PyKinematics)
    ///     list of kinematics
    /// scale_funcs : PyScales
    ///     `Scales` object
    #[new]
    #[must_use]
    pub fn new_grid(
        pid_basis: PyPidBasis,
        channels: Vec<PyRef<PyChannel>>,
        orders: Vec<PyRef<PyOrder>>,
        bin_limits: Vec<f64>,
        convolutions: Vec<PyRef<PyConv>>,
        interpolations: Vec<PyRef<PyInterp>>,
        kinematics: Vec<PyRef<PyKinematics>>,
        scale_funcs: PyRef<PyScales>,
    ) -> Self {
        Self {
            grid: Grid::new(
                pid_basis.into(),
                channels.into_iter().map(|pyc| pyc.entry.clone()).collect(),
                orders.into_iter().map(|pyo| pyo.order.clone()).collect(),
                bin_limits,
                convolutions
                    .into_iter()
                    .map(|pyx| pyx.conv.clone())
                    .collect(),
                interpolations
                    .into_iter()
                    .map(|pyi| pyi.interp.clone())
                    .collect(),
                kinematics
                    .into_iter()
                    .map(|pyk| pyk.clone().into())
                    .collect(),
                scale_funcs.scales.clone(),
            ),
        }
    }

    /// Add a point to the grid.
    ///
    /// Parameters
    /// ----------
    /// order : int
    ///     order index
    /// observable : float
    ///     reference point (to be binned)
    /// channel : int
    ///     channel index
    /// ntuple: list(float)
    ///     list containing information on kinematics
    /// weight : float
    ///     cross section weight
    pub fn fill(
        &mut self,
        order: usize,
        observable: f64,
        channel: usize,
        ntuple: Vec<f64>,
        weight: f64,
    ) {
        self.grid.fill(order, observable, channel, &ntuple, weight);
    }

    /// Retrieve a subgrid.
    #[must_use]
    pub fn subgrid(&self, order: usize, bin: usize, channel: usize) -> PySubgridEnum {
        PySubgridEnum {
            subgrid_enum: self.grid.subgrids()[[order, bin, channel]].clone(),
        }
    }

    /// Add an array to the grid.
    ///
    /// Useful to avoid multiple python calls, leading to performance improvement.
    ///
    /// Parameters
    /// ----------
    /// order : int
    ///     order index
    /// observables : list(float)
    ///     list of reference point (to be binned)
    /// channel : int
    ///     channel index
    /// ntuples: list(list(float))
    ///     list of `ntuple` kinematics
    /// weights : np.array(float)
    ///     cross section weight for all events
    pub fn fill_array(
        &mut self,
        order: usize,
        observables: Vec<f64>,
        channel: usize,
        ntuples: Vec<Vec<f64>>,
        weights: Vec<f64>,
    ) {
        for (ntuple, &observable, &weight) in
            izip!(ntuples.iter(), observables.iter(), weights.iter())
        {
            self.grid.fill(order, observable, channel, ntuple, weight);
        }
    }

    /// Add a point to the grid for all channels.
    ///
    /// Parameters
    /// ----------
    /// order : int
    ///     order index
    /// observable : float
    ///     reference point (to be binned)
    /// ntuple: list(float)
    ///     list containing information on kinematics
    /// weights : np.array(float)
    ///     cross section weights, one for each channels
    pub fn fill_all(&mut self, order: usize, observable: f64, ntuple: Vec<f64>, weights: Vec<f64>) {
        for (channel, &weight) in weights.iter().enumerate() {
            self.grid.fill(order, observable, channel, &ntuple, weight);
        }
    }

    /// Set a subgrid.
    ///
    /// Parameters
    /// ----------
    /// order : int
    ///     order index
    /// bin : int
    ///     bin index
    /// channel : int
    ///     channel index
    /// subgrid : PySubgridEnum
    ///     subgrid object
    pub fn set_subgrid(
        &mut self,
        order: usize,
        bin: usize,
        channel: usize,
        subgrid: PySubgridEnum,
    ) {
        self.grid.subgrids_mut()[[order, bin, channel]] = subgrid.subgrid_enum;
    }

    /// Set the bin normalizations.
    ///
    /// # Panics
    ///
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// remapper: BinRemapper
    ///     Remapper object
    pub fn set_remapper(&mut self, remapper: PyBinRemapper) {
        self.grid.set_remapper(remapper.bin_remapper).unwrap();
    }

    /// Convolve the grid with as many distributions.
    ///
    /// # Panics
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// pdg_convs : list(PyConv)
    ///     list containing the types of convolutions and PID
    /// xfxs : list(callable)
    ///     list of lhapdf-like callable with arguments `pid, x, Q2` returning x*pdf
    /// alphas : callable
    ///     lhapdf like callable with arguments `Q2` returning :math:`\alpha_s`
    /// order_mask : numpy.ndarray(bool)
    ///     Mask for selecting specific orders. The value `True` means the corresponding order
    ///     is included. An empty list corresponds to all orders being enabled.
    /// bin_indices : numpy.ndarray(int)
    ///     A list with the indices of the corresponding bins that should be calculated. An
    ///     empty list means that all bins should be calculated.
    /// channel_mask : numpy.ndarray(bool)
    ///     Mask for selecting specific channels. The value `True` means the
    ///     corresponding channel is included. An empty list corresponds to all channels being
    ///     enabled.
    /// xi : list((float, float))
    ///     A list with the scale variation factors that should be used to calculate
    ///     scale-varied results. The first entry of a tuple corresponds to the variation of
    ///     the renormalization scale, the second entry to the variation of the factorization
    ///     scale. If only results for the central scale are need the list should contain
    ///     `(1.0, 1.0)`.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray(float) :
    ///     cross sections for all bins, for each scale-variation tuple (first all bins, then
    ///     the scale variation)
    #[must_use]
    #[pyo3(signature = (pdg_convs, xfxs, alphas, order_mask = None, bin_indices = None, channel_mask = None, xi = None))]
    pub fn convolve<'py>(
        &self,
        pdg_convs: Vec<PyRef<PyConv>>,
        xfxs: Vec<PyObject>,
        alphas: PyObject,
        order_mask: Option<Vec<bool>>,
        bin_indices: Option<Vec<usize>>,
        channel_mask: Option<Vec<bool>>,
        xi: Option<Vec<(f64, f64, f64)>>,
        py: Python<'py>,
    ) -> Bound<'py, PyArray1<f64>> {
        let mut alphas = |q2: f64| {
            let result: f64 = alphas.call1(py, (q2,)).unwrap().extract(py).unwrap();
            result
        };

        let mut xfx_funcs: Vec<_> = xfxs
            .iter()
            .map(|xfx| {
                move |id: i32, x: f64, q2: f64| {
                    xfx.call1(py, (id, x, q2)).unwrap().extract(py).unwrap()
                }
            })
            .collect();

        let mut lumi_cache = ConvolutionCache::new(
            pdg_convs.into_iter().map(|pdg| pdg.conv.clone()).collect(),
            xfx_funcs
                .iter_mut()
                .map(|fx| fx as &mut dyn FnMut(i32, f64, f64) -> f64)
                .collect(),
            &mut alphas,
        );

        self.grid
            .convolve(
                &mut lumi_cache,
                &order_mask.unwrap_or_default(),
                &bin_indices.unwrap_or_default(),
                &channel_mask.unwrap_or_default(),
                &xi.unwrap_or_else(|| vec![(1.0, 1.0, 0.0)]),
            )
            .into_pyarray_bound(py)
    }

    /// Collect information for convolution with an evolution operator.
    ///
    /// # Panics
    ///
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// order_mask : numpy.ndarray(bool)
    ///     boolean mask to activate orders
    ///
    /// Returns
    /// -------
    /// PyEvolveInfo :
    ///     evolution informations
    #[must_use]
    pub fn evolve_info(&self, order_mask: Vec<bool>) -> PyEvolveInfo {
        PyEvolveInfo {
            evolve_info: self.grid.evolve_info(order_mask.as_slice()),
        }
    }

    /// Evolve the grid with as many EKOs as Convolutions.
    ///
    /// TODO: Expose `slices` to be a vector!!!
    ///
    /// # Panics
    /// TODO
    ///
    /// # Errors
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// slices : Iterable
    ///     list of (PyOperatorSliceInfo, 5D array) describing each convolution
    /// order_mask : numpy.ndarray(bool)
    ///     boolean mask to activate orders
    /// xi : (float, float)
    ///     factorization and renormalization variation
    /// ren1 : numpy.ndarray(float)
    ///     list of renormalization scales
    /// alphas : numpy.ndarray(float)
    ///     list with :math:`\alpha_s(Q2)` for the process scales
    ///
    /// Returns
    /// -------
    /// PyFkTable :
    ///     produced FK table
    pub fn evolve<'py>(
        &self,
        slices: &Bound<'py, PyIterator>,
        order_mask: Vec<bool>,
        xi: (f64, f64, f64),
        ren1: Vec<f64>,
        alphas: Vec<f64>,
    ) -> PyResult<PyFkTable> {
        Ok(self
            .grid
            .evolve(
                vec![slices.into_iter().map(|slice| {
                    let (info, op) = slice
                        .unwrap()
                        .extract::<(PyOperatorSliceInfo, PyReadonlyArray4<f64>)>()
                        .unwrap();
                    Ok::<_, std::io::Error>((
                        info.info,
                        // TODO: avoid copying
                        CowArray::from(op.as_array().to_owned()),
                    ))
                })],
                &order_mask,
                xi,
                &AlphasTable { ren1, alphas },
            )
            .map(|fk_table| PyFkTable { fk_table })
            // TODO: avoid unwrap and convert `Result` into `PyResult`
            .unwrap())
    }

    /// Load from file.
    ///
    /// # Panics
    ///
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// path : str
    ///     file path
    ///
    /// Returns
    /// -------
    /// PyGrid :
    ///     grid
    #[must_use]
    #[staticmethod]
    pub fn read(path: PathBuf) -> Self {
        Self {
            grid: Grid::read(BufReader::new(File::open(path).unwrap())).unwrap(),
        }
    }

    /// Write to file.
    ///
    /// # Panics
    ///
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// path : str
    ///     file path
    pub fn write(&self, path: PathBuf) {
        self.grid.write(File::create(path).unwrap()).unwrap();
    }

    /// Write to compressed file.
    ///
    /// # Panics
    ///
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// path : str
    ///     file path
    pub fn write_lz4(&self, path: PathBuf) {
        self.grid.write_lz4(File::create(path).unwrap()).unwrap();
    }

    /// Optimize content.
    pub fn optimize(&mut self) {
        self.grid.optimize();
    }

    /// Merge with another grid.
    ///
    /// # Panics
    ///
    /// TODO
    ///
    /// # Errors
    ///
    /// TODO
    pub fn merge(&mut self, other: Self) -> PyResult<()> {
        match self.grid.merge(other.grid) {
            Ok(()) => Ok(()),
            Err(x) => Err(PyValueError::new_err(format!("{x:?}"))),
        }
    }

    /// Extract the number of dimensions for bins.
    ///
    /// E.g.: two differential cross-sections will return 2.
    ///
    /// Returns
    /// -------
    /// int :
    ///     bin dimension
    #[must_use]
    pub fn bin_dimensions(&self) -> usize {
        self.grid.bin_info().dimensions()
    }

    /// Extract the normalizations for each bin.
    ///
    /// Returns
    /// -------
    /// np.ndarray
    ///     bin normalizations
    #[must_use]
    pub fn bin_normalizations<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.grid.bin_info().normalizations().into_pyarray_bound(py)
    }

    /// Extract the left edges of a specific bin dimension.
    ///
    /// Parameters
    /// ----------
    /// dimension : int
    ///     bin dimension
    ///
    /// Returns
    /// -------
    /// numpy.ndarray(float) :
    ///     left edges of bins
    #[must_use]
    pub fn bin_left<'py>(&self, dimension: usize, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.grid.bin_info().left(dimension).into_pyarray_bound(py)
    }

    /// Extract the right edges of a specific bin dimension.
    ///
    /// Parameters
    /// ----------
    /// dimension : int
    ///     bin dimension
    ///
    /// Returns
    /// -------
    /// numpy.ndarray(float) :
    ///     right edges of bins
    #[must_use]
    pub fn bin_right<'py>(&self, dimension: usize, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.grid.bin_info().right(dimension).into_pyarray_bound(py)
    }

    /// Return the number of bins.
    ///
    /// Returns
    /// -------
    /// int :
    ///     Number of bins
    #[must_use]
    pub fn bins(&self) -> usize {
        self.grid.bin_info().bins()
    }

    /// Extract the available perturbative orders and scale variations.
    ///
    /// Returns
    /// -------
    /// list(PyOrder) :
    ///     list with perturbative orders and scale variations
    #[must_use]
    pub fn orders(&self) -> Vec<PyOrder> {
        self.grid
            .orders()
            .iter()
            .map(|order| PyOrder {
                order: order.clone(),
            })
            .collect()
    }

    /// Get the type(s) of convolution(s) for the current Grid.
    ///
    /// Returns
    /// list(PyConv):
    ///     list of convolution type with the corresponding PIDs
    #[must_use]
    pub fn convolutions(&self) -> Vec<PyConv> {
        self.grid
            .convolutions()
            .iter()
            .map(|conv| PyConv { conv: conv.clone() })
            .collect()
    }

    /// Extract channels.
    ///
    /// Returns
    /// -------
    /// list(list(tuple(float,float,int))) :
    ///     channels as tuples (pid, pid, factor) (multiple tuples can be associated to the same
    ///     contribution)
    #[must_use]
    pub fn channels(&self) -> Vec<Vec<(Vec<i32>, f64)>> {
        self.grid
            .channels()
            .iter()
            .map(|entry| entry.entry().to_vec())
            .collect()
    }

    /// Scale all subgrids.
    ///
    /// Parameters
    /// ----------
    /// factor : float
    ///     scalar factor by which to scale
    pub fn scale(&mut self, factor: f64) {
        self.grid.scale(factor);
    }

    /// Scale subgrids bin by bin.
    ///
    /// # Panics
    ///
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// factors : numpy.ndarray[float]
    ///     bin-dependent factors by which to scale
    pub fn scale_by_bin(&mut self, factors: Vec<f64>) {
        self.grid.scale_by_bin(&factors);
    }

    /// Delete bins.
    ///
    /// # Panics
    ///
    /// TODO
    ///
    /// Repeated bins and those exceeding the length are ignored.
    ///
    /// Parameters
    /// ----------
    /// bin_indices : numpy.ndarray[int]
    ///     list of indices of bins to removed
    pub fn delete_bins(&mut self, bin_indices: Vec<usize>) {
        self.grid.delete_bins(&bin_indices);
    }
}

/// Register submodule in parent.
/// # Errors
///
/// Raises an error if (sub)module is not found.
pub fn register(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new_bound(parent_module.py(), "grid")?;
    m.setattr(pyo3::intern!(m.py(), "__doc__"), "Grid interface.")?;
    pyo3::py_run!(
        parent_module.py(),
        m,
        "import sys; sys.modules['pineappl.grid'] = m"
    );
    m.add_class::<PyGrid>()?;
    m.add_class::<PyOrder>()?;
    parent_module.add_submodule(&m)
}
