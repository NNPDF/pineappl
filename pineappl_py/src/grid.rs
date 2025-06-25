//! Grid interface.

use super::boc::{PyBin, PyBinsWithFillLimits, PyChannel, PyKinematics, PyOrder, PyScales};
use super::convolutions::PyConv;
use super::evolution::{PyEvolveInfo, PyOperatorSliceInfo};
use super::fk_table::PyFkTable;
use super::interpolation::PyInterp;
use super::pids::PyPidBasis;
use super::subgrid::PySubgridEnum;
use itertools::izip;
use ndarray::CowArray;
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray4};
use pineappl::boc::Kinematics;
use pineappl::convolutions::ConvolutionCache;
use pineappl::evolution::AlphasTable;
use pineappl::grid::Grid;
use pineappl::pids::PidBasis;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

/// PyO3 wrapper to :rustdoc:`pineappl::grid::Grid <grid/struct.Grid.html>`.
#[pyclass(name = "Grid", subclass)]
#[repr(transparent)]
#[derive(Clone)]
pub struct PyGrid {
    pub(crate) grid: Grid,
}

#[pymethods]
impl PyGrid {
    /// Constructor to instantiate a new PineAPPL Grid.
    ///
    /// # Panics
    ///
    /// Panics when the number of PIDs in `channels` is not equal to `convolutions.len()`, or
    /// `interps` and `kinematics` have different lengths or if `kinematics` are not compatible
    /// with `scales`.
    ///
    /// Parameters
    /// ----------
    /// pid_basis : PidBasis
    ///     choice of basis which can be `Evol` or `Pdg`
    /// channels : list(PyChannel)
    ///     channels
    /// orders : list(PyOrder)
    ///     orders
    /// bins : PyBinsWithFillLimits
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
        bins: PyRef<PyBinsWithFillLimits>,
        convolutions: Vec<PyRef<PyConv>>,
        interpolations: Vec<PyRef<PyInterp>>,
        kinematics: Vec<PyRef<PyKinematics>>,
        scale_funcs: PyRef<PyScales>,
    ) -> Self {
        Self {
            grid: Grid::new(
                bins.bins_fill_limits.clone(),
                orders.into_iter().map(|pyo| pyo.order.clone()).collect(),
                channels.into_iter().map(|pyc| pyc.entry.clone()).collect(),
                pid_basis.into(),
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
    pub fn fill_all_channels(
        &mut self,
        order: usize,
        observable: f64,
        ntuple: Vec<f64>,
        weights: Vec<f64>,
    ) {
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

    /// Get the bin specifications for this grid.
    ///
    /// Returns
    /// -------
    /// PyBinsWithFillLimits:
    ///     a `PyBinsWithFillLimits` object with containing the bin specifications
    #[must_use]
    pub fn bwfl(&self) -> PyBinsWithFillLimits {
        PyBinsWithFillLimits {
            bins_fill_limits: self.grid.bwfl().clone(),
        }
    }

    /// Set the bin specifications for this grid.
    ///
    /// # Errors
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// specs: PyBinsWithFillLimits
    ///     the object to define the bin specs
    pub fn set_bwfl(&mut self, specs: PyBinsWithFillLimits) -> PyResult<()> {
        match self.grid.set_bwfl(specs.bins_fill_limits) {
            Ok(()) => Ok(()),
            Err(msg) => Err(PyValueError::new_err(format!("{msg:?}"))),
        }
    }

    /// Return the number of bins.
    ///
    /// Returns
    /// -------
    /// int :
    ///     Number of bins
    #[must_use]
    pub fn bins(&self) -> usize {
        self.grid.bwfl().len()
    }

    /// Get the length/size of the bins.
    ///
    /// Returns
    /// -------
    /// int:
    ///     the size/length of the bins
    pub fn len(&mut self) -> usize {
        self.grid.bwfl().len()
    }

    /// Get the dimension of the bins that define the observable.
    ///
    /// Returns
    /// -------
    /// int:
    ///     the dimention of the bins
    pub fn bin_dimensions(&mut self) -> usize {
        self.grid.bwfl().dimensions()
    }

    /// Get the limits/edges of all the bins.
    ///
    /// Returns
    /// -------
    /// list(list(float)):
    ///     limits/edges of the bins with shape (n_bins, n_dimension, 2)
    #[must_use]
    pub fn bin_limits(&self) -> Vec<Vec<(f64, f64)>> {
        self.grid
            .bwfl()
            .bins()
            .iter()
            .map(|b| b.limits().to_vec())
            .collect()
    }

    /// Extract the normalizations for each bin.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    ///     bin normalizations
    #[must_use]
    pub fn bin_normalizations<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.grid.bwfl().normalizations().into_pyarray(py)
    }

    /// Get the bin slices for this grid.
    ///
    /// Returns
    /// -------
    /// list(list(int)):
    ///     a list of indices representing the slices
    pub fn bin_slices(&mut self) -> Vec<Vec<usize>> {
        self.grid
            .bwfl()
            .slices()
            .iter()
            .map(|slice| slice.clone().collect())
            .collect()
    }

    /// Get the removed bin using the index for this grid.
    ///
    /// Parameters
    /// ----------
    /// index: int
    ///     index of the bin to be removed
    /// Returns
    /// -------
    /// Bin:
    ///     a `Bin` object from the removed index
    pub fn removed_bin(&mut self, index: usize) -> PyBin {
        PyBin {
            bin: self.grid.bwfl().clone().remove(index),
        }
    }

    /// Set a metadata key-value pair in the grid.
    ///
    /// # Panics
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// key : str
    ///     key
    /// value : str
    ///     value
    pub fn set_metadata(&mut self, key: &str, value: &str) {
        self.grid
            .metadata_mut()
            .insert(key.to_owned(), value.to_owned());
    }

    /// Get metadata values stored in the grid.
    ///
    ///
    /// Returns
    /// -------
    /// dict :
    ///     key, value map
    #[getter]
    #[must_use]
    pub fn metadata(&self) -> BTreeMap<String, String> {
        self.grid.metadata().clone()
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

        let mut convolution_cache = ConvolutionCache::new(
            pdg_convs.into_iter().map(|pdg| pdg.conv.clone()).collect(),
            xfx_funcs
                .iter_mut()
                .map(|fx| fx as &mut dyn FnMut(i32, f64, f64) -> f64)
                .collect(),
            &mut alphas,
        );

        self.grid
            .convolve(
                &mut convolution_cache,
                &order_mask.unwrap_or_default(),
                &bin_indices.unwrap_or_default(),
                &channel_mask.unwrap_or_default(),
                &xi.unwrap_or_else(|| vec![(1.0, 1.0, 1.0)]),
            )
            .into_pyarray(py)
    }

    /// Collect information for convolution with an evolution operator.
    ///
    /// # Panics
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
    /// # Panics
    ///
    /// Panics when the operators returned by either slice have different dimensions than promised
    /// by the corresponding [`OperatorSliceInfo`].
    ///
    /// # Errors
    ///
    /// Raises error if either the `operator` or its `info` is incompatible with the Grid.
    /// Another error is raised if the iterator from `slices` themselves return an error.
    ///
    /// Parameters
    /// ----------
    /// slices : list(Generator(tuple(PyOperatorSliceInfo, PyReadOnlyArray4)))
    ///     list of EKOs where each element is in turn a list of (PyOperatorSliceInfo, 4D array)
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
    pub fn evolve(
        &self,
        slices: Vec<Bound<PyAny>>,
        order_mask: Vec<bool>,
        xi: (f64, f64, f64),
        ren1: Vec<f64>,
        alphas: Vec<f64>,
    ) -> PyResult<PyFkTable> {
        Ok(self
            .grid
            .evolve(
                slices
                    .into_iter()
                    .map(|subslice| {
                        // create lazy iterators from Python object
                        subslice.try_iter().unwrap().map(|item| {
                            let item = item.unwrap();
                            let op_tuple = item.downcast::<PyTuple>().unwrap();
                            let info: PyOperatorSliceInfo =
                                op_tuple.get_item(0).unwrap().extract().unwrap();
                            let op: PyReadonlyArray4<f64> =
                                op_tuple.get_item(1).unwrap().extract().unwrap();

                            Ok::<_, std::io::Error>((
                                info.info,
                                CowArray::from(op.as_array().to_owned()),
                            ))
                        })
                    })
                    .collect(),
                &order_mask,
                xi,
                &AlphasTable { ren1, alphas },
            )
            .map(|fk_table| PyFkTable { fk_table })
            .unwrap())
    }

    /// Load from file.
    ///
    /// # Panics
    ///
    /// Panics if the grid specified by the path is non-existent.
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
    /// Panics if the specified path to write the grid is non-existent or requires permission.
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
    /// Panics if the specified path to write the grid is non-existent or requires permission.
    ///
    /// Parameters
    /// ----------
    /// path : str
    ///     file path
    pub fn write_lz4(&self, path: PathBuf) {
        self.grid.write_lz4(File::create(path).unwrap()).unwrap();
    }

    /// Return the convention by which the channels' PIDS are encoded.
    #[getter]
    #[must_use]
    pub const fn pid_basis(&self) -> PyPidBasis {
        match self.grid.pid_basis() {
            PidBasis::Pdg => PyPidBasis::Pdg,
            PidBasis::Evol => PyPidBasis::Evol,
        }
    }

    /// Return the convention by which the Kinematics are encoded.
    #[getter]
    #[must_use]
    pub fn kinematics(&self) -> Vec<PyKinematics> {
        self.grid
            .kinematics()
            .iter()
            .map(|&kin| match kin {
                Kinematics::X(v) => PyKinematics::X(v),
                Kinematics::Scale(v) => PyKinematics::Scale(v),
            })
            .collect()
    }

    /// Return the convention by which the Scales are encoded.
    #[getter]
    #[must_use]
    pub fn scales(&self) -> PyScales {
        PyScales {
            scales: self.grid.scales().clone(),
        }
    }

    /// Optimize the contents of the Grid.
    pub fn optimize(&mut self) {
        self.grid.optimize();
    }

    /// Merge with another grid.
    ///
    /// # Panics
    /// TODO
    ///
    /// # Errors
    ///
    /// If the bin limits of `self` and `other` are different and if the bin limits of `other` can
    /// not be merged with `self` an error is returned.
    pub fn merge(&mut self, other: Self) -> PyResult<()> {
        match self.grid.merge(other.grid) {
            Ok(()) => Ok(()),
            Err(msg) => Err(PyValueError::new_err(format!("{msg:?}"))),
        }
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
    #[getter]
    #[must_use]
    pub fn convolutions(&self) -> Vec<PyConv> {
        self.grid
            .convolutions()
            .iter()
            .map(|conv| PyConv { conv: conv.clone() })
            .collect()
    }

    /// Get the interpolation specifications for the current grid.
    ///
    /// Returns
    /// list(PyInterp):
    ///     list of interpolation specifications
    #[getter]
    #[must_use]
    pub fn interpolations(&mut self) -> Vec<PyInterp> {
        self.grid
            .interpolations()
            .iter()
            .map(|interp| PyInterp {
                interp: interp.clone(),
            })
            .collect()
    }

    /// Extract channels.
    ///
    /// Returns
    /// -------
    /// list(list(tuple(list[float],int))) :
    ///     channels as tuples (List of PIDs, factor) (multiple tuples can be associated
    ///     to the same contribution)
    #[must_use]
    pub fn channels(&self) -> Vec<Vec<(Vec<i32>, f64)>> {
        self.grid
            .channels()
            .iter()
            .map(|entry| entry.entry().to_vec())
            .collect()
    }

    /// Extract the factors from all the channels.
    ///
    /// Returns
    /// -------
    /// list(float) :
    ///     list containing the factor values
    #[must_use]
    pub fn channels_factors(&self) -> Vec<f64> {
        self.grid
            .channels()
            .iter()
            .flat_map(|entry| entry.entry().iter().map(|(_, f)| *f))
            .collect()
    }

    /// Deduplicate channels
    ///
    /// Parameters
    /// ----------
    /// ulps: i64
    ///    value of the tolerance to be used
    pub fn dedup_channels(&mut self, ulps: i64) {
        self.grid.dedup_channels(ulps);
    }

    /// Rotate the Grid into the specified basis.
    ///
    /// Parameters
    /// ----------
    /// pid_basis: PyPidBasis
    ///     PID basis of the resulting Grid
    pub fn rotate_pid_basis(&mut self, pid_basis: PyPidBasis) {
        self.grid.rotate_pid_basis(pid_basis.into());
    }

    /// Merge the factors of all the channels.
    pub fn merge_channel_factors(&mut self) {
        self.grid.merge_channel_factors();
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
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// factors : list[float]
    ///     bin-dependent factors by which to scale
    pub fn scale_by_bin(&mut self, factors: Vec<f64>) {
        self.grid.scale_by_bin(&factors);
    }

    /// Scale subgrids by order.
    ///
    /// # Panics
    /// TODO
    ///
    /// Parameters
    /// ----------
    /// alphas : float
    ///     value of the strong coupling constant
    /// alpha : float
    ///     value of the electroweak constant
    /// logxir : float
    ///     value of the renormalization scale
    /// logxif : float
    ///     value of the factorization scale
    /// logxia : float
    ///     value of the fragmentation scale
    pub fn scale_by_order(
        &mut self,
        alphas: f64,
        alpha: f64,
        logxir: f64,
        logxif: f64,
        logxia: f64,
        global_factor: f64,
    ) {
        self.grid
            .scale_by_order(alphas, alpha, logxir, logxif, logxia, global_factor);
    }

    /// Delete orders with the corresponding `order_indices`. Repeated indices and indices larger
    /// or equal than the number of orders are ignored.
    ///
    /// Parameters
    /// ----------
    /// order_indices : list[int]
    ///     list of indices of orders to be removed
    pub fn delete_orders(&mut self, order_indices: Vec<usize>) {
        self.grid.delete_orders(&order_indices);
    }

    /// Delete bins.
    ///
    /// # Panics
    /// TODO
    ///
    /// Repeated bins and those exceeding the length are ignored.
    ///
    /// Parameters
    /// ----------
    /// bin_indices : list[int]
    ///     list of indices of bins to be removed
    pub fn delete_bins(&mut self, bin_indices: Vec<usize>) {
        self.grid.delete_bins(&bin_indices);
    }

    /// Deletes channels with the corresponding `channel_indices`. Repeated indices and indices
    /// larger or equal than the number of channels are ignored.
    ///
    /// Parameters
    /// ----------
    /// bin_indices : list[int]
    ///     list of indices of bins to be removed
    pub fn delete_channels(&mut self, channel_indices: Vec<usize>) {
        self.grid.delete_channels(&channel_indices);
    }

    /// Splits the grid such that each channel contains only a single tuple of PIDs.
    pub fn split_channels(&mut self) {
        self.grid.split_channels();
    }
}

/// Register submodule in parent.
///
/// # Errors
///
/// Raises an error if (sub)module is not found.
pub fn register(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "grid")?;
    m.setattr(pyo3::intern!(m.py(), "__doc__"), "Grid interface.")?;
    pyo3::py_run!(
        parent_module.py(),
        m,
        "import sys; sys.modules['pineappl.grid'] = m"
    );
    m.add_class::<PyGrid>()?;
    parent_module.add_submodule(&m)
}
