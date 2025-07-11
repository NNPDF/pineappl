# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.1.0] - 08/07/2025

### Added

- added a new `V3` metadata reader to the `pineappl evolve` CLI for EKOs
  generated with `v0.15.0` or higher
- C API: added new functions `pineappl_grid_evolve_info_shape`,
  `pineappl_grid_evolve_info`, and `pineappl_grid_evolve` to evolve grids
- C API: added `pineappl_fktable_optimize` to optimize FK Table-like objects
  given an optimization assumption
- added methods `Grid::merge_channel_factors` and `Channel::factor`

### Fixed

- fixed a bug that caused `pineappl export` to fail when called with grid
  having non-trivial factors in their channel definitions

## [1.0.0] - 10/06/2025

PineAPPL 1.0 is a major rewrite from the previous version, allowing grids to
have an arbitrary number of convolutions in the initial (PDFs with a
factorization scale) and in the final (FFs with a fragmentation scale) state.
This required a change in the file format that is used to write out grids, but
the old file format can still be read with this new version.

### Added

- added new method `Grid::delete_orders` and the corresponding switch
  `--delete-orders` in the subcommand `write` of the CLI
- added the switches `--xir` and `--xif`, which allow varying the
  renormalization and factorization scales with a custom factor in the
  subcommand `convolve`
- the CLI now allows the user to mark convolution functions as polarized
  by adding `+p` to its LHAPDF name, as a fragmentation function by adding
  `+f` and both by adding `+pf` or `+fp`
- added switches `--fk-table-fac0` and `--fk-table-frg0` to `pineappl read` to
  read out the squared factorization and fragmentation scales of FK-tables
- C API: added new functions `pineappl_channels_new` and `pineappl_channels_delete`
  to create and delete an instance of `Channels` object
- C API: added `pineappl_grid_channels` and `pineappl_channels_count` to get
  the channel objects for a given grid and their numbers
- C API: added a function `pineappl_channels_combinations` to retrieve the
  number of combinations of channels for a specified entry, and
  `pineappl_channels_entry` to retrieve the channel combination for a given
  entry
- C API: added a new function `pineappl_grid_new2` to create a grid with the
  new features introduced in `v1.0.0`; this includes the support for an
  arbitrary number of initial and final state hadrons
- C API: added new functions to fill grids with an arbitrary number of
  initial and final state hadrons; these include `pineappl_grid_fill2`,
  `pineappl_grid_fill_all2`, and `pineappl_grid_fill_array2`
- C API: added new functions to extract the various properties of a given
  grid; these include `pineappl_grid_conv_types`, `pineappl_grid_convolutions_len`,
  and `pineappl_grid_kinematics_len`
- C API: added a new function `pineappl_grid_convolve` to convolve grids
  with an arbitrary combination of initial and final state hadrons
- C API: added various functions to extract the subgrids of a given grid;
  these include `pineappl_grid_subgrid_shape`, `pineappl_grid_subgrid_node_values`,
  and `pineappl_grid_subgrid_array`

### Changed

- the macro `channel!` now accepts a channel specification that is of the
  format `factor * (pid, ..) + ...`
- Python API: dropped top-level Python interface layer
- Python API: renamed `lumi` to `channel` in PyO3 Python interface. This
  concerns 1) the argument names of `convolute_with_one` and similar functions;
  2) the module `pineappl.lumi` was moved to `pineappl.boc`; 3) the class
  `LumiEntry` was renamed to `Channel`
- Python API: `.into()` needs to be explicitly called on subgrids when calling
  `pineappl.grid.set_subgrid()`
- Python API: replaced `pineappl.grid.PyPidBasis` with
  `pineappl.evolution.PidBasis`
- Python API: replaced `pineappl.grid.PyOperatorSliceInfo` with
  `pineappl.evolution.OperatorSliceInfo`
- Python API: drop all `Py` prefixes, for instance `PyEvolveInfo` was renamed
  to `EvolveInfo`
- by default `pineappl plot` no longer shows a channel breakdown in the panel
  with absolute PDF predictions. However, this feature can be enabled with via
  a new array added at the start of the script
- raised MSRV to 1.80.1
- changed the order of elements in `Grid::fill` of the parameter `ntuple` to
  reflect the ordering of `kinematics` given to `Grid::new`
- renamed the following switches of `pineappl write`: `--remap` to
  `--set-bins`, `--remap-norm-ignore` to `--div-bin-norm-dims` and
  `--remap-norm` to `--mul-bin-norm`. These names should reflect the
  corresponding operations
- renamed the switch `--fktable` to `--fk-table` of `pineappl read`

### Removed

- Python API: removed `pineappl.grid.Grid.create()` and
  `pineappl.fk_table.FkTable.from_grid()` methods; use the constructors
  of the respective class instead
- removed the constructor `Grid::with_subgrid_type`
- removed `Grid::convolve_subgrid` and `--subgrid-pull` from `pineappl plot`
  that was using the method. The CLI subgrid-pull plotting routine only ever
  worked for grids with two convolutions

## [0.8.7] - 22/01/2025

### Added

- added support for Python 3.13 to the Python interface

## [0.8.6] - 18/10/2024

### Fixed

- fixed [Issue #318](https://github.com/NNPDF/pineappl/issues/318) that caused
  fastNLO tables with `NPDFDim = 2` to be incorrectly imported

## [0.8.5] - 07/10/2024

### Fixed

- fixed a bug in `pineappl_applgrid` that lead to linking problems with ROOT
  and `gfortran`

## [0.8.4] - 04/10/2024

### Fixed

- fixed a bug that lead to inconsistent convolution metadata
  (https://github.com/NNPDF/pineappl/issues/316)

## [0.8.3] - 30/08/2024

### Fixed

- fixed the scale-variation label in the plotting script produced by `pineappl
  plot`. Previously this would always show a 7-pt. variation irrespective of
  the parameter passed to `--scales`
- fixed a problem in the evolution when an EKO with 'similar' Q2 slices was
  used to evolve; this caused the Q2 slices of the grids to be evolved several
  times, leading to wrong results

## [0.8.2] - 22/07/2024

### Changed

- changed the implementation of evolution functions so that they are much
  faster and that NaNs are correctly propagated

## [0.8.1] - 18/07/2024

### Added

- added new method `Grid::evolve_with_slice_iter2` which is able to perform
  evolutions with two different EKOs

### Fixed

- fixed CI to build CAPI and CLI

### Changed

- PineAPPL's CAPI requires now GLIBC 2.29 since 0.8.0

## [0.8.0] - 05/07/2024

### Added

- added new type `Convolution`
- added new methods `Grid::convolutions` and `Grid::set_convolution`
- added the function `pineappl_grid_convolve_with_one` and
  `pineappl_grid_convolve_with_two` which replace the deprecated function
  similarly named with `convolute` in CAPI
- added `PidBasis::charge_conjugate` and `PidBasis::guess`
- added `Grid::set_pid_basis` method
- added `Grid::subgrids` and `Grid::subgrids_mut` methods
- added new switch `--conv-fun-uncert-from` to subcommand `plot` to allow
  choosing which convolution function uncertainty should be plotted

### Changed

- moved `Order` and `ParseOrderError` to the new module `boc`
- renamed switch `--split-lumi` of `pineappl write` to `--split-channels`. The
  old switch can still be used
- renamed switch `--lumis` of `pineappl read` to `--channels`. The old switch
  can still be used
- renamed switch `--ignore-lumis` of `pineappl diff` to `--ignore-channels`.
  The old switch can still be used
- renamed `Grid::lumi` to `Grid::channels`, `Grid::split_lumi` to
  `Grid::split_channels`, `Grid::rewrite_lumi` to `Grid::rewrite_channels` and
  `Grid::set_lumis` to `Grid::set_channels`. The term 'channel' is now used
  everywhere instead of 'lumi', 'luminosity function', etc.
- renamed the struct `LumiEntry` to `Channel` and `ParseLumiEntryError` to
  `ParseChannelError`. Both structures have been moved to the module `boc`
- renamed the macro `lumi_entry` to `channel`
- renamed `Grid::set_channels` to `Grid::channels_mut`
- renamed `TryFromGridError::InvalidLumi` to `TryFromGridError::InvalidChannel`
- changed member `lumi_id_types` of `OperatorInfo` and `OperatorSliceInfo` to
  `pid_basis`, which is now of type `PidBasis`
- renamed module `pineappl::lumi` to `pineappl::convolutions`
- renamed switch `--pdf` to `--conv-fun` in the subcommand `uncert`. This
  switch now optionally accepts a list of indices, which determines the
  corresponding convolution function (PDF/FF), for which the uncertainty should
  calculated
- renamed `no_pdf_unc` to `no_conv_fun_unc` in subcommand `plot`

### Removed

- removed support for Python 3.6
- removed the deprecated evolution methods `Grid::axes`, `Grid::convolute_eko`
  and the structs `EkoInfo` and `GridAxes`
- removed methods `Grid::has_pdf1`, `Grid::has_pdf2`, `Grid::initial_state_1`
  and `Grid::initial_state_2`
- removed `pids::charge_conjugate`; this function has been replaced with the
  new function `PidBasis::charge_conjugate`
- removed `pids::determine_lumi_id_types`; this function has been replaced with
  the new function `PidBasis::guess`
- removed `TryFromGridError::MetadataMissing`
- removed `Grid::subgrid` and `Grid::set_subgrid` methods; these functions have
  been replaced with `Grid::subgrids` and `Grid::subgrids_mut`
- removed the switch `--pdf-with-scale-cov` from `pineappl uncert`

## [0.7.4] - 23/05/2024

### Added

- added `Grid::evolve_with_slice_iter`, `AlphasTable` and `OperatorSliceInfo`,
  which define a new interface supporting very large evolution kernels that
  have been introduced in EKO v0.13. This interface will replace `Grid::evolve`
- added `--dont-sort` switch to `pineappl channels`, which displays the channel
  sizes ordered by channel index (instead of channel size)
- added `Grid::rotate_pid_basis` and `pineappl write --rotate-pid-basis`. This
  allows to change the meaning of the used particle IDs, and supported formats
  are PDG MC IDs and the evolution basis
- added `pineappl write --rewrite-order` that lets the user change the
  exponents of each order

### Changed

- changed the official name of the CLI subcommand `convolute` to `convolve`,
  because the latter is the proper verb of 'convolution'. The old name
  `convolute` is now an alias of `convolve`, which means both can be used. The
  methods `Grid::convolute*` are left unchanged and will be renamed in later
  version
- changed switch `--silence-lhapdf` to `--lhapdf-banner` and suppress LHAPDF's
  banners by default, unless `--lhapdf-banner` is given
- `Grid::evolve` has now been marked deprecated
- switched from `lhapdf` to `managed-lhapdf` crate which automatically
  downloads PDF sets when they are needed

### Fixed

- fixed yet another problem that prevent the Python interface for Python 3.6
  from being successfully installed
- fixed `Grid::delete_channels` and its CLI variant `pineappl write
  --delete-channels`. This command wasn't working properly before

## [0.7.3] - 23/02/2024

### Fixed

- fixed a problem that prevent the Python interface for Python 3.6 from being
  successfully installed

## [0.7.2] - 23/02/2024

### Fixed

- fixed problems in the CI that prevented the CAPI for Linux being uploaded

## [0.7.1] - 23/02/2024

### Fixed

- fixed problems in the CI that prevented the previous version from being
  released

## [0.7.0] - 23/02/2024

### Added

- added a new global option `--allow-extrapolation`, which, when present,
  allows LHAPDF to extrapolate PDF values outside its valid region in `x`.
  If this option is not present the PDFs are set to zero, which remains the
  default behavior

### Changed

- the Python interface for MacOS is now shipped separately for the two targets
  (aarch64 and x86_64) instead of a single universal wheel previously
- raised MSRV to 1.70.0
- when calling `BinRemapper::new`, the limits are now checked for overlaps, in
  which case a new error is returned
- changed the type `ParseBinRemapperError` to allow capturing errors from
  `BinRemapper::new`

### Fixed

- fixed the missing generation of CPython 3.7 to 3.10 wheels for MacOS targets

### Removed

- removed the type `SubGrid` that was exported in the CAPI as the type
  `pineappl_subgrid`. This type was not used anymore and was a left-over of the
  changes in the previous version.

## [0.6.3] - 12/12/2023

### Added

- added `Grid::dedup_channels`, the C function `pineappl_grid_dedup_channels`
  and the new switch `--dedup-channels` to the CLI that detect equal subgrids
  and reduce the space requirements of a grid
- added a new panel `double_ratio_pdf` to the `pineappl plot` tool
- added a Python binding for the `Grid::merge` method

### Removed

- removed the functions of `pineappl_grid_export_mu2_slice`,
  `pineappl_grid_nonzero_mu2_slices`, `pineappl_subgrid_delete`,
  `pineappl_grid_replace_and_delete`, `pineappl_subgrid_import_mu2_slice` and
  `pineappl_subgrid_new2`, which were only used in the C++ programs that were
  the predecessors of the `export` and `import` functionality now implemented
  in the CLI

### Fixed

- fixed a bug that caused channels with PID=0 gluons to not evolve when the
  metadata key-value pair `lumi_id_types=pdg_mc_ids` was not present. Now when
  this metadata is not present it is assumed PDG MC IDs are used

## [0.6.2] - 09/10/2023

### Added

- added support for Python 3.6 on Linux and PyPy 3.7, 3.8, 3.9 and 3.10 on
  MacOS and Windows for packages from <https://pypi.org/project/pineappl/>
- added support for scale uncertainties calculated with the covariance method,
  and for combined PDF and scale uncertainties
- added support for installing man pages. See installation instructions for
  more information
- added new method `Grid::optimize_using` and its corresponding C function
  `pineappl_grid_optimize_using`, which optimize a `Grid` like
  `Grid::optimize`, but allow for more fine-grained selection of optimization
  options, which are listed in the new bitflags `GridOptFlags`

### Changed

- `pineappl plot` now produces a more cleaned-up matplotlib script, in which
  the most important information, which a user may whish to change, is at the
  top of the file
- `pineappl pdfunc` was renamed into `pineappl uncert`, which is now also able
  to calculate scale uncertainties, which in turn have been removed from
  `pineappl convolute`
- `pineappl help` now relies on installed man pages

### Fixed

- fixed a bug in the calculation of asymmetries when multiple scales had been
  used
- fixed a bug causing `Grid::evolve` to require more x-grid values than needed;
  this result in errors similar to the following: 'no operator for x =
  0.0018585113621881083 found' where this x-grid value is a point that is used
  in a masked order

## [0.6.1] - 18/07/2023

### Added

- added switch `--no-pdf-unc` to `plot` subcommand to skip the time-consuming
  computation of PDF uncertainties
- added new function `pineappl_grid_merge_bins` to the CAPI. This function
  corresponds to `Grid::merge_bins` and merges a range of bins together into a
  new bin.
- added new function `Grid::split_lumi`, which splits the luminosity such that
  it contains a single combination per partonic channel. This function is
  available through the CAPI via `pineappl_grid_split_lumi` and via the CLI
  through `pineappl write --split-lumi`

### Fixed

- fixed panic when trying to plot DIS observables

## [0.6.0] - 01/06/2023

### Added

- added `--orders` switch to pull to allow the selection of a subset of
  perturbative orders
- added the subcommand `help` to show manpages of for each subcommand of
  `pineappl`
- the switch `--limit` in `pull` now allows the value `0` for faster pull
  computation
- added PDF ratio panel to plot script generated by `pineappl plot`
- added new subcommand `export`, which is able to convert *some* PineAPPL grids
  into the APPLgrid format
- added support for [EKO](https://eko.readthedocs.io/en/latest/)'s new file
  formats (old file format is still supported) in `pineappl evolve`
- `pineappl import` now also imports fastNLO's metadata:
  - `fastNLOTable::GetDimLabels()` is stored as the value of the key
     `x{1,2,...}_labels`
  - `fastNLOTable::GetXSDescr()` is stored as the value of the key `y_label`
  - `fastNLOTable::GetScDescr()` is stored as the value of the key
    `fastnlo_scenario`
- removed switch `--silence-libraries` which did not work properly and relied
  on non-portable code

### Changed

- raised MSRV to 1.64.0
- the switch `--force-positive` must be given at the same level as
  `--silence-lhapdf`
- the CLI subcommands `delete`, `optimize`, `set` and `update` were merged into
  the new `write` subcommand. Also the options `--merge-bins`, `--scale` and
  `--scale-by-order` were merged into `write`
- the CLI subcommand `remap` was merged into `write` as the option `--remap`,
  which expects as a single argument the remapping string. The option  `--norm`
  was renamed to `--remap-norm`, and `--ignore-obs-norm` was renamed to
  `--remap-norm-ignore`. Note that the latter now counts dimensions from zero
  onward, instead from one (old behavior)
- the CLI subcommand `info` and `obl` were merged into the new subcommand
  `read`

### Fixed

- fixed a bug introduced in v0.5.5 that caused `pineappl convolute` to not show
  results for more than two PDF sets
- fixed a bug that caused `pineappl convolute` to ignore a possible `--order`
  parameter for additionally given PDF sets
- fixed a bug that caused `pineappl plot` to show wrong pulls if the central
  PDF set had asymmetric uncertainties

### Removed

- removed `pineappl sum --integrated`. This subcommand did not always work as
  expected, use `pineappl obs --merge-bins` instead

## [0.5.9] - 02/01/2023

### Added

- added new function `pineappl_grid_scale_by_bin` to the CAPI that corresponds
  to `Grid::scale_by_bin`
- added new subcommand `analyze`, which performs various analyses. For time
  being only one analysis is available: `ckf`. This calculates the per-channel
  K factors which can be compared against the per-bin K factors
- added switch `--orders` to `evolve` to allow evolving a subset of orders of
  the full grid
- added new switch `--fktable` to the CLI subcommand `pineappl obl` to detect
  whether a grid is also an FK table
- added new method `Grid::evolve_info` which extracts the information need to
  generate EKOs. This function is less restrictive than `Grid::axes`, which it
  will replace
- `pineappl import` now converts also scale-variation log-grids from
  flexible-scale fastNLO tables

### Fixed

- fixed bug in the calculation of asymmetries when multiple scales were
  requested

## [0.5.8] - 21/10/2022

### Added

- added new method `Grid::evolve`, which will succeed `Grid::convolute_eko` and
  is more robust, faster and easier to understand

### Fixed

- fixed a bug in `LumiCache::setup`, which caused wrong results returned by
  `Grid::convolute`, if `lumi_cache` is reused for a second grid with different
  Q2/x1/x2 parameters

## [0.5.7] - 05/10/2022

### Fixed

- the importer did not depend on the right versions of `pineappl_applgrid` and
  `pineappl_fastnlo`

## [0.5.6] - 04/10/2022

### Added

- added `--dis-pid` switch to `import` subcommand
- added `--no-optimize` switch to `import` subcommand

### Changed

- improved documentation
- changed `--silence-fastnlo` to `--silence-libraries`, which now silences also
  APPLgrid
- improved `import` converter, which now is able to convert most APPLgrids and
  fastNLO tables on ploughshare

### Fixed

- fixed corner cases of the `REMAPPING` string that were handled incorrectly
- made most tests of `pineappl` integration tests, which resolves that problem
  that `cargo build` had to be called before
- fixed problem in `Grid::convolute_eko` that lead to wrong result when
  subgrids had different x-grids

## [0.5.5] - 25/08/2022

### Added

- added support for choosing a specific member instead of the averaging over
  all PDFs for the central results for the `pdfunc`, `plot` and `pull`
  subcommands. For instance, `pineappl pdfunc ... NNPDF40_nnlo_as_01180`
  calculates the central value using the average over all replicas, whereas
  `pineappl pdfunc ... NNPDF40_nnlo_as_01180/0` uses the zeroth member. The
  calculated uncertainties are the same for both
- added support for converting APPLgrids in `import`
- added support plotting asymmetries with `pineappl plot`
- added new function to the C API, `pineappl_grid_clone`, which clones a given
  grid

### Changed

- renamed the switch `--silence-fastnlo` to `--silence-libraries`, which also
  covers APPLgrid. For backwards compatibility the previous name still works

## [0.5.4] - 08/07/2022

### Added

- added new switch `--force-positive` to CLI to cut out negative PDF values
- exposed scale variation parameters to the Python API
- added the C functions `pineappl_grid_key_value` to read out metadata and
  `pineappl_string_delete` to delete string allocated from the previous
  function

## [0.5.3] - 22/06/2022

### Added

- added new switch `--ignore-bin-limits` to `diff` that ignores possibly
  different bin limits. However, the number of total bins must still be the
  same between both grids
- added support for
  [FK tables](https://docs.nnpdf.science/data/th-data-files.html) in the
  subcommand `import`
- added new switches `--digits-abs`, `--digits-rel` and `--digits` to various
  subcommand to influence the number of (fractional) digits printed
- PDF relabeling support has been added. This means that all subcommands
  understand PDF specifications of the type `pdfset=label`, where `pdfset` must
  be an LHAPDF identifier, and results using this set are displayed using
  `label` instead of `pdfset`.
- added bin-dependent rescaling method. To use it through the CLI use `pineappl
  ops --scale-by-bin=`

### Changed

- slightly changed the output of the CLI; indices of orders, bins and
  luminosities are now consistently abbreviated by `o`, `b` and `l`
- changed the output of the CLI to also print the units of the numbers

## [0.5.2] - 29/03/2022

### Added

- added new subcommand `import`, which at this stage converts fastNLO tables to
  PineAPPL grids. This program was previously an example and written in C++ and
  now has been removed.

## [0.5.1] - 01/03/2022

### Added

- added the `--bins` option to the CLI subcommand `sum`, which allows to sum
  bins together
- added the `--fk-table` option to the CLI subcommand `optimize`, which allows
  the optmization of FK tables based on assumptions of heavy-flavor PDFs
- added new subcommand `ops` which collects various modifying operations on
  grids. The switches `--cc1` and `--cc2` charge conjugate the first and second
  PDF, respectively, and charge conjugates the luminosity function
  correspondingly such that the convolutions are unchanged

### Changed

- when running `pineappl convolute ... -s 1` the scale-variation columns are no
  longer shown. The output would be zero, but this doesn't make sense. All
  other values of `-s` are unaffected.
- added a further optimization to `Grid::optimize` that strips empty orders

### Fixed

- fixed `pineappl obl --bins`, which had the wrong column headers if there were
  more than one observable

## [0.5.0] - 11/02/2022

### Added

- added support for DIS processes, which feature only a single hadron in the
  initial state
- added new C API function `pineappl_grid_set_remapper`
- added new subcommand `sum` to sum over bins of a grid
- added new subcommand `pull` to view where the differences between two PDF
  sets are coming from
- added an example program using the C API to convert fastNLO tables (of type
  `fastNLOCoeffAddFix`). Tables of type `fastNLOCoeffAddFlex` are not supported
  yet
- added a new switch `--subgrid-pull` to the `plot` subcommand; this generates
  several plots showing the where the pull between two PDFs comes from in `x1`
  and `x2`, and also in rapidity `y` and invariant mass `M`
- improved the `optimize` method such that it removes all entries of a
  luminosity function that are empty across orders and bins. This required a
  change in `merge`, which now allows the merge of grids with different
  luminosities and orders
- enabled relabeling PDF sets, which means that they can now specified as
  `LHAPDF-set-name=label`, which uses the set `LHAPDF-set-name` to load it from
  LHAPDF, but `label` is the name that appears in the plot legend
- added new options to the `subgrids` command: `--muf`, `--mur`, `--x1` and
  `--x2` to print the corresponding grid values and `--stats` to print
  information about the space requirements of the subgrid type
- added new subcommand `--delete`, which allows to delete bins from a grid
- added new subcommand `obl`, which stands for orders (o), bins (b) and lumis
  (l), that shows information about the all contained subgrids in a grid;
  `pineappl obl -l` will replace `pineappl lumis` in the future
- added options `--orders1`, `--orders2`, `--scale1` and `--scale2` to the
  subcommand `pineappl diff`, which allow the selection of orders and scaling
  of grids that should be compared

### Changed

- The command-line interface `pineappl` no longer prints both the differential
  and integrated cross sections. Instead it either prints the differential
  cross sections, or, if the switch `-i` or `--integrated` is given, the
  integrated cross sections (without bin limits and/or normalizations) are
  printed
- the C API functions `pineappl_subgrid_q2_slice`,
  `pineappl_subgrid_filled_q2_slices` and `pineappl_subgrid_replace_and_delete`
  have been replaced by `pineappl_grid_export_q2_slice`,
  `pineappl_grid_nonzero_q2_slices` and `pineappl_grid_replace_and_delete`,
  respectively.
- the C API function `pineappl_subgrid_fill_q2_slice` has been replaced by
  function `pineappl_subgrid_import_q2_slice`
- the C API function `pineappl_subgrid_new` has been replaced by a function
  with the similar name but different arguments
- changed the output of the `pineappl` subcommand `channels`, `convolute` and
  `pdf_uncertainty`. By default only the differential cross sections are shown
  (integrated numbers divided by bin widths), but the flag `-i` or
  `--integrated` can be given to switch to the integrated numbers, which are
  not divided by bin widths.
- slightly improved file sizes by introducing a type for empty subgrids:
  `EmptySubgridV1`
- vastly improved the output of the `plot` subcommand: bounding boxes are
  properly calculated now, added support for higher-dimensional distributions
- changed the `plot` subcommand such that the legend is put in between panels,
  thereby producing less overlapping elements
- replaced the C API functions `pineappl_subgrid_new` and
  `pineappl_subgrid_import_q2_slice` with `pineappl_subgrid_new2` and
  `pineappl_subgrid_import_mu2_slice`, respectively. The latter support
  independent renormalization and factorization scales
- changed the names of `pineappl_grid_export_q2_slice` and
  `pineappl_grid_nonzero_q2_slices` to `pineappl_grid_export_mu2_slice` and
  `pineappl_grid_nonzero_mu2_slices`, respectively
- replaced C API-based Python interface with PyO3-based one
- the subcommand `subgrids` now does not print empty grids; the old behavior
  can be restored with the new switch `--show-empty`
- `pineappl diff` now behaves differently whenever luminities are different and
  errors out when this is the case. This can be adjusted using
  `--ignore-lumis`

## [0.4.1] - 25/03/2021

### Fixed

- added fallback options to `git_version` that prevented uploading
  `pineappl_capi` and `pineappl_cli` to crates.io

## [0.4.0] - 25/03/2021

### Added

- added access to the contents of Lagrange-interpolation grids
- added more C functions previously missing: `pineappl_grid_bin_limits`,
  `pineappl_grid_lumi`, `pineappl_lumi_combinations`, `pineappl_lumi_count`,
  `pineappl_lumi_entry`
- added `remap` subcommand which allows to change the binning and dimensions of
  a distribution
- added global `--silence_lhapdf` switch to suppress the output of LHAPDF
- added new subgrid type: `LagrangeSparseGrid`, which can be generated using
  the new subcommand `optimize` from existing grids; this data structure uses
  less memory and is a bit faster
- added new subcommand `optimize`, which optimizes the size for each grid, both
  in memory and when written as a file
- added new subgrid type `LagrangeSubgridV2`, which allows for different `x1`
  and `x2` limits, for example for DIS processes
- added new switches to the subcommand `info`: `--get`, `--keys`, `--show`,
  which lets one read a single or all key-values pairs of each grid file.
- added new subcommand `set` which allows to modify the key-value storage
- added new C functions `pineappl_grid_set_key_value` and
  `pineappl_grid_optimize`
- added a new switch `--ignore_orders` to the diff subcommand, which sums over
  all orders and is therefore useful if two grids are compared that should have
  the same contents in principle, but in practice have different orders
- added new subcommand `plot`, which allows to plot the information contained
  in a grid

### Changed

- the order columns of the subcommand `diff` are now properly sorted and do not
  change randomly.
- the subcommand `diff` now shows the differential cross sections of both grids
  with the same number of digits as the subcommands `convolute` and similar.
- changed the default maximum value of Q from 1 TeV to 10 TeV and the number of
  grid points in Q^2 from 30 to 40
- Removed `Grid::bin_limit` and replaced it with `Grid::bin_info`
- in the C API the type `uintptr_t` has been changed to the more common type
  `size_t`
- changed the default `LagrangeSubgrid` type from `V1` to `V2`. This subgrid
  type supports DIS and allows to reconstruct static scales. If a static scale
  has been detected, `optimize` can remove the interpolation in the Q^2
  dimension.

### Fixed

- the subcommand `diff` did not show differences in per cent, although it
  printed a per cent sign. Now it shows relative differences, which is always
  useful, even when the differences are smaller than sub-per mille.

## [0.3.0] - 20/08/2020

### Added

- added the options `--absolute`, `--orders`, and `--lumis` to the subcommand
  `channels`
- added Python interface to the C API, the the folder `wrappers/python` and the
  example in `examples/python-dy-aa`
- added the option `--normalize` to the subcommand `orders`, which can be used
  to specify the orders that should be used to normalize to.

### Fixed

- added missing support for LHAIDs for `pdf_uncertainty`
- fixed a case where merging two grids fails, because the bin limits are
  numerically not exactly the same

## [0.2.0] - 02/08/2020

### Added

- the Lagrange-interpolation grid of PineAPPL now supports fully dynamic
  factorisation/renormalisation scales
- `pineappl` has a new subcommand `pdf_uncertainty` to calculate PDF
  uncertainties
- in `examples/capi-dy-aa` one can find an example how to use the C API
- added method `Grid::with_subgrid_type` that allows the selection of a custom
  Subgrid type
- added `--absolute` switch to the subcommand `convolute` of `pineappl`, which
  shows all values of the scale variation as absolute numbers
- added a PDF and alphas to `Grid::convolute` to speed up convolutions
- added `--orders` switch to the subcommand `pdf_uncertainty`, to caclulate the
  PDF uncertainty only with the specified perturbative orders
- added new subcommand `info` which shows the highest orders in alphas and
  alpha
- added new subcommand `diff` which compares two grids with each other

### Changed

- the parameter `reweight` of `SubgridParams` is set to `true` by default now
  (was `false` before)
- improved the output of `pineappl`, bin widths are printed now

### Fixed

- fixed subcommand `orders` of `pineappl` to normalise each leading order to
  the sum of all leading orders, not only the first
- fixed subcommand `pdf_uncertainty` of `pineappl` to show the correct
  integrated results (bin width was missing)
- fixed calculation of the scale variation

## [0.1.0] - 11/06/2020

- first release

[Unreleased]: https://github.com/NNPDF/pineappl/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/NNPDF/pineappl/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/NNPDF/pineappl/compare/v0.8.2...v1.0.0
[0.8.7]: https://github.com/NNPDF/pineappl/compare/v0.8.6...v0.8.7
[0.8.6]: https://github.com/NNPDF/pineappl/compare/v0.8.5...v0.8.6
[0.8.5]: https://github.com/NNPDF/pineappl/compare/v0.8.4...v0.8.5
[0.8.4]: https://github.com/NNPDF/pineappl/compare/v0.8.3...v0.8.4
[0.8.3]: https://github.com/NNPDF/pineappl/compare/v0.8.2...v0.8.3
[0.8.2]: https://github.com/NNPDF/pineappl/compare/v0.8.1...v0.8.2
[0.8.1]: https://github.com/NNPDF/pineappl/compare/v0.8.0...v0.8.1
[0.8.0]: https://github.com/NNPDF/pineappl/compare/v0.7.4...v0.8.0
[0.7.5]: https://github.com/NNPDF/pineappl/compare/v0.7.4...v0.7.5
[0.7.4]: https://github.com/NNPDF/pineappl/compare/v0.7.3...v0.7.4
[0.7.3]: https://github.com/NNPDF/pineappl/compare/v0.7.2...v0.7.3
[0.7.2]: https://github.com/NNPDF/pineappl/compare/v0.7.1...v0.7.2
[0.7.1]: https://github.com/NNPDF/pineappl/compare/v0.7.0...v0.7.1
[0.7.0]: https://github.com/NNPDF/pineappl/compare/v0.6.3...v0.7.0
[0.6.3]: https://github.com/NNPDF/pineappl/compare/v0.6.2...v0.6.3
[0.6.2]: https://github.com/NNPDF/pineappl/compare/v0.6.1...v0.6.2
[0.6.1]: https://github.com/NNPDF/pineappl/compare/v0.6.0...v0.6.1
[0.6.0]: https://github.com/NNPDF/pineappl/compare/v0.5.9...v0.6.0
[0.5.9]: https://github.com/NNPDF/pineappl/compare/v0.5.8...v0.5.9
[0.5.8]: https://github.com/NNPDF/pineappl/compare/v0.5.7...v0.5.8
[0.5.7]: https://github.com/NNPDF/pineappl/compare/v0.5.6...v0.5.7
[0.5.6]: https://github.com/NNPDF/pineappl/compare/v0.5.5...v0.5.6
[0.5.5]: https://github.com/NNPDF/pineappl/compare/v0.5.4...v0.5.5
[0.5.4]: https://github.com/NNPDF/pineappl/compare/v0.5.3...v0.5.4
[0.5.3]: https://github.com/NNPDF/pineappl/compare/v0.5.2...v0.5.3
[0.5.2]: https://github.com/NNPDF/pineappl/compare/v0.5.1...v0.5.2
[0.5.1]: https://github.com/NNPDF/pineappl/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/NNPDF/pineappl/compare/v0.4.1...v0.5.0
[0.4.1]: https://github.com/NNPDF/pineappl/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/NNPDF/pineappl/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/NNPDF/pineappl/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/NNPDF/pineappl/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/NNPDF/pineappl/compare/v0.0.0...v0.1.0
