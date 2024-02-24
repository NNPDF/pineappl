# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- fixed yet another problem that prevent the Python interface for Python 3.6
  from being successfully installed

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

[Unreleased]: https://github.com/NNPDF/pineappl/compare/v0.7.3...HEAD
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
