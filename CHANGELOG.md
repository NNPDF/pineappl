# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- The command-line interface `pineappl` no longer prints both the differential
  and integrated cross sections. Instead it either prints the differential
  cross sections, or, if the switch `-i` or `--integrated` is given, the
  integrated cross sections (without bin limits and/or normalizations) are
  printed

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

[Unreleased]: https://github.com/N3PDF/pineappl/compare/v0.4.1...HEAD
[0.4.1]: https://github.com/N3PDF/pineappl/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/N3PDF/pineappl/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/N3PDF/pineappl/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/N3PDF/pineappl/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/N3PDF/pineappl/compare/v0.0.0...v0.1.0
