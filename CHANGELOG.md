# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- added the options `--absolute`, `--orders`, and `--lumis` to the subcommand
  `channels`

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

[Unreleased]: https://github.com/N3PDF/pineappl/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/N3PDF/pineappl/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/N3PDF/pineappl/compare/v0.0.0...v0.1.0
