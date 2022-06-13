[![Rust](https://github.com/N3PDF/pineappl/workflows/Rust/badge.svg)](https://github.com/N3PDF/pineappl/actions?query=workflow%3ARust)
[![codecov](https://codecov.io/gh/N3PDF/pineappl/branch/master/graph/badge.svg)](https://codecov.io/gh/N3PDF/pineappl)
[![Documentation](https://docs.rs/pineappl/badge.svg)](https://docs.rs/pineappl)
[![crates.io](https://img.shields.io/crates/v/pineappl.svg)](https://crates.io/crates/pineappl)
[![Minimum cargo version](https://img.shields.io/badge/cargo-1.54+-lightgray.svg)](https://github.com/N3PDF/pineappl#installation)

# Introduction

This repository contains tools, libraries and interfaces to read and write
`PineAPPL` interpolation grids, which store theoretical predictions for
[high-energy collisions] independently from their [PDFs].

Similar projects are:

- [APPLgrid] and
- [fastNLO].

This repository hosts four main crates:

- [`pineappl`] is the crate containing the main functionality,
- [`pineappl_capi`] installs a library and a C header, to use PineAPPL from
  your C, C++ or Fortran programs,
- [`pineappl_cli`] installs the program `pineappl` to use PineAPPL from the
  command line and
- [`pineappl_py`] is the Python interface.

[APPLgrid]: https://applgrid.hepforge.org/
[fastNLO]: https://fastnlo.hepforge.org/
[high-energy collisions]: https://en.wikipedia.org/wiki/Particle_physics
[PDFs]: https://en.wikipedia.org/wiki/Parton_(particle_physics)#Parton_distribution_functions
[`pineappl`]: https://crates.io/crates/pineappl/
[`pineappl_capi`]: https://crates.io/crates/pineappl_capi/
[`pineappl_cli`]: https://crates.io/crates/pineappl_cli/
[`pineappl_py`]: https://pypi.org/project/pineappl/

# Documentation

Documentation is available [here](docs/README.md).

# Installation

Installation instructions are [here](docs/installation.md).

# Contributions

Before submitting a pull request please read the
[contribution guidelines](CONTRIBUTING.md).

# Citation

[![arXiv](https://img.shields.io/badge/arXiv-2008.12789-b31b1b?labelColor=222222)](https://arxiv.org/abs/2008.12789)
[![DOI](https://zenodo.org/badge/248306479.svg)](https://zenodo.org/badge/latestdoi/248306479)

If you use PineAPPL, please cite

1. the zenodo DOI above and
2. the following reference:

   ```
   @article{Carrazza:2020gss,
       author = "Carrazza, S. and Nocera, E. R. and Schwan, C. and Zaro, M.",
       title = "{PineAPPL: combining EW and QCD corrections for fast evaluation of LHC processes}",
       eprint = "2008.12789",
       archivePrefix = "arXiv",
       primaryClass = "hep-ph",
       doi = "10.1007/JHEP12(2020)108",
       journal = "JHEP",
       volume = "12",
       pages = "108",
       year = "2020"
   }
   ```
