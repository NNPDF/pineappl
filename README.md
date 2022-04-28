[![Rust](https://github.com/N3PDF/pineappl/workflows/Rust/badge.svg)](https://github.com/N3PDF/pineappl/actions?query=workflow%3ARust)
[![codecov](https://codecov.io/gh/N3PDF/pineappl/branch/master/graph/badge.svg)](https://codecov.io/gh/N3PDF/pineappl)
[![Documentation](https://docs.rs/pineappl/badge.svg)](https://docs.rs/pineappl)
[![crates.io](https://img.shields.io/crates/v/pineappl.svg)](https://crates.io/crates/pineappl)

# Introduction

This repository contains libraries, tools, and interfaces to read and write
`PineAPPL` interpolation grids.

There are four main crates in this repository:

- [`pineappl`](https://crates.io/crates/pineappl) is the crate containing the
  main functionality,
- [`pineappl_capi`](https://crates.io/crates/pineappl_capi/) installs a library
  and a C header, to use PineAPPL from your C, C++ or Fortran programs,
- [`pineappl_cli`](https://crates.io/crates/pineappl_cli/) installs the program
  `pineappl` to use PineAPPL from the command line and
- [`pineappl_py`](https://pypi.org/project/pineappl/) is the Python interface.

Similar projects are:

- [APPLgrid](https://applgrid.hepforge.org/) and
- [fastNLO](https://fastnlo.hepforge.org/).

# Installation

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pineappl/badges/installer/conda.svg)](https://anaconda.org/conda-forge/pineappl)
[![AUR](https://img.shields.io/aur/version/pineappl)](https://aur.archlinux.org/packages/pineappl)

`PineAPPL` is written in [`Rust`](https://www.rust-lang.org/) and therefore
needs the Rust compiler and its build system `cargo`. If `cargo` isn't
installed, use your favourite package manager to install it, or go to
<https://www.rust-lang.org/tools/install> and follow the instructions there.

Next, install the command-line interface (CLI) by choosing either the *release*
or *development version* below. In both cases the binary `pineappl` will be
installed user-wide, typically into `~/.cargo/bin`. You can use this binary to
perform all kinds of operations on PineAPPL grids.

For most users the release version is recommended, as we guarantee that all
grids generated with release versions will be supported in all future release
versions (backwards compatibility guarantee). The advantage of the development
version is that it typically supports more features.

## Release version (recommended)

Simply run

    cargo install pineappl_cli

anywhere and you are done; this will automatically download the most-recently
released version from [crates.io](https://crates.io).

## Development version (alternative)

To use the most recent version available run

    cargo install --git https://github.com/N3PDF/pineappl.git

Instead, if you plan to make changes to the source code it's better to checkout
this repository and run

    cargo install --path pineappl_cli

inside it.

## Optional: fastNLO converter

If you'd like to convert fastNLO tables to PineAPPL, make sure to install
[fastNLO](https://fastnlo.hepforge.org/) first and add the switch
`--features=fastnlo` during the CLI's installation, for instance for the
development version:

    cargo install --features=fastnlo --path pineappl_cli

## Optional: C interface

If you plan to use one of the supported Monte Carlo programs to *generate*
PineAPPL grids, or if you want to access the contents of grids from your own
program, you will likely need the C interface (unless you are using Python, see
below). In that case proceed by installing

- `cargo-c`, which is required for the next step:

      cargo install cargo-c

- Now install `pineappl_capi`, PineAPPL's C API:

      cd pineappl_capi
      cargo cinstall --release --prefix=${prefix}
      cd ..

  where `${prefix}` points to the desired installation directory.

- Finally, you need to set the environment variables `PKG_CONFIG_PATH` and
  `LD_LIBRARY_PATH` to the right directories. Adding

      export LD_LIBRARY_PATH=${prefix}/lib:$LD_LIBRARY_PATH
      export PKG_CONFIG_PATH=${prefix}/lib/pkgconfig:$PKG_CONFIG_PATH

  to your `~/.bashrc` should do the trick (remember to replace `${prefix}` with
  the correct directory). You can check `PKG_CONFIG_PATH` by running

      pkg-config pineappl_capi --libs

  which should print the library flags needed to link against the C API. If
  there's no output or an error, double-check that `PKG_CONFIG_PATH` is in the
  environment and that it points to a directory containing the
  `pineappl_capi.pc` file.

## Optional: Python interface

[![PyPI version](https://badge.fury.io/py/pineappl.svg)](https://badge.fury.io/py/pineappl)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pineappl/badges/installer/conda.svg)](https://anaconda.org/conda-forge/pineappl)
[![AUR](https://img.shields.io/aur/version/pineappl)](https://aur.archlinux.org/packages/pineappl)

To install the Python interface, run

    pip install pineappl

For more documentation and more information see its
[README](pineappl_py/README.md).

# Tutorials

To get to know PineAPPL, take the [tutorial](docs/cli-tutorial.md)!

# Contributions

Before submitting a pull request please read the
[contribution guidelines](CONTRIBUTING.md).

# Citation

[![arXiv](https://img.shields.io/badge/arXiv-2008.12789-b31b1b?labelColor=222222)](https://arxiv.org/abs/2008.12789)
[![DOI](https://zenodo.org/badge/248306479.svg)](https://zenodo.org/badge/latestdoi/248306479)

If you use PineAPPL, please cite the zenodo DOI above and the following reference:

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
