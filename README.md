[![Rust](https://github.com/N3PDF/pineappl/workflows/Rust/badge.svg)](https://github.com/N3PDF/pineappl/actions?query=workflow%3ARust)
[![codecov](https://codecov.io/gh/N3PDF/pineappl/branch/master/graph/badge.svg)](https://codecov.io/gh/N3PDF/pineappl)
[![Documentation](https://docs.rs/pineappl/badge.svg)](https://docs.rs/pineappl)
[![crates.io](https://img.shields.io/crates/v/pineappl.svg)](https://crates.io/crates/pineappl)
[![DOI](https://zenodo.org/badge/248306479.svg)](https://zenodo.org/badge/latestdoi/248306479)

# Introduction

This repository contains libraries, tools, and interfaces to read and write
`PineAPPL` grids.

There are three crates in this repository:

- [`pineappl`](https://crates.io/crates/pineappl) is the crate containing the
  main functionality
- [`pineappl_capi`](https://crates.io/crates/pineappl) installs a library and a
  C header, to use PineAPPL from your C, C++, Fortran, or Python programs
- [`pineappl_cli`](https://crates.io/crates/pineappl) installs the program
  `pineappl` to use PineAPPL from the command line

# Installation

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pineappl/badges/installer/conda.svg)](https://anaconda.org/conda-forge/pineappl)
[![AUR](https://img.shields.io/aur/version/pineappl)](https://aur.archlinux.org/packages/pineappl)

- `PineAPPL` is written in [`Rust`](https://www.rust-lang.org/) and therefore
  needs the Rust compiler and build system `cargo`. If `cargo` isn't installed,
  use your favourite package manager to install it, or go to
  <https://www.rust-lang.org/tools/install> and follow the instructions there.

- Next install the command-line interface:

      cargo install --path pineappl_cli

  This will install the binary `pineappl` user-wide, typically into
  `~/.cargo/bin`. You can use this binary to perform all kinds of operations
  on PineAPPL grid files.

- Proceed by installing `cargo-c`, which is required for the next step:

      cargo install cargo-c

- Install `pineappl_capi` (the C API, needed if you plan to use PineAPPL in
  your C, C++, Fortran, or Python program):

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

- For the python interface (optional) please look into the subfolder `./pineappl_py`

# Citation

If you use PineAPPL, please cite all of the following references:

[![arXiv](https://img.shields.io/badge/arXiv-2008.12789-b31b1b?labelColor=222222)](https://arxiv.org/abs/2008.12789)
[![DOI](https://zenodo.org/badge/248306479.svg)](https://zenodo.org/badge/latestdoi/248306479)
