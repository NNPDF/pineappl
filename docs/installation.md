# Installation

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pineappl/badges/version.svg)](https://anaconda.org/conda-forge/pineappl)
[![AUR](https://img.shields.io/aur/version/pineappl)](https://aur.archlinux.org/packages/pineappl)

`PineAPPL` is written in [`Rust`](https://www.rust-lang.org/), but besides Rust
it also offers interfaces for the most popular programming languages: C and
C++, Fortran and Python. Furthermore the program `pineappl` can be installed
that will allow you to perform many operations on grids in your favorite shell:
the command-line interface (CLI).

Simply pick the interface(s) you're planning to use and follow the
corresponding instructions below. If you don't know which interfaces you'll
likely use, here are a few guidelines:

- if you're planning to use PineAPPL within the NNPDF fitting framework, you'll
  need the [Python interface](#python);
- if you want to run a Monte Carlo to *generate* PineAPPL grids, you'll need
  the [C interface](#c-c-and-fortran-the-capi);
- if you want to quickly produce predictions, plots and small analyses install
  the [CLI](#cli-pineappl-for-your-shell).

## C, C++ and Fortran: the CAPI

### Using pre-built binaries

The fastest way to install the CAPI is to download the pre-built binaries:

    curl --proto '=https' --tlsv1.2 -sSf https://nnpdf.github.io/pineappl/install-capi.sh | sh

This may not work for the system you're working with. If that's the case,
please open an [Issue](https://github.com/NNPDF/pineappl/issues/new) for that.

### From source

If you want to build the CAPI from source instead, you first need to install
Rust and `cargo`, see the [instructions](#rust) below.

1. Then install `cargo-c`, which is required for the next step:

       cargo install cargo-c

   It is possible that the installation fails if your Rust compiler is too old.
   In that case update Rust or try installing an older version of `cargo-c`:

       cargo install cargo-c --version 0.9.14+cargo-0.67

2. Now install `pineappl_capi`, PineAPPL's C API:

       cd pineappl_capi
       cargo cinstall --release --prefix=${prefix}
       cd ..

   where `${prefix}` points to the desired installation directory.

3. Finally, you need to set the environment variables `PKG_CONFIG_PATH` and
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

## CLI: `pineappl` for your shell

You need to install [Rust](#rust) first (see below). Then simply run

    cargo install pineappl_cli

anywhere and you are done; this will automatically download the most-recently
released version from [crates.io](https://crates.io).

More functionality can be added by adding the `--feature=feature1,feature2,...`
flag, see below.

### Optional: APPLgrid converters

If you'd like to convert APPLgrids to PineAPPL grids, or vice versa, make sure
to

1. install [APPLgrid](https://applgrid.hepforge.org/) (you need at least
   version 1.6.27),
2. set the environment variable `APPL_IGRID_DIR` to the `src` directory of
   APPLgrid and
3. add the switch `--features=applgrid` during the CLI's installation, for
   instance:

       APPL_IGRID_DIR=/tmp/applgrid-1.6.27/src cargo install --features=applgrid pineappl_cli

### Optional: Evolution/EKO support

If you'd like to convert PineAPPL grids into FK tables using evolution kernel
operators (EKO), add the switch `--features=evolve` during the CLI's
installation, for instance:

    cargo install --features=evolve pineappl_cli

### Optional: fastNLO converter

If you'd like to convert fastNLO tables to PineAPPL grids, make sure to install
[fastNLO](https://fastnlo.hepforge.org/) first and add the switch
`--features=fastnlo` during the CLI's installation, for instance:

    cargo install --features=fastnlo pineappl_cli

### Optional: FK table converter

If you'd like to convert NNPDF's legacy FK tables to PineAPPL grids, add the switch
`--features=fktable` during the CLI's installation, for instance:

    cargo install --features=fktable pineappl_cli

### Alternative: development version

To use the most recent version available run

    cargo install --git https://github.com/NNPDF/pineappl.git

Instead, if you plan to make changes to the source code it's better to checkout
this repository and run

    cargo install --path pineappl_cli

inside it.

## Python

[![PyPI version](https://badge.fury.io/py/pineappl.svg)](https://badge.fury.io/py/pineappl)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pineappl/badges/version.svg)](https://anaconda.org/conda-forge/pineappl)
[![AUR](https://img.shields.io/aur/version/pineappl)](https://aur.archlinux.org/packages/pineappl)

To install the Python interface, run

    pip install pineappl

For more documentation and more information see its
[README](../pineappl_py/README.md).

## Rust

You will need the Rust compiler and its build system `cargo`. If `cargo` is
already installed, make sure it is recent enough:

    cargo --version

This should show a version that is at least 1.64.0. If you do not have `cargo`
or if it is too old, go to <https://www.rust-lang.org/tools/install> and follow
the instructions there.

The Rust crate is called [pineappl](https://docs.rs/pineappl/latest/pineappl/).
