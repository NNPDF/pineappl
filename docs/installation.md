# Installation

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pineappl/badges/version.svg)](https://anaconda.org/conda-forge/pineappl)
[![AUR](https://img.shields.io/aur/version/pineappl)](https://aur.archlinux.org/packages/pineappl)

`PineAPPL` and its interface is written in [Rust](https://www.rust-lang.org/),
but besides Rust PineAPPL also offers interfaces for the most popular
programming languages: C, C++ and Fortran via the C API and Python through
[PyO3](https://pyo3.rs/).

Furthermore the program `pineappl` can be installed that will allow you to
perform many operations on grids in your favorite shell: the command-line
interface (CLI).

Pick the interface(s) you're planning to use and follow the corresponding
instructions below. If you don't know which interfaces you'll likely use, here
are a few guidelines:

- if you're planning to use PineAPPL with Python and/or within the NNPDF
  fitting framework, you'll need the [Python interface](#python);
- if you want to run a Monte Carlo to *generate* PineAPPL grids, you'll likely
  need the [CAPI](#capi-interface-for-c-c-and-fortran);
- if you want to quickly produce predictions, plots and small analyses install
  the [CLI](#cli-pineappl-for-your-shell).

## CAPI: Interface for C, C++ and Fortran

You can install the CAPI in two different ways:

- [using pre-built libraries](#using-pre-built-libraries), which doesn't
  require installing any other dependencies, or
- alternatively you can [build it from source](#building-the-capi-from-source),
  but this requires Rust and Cargo.

### Using pre-built libraries

The fastest way to install the CAPI is to download the pre-built libraries:

    curl --proto '=https' --tlsv1.2 -sSf https://nnpdf.github.io/pineappl/install-capi.sh | sh

The installation script will prompt you for the installation directory where
the files shall be installed to. If you want to pass this directory on the
command line, change the arguments to the shell:

    .. | sh -s -- --prefix /my/custom/installation/path

By default `install-capi.sh` will download the latest stable release. If you'd
like a specific version, pass the version along with `--version`:

    .. | sh -s -- --version 0.6.0-alpha.18

After the installation the script will check whether certain environment
variables are set correctly and will complain if this is not the case. Read and
follow the instructions of the script.

### Building the CAPI from source

If you want to build the CAPI from its sources instead, you first need to

0. Install Rust, see the [instructions](#rust) below.
1. Then install `cargo-c`, which is required for the next step:

       cargo install cargo-c

   It is possible that the installation fails if your Rust compiler is too old.
   In that case update Rust or try installing an older version of `cargo-c`:

       cargo install cargo-c --version 0.9.14+cargo-0.67

2. Next, check out PineAPPL's repository, for instance with

       git clone https://github.com/NNPDF/pineappl.git
       cd pineappl

3. Now install `pineappl_capi`, PineAPPL's C API:

       cd pineappl_capi
       cargo cinstall --release --prefix=${prefix} --libdir=${prefix}/lib
       cd ..

   where `${prefix}` points to the desired installation directory.
4. Finally, if you installed the library into a non-standard prefix you need to
   set the environment variables `PKG_CONFIG_PATH` and `LD_LIBRARY_PATH` to the
   right directories. Adding

       export LD_LIBRARY_PATH=${prefix}/lib:$LD_LIBRARY_PATH
       export PKG_CONFIG_PATH=${prefix}/lib/pkgconfig:$PKG_CONFIG_PATH

   to your `~/.bashrc` should do the trick (remember to replace `${prefix}` with
   the correct directory). You can check `PKG_CONFIG_PATH` by running

       pkg-config pineappl_capi --libs

   which should print the library flags needed to link against the C API. If
   there's no output or an error, double-check that `PKG_CONFIG_PATH` is in the
   environment and that it points to a directory containing the
   `pineappl_capi.pc` file.

If you want to update the CAPI, first update the PineAPPL repository with `git
pull` and then repeat step 3.

## CLI: `pineappl` for your shell

You can install the CLI in two different ways:

- [using pre-built binaries](#using-pre-built-binaries), which doesn't require
  installing any other dependencies, or
- alternatively you can [build it from source](#building-the-cli-from-source),
  but this requires Rust and Cargo.

### Using pre-built binaries

The fastest way to install the CLI is to download the pre-built binaries, for
instance using `pip`:

    pip install pineappl-cli

or using the installation script:

    curl --proto '=https' --tlsv1.2 -sSf https://nnpdf.github.io/pineappl/install-cli.sh | sh

The installation script will prompt you for the installation directory where
the files shall be installed to. If you want to pass this directory on the
command line, change the arguments to the shell:

    .. | sh -s -- --prefix /my/custom/installation/path

By default `install-cli.sh` will download the latest stable release. If you'd
like a specific version, pass the version along with `--version`:

    .. | sh -s -- --version 0.6.0-alpha.18

After the installation the script will check whether certain environment
variables are set correctly and will complain if this is not the case. Read and
follow the instructions of the script.

### Building the CLI from source

You need to install [Rust](#rust) first (see below). Then simply run

    cargo install --locked pineappl_cli

anywhere and you are done; this will automatically download the most-recently
released version from [crates.io](https://crates.io).

More functionality can be added by appending the `--feature=feature1,...` flag,
see below.

To update the CLI simply repeat the `cargo install` step; if a newer version is
available it will be automatically downloaded, otherwise `cargo` will exit with
a message saying that the most recent version is already installed.

If `LHAPDF` is not found, you need to update your `PKG_CONFIG_PATH` variable,
which must point to the directory where LHAPDF's `.pc` file is installed:

    export PKG_CONFIG_PATH=${lhapdf_prefix}/lib/pkgconfig:$PKG_CONFIG_PATH

(replace `${lhapdf_prefix}` with your LHAPDF installation root).

#### Optional: APPLgrid exporter/importer

If you'd like to convert APPLgrids to PineAPPL grids, or vice versa, make sure
to

1. install [APPLgrid](https://applgrid.hepforge.org/) (you need at least
   version 1.6.27),
2. set the environment variable `APPL_IGRID_DIR` to the `src` directory of
   APPLgrid and
3. add the switch `--features=applgrid` during the CLI's installation, for
   instance:

       APPL_IGRID_DIR=/tmp/applgrid-1.6.27/src cargo install --locked --features=applgrid pineappl_cli

#### Optional: Evolution/EKO support

If you'd like to convert PineAPPL grids into FK tables using [evolution kernel
operators (EKO)](https://eko.readthedocs.io/), add the switch
`--features=evolve` during the CLI's installation, for instance:

    cargo install --locked --features=evolve pineappl_cli

### Optional: fastNLO importer

If you'd like to convert fastNLO tables to PineAPPL grids, make sure to install
[fastNLO](https://fastnlo.hepforge.org/) first and add the switch
`--features=fastnlo` during the CLI's installation, for instance:

    cargo install --locked --features=fastnlo pineappl_cli

#### Optional: FK table converter

If you'd like to convert NNPDF's legacy FK tables to PineAPPL grids, add the switch
`--features=fktable` during the CLI's installation, for instance:

    cargo install --locked --features=fktable pineappl_cli

#### Alternative: development version

To use the most recent version available run

    cargo install --locked --git https://github.com/NNPDF/pineappl.git pineappl_cli

Instead, if you plan to make changes to the source code it's better to checkout
this repository and run

    cargo install --locked --path pineappl_cli

inside it.

#### Man pages

Most of PineAPPL's help is available as `pineappl --help` and as `pineappl
convolve --help`, for example. However, the same information can also be
installed as man pages, which can then be used as either

    man pineappl
    man pineappl-convolve

or

    pineappl help
    pineappl help convolve

For this to work first install the [CLI](#cli-pineappl-for-your-shell). Next,
go to your local copy of PineAPPL's repository and create the man page
directory relative to the CLI,

    mkdir -p $(dirname $(which pineappl))/../share/man/man1/

and finally install the manpages there:

    cargo xtask install-manpages $(dirname $(which pineappl))/../share/man/man1/

## Python

[![PyPI version](https://badge.fury.io/py/pineappl.svg)](https://badge.fury.io/py/pineappl)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pineappl/badges/version.svg)](https://anaconda.org/conda-forge/pineappl)
[![AUR](https://img.shields.io/aur/version/pineappl)](https://aur.archlinux.org/packages/pineappl)

To install the Python interface, run

    pip install pineappl

This will not require any previous installation of Rust. For more documentation
and more information see its [README](../pineappl_py/README.md).

## Rust

You will need the Rust compiler and its build system `cargo`. If `cargo` is
already installed, make sure it is recent enough:

    cargo --version

This should show a version that is at least 1.80.1. If you do not have `cargo`
or if it is too old, go to <https://www.rust-lang.org/tools/install> and follow
the instructions there.

The Rust crate is called [pineappl](https://docs.rs/pineappl/latest/pineappl/).
