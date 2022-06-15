# Installation

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pineappl/badges/installer/conda.svg)](https://anaconda.org/conda-forge/pineappl)
[![AUR](https://img.shields.io/aur/version/pineappl)](https://aur.archlinux.org/packages/pineappl)

`PineAPPL` is written in [`Rust`](https://www.rust-lang.org/) and therefore
needs the Rust compiler and its build system `cargo`. If `cargo` is already
installed, make sure it is recent enough:

    cargo --version

This should show a version that 1.54 or newer. If you do not have `cargo` or it
is too old, go to <https://www.rust-lang.org/tools/install> and follow the
instructions there.

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

  It is possible that the installation fails if your Rust compiler is too old.
  In that case update Rust or try installing an older version of `cargo-c`:

      cargo install cargo-c --version 0.7.3

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
[README](../pineappl_py/README.md).