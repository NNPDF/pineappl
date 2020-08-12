[![Rust](https://github.com/N3PDF/pineappl/workflows/Rust/badge.svg)](https://github.com/N3PDF/pineappl/actions?query=workflow%3ARust)
[![codecov](https://codecov.io/gh/N3PDF/pineappl/branch/master/graph/badge.svg)](https://codecov.io/gh/N3PDF/pineappl)
[![Documentation](https://docs.rs/pineappl/badge.svg)](https://docs.rs/pineappl)
[![crates.io](https://img.shields.io/crates/v/pineappl.svg)](https://crates.io/crates/pineappl)
[![DOI](https://zenodo.org/badge/248306479.svg)](https://zenodo.org/badge/latestdoi/248306479)
[![Documentation Status](https://readthedocs.org/projects/pineappl/badge/?version=latest)](https://pineappl.readthedocs.io/en/latest/?badge=latest)

# Introduction

This repository contains libraries, tools, and interfaces to read and write
`PineAPPL` grids.

There are three crates in this repository:

- [`pineappl`](https://crates.io/crates/pineappl) is the crate containing the
  main functionality
- [`pineappl_capi`](https://crates.io/crates/pineappl) installs a library and a
  C header, to use PineAPPL inside a C program
- [`pineappl_cli`](https://crates.io/crates/pineappl) installs a program to use
  PineAPPL from the command line

# Installation

`PineAPPL` depends on [`Rust`](https://www.rust-lang.org/). If it's already
installed make sure that you have a recent version, otherwise the following
steps might break during compilations. If it's not installed yet, use your
favourite package manager to install it, or go to
<https://www.rust-lang.org/tools/install> and follow the instructions there.

Proceed by installing `cargo-c`, which is required by `pineappl_capi`:

    cargo install cargo-c

Next, install `pineappl_capi`:

    cd pineappl_capi
    cargo cinstall --release --prefix=${prefix}
    cd ..

and finally the command-line program:

    cargo install --path pineappl_cli

Make sure that all the required environment variables are set. See the
`README.md` of `pineappl_capi` for further instructions.

For the python interface please refer to the documentation [here](https://pineappl.readthedocs.io/).