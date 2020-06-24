[![Rust](https://github.com/N3PDF/pineappl/workflows/Rust/badge.svg)](https://github.com/N3PDF/pineappl/actions?query=workflow%3ARust)
[![Documentation](https://docs.rs/pineappl/badge.svg)](https://docs.rs/pineappl_capi)
[![crates.io](https://img.shields.io/crates/v/pineappl.svg)](https://crates.io/crates/pineappl_capi)

# C API to the PineAPPL library

To use [`pineappl`](https://crates.io/crates/pineappl) via the `C` API, you
first need [`cargo-c`](https://crates.io/crates/cargo-c) to generate the C
header. First run `cargo install cargo-c` to install it, and then install the C
API:

    cargo cinstall --release --prefix=${prefix}

Make sure to replace `${prefix}` with the directory you want it installed to.
This crate installs a header file, a library, and a pkg-config file, so make
sure to set the necessary environment variables.

On Linux, you need to set at least `PKG_CONFIG_PATH` to the directory where the
`pineappl_capi.pc` file is. It usually is in `${prefix}/lib/pkgconfig`. If
you've set it to right value the following command

    pkg-config pineappl_capi --libs

should print the library flags needed to link against the library. If there's
no output, double-check your installation and environment variables. Finally,
it's probably necessary to set `LD_LIBRARY_PATH` to the directory where the
PineAPPL shared/static library was installed to.
