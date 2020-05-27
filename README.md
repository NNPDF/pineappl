# Introduction

This repository contains libraries, tools, and interfaces to read and write
`PineAPPL` grids.

# Installation

`PineAPPL` needs `Rust`. If it's not already available on your system, go to
<https://www.rust-lang.org/tools/install> and follow the instructions.

Next, run `cargo build --release` in the top-level directory of this
repository. To install the program `pineappl` run `cargo install --path
pineappl_cli`.

## C API

To use `PineAPPL` via the `C` API, you need `cbindgen` to generate the C
header. First run `cargo install cbindgen` to install it, then run

    cbindgen -c pineappl_capi/cbindgen.toml pineappl_capi/ > pineappl_capi.h

This will create the header in the top-level directory of the repository. The
corresponding library can be found after building the project in

    target/release/pineappl_capi.{so,dylib,dll}
