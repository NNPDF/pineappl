[![Rust](https://github.com/NNPDF/pineappl/workflows/Rust/badge.svg)](https://github.com/NNPDF/pineappl/actions?query=workflow%3ARust)
[![crates.io](https://img.shields.io/crates/v/pineappl.svg)](https://crates.io/crates/pineappl_cli)

# PineAPPL CLI

This is the command-line interface (CLI) to [`pineappl`](https://crates.io/crates/pineappl).

# Installation

The easiest way to install `pineappl` CLI is to download the pre-built binaries, for
example using `pip`:

```sh
pip install pineappl-cli
```

or using the installation script (see this [guide](https://nnpdf.github.io/pineappl/docs/installation.html)
for more details):

```sh
curl --proto '=https' --tlsv1.2 -sSf https://nnpdf.github.io/pineappl/install-cli.sh | sh
```

Alternatively, the CLI can be installed in a development mode by running the following
command in the root of the repository:

```sh
cargo install pineappl_cli
```

# Documentation

A good starting point to learn what `pineappl` can do is the [CLI tutorial].
For more information, see:

- [CLI tutorial](https://nnpdf.github.io/pineappl/docs/cli-tutorial.html): a tutorial
  on how to use the `pineappl` CLI interface.
- [CLI reference](https://nnpdf.github.io/pineappl/docs/cli-reference.html): a reference
  for the different parameters of the CLI.
