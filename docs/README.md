# PineAPPL documentation

A good starting point to learn what PineAPPL and its command-line program
`pineappl` can do is the [tutorial](cli-tutorial.md)!

## Documentation

Here is a list of all the documentation is this directory:

- [CLI tutorial](cli-tutorial.md): a tutorial for the command-line interface
  (CLI).
- [CLI reference](cli-reference.md): a reference for all parameters of the CLI.
- [Changelog](../CHANGELOG.md): a list of additions and changes for all
  released and unreleased versions of PineAPPL.
- [Contribution guidelines](../CONTRIBUTING.md): technical guidelines for how
  to contribute to PineAPPL's development.
- [Grid repository](): pre-computed grids for specific experimental setups.
- [Installation](installation.md): installation instructions.
- [Madgraph5_aMC@NLO](mg5_aMC.md): how to create PineAPPL grids with
  [Madgraph5_aMC@NLO](https://launchpad.net/mg5amcnlo/).
- [Metadata](metadata.md): a list of all recognized keys and values of the
  metadata.

## API documentation

- [C](https://docs.rs/pineappl_capi/)
- for Fortran there's no dedicated documentation available, because it's a
  [wrapper](../examples/fortran/pineappl.f90) of the C API; simply read the C
  API documentation, the function names are exactly the same
- [Python](https://pineappl.readthedocs.io/)
- [Rust](https://docs.rs/pineappl/)

### Code examples

Another way to learn using the APIs is to have a look/modify the
[examples](../examples/).
