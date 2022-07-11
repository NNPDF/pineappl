# PineAPPL documentation

A good starting point to learn what PineAPPL and its command-line program
`pineappl` can do is the [tutorial](cli-tutorial.md)!

## Documentation

Here is an alphabetically ordered list of all documentation:

- [CLI tutorial](cli-tutorial.md): a tutorial for the command-line interface
  (CLI) `pineappl`.
- [CLI reference](cli-reference.md): a reference for all parameters of the CLI.
- [Changelog](../CHANGELOG.md): a list of additions and changes for all
  released and unreleased versions of PineAPPL.
- [Grid repository](https://github.com/NNPDF/pineapplgrids/): pre-computed grids
  for specific experimental setups (currently private).
- [Installation](installation.md): installation instructions.
- [Madgraph5_aMC@NLO](mg5_aMC.md): how to create PineAPPL grids with
  [Madgraph5_aMC@NLO](https://launchpad.net/mg5amcnlo/).
- [Metadata](metadata.md): a list of all recognized keys and values of the
  metadata.

## API documentation

- [C and C++](https://docs.rs/pineappl_capi/)
- for Fortran there is no dedicated documentation available, because it is a
  [wrapper](../examples/fortran/pineappl.f90) of the C API; simply read the C
  API documentation, the function names are exactly the same
- [Python](https://pineappl.readthedocs.io/)
- [Rust](https://docs.rs/pineappl/)

### Code examples

Another way to learn using the APIs is to have a look/modify the
[examples](../examples/).

## Developer documentation

- [Contribution guidelines](../CONTRIBUTING.md): technical guidelines for how
  to contribute to PineAPPL's development.
