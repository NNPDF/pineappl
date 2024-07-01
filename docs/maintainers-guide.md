# Repository structure

- `.cargo/config.toml`: configuration file that adds support for [cargo-xtask]
- `.github/workflows/`: files that steer the CI, see below
  - `capi.yml`: compiles, installs and tests the CAPI by running the the C, C++
    and Fortran examples. It also creates code coverage of the CAPI
  - `container.yml`: creates the container, which is used in many CI jobs. This
    github action must be run manually if a new container should be created
  - `msrv.yml`: checks whether the Rust code conforms to the MSRV specified in
    the top-most `Cargo.toml` file
  - `python.yml`: checks the Python bindings
  - `release.yml`: this action creates a new release whenever a certain git tag
    is pushed to the repository. It then uploads pre-compiled libraries of the
    CAPI for Linux, MacOS and Windows to Github, publishes the crates to
    <https://crates.io/> and the Python wheels to <https://pypi.org/>. Finally,
    it creates a new release on PineAPPL's Github [release page]
  - `rust.yml`: checks whether the Rust code succesfully compiles and generates
    code coverage
- `.github/codecov.yml`: configuration file for PineAPPL's [codecov page]
- `_includes/head-custom.html`: this file is part of PineAPPL's Github-pages
  website and includes KaTeX to support rendering formulas
- `docs`: houses the markdown files and images used to render PineAPPL's
  Github-pages website
- `examples`: contains examples programs to learn how to use PineAPPL's C, C++,
  Fortran and Python APIs
- `maintainer`: contains [maintainer-specific tools]
- `pineappl`: the main Rust crate, documentation at <https://docs.rs/pineappl>
- `pineappl_applgrid`: interface to [APPLgrid], which the CLI uses to convert
  APPLgrids to PineAPPL grids (and vice versa)
- `pineappl_capi`: the crate that builds PineAPPL's CAPI, documentation at
  <https://docs.rs/pineappl_capi>
- `pineappl_cli`: the crate that builds PineAPPL's CLI
- `pineappl_fastnlo`: interface to [fastNLO], which the CLI uses to convert
  fastNLO tables to PineAPPL grids
- `pineappl_py`: the crate that builds PineAPPL's Python interface,
  documentation at <https://pineappl.readthedocs.io>
- `xtask`: crate for [cargo-xtask] commands
- `.envrc`: [direnv] hooks (mainly enable Nix shell)
- `.gitignore`: PineAPPL's Git ignore rules
- `.readthedocs.yml`: configuration for PineAPPL's [Read-the-Docs] Python
  interface documentation
- `CHANGELOG.md`: change log for the code part of this repository
- `CONTRIBUTING.md`: rules for how to contribute to this repository
- `Cargo.lock`: Cargo [lock file]
- `Cargo.toml`: Cargo configuration file
- `LICENSE`: GNU General public [license v3]
- `README.md`: PineAPPL's repository README file. This file also serves as the
  homepage of [PineAPPL's website]
- `_config.yml`: configuration file for PineAPPL's Github-pages website
- `flake.lock`: Nix [flake lock file]
- `flake.nix`: Nix [flake], which tracks cross-language dependencies to define a
  reproducible build for the PineAPPL packages, and a suitable development shell
- `install-capi.sh`: POSIX-compliant shell script to download and install
  PineAPPL's pre-built CAPI

[cargo-xtask]: https://github.com/matklad/cargo-xtask
[release page]: https://github.com/NNPDF/pineappl/releases
[codecov page]: https://app.codecov.io/gh/NNPDF/pineappl
[maintainer-specific tools]: ../maintainer/README.md
[APPLgrid]: https://applgrid.hepforge.org/
[fastNLO]: https://fastnlo.hepforge.org/
[Read-the-Docs]: https://pineappl.readthedocs.io/
[lock file]: https://doc.rust-lang.org/cargo/guide/cargo-toml-vs-cargo-lock.html
[license v3]: https://www.gnu.org/licenses/gpl-3.0.en.html
[PineAPPL's website]: https://nnpdf.github.io/pineappl/
[flake]: https://nixos.org/manual/nix/stable/command-ref/new-cli/nix3-flake.html#flake-format
[flake lock file]: https://nixos.org/manual/nix/stable/command-ref/new-cli/nix3-flake.html#lock-files
[direnv]: https://direnv.net/

## Abbreviations and terminology

- CAPI: C application programmer's interface, a bunch of C functions that allow
  to PineAPPL without Rust. The CAPI allows C, C++ and Fortran programmers to
  use PineAPPL
- Cargo: Rust's [package manager]
- CI: continuous integration, the Github actions that run for every commit
  pushed to the PineAPPL repository
- CLI: command-line interface, the program `pineappl` that allow quick
  operations on PineAPPL grids
- code coverage: this tests which percentage of lines of code are run in the CI
  through tests. The lower this number, the more likely it is that bugs and
  regressions are introduced without noticing them
- container: bundles software with its dependencies in a way that it can be run
  reliably on different computers
- MSRV: minimally supported Rust version, the lowest version of the Rust
  compiler toolchain with which PineAPPL can be built. Run `cargo --version` to
  figure which version you've got

[package manager]: https://doc.rust-lang.org/cargo/index.html

# Connected accounts

- Codecov: <https://app.codecov.io/gh/NNPDF/pineappl>
- Conda: <https://github.com/conda-forge/pineappl-feedstock>
- Crates.io: <https://crates.io/crates/pineappl>
- PyPI: <https://pypi.org/project/pineappl>
- ReadTheDocs: <https://pineappl.readthedocs.io>
- Zenodo: <https://zenodo.org/records/10700087>

We also store testing data in <https://data.nnpdf.science/pineappl/test-data>.
