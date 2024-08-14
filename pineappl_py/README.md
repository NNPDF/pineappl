[![PyPI version](https://badge.fury.io/py/pineappl.svg)](https://badge.fury.io/py/pineappl)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/pineappl/badges/installer/conda.svg)](https://anaconda.org/conda-forge/pineappl)
[![AUR](https://img.shields.io/aur/version/pineappl)](https://aur.archlinux.org/packages/pineappl)
[![Documentation Status](https://readthedocs.org/projects/pineappl/badge/?version=latest)](https://pineappl.readthedocs.io/en/latest/?badge=latest)

# Python bindings for PineAPPL

This crate uses [PyO3] to provide Python bindings to PineAPPL's [Rust API].

For installation instructions see the [documentation].

[PyO3]: https://pyo3.rs
[Rust API]: https://docs.rs/pineappl
[documentation]: https://pineappl.readthedocs.io/en/latest/installation.html

## Development

Run

```shell
python -m venv env && . env/bin/activate
```

to setup a new environment and check that `pip --version` returns at least `pip
22.0 from ...`. If not, upgrade `pip` via

```shell
pip install -U pip
```

Next, install `maturin`:

```shell
pip install maturin
```

Run

```shell
maturin develop
```

to build the project, which also installs it into the environment so that it
can be used in Python projects that use the same environment.

### Documentation

Run the following once to install the documentation's dependencies:

```shell
pip install '.[docs]'
```

Then run

```shell
( cd docs && make clean html )
```

to generate the documentation.
