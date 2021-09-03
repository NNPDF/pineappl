# Python bindings for PineAPPL

We're using [PyO3](https://pyo3.rs) to provide python binding to `pineappl` written in Rust.
This package will superseed the former Python wrappers written with
[ctypes](https://docs.python.org/3/library/ctypes.html).

The layout consist of 2 layers:
- a Rust layer (under `./src/`) that provides and declares the bridge from the original
  `pineappl` library to PyO3 (e.g. casting Python-understable types to Rust types)
- a Python layer (under `./pineappl/`) that hides the PyO3 objects with pure
  Python objects and provides some convenience methods

## How to compile (for development)

1. Make a virtual environment in your favorite way (suggested: `virtualenv`)

```sh
virtualenv env # --system-site-packages
```

2. Activate the environment and install `maturin` via `pip`

```sh
. env/bin/activate
pip install maturin
```

3. Run `maturin` to compile and install the library as a python package in the
   current environment

```sh
maturin develop
```

**NOTE:** the packaged built and installed with `maturin` will contain not
only the compiled of `src/*.rs`, but also the python package in `pineappl/`
(i.e. `<root>/pineappl_py/pineappl`, to disambiguate a little the
name-clashing).

## And for deploying?

The very same way, `maturin` is also able to deploy to [PyPI](https://pypi.org/):

```sh
maturin publish
```
