# PyO3 bindings for PineAPPL

## How to..

compile the bindings in a very easy way and immediately use from python:

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

The very same way, `maturin` it's also able to deploy to PyPI:

```sh
maturin publish
```

## Getting Started (`->` carbon copying...)

I'm using as an example the [`pypolars`](https://github.com/ritchie46/polars/blob/master/py-polars/) package.

The reason is that it is very close to this use case, even if far more
complicated (hopefully). In particular:

- the main library it's developed in `rust` for `rust`
- then python bindings are provided on top
- it's using `pyo3` and `maturin`
- it contains both functions and objects
- it includes python tests with `pytest`

I chose the following as prototypes for implementing `pineappl` ingredients:

- **function**: the default `pyo3` example it's sufficient at the moment
- **class**: I chose to follow `pypolars.Series` as a prototype
  - binding definition: `src/series.rs`
  - binding export: `src/lib.rs`
  - python packaging: `pypolars/series.py` (and finally `pypolars/__init__.py`)
  - usage: `example/small_intro.ipynb`
- **packaging**: how to collect all the ingredients for the python package
  - the basic workflows are in `tasks.sh`
  - also `Dockerfile`s are available, for building with multiple python versions
    in isolation
