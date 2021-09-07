# Python bindings for PineAPPL

This crate uses [PyO3] to provide Python bindings to PineAPPL's [Rust API]. It
will supersede the [Python wrapper] written with [ctypes].

# Installation

TODO

# Compilation (for development)

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

[PyO3]: https://pyo3.rs
[Rust API]: https://docs.rs/pineappl
[Python wrapper]: ../wrappers/python/README.md
[ctypes]: https://docs.python.org/3/library/ctypes.html
