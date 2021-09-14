# Python bindings for PineAPPL

This crate uses [PyO3] to provide Python bindings to PineAPPL's [Rust API].

# Installation

```sh
pip install pineappl
```

# Compilation (for development)

1. Make a virtual environment in your favorite way (suggested: `virtualenv`)

```sh
virtualenv env # --system-site-packages
```

2. Activate the environment and install `maturin` via `pip`

```sh
. env/bin/activate
pip install -r dev.requirements.txt
```

3. Run `maturin` to compile and install the library as a python package in the
   current environment

```sh
maturin develop
```

[PyO3]: https://pyo3.rs
[Rust API]: https://docs.rs/pineappl
