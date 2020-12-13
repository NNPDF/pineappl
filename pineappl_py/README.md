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

## And for deploying?
The very same way, `maturin` it's also able to deploy to PyPI:

```sh
maturin publish
```
