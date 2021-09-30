Installation
============

Installing with pip
-------------------
Thanks to the `container build <https://github.com/N3PDF/pineappl/blob/master/pineappl_py/package/README.md>`_, installing should be as easy as

.. code:: sh

    pip install pineappl


Development installation
------------------------

1. Make a virtual environment in your favorite way (suggested: `virtualenv`)

.. code:: sh

    virtualenv env # --system-site-packages


2. Activate the environment and install `maturin` via `pip`

.. code:: sh

    . env/bin/activate
    pip install -r dev.requirements.txt


3. Run `maturin` to compile and install the library as a python package in the
   current environment

.. code:: sh

    maturin develop
