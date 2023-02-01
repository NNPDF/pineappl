Installation
============

Installing with pip
-------------------
Thanks to the `container build <https://github.com/NNPDF/pineappl/blob/master/pineappl_py/package/README.md>`_, installing should be as easy as

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
    pip install maturin


3. Run `maturin` to compile and install the library as a python package in the
   current environment

.. code:: sh

    maturin develop

Notebooks
~~~~~~~~~

To be able to edit the documentation's notebooks in your virtual environment,
you need to install the relative kernel:

.. code:: sh

    # enter the environment
    pip install ipykernel
    python -m ipykernel install --user --name=$(basename $PWD)

And select in the notebook the kernel with the given name (notice that
``$(basename $PWD)`` is only a suggestion for the name, you can choose whatever
string you prefer).
