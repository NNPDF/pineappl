Installing PineAPPL
===================

The PineAPPL package comes with the following modules:

* :ref:`installing-with-conda`
* :ref:`installing-with-pip`
* :ref:`installing-from-source`

_______________________

. _installing-with-conda:

Installing with conda
---------------------

The installation using ``conda`` is the recommended approach to install PineAPPL.

Make sure you have a working conda environment and then use the `conda` program to install
``pineappl`` with:

.. code-block:: bash

      conda install -c conda-forge pineappl

The ``conda`` program will download and install all the required
dependencies for the PineAPPL (precompiled library) and the python interface.

.. _installing-with-pip:

Installing with pip
-------------------

Make sure you have Python 3.6 or greater, then use ``pip`` to install
``pineappl`` with:

.. code-block:: bash

      pip install pineappl

The ``pip`` program will download and install all the required
dependencies for the PineAPPL python wrapper.

.. note::
    The ``pip`` packages install just the python wrapper for the Rust PineAPPL
    library. Before using the python wrapper make sure you have installed the
    Rust PineAPPL library properly.

.. _installing-from-source:

Installing from source
----------------------

In order to install PineAPPL from source, you can simply clone the GitHub
repository with

.. code-block::

      git clone https://github.com/N3PDF/pineappl.git
      cd wrappers/python

then proceed with the installation of requirements with:

.. code-block::

      pip install .

If you prefer to keep changes always synchronized with the code then install using the develop option:

.. code-block::

      pip install -e .
      # or
      python setup.py develop

