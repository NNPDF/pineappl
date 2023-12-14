Welcome to PineAPPL
===================

This is the Python wrapper for the `Rust PineAPPL library <https://nnpdf.github.io/pineappl/>`_.

PineAPPL is a computer library that makes it possible to produce fast-interpolation grids for fitting parton distribution functions (PDFs) including corrections of strong and electroweak origin.

The :doc:`installation` instructions are given :doc:`here <installation>`.

A practical example can be found in the ``example/`` subfolder of the `repository <https://github.com/NNPDF/pineappl/>`_.
The Python wrapper is also used in :yadism:`\ ` and :pineko:`\ `. We also list some common :doc:`recipes` here.

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Contents:

   installation
   recipes
   implementation
   API <modules/pineappl/pineappl>
   indices

.. important::

   If you are looking for the methods of a specific class, be aware that part of
   them are just passed to the underlying Rust object, whose class is the same
   of the user-facing one, but prefixed with a ``Py``, e.g.:
   :class:`pineappl.grid.Grid` and :class:`pineappl.pineappl.grid.PyGrid`.

   You will find the documentation of the unwrapped method in the raw ``Py``
   class, while part of the methods are wrapped and thus even documented in the
   user-facing class.
