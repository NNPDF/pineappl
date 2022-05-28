Recipes
=======

Below we list some common use cases with their solutions.


How can I convolute a given PineAPPL grid with my PDF?
------------------------------------------------------

.. code:: python

    import pineappl
    import lhapdf
    g = pineappl.grid.Grid.read("path/to/grid.pineappl")
    pdf = lhapdf.mkPDF("YourPDF", 0)
    bins = g.convolute(pdf.xfxQ2, pdf.xfxQ2, pdf.alphasQ2)

If the grid is actually an FkTable just replace

.. code:: python

    g = pineappl.fk_table.FKTable.read("path/to/grid.pineappl")


How can I get the bin configurations from a given PineAPPL grid?
----------------------------------------------------------------

.. code:: python

    import pineappl
    g = pineappl.grid.Grid.read("path/to/grid.pineappl")
    # first get the number of dimensions
    bin_dims = g.bin_dimensions()
    # now you can get each of them
    for bin_dim in range(bin_dims):
        bin_left = g.bin_left(bin_dim)
        bin_right = g.bin_right(bin_dim)
