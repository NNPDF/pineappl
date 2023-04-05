import numpy as np

from .pineappl import PyImportOnlySubgridV1
from .pineappl import PyImportOnlySubgridV2
from .utils import PyWrapper


class ImportOnlySubgridV1(PyWrapper):
    """
    Python wrapper object to :class:`~pineappl.pineappl.PyImportOnlySubgridV1`.

    Parameters
    ----------
        array : numpy.ndarray(float, dim=3)
            3-dimensional subgrid content
        q2_grid : sequence(float)
            scale grid
        x1_grid : sequence(float)
            interpolation grid for :math:`x_1`
        x2_grid : sequence(float)
            interpolation grid for :math:`x_2`
    """

    def __init__(self, array, q2_grid, x1_grid, x2_grid):
        self._raw = PyImportOnlySubgridV1(
            np.array(array), np.array(q2_grid), np.array(x1_grid), np.array(x2_grid)
        )

class ImportOnlySubgridV2(PyWrapper):
    """
    Python wrapper object to :class:`~pineappl.pineappl.PyImportOnlySubgridV2`.

    Parameters
    ----------
        array : numpy.ndarray(float, dim=3)
            3-dimensional subgrid content
        mu2_grid : sequence(float)
            scale grid
        x1_grid : sequence(float)
            interpolation grid for :math:`x_1`
        x2_grid : sequence(float)
            interpolation grid for :math:`x_2`
    """

    def __init__(self, array, mu2_grid, x1_grid, x2_grid):
        self._raw = PyImportOnlySubgridV2(
            np.array(array), mu2_grid, np.array(x1_grid), np.array(x2_grid)
        )
