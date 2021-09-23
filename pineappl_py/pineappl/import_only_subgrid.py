from .pineappl import PyImportOnlySubgridV1
from .utils import PyWrapper


class ImportOnlySubgridV1(PyWrapper):
    """
    Python wrapper object to :class:`~pineappl.pineappl.PyImportOnlySubgridV1`.

    Parameters
    ----------
        array : numpy.ndarray
            3-dimensional subgrid content
        q2_grid : list(float)
            scale grid
        x1_grid : list(float)
            interpolation grid for :math:`x_1`
        x2_grid : list(float)
            interpolation grid for :math:`x_2`
    """

    def __init__(self, array, q2_grid, x1_grid, x2_grid):
        self._raw = PyImportOnlySubgridV1(array, q2_grid, x1_grid, x2_grid)
