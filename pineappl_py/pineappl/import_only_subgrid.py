try:
    from .pineappl import PyImportOnlySubgridV1
except:
    import warnings

    warnings.warn("binary files missing")

from .utils import PyWrapper


class ImportOnlySubgridV1(PyWrapper):
    def __init__(self, array, q2_grid, x1_grid, x2_grid):
        self._raw = PyImportOnlySubgridV1(array, q2_grid, x1_grid, x2_grid)
