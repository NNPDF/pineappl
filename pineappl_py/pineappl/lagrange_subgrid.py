try:
    from .pineappl import PyLagrangeSubgridV2
except:
    import warnings

    warnings.warn("binary files missing")

from .utils import PyWrapper


class LagrangeSubgridV2(PyWrapper):
    def __init__(self, subgrid_params, extra_params):
        self._raw = PyLagrangeSubgridV2(subgrid_params.raw, extra_params.raw)
