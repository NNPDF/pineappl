try:
    from .pineappl import PyExtraSubgridParams, PySubgridParams
except:
    import warnings

    warnings.warn("binary files missing")

from .utils import PyWrapper


class SubgridParams(PyWrapper):
    def __init__(self):
        self._raw = PySubgridParams()

    # def set_x_bins(self, x_bins):
    # self._raw.set_x_bins(x_bins)

    # def set_x_bins(self, x_bins):
    # self._raw.set_x_bins(x_bins)


class ExtraSubgridParams(PyWrapper):
    def __init__(self):
        self._raw = PyExtraSubgridParams()
