from .pineappl import PyChannel
from .utils import PyWrapper


class Channel(PyWrapper):
    """
    Python wrapper object to :class:`~pineappl.pineappl.PyChannel`.

    Parameters
    ----------
    entry : list(tuple(int,int,float))
        sequence describing a channel combination.
    """

    def __init__(self, entry):
        self._raw = PyChannel(entry)
