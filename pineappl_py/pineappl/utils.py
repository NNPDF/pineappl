class PyWrapper:
    """
    Python wrapper helper to delegate function calls to the underlying
    raw object.
    """

    _raw = None

    @property
    def raw(self):
        """Raw PyO3 object"""
        return self._raw

    def __getattr__(self, name):
        """Delegate function calls down."""
        if name[0] != "_":
            return self._raw.__getattribute__(name)
        else:
            raise AttributeError
