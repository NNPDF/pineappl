class PyWrapper:
    """
    Python wrapper helper to delegate function calls to the underlying
    raw object.
    """
    _raw = None

    @property
    def raw(self):
        return self._raw

    def __getattr__(self, name):
        if name[0] != "_":
            return self._raw.__getattribute__(name)
        else:
            raise AttributeError
