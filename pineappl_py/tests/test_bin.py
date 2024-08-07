import numpy as np
import pineappl
import pytest


class TestBinRemapper:
    def test_init(self):
        br = pineappl.bin.BinRemapper(np.array([1.0]), [(2, 3)])

        assert isinstance(br, pineappl.bin.BinRemapper)

        with pytest.raises(AttributeError):
            br._bla()
