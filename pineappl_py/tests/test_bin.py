import numpy as np
import pytest
from pineappl.bin import BinRemapper


class TestBinRemapper:
    def test_init(self):
        br = BinRemapper(np.array([1.0]), [(2, 3)])

        assert isinstance(br, BinRemapper)

        with pytest.raises(AttributeError):
            br._bla()
