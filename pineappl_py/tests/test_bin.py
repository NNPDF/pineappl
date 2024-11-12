import numpy as np
import pytest

from pineappl.boc import Channel, Order
from pineappl.bin import BinRemapper
from pineappl.convolutions import Conv, ConvType


class TestBinRemapper:
    def test_init(self):
        br = BinRemapper(np.array([1.0]), [(2, 3)])

        assert isinstance(br, BinRemapper)

        with pytest.raises(AttributeError):
            br._bla()

    def test_binremapper(self, fake_grids):
        h = ConvType(polarized=True, time_like=False)
        h_conv = Conv(conv_type=h, pid=2212)
        convolutions = [h_conv]

        down_channel = [([1], 1.0)]  # DIS-case
        up_channel = [([2], 1.0)]  # DIS-case
        channels = [Channel(down_channel), Channel(up_channel)]

        orders = [Order(3, 0, 0, 0, 0)]
        g = fake_grids.grid_with_generic_convolution(
            nb_convolutions=1,
            channels=channels,
            orders=orders,
            convolutions=convolutions,
            bins=np.linspace(1e-2, 1, num=20),
        )

        # Extract the left & right bin limits and redefine the normalization
        bin_dims = g.bin_dimensions()
        bin_limits = [
            (left, right)
            for left, right in zip(
                g.bin_left(bin_dims - 1), g.bin_right(bin_dims - 1)
            )
        ]
        normalizations = [10.0 for _ in g.bin_normalizations()]

        remapper = BinRemapper(np.array(normalizations), bin_limits)
        # Modify the bin normalization
        g.set_remapper(remapper)
        new_normalizations = g.bin_normalizations()

        # Check that the bin normalizations have been updated
        np.testing.assert_allclose(new_normalizations, normalizations)
