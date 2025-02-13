import pytest
from pineappl.convolutions import Conv, ConvType


class TestConvolutions:
    """Test that the getter methods are returning the exptected values.
    For more realistic tests, see `test_grid`.
    """

    @pytest.mark.parametrize(
        "polarized, time_like",
        [
            (False, False),
            (True, True),
            (True, False),
            (False, True),
        ],
    )
    def test_init(self, polarized: bool, time_like: bool):
        convolution_types = ConvType(polarized=polarized, time_like=time_like)
        convolutions = Conv(convolution_types=convolution_types, pid=2212)

        assert convolution_types.polarized == polarized
        assert convolution_types.time_like == time_like
        assert convolutions.convolution_types.polarized == polarized
        assert convolutions.convolution_types.time_like == time_like
