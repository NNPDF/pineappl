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
        conv_type = ConvType(polarized=polarized, time_like=time_like)
        convolutions = Conv(conv_type=conv_type, pid=2212)

        assert conv_type.polarized == polarized
        assert conv_type.time_like == time_like
        assert convolutions.conv_type.polarized == polarized
        assert convolutions.conv_type.time_like == time_like
