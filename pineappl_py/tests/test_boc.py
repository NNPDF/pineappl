import pineappl


class TestChannel:
    def test_init(self):
        le = pineappl.boc.Channel([(2, 2, 0.5)])
        assert isinstance(le, pineappl.boc.Channel)
