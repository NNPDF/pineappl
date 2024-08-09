import pineappl


class TestLumiEntry:
    def test_init(self):
        le = pineappl.channel.Channel([(2, 2, 0.5)])
        assert isinstance(le, pineappl.channel.Channel)
