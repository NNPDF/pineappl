import pineappl


class TestLumiEntry:
    def test_init(self):
        le = pineappl.LumiEntry([(2, 2, 0.5)])

        assert isinstance(le, pineappl.LumiEntry)
