import pineappl


class TestLumiEntry:
    def test_init(self):
        le = pineappl.lumi.LumiEntry([(2, 2, 0.5)])

        assert isinstance(le, pineappl.lumi.LumiEntry)
        assert isinstance(le.raw, pineappl.pineappl.PyLumiEntry)
