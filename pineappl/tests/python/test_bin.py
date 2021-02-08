import pineappl


class TestBinRemapper:
    def test_init(self):
        br = pineappl.bin.BinRemapper([1], [(2, 3)])

        assert isinstance(br, pineappl.bin.BinRemapper)
        assert isinstance(br.raw, pineappl.pineappl.PyBinRemapper)
