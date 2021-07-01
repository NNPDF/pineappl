import pineappl


class TestSubgridParams:
    def test_init(self):
        sp = pineappl.subgrid.SubgridParams()

        assert isinstance(sp, pineappl.subgrid.SubgridParams)
        assert isinstance(sp.raw, pineappl.pineappl.PySubgridParams)
