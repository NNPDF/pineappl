import pineappl


class TestSubgridParams:
    def test_init(self):
        sp = pineappl.subgrid.SubgridParams()

        assert isinstance(sp, pineappl.subgrid.SubgridParams)
        assert isinstance(sp.raw, pineappl.pineappl.PySubgridParams)


class TestExtraSubgridParams:
    def test_init(self):
        esp = pineappl.subgrid.ExtraSubgridParams()

        assert isinstance(esp, pineappl.subgrid.ExtraSubgridParams)
        assert isinstance(esp.raw, pineappl.pineappl.PyExtraSubgridParams)
