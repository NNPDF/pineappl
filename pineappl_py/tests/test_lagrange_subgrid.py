import pineappl


class TestLagrangeSubgridV2:
    def test_init(self):
        lsv2 = pineappl.lagrange_subgrid.LagrangeSubgridV2(
            pineappl.subgrid.SubgridParams(), pineappl.subgrid.ExtraSubgridParams()
        )

        assert isinstance(lsv2, pineappl.lagrange_subgrid.LagrangeSubgridV2)
        assert isinstance(lsv2.raw, pineappl.pineappl.PyLagrangeSubgridV2)
