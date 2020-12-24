import pineappl


class TestOrder:
    def test_init(self):
        o = pineappl.grid.Order(2, 1, 0, 1)

        assert isinstance(o, pineappl.grid.Order)
        assert isinstance(o.raw, pineappl.pineappl.PyOrder)


class TestGrid:
    def test_init(self):
        lumi = []
        orders = [pineappl.grid.Order(3, 0, 1, 0)]
        bin_limits = [1e-7, 1e-3, 1]
        subgrid_params = pineappl.subgrid.SubgridParams()
        g = pineappl.grid.Grid(None, orders, bin_limits, subgrid_params)

        assert isinstance(g, pineappl.grid.Grid)
        assert isinstance(g.raw, pineappl.pineappl.PyGrid)
