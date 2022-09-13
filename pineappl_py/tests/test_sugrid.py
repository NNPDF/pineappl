import pineappl
import pytest


class TestSubgridParams:
    def test_init(self):
        sp = pineappl.subgrid.SubgridParams()

        assert isinstance(sp, pineappl.subgrid.SubgridParams)
        assert isinstance(sp.raw, pineappl.pineappl.PySubgridParams)


def test_issue_164(pdf):
    luminosities = [pineappl.lumi.LumiEntry([(1, 2, 1.0)])]
    orders = [pineappl.grid.Order(0, 0, 0, 0)]
    params = pineappl.subgrid.SubgridParams()

    def convolute_grid():
        grid = pineappl.grid.Grid.create(luminosities, orders, [0.0, 1.0], params)
        grid.fill(0.2, 0.2, 10, 0, 0.5, 0, 0.5)
        return grid.convolute_with_one(2212, pdf.xfxQ, pdf.alphasQ)

    # default minimum is q2=100
    res = convolute_grid()
    assert res == 0.0

    # lower minimum to q2=1
    params.set_q2_min(1.0)
    res = convolute_grid()
    assert pytest.approx(res) != 0.0
