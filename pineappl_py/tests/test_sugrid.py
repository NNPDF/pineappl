import pineappl
import pytest

import numpy as np


class TestSubgridParams:
    def test_init(self):
        sp = pineappl.subgrid.SubgridParams()

        assert isinstance(sp, pineappl.subgrid.SubgridParams)
        assert isinstance(sp.raw, pineappl.pineappl.PySubgridParams)


def test_issue_164(pdf):
    luminosities = [pineappl.lumi.LumiEntry([(1, 2, 1.0)])]
    orders = [pineappl.grid.Order(0, 0, 0, 0)]
    params = pineappl.subgrid.SubgridParams()

    def convolve_grid():
        grid = pineappl.grid.Grid.create(luminosities, orders, [0.0, 1.0], params)
        grid.fill(0.2, 0.2, 10, 0, 0.5, 0, 0.5)
        return grid.convolve_with_one(2212, pdf.xfxQ, pdf.alphasQ)

    # default minimum is q2=100
    res = convolve_grid()
    assert res == 0.0

    # lower minimum to q2=1
    params.set_q2_min(1.0)
    res = convolve_grid()
    assert pytest.approx(res) != 0.0

class TestSubgrid:
    def fake_grid(self):
        luminosities = [pineappl.lumi.LumiEntry([(1, 2, 1.0)])]
        orders = [pineappl.grid.Order(0, 0, 0, 0)]
        params = pineappl.subgrid.SubgridParams()
        grid = pineappl.grid.Grid.create(luminosities, orders, [0.0, 1.0], params)
        return grid
    
    def fake_importonlysubgrid(self):
        x1s = np.linspace(0.1, 1, 2)
        x2s = np.linspace(0.5, 1, 2)
        Q2s = np.linspace(10, 20, 2)
        mu2s = [tuple([q2, q2]) for q2 in Q2s]
        array = np.random.rand(len(Q2s), len(x1s), len(x2s))
        subgrid = pineappl.import_only_subgrid.ImportOnlySubgridV2(array, mu2s , x1s, x2s)
        return subgrid, [x1s, x2s, mu2s, array]

    def test_subgrid_methods(self):
        grid = self.fake_grid()
        test_subgrid, infos = self.fake_importonlysubgrid()
        x1s, x2s, mu2s, _ = (obj for obj in infos)
        grid.set_subgrid(0,0,0, test_subgrid)
        extr_subgrid = grid.subgrid(0,0,0)
        facgrid = np.array([mu2.fac for mu2 in extr_subgrid.mu2_grid()])
        rengrid = np.array([mu2.ren for mu2 in extr_subgrid.mu2_grid()])
        np.testing.assert_allclose([mu2[0] for mu2 in mu2s], rengrid)
        np.testing.assert_allclose([mu2[1] for mu2 in mu2s], facgrid)
        np.testing.assert_allclose(extr_subgrid.x1_grid(), x1s)
        np.testing.assert_allclose(extr_subgrid.x2_grid(), x2s)
    
    def test_to_array3(self):
        grid = self.fake_grid()
        test_subgrid, infos = self.fake_importonlysubgrid()
        _, _, _, array = (obj for obj in infos)
        grid.set_subgrid(0,0,0, test_subgrid)
        extr_subgrid = grid.subgrid(0,0,0)
        test_array = extr_subgrid.to_array3()
        print(test_array)
        print(array)
        np.testing.assert_allclose(test_array, array)
