import numpy as np
import pytest
from pineappl.boc import (
    Bin,
    BinsWithFillLimits,
    Channel,
    Kinematics,
    Order,
    ScaleFuncForm,
    Scales,
)


class TestChannel:
    def test_init(self):
        le = Channel([([2, 2], 0.5)])
        assert isinstance(le, Channel)
        assert le.into_array() == [([2, 2], 0.5)]


class TestKinematics:
    @pytest.mark.parametrize(
        "kintype, argument",
        [
            ("Scale", 0),
            ("Scale", 1),
            ("Scale", 2),
            ("X", 0),
            ("X", 1),
            ("X", 2),
        ],
    )
    def test_init(self, kintype: str, argument: int):
        kin_method = getattr(Kinematics, kintype)
        result = kin_method(argument)
        assert isinstance(result, Kinematics)


class TestScaleFuncForm:
    @pytest.mark.parametrize(
        "scaletype, argument",
        [
            ("NoScale", [0]),
            ("Scale", [0]),
            ("QuadraticSum", [0, 1]),
            ("QuadraticMean", [0, 1]),
            ("QuadraticSumOver4", [0, 1]),
            ("LinearMean", [0, 1]),
            ("LinearSum", [0, 1]),
            ("ScaleMax", [0, 1]),
            ("ScaleMin", [0, 1]),
            ("Prod", [0, 1]),
            ("S2plusS1half", [0, 1]),
            ("Pow4Sum", [0, 1]),
            ("WgtAvg", [0, 1]),
            ("S2plusS1fourth", [0, 1]),
            ("ExpProd2", [0, 1]),
        ],
    )
    def test_init(self, scaletype: ScaleFuncForm, argument: list):
        scale_method = getattr(ScaleFuncForm, scaletype)
        result = scale_method(*argument)
        scale_funcs = Scales(ren=result, fac=result, frg=result)
        assert isinstance(result, ScaleFuncForm)
        assert isinstance(scale_funcs, Scales)


class TestOrder:
    def create_order(self, args: tuple = (2, 1, 0, 1, 0)) -> Order:
        return Order(*args)

    def test_init(self):
        args = (2, 1, 0, 1, 0)
        o = self.create_order(args=args)

        assert isinstance(o, Order)
        assert o.as_tuple() == args

    def test_mask(self):
        o = self.create_order()
        mask = o.create_mask(orders=[o], max_as=2, max_al=1, logs=True)
        assert np.all(mask)


class TestBin:
    @pytest.mark.parametrize(
        "bin_limits, normalization, dimensions, norm, bin_edges",
        [
            ([(1.0, 2.0)], 1.0, 1, 1.0, [(1.0, 2.0)]),
            ([(1.0, 2.0)], 2.5, 1, 2.5, [(1.0, 2.0)]),
            ([(1.0, 2.0), (125, 125)], 2.5, 2, 2.5, [(1.0, 2.0), (125, 125)]),
            (
                [(1.0, 2.0), (125, 125), (0.25, 0.50)],
                2.5,
                3,
                2.5,
                [(1.0, 2.0), (125, 125), (0.25, 0.50)],
            ),
        ],
    )
    def test_bin(
        self,
        bin_limits: list,
        normalization: float,
        dimensions: float,
        norm: int,
        bin_edges: np.ndarray,
    ):
        bin_obj = Bin(bin_limits=bin_limits, normalization=normalization)
        assert bin_obj.dimensions == dimensions
        assert bin_obj.normalization == norm
        np.testing.assert_allclose(bin_obj.bin_limits, bin_edges)


class TestBinsWithFillLimits:
    def test_from_fill_limits(self):
        fill_limits = [float(i) for i in range(6)]
        fill_edges = [
            [(fill_limits[i], fill_limits[i + 1])] for i in range(len(fill_limits) - 1)
        ]
        bin_config = BinsWithFillLimits.from_fill_limits(fill_limits)

        assert isinstance(bin_config.bins()[0], Bin)
        assert bin_config.len() == len(fill_limits) - 1
        assert bin_config.dimensions() == 1
        np.testing.assert_allclose(bin_config.bin_limits(), fill_edges)

    @pytest.mark.parametrize(
        "limits, normalizations, dimensions, slices",
        [
            (
                [[(i, i + 1)] for i in range(6)],
                [(i + 1) / 10 for i in range(6)],
                1,
                [[i for i in range(6)]],
            ),
            (
                [[(i + j, i + 2 * j) for j in range(4)] for i in range(6)],
                [(i + 1) / 10 for i in range(6)],
                4,
                [[i] for i in range(6)],
            ),
        ],
    )
    def test_from_limits_and_normalizations(
        self, limits: list, normalizations: list, dimensions: int, slices: list
    ):
        bin_config = BinsWithFillLimits.from_limits_and_normalizations(
            limits, normalizations
        )

        assert isinstance(bin_config.bins()[0], Bin)
        assert bin_config.len() == len(limits)
        assert bin_config.dimensions() == dimensions
        assert bin_config.slices() == slices
        np.testing.assert_allclose(bin_config.bin_normalizations(), normalizations)
        np.testing.assert_allclose(bin_config.bin_limits(), limits)
