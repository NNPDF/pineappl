from dataclasses import dataclass, field
import numpy as np
import pytest
from typing import List, Optional

from pineappl.boc import (
    Bin,
    BinsWithFillLimits,
    Channel,
    Kinematics,
    Order,
    ScaleFuncForm,
    Scales,
)


@dataclass
class BinFields:
    bin_limits: List
    normalization: float
    dimensions: int


@dataclass
class BwflFields:
    bins: Optional[Bin] = Bin([(1, 2)], 1.0)
    fill_limits: Optional[List] = field(default_factory=list)
    limits: Optional[List] = field(default_factory=list)
    normalizations: Optional[List] = field(default_factory=list)
    dimensions: Optional[int] = 1


def _generate_bin_fields(n_dimensions: int, normalization: float) -> BinFields:
    bin_limits = [(i, i + 1) for i in range(n_dimensions)]
    return BinFields(
        bin_limits=bin_limits, normalization=normalization, dimensions=n_dimensions
    )


def _generate_new_bwfl(
    n_bins: int, n_dimensions: int, normalization: int
) -> BwflFields:
    limits = [[(i + j, i + 2 * j) for j in range(n_dimensions)] for i in range(n_bins)]
    bins = [Bin(lim, normalization) for lim in limits]
    fill_limits = [i for i in range(n_bins + 1)]
    normalizations = [normalization for _ in range(n_bins)]
    return BwflFields(
        bins=bins, fill_limits=fill_limits, normalizations=normalizations, limits=limits
    )


def _generated_bwfl_fields(n_bins: int, n_dimensions: int) -> BwflFields:
    limits = [[(i + j, i + 2 * j) for j in range(n_dimensions)] for i in range(n_bins)]
    normalizations = [(i + 1) / 10 for i in range(n_bins)]
    return BwflFields(
        limits=limits, normalizations=normalizations, dimensions=n_dimensions
    )


class TestChannel:
    def test_init(self):
        channel = Channel([([2, -2], 0.5)])
        assert isinstance(channel, Channel)
        assert channel.into_array() == [([2, -2], 0.5)]

        channels = Channel([([2, -2], 0.5), ([3, -3], 1.5)])
        assert isinstance(channels, Channel)
        assert channels.into_array() == [([2, -2], 0.5), ([3, -3], 1.5)]


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
    @pytest.mark.parametrize("n_dimensions, normalization", [(1, 1.0), (4, 0.25)])
    def test_bin(self, n_dimensions: int, normalization: float):
        bin_field = _generate_bin_fields(n_dimensions, normalization)
        bin_obj = Bin(bin_field.bin_limits, bin_field.normalization)
        assert bin_obj.dimensions == n_dimensions
        assert bin_obj.normalization == bin_field.normalization
        np.testing.assert_allclose(bin_obj.bin_limits, bin_field.bin_limits)


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

        index_to_remove = 0
        removed_bin = bin_config.removed_index(index=0)
        np.testing.assert_allclose(
            removed_bin.bin_limits, bin_config.bin_limits()[index_to_remove]
        )

    @pytest.mark.parametrize(
        "n_bins, n_dimensions, slices",
        [(6, 1, [[i for i in range(6)]]), (6, 4, [[i] for i in range(6)])],
    )
    def test_from_limits_and_normalizations(
        self, n_bins: int, n_dimensions: int, slices: List
    ):
        bwfl_field = _generated_bwfl_fields(n_bins, n_dimensions)
        bin_config = BinsWithFillLimits.from_limits_and_normalizations(
            bwfl_field.limits, bwfl_field.normalizations
        )

        assert isinstance(bin_config.bins()[0], Bin)
        assert bin_config.len() == n_bins
        assert bin_config.dimensions() == n_dimensions
        assert bin_config.slices() == slices
        np.testing.assert_allclose(
            bin_config.bin_normalizations(), bwfl_field.normalizations
        )
        np.testing.assert_allclose(bin_config.bin_limits(), bwfl_field.limits)

    @pytest.mark.parametrize("n_bins, n_dimensions, normalization", [(4, 1, 10.0)])
    def test_new_bwfl(self, n_bins: int, n_dimensions: int, normalization: int):
        bwfl_field = _generate_new_bwfl(n_bins, n_dimensions, normalization)
        bin_config = BinsWithFillLimits(bwfl_field.bins, bwfl_field.fill_limits)

        assert isinstance(bin_config.bins()[0], Bin)
        assert bin_config.len() == n_bins
        assert bin_config.dimensions() == n_dimensions
        assert bin_config.slices() == [[i for i in range(n_bins)]]
        np.testing.assert_allclose(
            bin_config.bin_normalizations(), bwfl_field.normalizations
        )
        np.testing.assert_allclose(bin_config.bin_limits(), bwfl_field.limits)
