import numpy as np
import pytest
from pineappl.boc import Channel, Kinematics, Order, ScaleFuncForm


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
        assert isinstance(result, ScaleFuncForm)


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
