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
        "kintype, argument, expected_type",
        [
            ("Scale", 0, Kinematics),
            ("Scale", 1, Kinematics),
            ("Scale", 2, Kinematics),
            ("X", 0, Kinematics),
            ("X", 1, Kinematics),
            ("X", 2, Kinematics),
        ],
    )
    def test_init(self, kintype: str, argument: int, expected_type: Kinematics):
        kin_method = getattr(Kinematics, kintype)
        result = kin_method(argument)
        assert isinstance(result, expected_type)


class TestScaleFuncForm:
    @pytest.mark.parametrize(
        "scaletype, argument, expected_type",
        [
            ("NoScale", [0], ScaleFuncForm),
            ("Scale", [0], ScaleFuncForm),
            ("QuadraticSum", [0, 1], ScaleFuncForm),
            ("QuadraticMean", [0, 1], ScaleFuncForm),
            ("QuadraticSumOver4", [0, 1], ScaleFuncForm),
            ("LinearMean", [0, 1], ScaleFuncForm),
            ("LinearSum", [0, 1], ScaleFuncForm),
            ("ScaleMax", [0, 1], ScaleFuncForm),
            ("ScaleMin", [0, 1], ScaleFuncForm),
            ("Prod", [0, 1], ScaleFuncForm),
            ("S2plusS1half", [0, 1], ScaleFuncForm),
            ("Pow4Sum", [0, 1], ScaleFuncForm),
            ("WgtAvg", [0, 1], ScaleFuncForm),
            ("S2plusS1fourth", [0, 1], ScaleFuncForm),
            ("ExpProd2", [0, 1], ScaleFuncForm),
        ],
    )
    def test_init(
        self,
        scaletype: ScaleFuncForm,
        argument: list,
        expected_type: ScaleFuncForm,
    ):
        scale_method = getattr(ScaleFuncForm, scaletype)
        result = scale_method(*argument)
        assert isinstance(result, expected_type)


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
