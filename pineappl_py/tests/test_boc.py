import numpy as np
from pineappl.boc import Channel, Kinematics, Order


class TestChannel:
    def test_init(self):
        le = Channel([([2, 2], 0.5)])
        assert isinstance(le, Channel)
        assert le.into_array() == [([2, 2], 0.5)]


class TestKinematics:
    def test_init(self):
        kin = Kinematics(0)
        assert isinstance(kin, Kinematics)


class TestOrder:
    def create_order(self, args=(2, 1, 0, 1, 0)):
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
