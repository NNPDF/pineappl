import pytest


class PDF:
    def xfxQ(self, pid, x, q):
        return self.xfxQ2(pid, x, q**2)

    def xfxQ2(self, pid, x, q2):
        if pid in range(-6, 6):
            return x * (1 - x)
        else:
            return 0.0

    def alphasQ(self, q):
        return 1.0


@pytest.fixture
def pdf():
    return PDF()
