# -*- coding: utf-8 -*-
#  Copyright (c) 2024 Kumagai group.

import numpy as np
import pytest

from pydefect_ccd.fitting_curve import QuadraticFittingCurve
from pydefect_ccd.potential_curve import SinglePoints, SinglePoint


class FakeSinglePoints:
    def __init__(self, a, dE):
        self._a = a
        self._dE = dE

    def fitting(self, func):
        # emulate SinglePoints.fitting returning (a, dE)
        return self._a, self._dE


def test_static_f_scalar():
    assert QuadraticFittingCurve.fitting_func(2.0, 0.1, 3.0) == pytest.approx(12.1)


def test_static_f_array():
    arr = np.array([0.0, 1.0, 2.0])
    res = QuadraticFittingCurve.fitting_func(arr, 0.0, 2.0)
    assert np.allclose(res, 2.0 * arr**2)


def test_call_and_omega_and_shift():
    q = QuadraticFittingCurve(a=4.0, Q0=1.0, dE=0.5)
    # __call__
    assert q(1.5) == pytest.approx(4.0 * (0.5)**2 + 0.5)
    # omega property
    assert q.omega == pytest.approx(0.5 * 4.0)
    # shift
    new = q.shift(shift_Q=0.2, shift_energy=-0.1)
    assert new.Q0 == pytest.approx(1.2)
    assert new.dE == pytest.approx(0.4)
    assert new.a == pytest.approx(q.a)

def test_from_single_points():
    class FakeSinglePoint:
        def __init__(self, Q, energy):
            self.Q = Q
            self.ccd_corrected_energy = energy

    sp = SinglePoints([FakeSinglePoint(Q=-1.0, energy=3.0 * 1.0 + 0.25),
                       FakeSinglePoint(Q=0.0, energy=0.25),
                       FakeSinglePoint(Q=1.0, energy=3.0 * 1.0 + 0.25)])
    q = QuadraticFittingCurve.from_single_points(sp)
    assert q.a == pytest.approx(3.0)
    assert q.Q0 == pytest.approx(0.0)
    assert q.dE == pytest.approx(0.25)


