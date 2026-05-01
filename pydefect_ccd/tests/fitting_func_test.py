# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.

import pytest

from pydefect_ccd.fitting_func import QuadraticFittingFunc, QuarticFittingFunc, \
    intersections
from pydefect_ccd.potential_curve import SinglePoints


class FakeSinglePoints:
    def __init__(self, a, dE):
        self._a = a
        self._dE = dE


def test_static_f_scalar():
    assert QuadraticFittingFunc.fitting_func(2.0, 0.1, 3.0) == pytest.approx(12.1)


def test_call_and_omega_and_str_and_shift_for_quadratic():
    q = QuadraticFittingFunc(a=4.0, Q0=1.0, E0=0.5)
    assert q(1.5) == pytest.approx(4.0 * (0.5)**2 + 0.5) # __call__
    assert q.omega == pytest.approx(0.5 * 4.0)
    assert str(q) == "Quadratic Curve: omega=0.129 (eV), Q0=1.000 (amu**0.5*Å), Emin=0.500 (eV)"
    # shift
    new = q.shift(shift_Q=0.2, shift_energy=-0.1)
    assert new.Q0 == pytest.approx(1.2)
    assert new.E0 == pytest.approx(0.4)
    assert new.a == pytest.approx(q.a)


def test_call_and_omega_and_str_and_shift_for_quartic():
    q = QuarticFittingFunc(a=1.0, b=1.0, c=1.0, Q0=1.0, E0=0.5)
    assert q(0.) == pytest.approx(1-1+1+0.5) # __call__
    assert q.omega == pytest.approx(0.5)
    assert str(q) == "QuarticCurve: 1.0*(Q-Q0)^4 + 1.0*(Q-Q0)^3 + 1.0*(Q-Q0)^2 + 0.5 (eV), Q0=1.000 (amu**0.5*Å)"
    # shift
    new = q.shift(shift_Q=-1.0, shift_energy=-0.5, revert=True)
    assert new.Q0 == pytest.approx(0.0)
    assert new.E0 == pytest.approx(0.0)
    assert new.b == pytest.approx(-1.0)


def test_single_intersection():
    q1 = QuadraticFittingFunc(Q0=0.0, E0=0.0, a=1.0)
    q2 = QuadraticFittingFunc(Q0=1.0, E0=0.0, a=1.0)
    res = intersections(q1, q2, [-1.0, 1.0])
    assert len(res) == 1
    q, e = res[0]
    assert q == pytest.approx(0.5)
    assert e == pytest.approx(0.25)


