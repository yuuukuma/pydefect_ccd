# -*- coding: utf-8 -*-
#  Copyright (c) 2024 Kumagai group.

import numpy as np
import pytest

from pydefect_ccd.fitting_curve import QuadraticFittingCurve, QuarticFittingCurve, \
    intersections
from pydefect_ccd.potential_curve import SinglePoints


class FakeSinglePoints:
    def __init__(self, a, dE):
        self._a = a
        self._dE = dE


def test_static_f_scalar():
    assert QuadraticFittingCurve.fitting_func(2.0, 0.1, 3.0) == pytest.approx(12.1)


def test_call_and_omega_and_str_and_shift_for_quadratic():
    q = QuadraticFittingCurve(a=4.0, Q0=1.0, E0=0.5)
    assert q(1.5) == pytest.approx(4.0 * (0.5)**2 + 0.5) # __call__
    assert q.omega == pytest.approx(0.5 * 4.0)
    assert str(q) == "Quadratic Curve: omega=0.129 (eV), Q0=1.000 (amu**0.5*Å), Emin=0.500 (eV)"
    # shift
    new = q.shift(shift_Q=0.2, shift_energy=-0.1)
    assert new.Q0 == pytest.approx(1.2)
    assert new.E0 == pytest.approx(0.4)
    assert new.a == pytest.approx(q.a)


def test_call_and_omega_and_str_and_shift_for_quartic():
    q = QuarticFittingCurve(a=1.0, b=1.0, c=1.0, Q0=1.0, E0=0.5)
    assert q(0.) == pytest.approx(1-1+1+0.5) # __call__
    assert q.omega == pytest.approx(0.5)
    assert str(q) == "QuarticCurve: 1.0*(Q-Q0)^4 + 1.0*(Q-Q0)^3 + 1.0*(Q-Q0)^2 + 0.5 (eV), Q0=1.000 (amu**0.5*Å)"
    # shift
    new = q.shift(shift_Q=-1.0, shift_energy=-0.5, revert=True)
    assert new.Q0 == pytest.approx(0.0)
    assert new.E0 == pytest.approx(0.0)
    assert new.b == pytest.approx(-1.0)


class FakeSinglePoint:
    def __init__(self, Q, energy):
        self.Q = Q
        self.ccd_corrected_energy = energy


def test_quadratic_from_single_points():
    sp = SinglePoints([FakeSinglePoint(Q=-1.0, energy=3.0 * 1.0 + 0.25),
                       FakeSinglePoint(Q=0.0, energy=0.25),
                       FakeSinglePoint(Q=1.0, energy=3.0 * 1.0 + 0.25)])
    q = QuadraticFittingCurve.from_single_points(sp)
    assert q.a == pytest.approx(3.0)
    assert q.Q0 == pytest.approx(0.0)
    assert q.E0 == pytest.approx(0.25)


def test_quartic_from_single_points():
    sp = SinglePoints([FakeSinglePoint(Q=-1.0, energy=1.0 - 1.0 + 1.0 + 0.25),
                       FakeSinglePoint(Q=0.0, energy=0.25),
                       FakeSinglePoint(Q=1.0, energy=1.0 + 1.0 + 1.0 + 0.25),
                       FakeSinglePoint(Q=2.0, energy=16.0 + 8.0 + 4.0  + 0.25),
                       ])
    q = QuarticFittingCurve.from_single_points(sp)
    assert q.a == pytest.approx(1.0)
    assert q.b == pytest.approx(1.0)
    assert q.c == pytest.approx(1.0)
    assert q.Q0 == pytest.approx(0.0)
    assert q.E0 == pytest.approx(0.25)


def test_single_intersection():
    q1 = QuadraticFittingCurve(Q0=0.0, E0=0.0, a=1.0)
    q2 = QuadraticFittingCurve(Q0=1.0, E0=0.0, a=1.0)
    res = intersections(q1, q2, [-1.0, 1.0])
    assert len(res) == 1
    q, e = res[0]
    assert q == pytest.approx(0.5)
    assert e == pytest.approx(0.25)


