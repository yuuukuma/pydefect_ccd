# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.

import numpy as np
import pytest

from pydefect_ccd.fitting_func import QuadraticFittingFunc
from pydefect_ccd.potential_curve import SinglePoints


def test_from_single_points():
    class FakeSinglePoint:
        def __init__(self, Q, energy):
            self.Q = Q
            self.ccd_corrected_energy = energy

    sp = SinglePoints([FakeSinglePoint(Q=-1.0, energy=3.0 * 1.0 + 0.25),
                       FakeSinglePoint(Q=0.0, energy=0.25),
                       FakeSinglePoint(Q=1.0, energy=3.0 * 1.0 + 0.25)])
    q = sp.fitting(QuadraticFittingFunc.fitting_func)
    assert q.a == pytest.approx(3.0)
    assert q.Q0 == pytest.approx(0.0)
    assert q.E0 == pytest.approx(0.25)


class FakeSinglePoint:
    def __init__(self, Q, energy):
        self.Q = Q
        self.ccd_corrected_energy = energy


def test_quadratic_from_single_points():
    sp = SinglePoints([FakeSinglePoint(Q=-1.0, energy=3.0 * 1.0 + 0.25),
                       FakeSinglePoint(Q=0.0, energy=0.25),
                       FakeSinglePoint(Q=1.0, energy=3.0 * 1.0 + 0.25)])
    q = QuadraticFittingFunc.from_single_points(sp)
    assert q.a == pytest.approx(3.0)
    assert q.Q0 == pytest.approx(0.0)
    assert q.E0 == pytest.approx(0.25)


def test_quartic_from_single_points():
    sp = SinglePoints([FakeSinglePoint(Q=-1.0, energy=1.0 - 1.0 + 1.0 + 0.25),
                       FakeSinglePoint(Q=0.0, energy=0.25),
                       FakeSinglePoint(Q=1.0, energy=1.0 + 1.0 + 1.0 + 0.25),
                       FakeSinglePoint(Q=2.0, energy=16.0 + 8.0 + 4.0  + 0.25),
                       ])
    q = QuarticFittingFunc.from_single_points(sp)
    assert q.a == pytest.approx(1.0)
    assert q.b == pytest.approx(1.0)
    assert q.c == pytest.approx(1.0)
    assert q.Q0 == pytest.approx(0.0)
    assert q.E0 == pytest.approx(0.25)
