# -*- coding: utf-8 -*-
#  Copyright (c) 2024 Kumagai group.

import numpy as np
import pytest

from pydefect_ccd.fitting_curve import QuadraticFittingCurve
from pydefect_ccd.potential_curve import SinglePoints


def test_from_single_points():
    class FakeSinglePoint:
        def __init__(self, Q, energy):
            self.Q = Q
            self.ccd_corrected_energy = energy

    sp = SinglePoints([FakeSinglePoint(Q=-1.0, energy=3.0 * 1.0 + 0.25),
                       FakeSinglePoint(Q=0.0, energy=0.25),
                       FakeSinglePoint(Q=1.0, energy=3.0 * 1.0 + 0.25)])
    q = sp.fitting(QuadraticFittingCurve.fitting_func)
    assert q.a == pytest.approx(3.0)
    assert q.Q0 == pytest.approx(0.0)
    assert q.dE == pytest.approx(0.25)
