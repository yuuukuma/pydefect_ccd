# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.

import pytest
import numpy as np

from pydefect_ccd.local_enum import Carrier
from pydefect_ccd.sommerfeld_scaling import SommerfeldScaling


@pytest.fixture
def sommerfeld_scaling():
    return SommerfeldScaling(epsilon0=8.9,
                             electron_effective_mass=0.5,
                             hole_effective_mass=0.18,
                             Ts=[100, 200, 300, 400, 500])


def test_returns_correct_scaling_for_electrons(sommerfeld_scaling, monkeypatch):
    expected = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    calls = []

    def fake_sommerfeld_parameter(Ts, Z, mass, epsilon0, method):
        calls.append((Ts, Z, mass, epsilon0, method))
        return expected

    monkeypatch.setattr(
        "pydefect_ccd.sommerfeld_scaling.sommerfeld_parameter",
        fake_sommerfeld_parameter)

    actual = sommerfeld_scaling.scaling(Carrier.e, defect_charge=-1)

    assert actual is expected
    np.testing.assert_array_equal(calls[0][0],
                                  np.array(sommerfeld_scaling.Ts))
    assert calls[0][1:] == (1, 0.5, 8.9, "Integrate")
    assert sommerfeld_scaling.scaling(Carrier.e, defect_charge=-1) is expected
    assert len(calls) == 1
