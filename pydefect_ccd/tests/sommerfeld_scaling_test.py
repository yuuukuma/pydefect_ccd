# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.

import pytest
from matplotlib import pyplot as plt

from pydefect_ccd.sommerfeld_scaling import SommerfeldScaling


@pytest.fixture
def sommerfeld_scaling():
    return SommerfeldScaling(dielectric_constant=8.9, electron_effective_mass=0.5, hole_effective_mass=0.18, max_defect_charge=1)

def test_returns_correct_scaling_for_electrons(sommerfeld_scaling):
    result = sommerfeld_scaling.scaling
    print(result)
    # assert np.array_equal(result, np.array([1.0, 2.0]))
    ax = plt.gca()
    sommerfeld_scaling.plot(ax, keys=[("h", 1)])
    plt.show()