# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.

import pytest
from matplotlib import pyplot as plt

from pydefect_ccd.sommerfeld_scaling import SommerfeldScaling


@pytest.fixture
def sommerfeld_scaling():
    return SommerfeldScaling(dielectric_constant=8.9, electron_effective_mass=0.5,
                             hole_effective_mass=0.18)

def test_returns_correct_scaling_for_electrons(sommerfeld_scaling):
    result = sommerfeld_scaling.scaling
    print(result)
    # assert np.array_equal(result, np.array([1.0, 2.0]))
    ax = plt.gca()
    sommerfeld_scaling.plot(ax, keys=[("h", -1)])
    plt.show()


def test_s():
    ax = plt.gca()
    for mass in  [0.1, 1.0]:
        # for mass in  [0.1, 0.2, 0.5, 1.0, 2.0]:
        for diele, ls in zip([10, 20, 30], ["-", "--", ":"]):
            ss = SommerfeldScaling(dielectric_constant=diele, electron_effective_mass=0.5,
                             hole_effective_mass=mass)
            ss.plot(ax, key=("h", 1), ls=ls)
    plt.show()
