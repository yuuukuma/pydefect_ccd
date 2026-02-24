# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.

import pytest
from matplotlib import pyplot as plt

from pydefect_ccd.enum import Carrier
from pydefect_ccd.sommerfeld_scaling import SommerfeldScaling


@pytest.fixture
def sommerfeld_scaling():
    return SommerfeldScaling(epsilon0=8.9,
                             electron_effective_mass=0.5,
                             hole_effective_mass=0.18,
                             Ts=[100, 200, 300, 400, 500])

def test_returns_correct_scaling_for_electrons(sommerfeld_scaling):
    result = sommerfeld_scaling.scaling
    print(result)
    # assert np.array_equal(result, np.array([1.0, 2.0]))
    ax = plt.gca()
    sommerfeld_scaling.plot(ax, keys=[("h", -1)])
    plt.show()


def test_s(sommerfeld_scaling):
    ax = plt.gca()
    for mass in  [0.1, 1.0]:
        # for mass in  [0.1, 0.2, 0.5, 1.0, 2.0]:
        for diele, ls in zip([10, 20, 30], ["-", "--", ":"]):
            for dc in [-1, 0, 1]:
                ss = SommerfeldScaling(diele, mass, 1.0, Ts=[100, 200, 300, 400, 500])
                ss.add_to_ax(ax, carrier_type=Carrier.e, defect_charge=dc, ls=ls)
        ss.set_label(ax)
    plt.show()
