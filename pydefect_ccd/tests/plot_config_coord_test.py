# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.

import pytest

from pydefect_ccd.ccd import (Ccd, SinglePoint, CcdPlotter,
                              PotentialCurve)
from pydefect_ccd.local_enum import Carrier


@pytest.fixture
def single_ccd():
    return PotentialCurve(
        id_=CcdId("from_0_to_1", carriers=[Carrier.h, Carrier.e]),
        charge=0,
        single_points=[
            SinglePoint(2., 1.0, 3.3, False, used_for_fitting=True),
            SinglePoint(1., 0.5, 2.2, False, used_for_fitting=True),
            SinglePoint(0., 0.0, 3.4, False, used_for_fitting=True),
            SinglePoint(3., 1.5, 3.5, False, used_for_fitting=False)])


@pytest.fixture
def ccd(single_ccd):
    excited_state = PotentialCurve(
        id_=CcdId("from_1_to_0", carriers=[Carrier.h]),
        charge=1,
        single_points=[SinglePoint(3., 1.0, 10.1, False),
                       SinglePoint(2., 0.9, 10.2, False),
                       SinglePoint(1., 0.8, 10.3, False)])
    return Ccd(name="Va_O1", potential_curve_results=[single_ccd, excited_state])


def test_plot_ccd(ccd):
    plotter = CcdPlotter(ccd)
    plotter.construct_plot()
    plotter.plt.show()
