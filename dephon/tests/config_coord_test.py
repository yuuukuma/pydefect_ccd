# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from copy import deepcopy

import numpy as np
import pytest
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from vise.tests.helpers.assertion import assert_dataclass_almost_equal

from dephon.config_coord import ConfigCoordDiagram, SinglePoint, CcdPlotter, \
    PotentialCurve, spline3, captured_carrier, CarrierDiffError
from dephon.enum import CorrectionType, Carrier


@pytest.fixture
def single_point_min_info():
    return SinglePoint(dQ=1.0, disp_ratio=0.1)


@pytest.fixture
def single_point_max_info():
    orb_info = LocalizedOrbital(band_idx=2,
                                ave_energy=2.0,
                                occupation=1.0,
                                orbitals={"O": [0.0, 1.0, 0.0]})

    return SinglePoint(dQ=1.0,
                       disp_ratio=0.1,
                       corrected_energy=5.0,
                       magnetization=1.0,
                       localized_orbitals=[[], [orb_info]],
                       is_shallow=True,
                       correction_method=CorrectionType.extended_FNV,
                       used_for_fitting=True,
                       base_energy=1.0)


def test_single_point_info_corrected_energy(single_point_min_info,
                                            single_point_max_info):
    assert single_point_min_info.relative_energy is None
    assert single_point_max_info.relative_energy == 5.0 - 1.0


def test_image_structure_str(single_point_min_info,
                             single_point_max_info):
    actual = single_point_min_info.__str__()
    expected = """   dQ    disp ratio  corr. energy    relative energy    used for fitting?    is shallow?    localized orb
1.000         0.100"""
    assert actual == expected

    actual = single_point_max_info.__str__()
    expected = """   dQ    disp ratio    corr. energy    relative energy  used for fitting?    is shallow?    localized orb
1.000         0.100           5.000              4.000  True                 True           down-2(1.0)"""
    assert actual == expected


@pytest.fixture
def potential_curve():
    return PotentialCurve(
        name="from_0_to_1",
        carriers=[Carrier.h, Carrier.e],
        charge=0,
        points=[
            SinglePoint(2., 1.0, 3.3, is_shallow=False, used_for_fitting=True),
            SinglePoint(1., 0.5, 2.2, is_shallow=False, used_for_fitting=True),
            SinglePoint(0., 0.0, 3.4, is_shallow=False, used_for_fitting=True),
            SinglePoint(3., 1.5, 3.5, is_shallow=False, used_for_fitting=False)])


def test_single_ccd_sort_single_point_infos(potential_curve):
    actual = potential_curve.points[0]
    expected = SinglePoint(0., 0.0, 3.4, is_shallow=False,
                           used_for_fitting=True, base_energy=2.2)
    assert actual == expected


def test_single_ccd_set_quadratic_fitting_range(potential_curve):
    potential_curve.set_quadratic_fitting_range(q_range=[-0.1, 0.1])
    actual = [p_info.used_for_fitting for p_info in potential_curve.points]
    expected = [True, False, False, False]
    assert actual == expected

    potential_curve.set_quadratic_fitting_range()
    actual = [p_info.used_for_fitting for p_info in potential_curve.points]
    expected = [True, True, True, True]
    assert actual == expected


def test_single_ccd_dQs_and_energies(potential_curve):
    potential_curve.set_base_energy(0.1)

    actual = potential_curve.dQs_and_energies(False)
    dQs = [0.0, 1.0, 2.0, 3.0]
    energies = [3.3, 2.1, 3.2, 3.4]
    expected = (dQs, energies)
    assert np.array(actual) == pytest.approx(np.array(expected))

    potential_curve.set_base_energy()

    actual = potential_curve.dQs_and_energies(only_used_for_fitting=True)
    dQs = [0.0, 1.0, 2.0]
    energies = [0.0, -1.2, -0.1]
    expected = (dQs, energies)
    assert np.array(actual) == pytest.approx(np.array(expected))


def test_single_ccd_omega(potential_curve):
    assert potential_curve.omega() == pytest.approx(0.09805287531186566)


def test_energy_shifted_single_ccd(potential_curve):
    actual = deepcopy(potential_curve)
    actual.shift_energy(energy=1.0)
    assert actual.points[0].corrected_energy == 4.4


def test_dQ_reverted_single_ccd(potential_curve):
    actual = potential_curve.dQ_reverted_single_ccd()
    expected = deepcopy(potential_curve)
    expected.points = [
        SinglePoint(-1.0, 1.5, 3.5, is_shallow=False, used_for_fitting=False),
        SinglePoint(0.0, 1.0, 3.3, is_shallow=False, used_for_fitting=True),
        SinglePoint(1.0, 0.5, 2.2, is_shallow=False, used_for_fitting=True),
        SinglePoint(2.0, 0.0, 3.4, is_shallow=False, used_for_fitting=True)]
    assert_dataclass_almost_equal(actual, expected)


def test_single_ccd_str(potential_curve):
    potential_curve.set_base_energy()
    actual = potential_curve.__str__()
    expected = """name: from_0_to_1
charge: 0
omega: 0.098
carriers: h e
   dQ    disp ratio    corr. energy    relative energy  used for fitting?    is shallow?    localized orb
0.000         0.000           3.400              0.000  True                 False
1.000         0.500           2.200             -1.200  True                 False
2.000         1.000           3.300             -0.100  True                 False
3.000         1.500           3.500              0.100  False                False"""
    assert actual == expected


@pytest.fixture
def excited_ccd():
    return PotentialCurve(
        CcdId(name="from_1_to_0", carriers=[Carrier.e]),
        charge=1,
        points=[SinglePoint(3., 1.0, 10.1, is_shallow=False),
                SinglePoint(2., 0.9, 10.2, is_shallow=False),
                SinglePoint(1., 0.8, 10.3, is_shallow=False)])


@pytest.fixture
def ccd(potential_curve, excited_ccd):
    return ConfigCoordDiagram(name="Va_O1_1 ⇆ Va_O1_0", potential_curves=[potential_curve, excited_ccd])


def test_captured_carrier(potential_curve, excited_ccd):
    assert captured_carrier(potential_curve, excited_ccd) == Carrier.h
    with pytest.raises(CarrierDiffError):
        captured_carrier(excited_ccd, potential_curve)


def test_ccd_single_ccd(ccd, potential_curve):
    actual = ccd.single_ccd("from_0_to_1 + h + e")
    expected = potential_curve
    assert actual == expected

    with pytest.raises(ValueError):
        ccd.single_ccd("No existing name")


def test_ccd_initial_and_final_ccd_from_captured_carrier(ccd, potential_curve, excited_ccd):
    actual = ccd.initial_and_final_ccd_from_captured_carrier(Carrier.h)
    expected = (potential_curve, excited_ccd)
    assert actual == expected


def test_ccd_str(ccd):
    print(ccd)


def test_spline3():
    x = [0.0, 1.0, 2.0, 3.0]
    y = [0.0, 1.0, 8.0, 28.0]
    actual = spline3(x, y, num_points=11)
    assert actual[1][1] == pytest.approx(2.17924820)

    actual = spline3(x, y, num_points=11, xrange=[-1.0, 4.0])
    # assert actual[0][0] == pytest.approx(-1.0)


def test_plot_ccd(ccd):
    plotter = CcdPlotter(ccd)
    plotter.construct_plot()
    plotter.plt.show()


def test_plot_ccd_q_range(ccd):
    plotter = CcdPlotter(ccd, q_range=[-2, 4])
    plotter.construct_plot()
    plotter.plt.show()
