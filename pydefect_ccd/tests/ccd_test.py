# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from copy import deepcopy

import numpy as np
import pytest
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from vise.tests.helpers.assertion import assert_dataclass_almost_equal

from pydefect_ccd.ccd import Ccd, SinglePoint, PotentialCurve, spline3, \
    SinglePointSpec, PotentialCurveSpec, dQ_revert, calc_omega_and_Q0
from pydefect_ccd.enum import Carrier


@pytest.fixture
def single_point():
    orb_info = LocalizedOrbital(band_idx=2,
                                ave_energy=2.0,
                                occupation=1.0,
                                orbitals={"O": [0.0, 1.0, 0.0]})

    return SinglePoint(SinglePointSpec(dQ=1.0, disp_ratio=0.1),
                       energy=10.0,
                       magnetization=1.0,
                       localized_orbitals=[[orb_info]],
                       valence_bands=None,
                       ccd_correction_energy=1.0,
                       is_shallow=True)


def test_single_point_info_corrected_energy(single_point):
    assert single_point.ccd_corrected_energy == 10.0 + 1.0


def test_image_structure_str(single_point):
    actual = single_point.__str__()
    expected = """  corr. energy  is shallow?    localized orb
        11.000  True           up-2(1.0)"""
    assert actual == expected


@pytest.fixture
def potential_curve(single_point):
    spec = PotentialCurveSpec(charge=0,
                              correction_energy=1.0,
                              counter_charge=1,
                              Q_diff=10.0)

    return PotentialCurve(spec=spec,
                          single_points=[single_point],
                          shifted_energy=3.0)


def test_potential_curve_dQs_and_energies(potential_curve):
    actual = potential_curve.dQs_and_energies()
    dQs = [1.0]
    energies = [10.0+1.0+1.0+3.0]
    expected = (dQs, energies)
    assert np.array(actual) == pytest.approx(np.array(expected))


def test_potential_curve_dQs_and_energies_w_range(potential_curve):
    actual = potential_curve.dQs_and_energies((-1.0, 0.0))
    expected = ([], [])
    assert np.array(actual) == pytest.approx(np.array(expected))


def test_dQ_reverted_single_ccd(potential_curve):
    actual = dQ_revert(potential_curve)
    expected = deepcopy(potential_curve)
    expected.single_points[0].spec = SinglePointSpec(dQ=9.0, disp_ratio=0.9)
    assert_dataclass_almost_equal(actual, expected)


@pytest.fixture
def potential_curve_final():
    spec = PotentialCurveSpec(charge=1,
                              correction_energy=1.0,
                              counter_charge=0,
                              Q_diff=10.0)
    return PotentialCurve(spec=spec,
                          single_points=[],
                          shifted_energy=0.0)


@pytest.fixture
def ccd(potential_curve, potential_curve_final):
    return Ccd(name="Va_O1_1_Va_O1_0",
               ground_curve=potential_curve_final,
               excited_curve=potential_curve)


def test_ccd_captured_carrier(ccd):
    assert ccd.captured_carrier == Carrier.h


def test_ccd_str(ccd):
    print(ccd)


def test_spline3():
    x = [0.0, 1.0, 2.0, 3.0]
    y = [0.0, 1.0, 8.0, 28.0]
    actual = spline3(x, y, num_points=11)
    assert actual[1][1] == pytest.approx(2.17924820)


def _make_energies(omega, Q0, dE, Qs):
    return 0.5 * omega * (Qs - Q0)**2 + dE


def test_calc_omega_and_Q0_variable_Q0():
    omega_true = 2.0
    Q0_true = 0.5
    dE_true = -0.1
    Qs = np.array([-1.0, 0.0, 1.0, 2.0])
    energies = _make_energies(omega_true, Q0_true, dE_true, Qs)

    omega, Q0, dE = calc_omega_and_Q0(list(Qs), list(energies), None)

    assert pytest.approx(omega_true, rel=1e-6) == omega
    assert pytest.approx(Q0_true, rel=1e-6) == Q0
    assert pytest.approx(dE_true, rel=1e-6) == dE


# def test_plot_ccd(ccd):
#     spec = PotentialCurveSpec(charge=0,
#                               correction_energy=1.0,
#                               counter_charge=1,
#                               Q_diff=10.0)

#     potential_curve = PotentialCurve(spec=spec,
#                           single_points=[single_point],
#                           shifted_energy=3.0)

#     plotter = CcdPlotter(ccd)
#     plotter.construct_plot()
#     plotter.plt.show()


# def test_plot_ccd_q_range(ccd):
#     plotter = CcdPlotter(ccd, q_range=[-2, 4])
#     plotter.construct_plot()
#     plotter.plt.show()
