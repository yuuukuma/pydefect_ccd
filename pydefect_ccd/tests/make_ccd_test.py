# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.

import pytest

from pydefect_ccd.fitting_func import QuadraticFittingFunc
from pydefect_ccd.potential_curve import SinglePointSpec, SinglePoint, \
    CurveTransform, PotentialCurveSpec, SinglePoints
from pydefect_ccd.make_ccd import MakeCcd

band_edges = dict(vbm=1.0, cbm=3.0, supercell_vbm=1.1, supercell_cbm=2.9)

common = dict(is_shallow=False, used_for_fitting=True)


def test_make_ccd(excited_structure, ground_structure):
    ground_spec = PotentialCurveSpec(charge=0, correction_energy=0.0,
                                     counter_charge=1, Q_diff=5.0)
    excited_spec = PotentialCurveSpec(charge=1, correction_energy=1.0,
                                      counter_charge=0, Q_diff=5.0)

    ground_single_points = SinglePoints([
        SinglePoint(spec=SinglePointSpec(Q=0.0, disp_ratio=0.0),
                    energy=0.0,
                    magnetization=0.0,
                    ccd_correction_energy=0.0,
                    used_for_fitting=True),
        SinglePoint(spec=SinglePointSpec(Q=2.5, disp_ratio=0.5),
                    energy=10.0,
                    magnetization=0.0,
                    ccd_correction_energy=10.0,
                    used_for_fitting=True),
        SinglePoint(spec=SinglePointSpec(Q=5.0, disp_ratio=1.0),
                    energy=40.0,
                    magnetization=0.0,
                    ccd_correction_energy=40.0,
                    used_for_fitting=True)
    ])
    excited_single_points = SinglePoints([
        SinglePoint(spec=SinglePointSpec(Q=0.0, disp_ratio=0.0),
                    energy=3.0,
                    magnetization=1.0,
                    ccd_correction_energy=0.0,
                    used_for_fitting=True),
        SinglePoint(spec=SinglePointSpec(Q=2.5, disp_ratio=0.5),
                    energy=13.0,
                    magnetization=1.0,
                    ccd_correction_energy=10.0,
                    used_for_fitting=True),
        SinglePoint(spec=SinglePointSpec(Q=5.0, disp_ratio=1.0),
                    energy=43.0,
                    magnetization=1.0,
                    ccd_correction_energy=40.0,
                    used_for_fitting=True)
    ])

    actual = MakeCcd(ground_single_points=ground_single_points,
                     ground_pot_curve_spec=ground_spec,
                     ground_fitting_func=QuadraticFittingFunc,
                     excited_single_points=excited_single_points,
                     excited_pot_curve_spec=excited_spec,
                     excited_fitting_func=QuadraticFittingFunc,
                     vbm=100.0,
                     cbm=200.0,
                     name="test").ccd

    assert actual.name == "test"
    assert actual.ground_curve.curve_transform == CurveTransform(0.0, False)
    assert actual.excited_curve.curve_transform == CurveTransform(196.0, True)
    assert actual.excited_curve.single_points.single_points[0].spec \
        == SinglePointSpec(Q=5.0, disp_ratio=1.0)
    assert actual.excited_curve.single_points.lowest_energy == pytest.approx(200.0)
    assert actual.excited_curve.lowest_energy == pytest.approx(201.0)
    assert actual.ground_curve.fitting_curve
    assert actual.excited_curve.fitting_curve
