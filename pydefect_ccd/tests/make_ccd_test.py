# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.

from pydefect_ccd.ccd import Ccd
from pydefect_ccd.potential_curve import SinglePointSpec, SinglePoint, \
    PotentialCurveSpec, PotentialCurve
from pydefect_ccd.make_ccd import MakeCcd

band_edges = dict(vbm=1.0, cbm=3.0, supercell_vbm=1.1, supercell_cbm=2.9)

common = dict(is_shallow=False, used_for_fitting=True)


def test_make_ccd(excited_structure, ground_structure):
    ground_spec = PotentialCurveSpec(charge=0, correction_energy=0.0,
                                     counter_charge=1, Q_diff=5.0)
    excited_spec = PotentialCurveSpec(charge=1, correction_energy=1.0,
                                      counter_charge=0, Q_diff=5.0)

    ground_single_points = [
        SinglePoint(SinglePointSpec(Q=0.0, disp_ratio=0.0),
                    energy=0.0,
                    magnetization=0.0,
                    ccd_correction_energy=0.0),
        SinglePoint(SinglePointSpec(Q=2.5, disp_ratio=0.5),
                    energy=10.0,
                    magnetization=0.0,
                    ccd_correction_energy=10.0)
    ]
    excited_single_points = [
        SinglePoint(SinglePointSpec(Q=0.0, disp_ratio=0.0),
                    energy=3.0,
                    magnetization=1.0,
                    ccd_correction_energy=0.0),
        SinglePoint(SinglePointSpec(Q=2.5, disp_ratio=0.5),
                    energy=13.0,
                    magnetization=1.0,
                    ccd_correction_energy=10.0)
    ]

    ground = PotentialCurve(ground_spec, ground_single_points)
    excited = PotentialCurve(excited_spec, excited_single_points)

    actual = MakeCcd(ground, excited, vbm=100.0, cbm=200.0, name="test").ccd

    # ----------------
    expected_ground_curve = PotentialCurve(
        ground_spec,
        [
            SinglePoint(SinglePointSpec(Q=0.0, disp_ratio=0.0),
                        energy=0.0,
                        magnetization=0.0,
                        ccd_correction_energy=0.0),
            SinglePoint(SinglePointSpec(Q=2.5, disp_ratio=0.5),
                        energy=10.0,
                        magnetization=0.0,
                        ccd_correction_energy=10.0)
        ]
    )

    expected_excited_curve = PotentialCurve(
        excited_spec,
        [
            SinglePoint(SinglePointSpec(Q=5.0, disp_ratio=1.0),
                        energy=3.0,
                        magnetization=1.0,
                        ccd_correction_energy=0.0),
            SinglePoint(SinglePointSpec(Q=2.5, disp_ratio=0.5),
                        energy=13.0,
                        magnetization=1.0,
                        ccd_correction_energy=10.0)
        ],
        shifted_energy=200.0
    )

    expected = Ccd(name="test",
                   ground_curve=expected_ground_curve,
                   excited_curve=expected_excited_curve)

    assert actual == expected

