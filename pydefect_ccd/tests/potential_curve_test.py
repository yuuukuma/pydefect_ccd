# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.

import pytest

from pydefect_ccd.fitting_func import QuadraticFittingFunc, QuarticFittingFunc
from pydefect_ccd.potential_curve import make_fitting_func, make_shifter, \
    PotentialCurve, PotentialCurveSpec, Shifter, SinglePoint, SinglePointSpec, \
    SinglePoints


def make_single_point(Q: float,
                      disp_ratio: float,
                      energy: float,
                      ccd_correction_energy: float = 0.0) -> SinglePoint:
    return SinglePoint(energy=energy,
                       spec=SinglePointSpec(Q=Q, disp_ratio=disp_ratio),
                       ccd_correction_energy=ccd_correction_energy,
                       used_for_fitting=True,
                       magnetization=0.0,
                       localized_orbitals=[[], []],
                       valence_bands=[[], []],
                       conduction_bands=[[], []])


def test_single_point_spec_flip():
    spec = SinglePointSpec(Q=1.25, disp_ratio=0.2)

    actual = spec.flip(Q_diff=4.0)

    assert actual == SinglePointSpec(Q=2.75, disp_ratio=0.8)
    assert spec == SinglePointSpec(Q=1.25, disp_ratio=0.2)


def test_single_points_flip_shifts_energy_and_flips_coordinate():
    single_point = make_single_point(Q=1.0,
                                     disp_ratio=0.25,
                                     energy=10.0,
                                     ccd_correction_energy=0.5)
    single_points = SinglePoints([single_point])

    actual = single_points.flip(Shifter(shift_energy=2.0, flip=True),
                                Q_diff=5.0,
                                total_energy_correction=0.3)
    shifted_point = actual.single_points[0]

    assert shifted_point is not single_point
    assert shifted_point.spec == SinglePointSpec(Q=4.0, disp_ratio=0.75)
    assert shifted_point.energy == pytest.approx(12.3)
    assert shifted_point.ccd_corrected_energy == pytest.approx(12.8)
    assert single_point.spec == SinglePointSpec(Q=1.0, disp_ratio=0.25)
    assert single_point.energy == pytest.approx(10.0)


def test_make_shifter_places_lowest_energy_at_offset():
    single_points = SinglePoints([
        make_single_point(Q=0.0, disp_ratio=0.0, energy=5.0,
                          ccd_correction_energy=0.5),
        make_single_point(Q=1.0, disp_ratio=0.5, energy=3.0,
                          ccd_correction_energy=0.25)
    ])
    spec = PotentialCurveSpec(charge=1,
                              correction_energy=1.0,
                              counter_charge=0,
                              Q_diff=5.0)

    actual = make_shifter(spec, single_points, offset=2.0, flip=True)

    assert actual == Shifter(shift_energy=pytest.approx(-2.25), flip=True)


def test_potential_curve_applies_shifter_to_single_points():
    single_point = make_single_point(Q=1.0,
                                     disp_ratio=0.25,
                                     energy=4.0,
                                     ccd_correction_energy=0.2)
    spec = PotentialCurveSpec(charge=1,
                              correction_energy=1.5,
                              counter_charge=0,
                              Q_diff=5.0)
    curve = PotentialCurve(spec=spec,
                           original_single_points=SinglePoints([single_point]),
                           shifter=Shifter(shift_energy=-2.0, flip=True))

    shifted_point = curve.single_points.single_points[0]

    assert shifted_point.spec == SinglePointSpec(Q=4.0, disp_ratio=0.75)
    assert shifted_point.energy == pytest.approx(3.5)
    assert shifted_point.ccd_corrected_energy == pytest.approx(3.7)
    assert curve.Qs_and_energies == pytest.approx([(4.0, 3.7)])
    assert single_point.spec == SinglePointSpec(Q=1.0, disp_ratio=0.25)
    assert single_point.energy == pytest.approx(4.0)


def test_potential_curve_lowest_energy_adds_spec_correction_to_shifted_minimum():
    single_points = SinglePoints([
        make_single_point(Q=0.0, disp_ratio=0.0, energy=4.0,
                          ccd_correction_energy=0.2),
        make_single_point(Q=1.0, disp_ratio=0.5, energy=6.0,
                          ccd_correction_energy=0.3)
    ])
    spec = PotentialCurveSpec(charge=1,
                              correction_energy=1.5,
                              counter_charge=0,
                              Q_diff=5.0)
    curve = PotentialCurve(spec=spec,
                           original_single_points=single_points,
                           shifter=make_shifter(spec, single_points,
                                                offset=2.0))

    assert curve.single_points.lowest_energy == pytest.approx(2.0)
    assert curve.lowest_energy == pytest.approx(3.5)


def test_single_point_from_disp_uses_close_match():
    first = make_single_point(Q=0.0, disp_ratio=0.0, energy=1.0)
    second = make_single_point(Q=1.0, disp_ratio=0.3, energy=2.0)
    single_points = SinglePoints([first, second])

    assert single_points.single_point_from_disp(0.3 + 1e-10) is second


def test_single_point_from_disp_raises_when_missing():
    single_points = SinglePoints([
        make_single_point(Q=0.0, disp_ratio=0.0, energy=1.0)
    ])

    with pytest.raises(ValueError, match="disp_ratio=0.5"):
        single_points.single_point_from_disp(0.5)


def test_potential_curve_spec_effective_charge():
    spec = PotentialCurveSpec(charge=1,
                              correction_energy=0.0,
                              counter_charge=-1,
                              Q_diff=5.0)

    assert spec.effective_charge == 2


def test_from_single_points():
    class FakeSinglePoint:
        def __init__(self, Q, energy):
            self.Q = Q
            self.ccd_corrected_energy = energy

    sp = SinglePoints([FakeSinglePoint(Q=-1.0, energy=3.0 * 1.0 + 0.25),
                       FakeSinglePoint(Q=0.0, energy=0.25),
                       FakeSinglePoint(Q=1.0, energy=3.0 * 1.0 + 0.25)])
    q = make_fitting_func(QuadraticFittingFunc, sp)
    assert q.a == pytest.approx(3.0)
    assert q.Q0 == pytest.approx(0.0)
    assert q.E0 == pytest.approx(0.25)


class FakeSinglePoint:
    def __init__(self, Q, energy):
        self.Q = Q
        self.ccd_corrected_energy = energy


def test_quadratic_from_single_points():
    sp = SinglePoints([FakeSinglePoint(Q=-1.0, energy=3.0 * 1.0 + 0.25),
                       FakeSinglePoint(Q=0.0, energy=0.25),
                       FakeSinglePoint(Q=1.0, energy=3.0 * 1.0 + 0.25)])
    q = make_fitting_func(QuadraticFittingFunc, sp)
    assert q.a == pytest.approx(3.0)
    assert q.Q0 == pytest.approx(0.0)
    assert q.E0 == pytest.approx(0.25)


def test_quartic_from_single_points():
    sp = SinglePoints([FakeSinglePoint(Q=-2.0, energy=16.0 - 8.0 + 4.0 + 0.25),
                       FakeSinglePoint(Q=-1.0, energy=1.0 - 1.0 + 1.0 + 0.25),
                       FakeSinglePoint(Q=0.0, energy=0.25),
                       FakeSinglePoint(Q=1.0, energy=1.0 + 1.0 + 1.0 + 0.25),
                       FakeSinglePoint(Q=2.0, energy=16.0 + 8.0 + 4.0  + 0.25),
                       ])
    q = make_fitting_func(QuarticFittingFunc, sp)
    assert q.a == pytest.approx(1.0)
    assert q.b == pytest.approx(1.0)
    assert q.c == pytest.approx(1.0)
    assert q.Q0 == pytest.approx(0.0)
    assert q.E0 == pytest.approx(0.25)
