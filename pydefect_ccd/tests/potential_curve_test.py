# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.

import matplotlib
import pytest

matplotlib.use("Agg")
from matplotlib import pyplot as plt

from pydefect_ccd.fitting_func import QuadraticFittingFunc, QuarticFittingFunc
from pydefect_ccd.potential_curve import CurveTransform, make_curve_transform, \
    make_fitting_func, PotentialCurve, PotentialCurveSpec, SinglePoint, \
    SinglePointSpec, SinglePoints


def make_single_point(Q: float,
                      disp_ratio: float,
                      energy: float,
                      ccd_correction_energy: float = 0.0,
                      used_for_fitting: bool = True,
                      is_shallow: bool = None) -> SinglePoint:
    return SinglePoint(energy=energy,
                       spec=SinglePointSpec(Q=Q, disp_ratio=disp_ratio),
                       ccd_correction_energy=ccd_correction_energy,
                       used_for_fitting=used_for_fitting,
                       is_shallow=is_shallow,
                       magnetization=0.0,
                       localized_orbitals=[[], []],
                       valence_bands=[[], []],
                       conduction_bands=[[], []])


def test_single_point_spec_flip():
    spec = SinglePointSpec(Q=0.8, disp_ratio=0.2)
    actual = spec.flip(Q_diff=4.0)
    assert actual == SinglePointSpec(Q=3.2, disp_ratio=0.8)
    assert spec == SinglePointSpec(Q=0.8, disp_ratio=0.2)


def test_single_point_spec_flip_raises_for_inconsistent_q_diff():
    spec = SinglePointSpec(Q=1.0, disp_ratio=0.25)

    with pytest.raises(ValueError, match="Q_diff=5.0"):
        spec.flip(Q_diff=5.0)


def test_single_points_transform_shifts_energy_and_flips_coordinate():
    single_point = make_single_point(Q=1.0,
                                     disp_ratio=0.2,
                                     energy=10.0,
                                     ccd_correction_energy=0.5)
    single_points = SinglePoints([single_point])

    actual = single_points.transform(CurveTransform(shift_energy=2.0, flip=True),
                                     Q_diff=5.0,
                                     total_energy_correction=0.3)
    shifted_point = actual.single_points[0]

    assert shifted_point is not single_point
    assert shifted_point.spec == SinglePointSpec(Q=4.0, disp_ratio=0.8)
    assert shifted_point.energy == pytest.approx(12.3)
    assert shifted_point.ccd_corrected_energy == pytest.approx(12.8)
    assert single_point.spec == SinglePointSpec(Q=1.0, disp_ratio=0.2)
    assert single_point.energy == pytest.approx(10.0)


def test_make_curve_transform_places_lowest_energy_at_offset():
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

    actual = make_curve_transform(spec, single_points, offset=2.0, flip=True)

    assert actual == CurveTransform(shift_energy=pytest.approx(-2.25), flip=True)


def test_potential_curve_applies_curve_transform_to_single_points():
    single_point = make_single_point(Q=1.0,
                                     disp_ratio=0.2,
                                     energy=4.0,
                                     ccd_correction_energy=0.2)
    spec = PotentialCurveSpec(charge=1,
                              correction_energy=1.5,
                              counter_charge=0,
                              Q_diff=5.0)
    curve = PotentialCurve(spec=spec,
                           original_single_points=SinglePoints([single_point]),
                           curve_transform=CurveTransform(shift_energy=-2.0,
                                                          flip=True))

    shifted_point = curve.single_points.single_points[0]

    assert shifted_point.spec == SinglePointSpec(Q=4.0, disp_ratio=0.8)
    assert shifted_point.energy == pytest.approx(3.5)
    assert shifted_point.ccd_corrected_energy == pytest.approx(3.7)
    assert curve.single_points.Qs == pytest.approx([4.0])
    assert curve.single_points.corrected_energies == pytest.approx([3.7])
    assert single_point.spec == SinglePointSpec(Q=1.0, disp_ratio=0.2)
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
                           curve_transform=make_curve_transform(spec,
                                                                single_points,
                                                                offset=2.0))

    assert curve.single_points.lowest_energy == pytest.approx(2.0)
    assert curve.lowest_energy == pytest.approx(2.0)


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


def test_make_fitting_func_uses_only_points_marked_for_fitting():
    sp = SinglePoints([
        make_single_point(Q=-1.0, disp_ratio=-1.0, energy=3.25),
        make_single_point(Q=0.0, disp_ratio=0.0, energy=0.25),
        make_single_point(Q=1.0, disp_ratio=1.0, energy=3.25),
        make_single_point(Q=2.0, disp_ratio=2.0, energy=100.0,
                          used_for_fitting=False)
    ])

    q = make_fitting_func(QuadraticFittingFunc, sp)

    assert q.a == pytest.approx(3.0)
    assert q.E0 == pytest.approx(0.25)


def test_potential_curve_set_fitting_curve_uses_fitting_points():
    single_points = SinglePoints([
        make_single_point(Q=-1.0, disp_ratio=-1.0, energy=3.25),
        make_single_point(Q=0.0, disp_ratio=0.0, energy=0.25),
        make_single_point(Q=1.0, disp_ratio=1.0, energy=3.25),
        make_single_point(Q=2.0, disp_ratio=2.0, energy=100.0,
                          used_for_fitting=False)
    ])
    spec = PotentialCurveSpec(charge=0,
                              correction_energy=0.0,
                              counter_charge=1,
                              Q_diff=1.0)
    curve = PotentialCurve(spec=spec,
                           original_single_points=single_points,
                           curve_transform=CurveTransform(0.0, False))

    curve.set_fitting_curve(QuadraticFittingFunc)

    assert curve.fitting_func.a == pytest.approx(3.0)
    assert curve.fitting_func.E0 == pytest.approx(0.25)


def test_potential_curve_set_fitting_curve_uses_transformed_points():
    single_points = SinglePoints([
        make_single_point(Q=-1.0, disp_ratio=-1.0, energy=3.25),
        make_single_point(Q=0.0, disp_ratio=0.0, energy=0.25),
        make_single_point(Q=1.0, disp_ratio=1.0, energy=3.25)
    ])
    spec = PotentialCurveSpec(charge=0,
                              correction_energy=0.5,
                              counter_charge=1,
                              Q_diff=1.0)
    curve = PotentialCurve(spec=spec,
                           original_single_points=single_points,
                           curve_transform=CurveTransform(2.0, False))

    curve.set_fitting_curve(QuadraticFittingFunc)

    assert curve.fitting_func.a == pytest.approx(3.0)
    assert curve.fitting_func.E0 == pytest.approx(2.75)


def test_potential_curve_add_plot(mocker):
    single_points = SinglePoints([
        make_single_point(Q=0.0, disp_ratio=0.0, energy=0.25,
                          is_shallow=False),
        make_single_point(Q=1.0, disp_ratio=1.0, energy=3.25,
                          is_shallow=True)
    ])
    spec = PotentialCurveSpec(charge=0,
                              correction_energy=0.0,
                              counter_charge=1,
                              Q_diff=1.0)
    curve = PotentialCurve(spec=spec,
                           original_single_points=single_points,
                           curve_transform=CurveTransform(0.0, False),
                           fitting_func=QuadraticFittingFunc(Q0=0.0,
                                                             E0=0.25,
                                                             a=3.0))
    ax = mocker.Mock()

    curve.add_plot(ax, "red")

    ax.plot.assert_called_once()
    assert ax.scatter.call_args_list == [
        mocker.call([0.0], [0.25],
                    facecolors="red",
                    edgecolors="red",
                    label=None),
        mocker.call([1.0], [3.25],
                    facecolors="none",
                    edgecolors="red",
                    label=None)
    ]


def test_potential_curve_add_plot_without_fitting_curve_plots_points_only(mocker):
    single_points = SinglePoints([
        make_single_point(Q=0.0, disp_ratio=0.0, energy=0.25,
                          is_shallow=False),
        make_single_point(Q=1.0, disp_ratio=1.0, energy=3.25,
                          is_shallow=True)
    ])
    spec = PotentialCurveSpec(charge=0,
                              correction_energy=0.0,
                              counter_charge=1,
                              Q_diff=1.0)
    curve = PotentialCurve(spec=spec,
                           original_single_points=single_points,
                           curve_transform=CurveTransform(0.0, False))
    ax = mocker.Mock()

    curve.add_plot(ax, "red")

    ax.plot.assert_not_called()
    assert ax.scatter.call_args_list == [
        mocker.call([0.0], [0.25],
                    facecolors="red",
                    edgecolors="red",
                    label="q=0"),
        mocker.call([1.0], [3.25],
                    facecolors="none",
                    edgecolors="red",
                    label="q=0 shallow")
    ]


def test_potential_curve_add_plot_renders_png(tmp_path):
    single_points = SinglePoints([
        make_single_point(Q=0.0, disp_ratio=0.0, energy=0.25,
                          is_shallow=False),
        make_single_point(Q=1.0, disp_ratio=1.0, energy=3.25,
                          is_shallow=True)
    ])
    spec = PotentialCurveSpec(charge=0,
                              correction_energy=0.0,
                              counter_charge=1,
                              Q_diff=1.0)
    curve = PotentialCurve(spec=spec,
                           original_single_points=single_points,
                           curve_transform=CurveTransform(0.0, False),
                           fitting_func=QuadraticFittingFunc(Q0=0.0,
                                                             E0=0.25,
                                                             a=3.0))
    fig, ax = plt.subplots()

    try:
        curve.add_plot(ax, "red")
        fig_path = tmp_path / "potential_curve.png"
        fig.savefig(fig_path)
        print(tmp_path)

        assert fig_path.exists()
        assert fig_path.stat().st_size > 0
        assert len(ax.lines) == 1
        assert len(ax.collections) == 2
    finally:
        plt.close(fig)


def test_potential_curve_add_plot_without_fitting_curve_renders_points_only_png(tmp_path):
    single_points = SinglePoints([
        make_single_point(Q=0.0, disp_ratio=0.0, energy=0.25,
                          is_shallow=False),
        make_single_point(Q=1.0, disp_ratio=1.0, energy=3.25,
                          is_shallow=True)
    ])
    spec = PotentialCurveSpec(charge=0,
                              correction_energy=0.0,
                              counter_charge=1,
                              Q_diff=1.0)
    curve = PotentialCurve(spec=spec,
                           original_single_points=single_points,
                           curve_transform=CurveTransform(0.0, False))
    fig, ax = plt.subplots()

    try:
        curve.add_plot(ax, "red")
        fig_path = tmp_path / "potential_curve_points_only.png"
        fig.savefig(fig_path)
        print(tmp_path)

        assert fig_path.exists()
        assert fig_path.stat().st_size > 0
        assert len(ax.lines) == 0
        assert len(ax.collections) == 2
    finally:
        plt.close(fig)


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
