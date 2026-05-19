# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
from argparse import Namespace
from pathlib import Path

import matplotlib
import pytest
from monty.serialization import loadfn
from pydefect.analyzer.unitcell import Unitcell
from pymatgen.core import Structure
from pymatgen.electronic_structure.core import Spin
from vise.analyzer.effective_mass import EffectiveMass
from vise.input_set.incar import ViseIncar
from vise.input_set.prior_info import PriorInfo

matplotlib.use("Agg")

from pydefect_ccd.fitting_func import FittingFuncType, QuadraticFittingFunc
from pydefect_ccd.potential_curve import SinglePointSpec, PotentialCurveSpec
from pydefect_ccd.ccd_init import CcdInit
from pydefect_ccd.cli.main_function import make_ccd_init, make_ccd, plot_ccd, \
    make_ccd_dirs, make_wswq_dirs, plot_eigenvalues, \
    main_make_e_p_matrix_element, make_sommerfeld_scaling, \
    make_ccd_corrections, make_potential_curve, make_single_points, _fit_curve, \
    _plot_sommerfeld_scaling
from pydefect_ccd.local_enum import Carrier
from pydefect_ccd.relaxed_point import NearEdgeState, RelaxedPoint


def test_make_sommerfeld_scaling(test_files, tmpdir):
    tmpdir.chdir()
    args = Namespace(epsilon0=10.0,
                     unitcell=None,
                     electron_and_hole_effective_mass=[2.0, 5.0],
                     effective_mass_file=None,
                     temperatures=[100, 200])
    make_sommerfeld_scaling(args)
    assert str(loadfn("sommerfeld_scaling.json"))
    assert Path("sommerfeld_scaling_attractive.pdf").exists()
    assert Path("sommerfeld_scaling_repulsive.pdf").exists()


def test_make_sommerfeld_scaling_from_effective_mass_file(tmpdir):
    tmpdir.chdir()
    effective_mass = EffectiveMass(
        p=[[[3.0, 0.0, 0.0],
            [0.0, 6.0, 0.0],
            [0.0, 0.0, 9.0]]],
        n=[[[2.0, 0.0, 0.0],
            [0.0, 4.0, 0.0],
            [0.0, 0.0, 6.0]]],
        temperature=300.0,
        concentrations=[1e18])
    args = Namespace(epsilon0=10.0,
                     unitcell=None,
                     electron_and_hole_effective_mass=None,
                     effective_mass_file=effective_mass,
                     temperatures=[100, 200])

    make_sommerfeld_scaling(args)

    actual = loadfn("sommerfeld_scaling.json")
    assert actual.electron_effective_mass == pytest.approx(4.0)
    assert actual.hole_effective_mass == pytest.approx(6.0)


def test_make_sommerfeld_scaling_plots_attractive_and_repulsive(tmpdir, mocker):
    tmpdir.chdir()
    mock_plot = mocker.patch(
        "pydefect_ccd.cli.main_function._plot_sommerfeld_scaling")
    args = Namespace(epsilon0=10.0,
                     unitcell=None,
                     electron_and_hole_effective_mass=[2.0, 5.0],
                     effective_mass_file=None,
                     temperatures=[100, 200])

    make_sommerfeld_scaling(args)

    scaling = mock_plot.call_args_list[0].args[0]
    assert mock_plot.call_args_list == [
        mocker.call(scaling,
                    title="attractive",
                    carrier_charges=[(Carrier.e, 1),
                                     (Carrier.h, -1)],
                    filename="sommerfeld_scaling_attractive.pdf"),
        mocker.call(scaling,
                    title="repulsive",
                    carrier_charges=[(Carrier.e, -1),
                                     (Carrier.h, 1)],
                    filename="sommerfeld_scaling_repulsive.pdf")]


def test_plot_sommerfeld_scaling_adds_electron_and_hole_curves(mocker):
    scaling = mocker.Mock()
    fig = mocker.Mock()
    ax = mocker.Mock()
    ax.lines = []

    def add_line(**_kwargs):
        ax.lines.append(mocker.Mock())

    scaling.add_to_ax.side_effect = add_line
    mocker.patch("pydefect_ccd.cli.main_function.plt.subplots",
                 return_value=(fig, ax))
    mock_close = mocker.patch("pydefect_ccd.cli.main_function.plt.close")

    _plot_sommerfeld_scaling(scaling,
                             title="attractive",
                             carrier_charges=[(Carrier.e, 1),
                                              (Carrier.h, -1)],
                             filename="sommerfeld_scaling_attractive.pdf")

    assert scaling.add_to_ax.call_args_list == [
        mocker.call(ax=ax, carrier_type=Carrier.e, defect_charge=1),
        mocker.call(ax=ax, carrier_type=Carrier.h, defect_charge=-1)]
    ax.lines[0].set_label.assert_called_once_with(str(Carrier.e))
    ax.lines[1].set_label.assert_called_once_with(str(Carrier.h))
    ax.set_title.assert_called_once_with("attractive")
    scaling.set_label.assert_called_once_with(ax)
    fig.tight_layout.assert_called_once_with()
    fig.savefig.assert_called_once_with("sommerfeld_scaling_attractive.pdf")
    mock_close.assert_called_once_with(fig)


def test_make_ccd_init(test_files, tmpdir):
    print(tmpdir)
    tmpdir.chdir()
    dir_ = test_files / "GaN_C_N"

    args = Namespace(first_dir=dir_ / "relaxed" / "q_0",
                     second_dir=dir_ / "relaxed" / "q_-1",
                     unitcell=Unitcell.from_yaml(dir_ / "unitcell.yaml"),
                     p_state=loadfn(dir_ / "perfect_band_edge_state.json"))
    make_ccd_init(args)
    print(loadfn("C_N1_0_-1/ccd_init.json"))


def test_make_ccd_dirs(tmpdir, ground_structure, excited_structure,
                       intermediate_structure):
    print(tmpdir)
    tmpdir.chdir()
    ccd_init = CcdInit(
        first_relaxed_point=RelaxedPoint(name="test",
                                         charge=1,
                                         structure=ground_structure,
                                         formation_energy=10.0,
                                         correction_energy=200.0,
                                         magnetization=0.0,
                                         localized_orbitals=[[]],
                                         initial_site_symmetry="",
                                         final_site_symmetry="",
                                         parsed_dir="",
                                         valence_bands=[
                                             [NearEdgeState(band_index=10,
                                                            kpt_coord=[0.0] * 3,
                                                            kpt_index=1,
                                                            eigenvalue=1.0,
                                                            occupation=1.0)]],
                                         conduction_bands=[
                                             [NearEdgeState(band_index=11,
                                                            kpt_coord=[0.0] * 3,
                                                            kpt_index=1,
                                                            eigenvalue=2.0,
                                                            occupation=0.0)]]),
        second_relaxed_point=RelaxedPoint(name="test",
                                          charge=0,
                                          structure=excited_structure,
                                          formation_energy=10.0,
                                          correction_energy=100.0,
                                          magnetization=1.0,
                                          localized_orbitals=[[]],
                                          initial_site_symmetry="",
                                          final_site_symmetry="",
                                          parsed_dir="",
                                          valence_bands=[
                                              [NearEdgeState(band_index=10,
                                                             kpt_coord=[0.0] * 3,
                                                             kpt_index=1,
                                                             eigenvalue=1.0,
                                                             occupation=1.0)]],
                                          conduction_bands=[
                                              [NearEdgeState(band_index=11,
                                                             kpt_coord=[0.0] * 3,
                                                             kpt_index=1,
                                                             eigenvalue=2.0,
                                                             occupation=0.0)]],
                                          ),
        vbm=-100.0, cbm=100.0, supercell_vbm=-100.0, supercell_cbm=100.0)

    Path("test").mkdir()
    args = Namespace(ccd_init=ccd_init,
                     first_to_second_div_ratios=[0.5, 1.0],
                     second_to_first_div_ratios=[0.0, 1.0],
                     calc_dir=Path("test"))
    make_ccd_dirs(args)

    actual = Structure.from_file("q_1/disp_0.5/POSCAR")
    assert actual == intermediate_structure

    actual = PriorInfo.load_yaml("q_1/disp_0.5/prior_info.yaml")
    expected = PriorInfo(charge=1)
    assert actual == expected

    # dQ = sqrt((0.1*10)**2*6 * Element.H.atomic_mass)
    # dQ / 2 =1.2295974951178128
    actual = loadfn("q_1/disp_0.5/single_point_spec.json")
    dQ = 2.4591949902356256
    expected = SinglePointSpec(Q=dQ / 2, disp_ratio=0.5)
    assert actual == expected

    actual = loadfn("q_0/disp_0.0/single_point_spec.json")
    expected = SinglePointSpec(Q=0.0, disp_ratio=0.0)
    assert actual == expected

    actual = Structure.from_file("q_1/disp_1.0/POSCAR")
    assert actual == excited_structure

    actual = loadfn("q_0/potential_curve_spec.json")
    expected = PotentialCurveSpec(charge=0, correction_energy=100.0,
                                  counter_charge=1, Q_diff=dQ)
    assert actual == expected


def test_make_ccd_corrections(tmpdir, mocker):
    tmpdir.chdir()
    dir_ = Path("disp_0.5")
    dir_.mkdir()
    SinglePointSpec(Q=1.0, disp_ratio=0.5).to_json_file(
        dir_ / "single_point_spec.json")

    calc_results = mocker.Mock()
    no_disp_calc_results = mocker.Mock()
    ccd_correction = mocker.Mock()
    plotter = mocker.Mock()
    plotter.plt = mocker.Mock()

    mock_get_calc_results = mocker.patch(
        "pydefect_ccd.cli.main_function.get_calc_results",
        return_value=calc_results)
    mock_make_efnv_correction = mocker.patch(
        "pydefect_ccd.cli.main_function.make_efnv_correction",
        return_value=ccd_correction)
    mock_plotter_cls = mocker.patch(
        "pydefect_ccd.cli.main_function.SitePotentialMplPlotter")
    mock_plotter_cls.from_efnv_corr.return_value = plotter

    args = Namespace(
        dirs=[dir_],
        potential_curve_spec=PotentialCurveSpec(charge=1,
                                                correction_energy=0.0,
                                                counter_charge=-2,
                                                Q_diff=1.0),
        no_disp_calc_results=no_disp_calc_results,
        unitcell=Namespace(effective_ionic_diele_const=12.3),
        no_disp_defect_entry=Namespace(defect_center=[0.1, 0.2, 0.3]))

    make_ccd_corrections(args)

    mock_get_calc_results.assert_called_once_with(dir_, False)
    mock_make_efnv_correction.assert_called_once_with(
        1.5,
        calc_results,
        no_disp_calc_results,
        12.3,
        [0.1, 0.2, 0.3])
    ccd_correction.to_json_file.assert_called_once_with(
        dir_ / "ccd_correction.json")
    mock_plotter_cls.from_efnv_corr.assert_called_once_with(
        title="ccd_correction", efnv_correction=ccd_correction)
    plotter.construct_plot.assert_called_once_with()
    plotter.plt.savefig.assert_called_once_with(
        fname=dir_ / "ccd_correction.pdf")
    plotter.plt.clf.assert_called_once_with()


def test_make_single_points(tmpdir, mocker):
    tmpdir.chdir()
    dir_ = Path("disp_0.5")
    dir_.mkdir()

    calc_results = Namespace(energy=-10.5, magnetization=2.0)
    ccd_correction = Namespace(correction_energy=0.7)
    band_edge_states = Namespace(is_shallow=True)
    band_edge_orbital_infos = mocker.Mock()
    spec = SinglePointSpec(Q=1.0, disp_ratio=0.5)
    localized_orbitals = [["lo"]]
    valence_bands = [["vb"]]
    conduction_bands = [["cb"]]
    single_point = mocker.Mock()

    mock_get_calc_results = mocker.patch(
        "pydefect_ccd.cli.main_function.get_calc_results",
        return_value=calc_results)
    mock_loadfn = mocker.patch(
        "pydefect_ccd.cli.main_function.loadfn",
        side_effect=[
            ccd_correction,
            band_edge_states,
            band_edge_orbital_infos,
            spec])
    mock_get_band_edge_info = mocker.patch(
        "pydefect_ccd.cli.main_function._get_band_edge_info",
        return_value=(localized_orbitals, valence_bands, conduction_bands))
    mock_single_point = mocker.patch(
        "pydefect_ccd.cli.main_function.SinglePoint",
        return_value=single_point)

    args = Namespace(dirs=[dir_], parse_ccd_correction=True)

    make_single_points(args)

    mock_get_calc_results.assert_called_once_with(dir_, False)
    assert mock_loadfn.call_args_list == [
        mocker.call(dir_ / "ccd_correction.json"),
        mocker.call(dir_ / "band_edge_states.json"),
        mocker.call(dir_ / "band_edge_orbital_infos.json"),
        mocker.call(dir_ / "single_point_spec.json")]
    mock_get_band_edge_info.assert_called_once_with(
        band_edge_orbital_infos, band_edge_states)
    mock_single_point.assert_called_once_with(
        spec=spec,
        energy=-10.5,
        ccd_correction_energy=0.7,
        magnetization=2.0,
        localized_orbitals=localized_orbitals,
        valence_bands=valence_bands,
        conduction_bands=conduction_bands,
        is_shallow=True)
    single_point.to_json_file.assert_called_once_with(
        dir_ / "single_point.json")


def test_fit_curve_accepts_fitting_curve_type(mocker):
    potential_curve = mocker.Mock()

    actual = _fit_curve(potential_curve, FittingFuncType.quadratic)

    potential_curve.set_fitting_curve.assert_called_once_with(QuadraticFittingFunc)
    assert actual is potential_curve.fitting_func


def test_make_potential_curve_sets_plot_labels(mocker):
    potential_curve = mocker.Mock()
    ax = mocker.Mock()
    mocker.patch("pydefect_ccd.cli.main_function.parse_dirs",
                 return_value=[mocker.Mock()])
    mocker.patch("pydefect_ccd.cli.main_function.make_curve_transform",
                 return_value=mocker.Mock())
    mocker.patch("pydefect_ccd.cli.main_function.PotentialCurve",
                 return_value=potential_curve)
    mocker.patch("pydefect_ccd.cli.main_function.plt.gca", return_value=ax)
    mock_savefig = mocker.patch("pydefect_ccd.cli.main_function.plt.savefig")
    mock_tight_layout = mocker.patch(
        "pydefect_ccd.cli.main_function.plt.tight_layout")
    args = Namespace(dirs=[Path("disp_0.0")],
                     potential_curve_spec=mocker.Mock(),
                     fitting_func=None)

    make_potential_curve(args)

    potential_curve.add_plot.assert_called_once_with(ax)
    ax.set_xlabel.assert_called_once_with("Q (amu$^{1/2}$ Å)")
    ax.set_ylabel.assert_called_once_with("dE (eV)")
    ax.legend.assert_called_once_with()
    mock_tight_layout.assert_called_once_with()
    mock_savefig.assert_called_once_with("potential_curve.pdf")


# @pytest.fixture
# def ccd():
#     return Ccd(name="test",
#                image_infos=[SinglePointInfo(-0.2, -1.2, 12.0),
#                                     SinglePointInfo(0.0, 9.0, 10.0),
#                                     SinglePointInfo(0.2, -0.8, 8.0)],
#                ground_image_infos=[SinglePointInfo(-0.2, -2.2, -2.0),
#                                    SinglePointInfo(0.0, 18.0, 0.0),
#                                    SinglePointInfo(0.2, -1.8, 2.0)])
#
#


def test_make_ccd(mocker):
    ccd = mocker.Mock()
    mock_make_ccd = mocker.patch("pydefect_ccd.cli.main_function.MakeCcd")
    mock_make_ccd.return_value.ccd = ccd
    ground_curve = mocker.Mock()
    excited_curve = mocker.Mock()
    ccd_init = Namespace(vbm=1.0, cbm=2.0, name="test")
    args = Namespace(ground_pot_curve=ground_curve,
                     excited_pot_curve=excited_curve,
                     ccd_init=ccd_init)

    make_ccd(args)

    mock_make_ccd.assert_called_once_with(ground_pot_curve=ground_curve,
                                          excited_pot_curve=excited_curve,
                                          vbm=1.0,
                                          cbm=2.0,
                                          name="test")
    ccd.to_json_file.assert_called_once_with()


def test_plot_ccd(tmpdir, mocker):
    plotter = mocker.Mock()
    mock_plotter = mocker.patch("pydefect_ccd.cli.main_function.CcdPlotter",
                                return_value=plotter)
    mock_savefig = mocker.patch("pydefect_ccd.cli.main_function.plt.savefig")
    mock_show = mocker.patch("pydefect_ccd.cli.main_function.plt.show")
    ccd = mocker.Mock()
    args = Namespace(ccd=ccd,
                     fig_name=tmpdir / "ccd.pdf",
                     ground_q_range=[-1.0, 1.0],
                     excited_q_range=[0.0, 2.0])

    plot_ccd(args)

    mock_plotter.assert_called_once_with(ccd,
                                         mocker.ANY,
                                         ground_q_range=[-1.0, 1.0],
                                         excited_q_range=[0.0, 2.0])
    plotter.construct_plot.assert_called_once_with()
    mock_savefig.assert_called_once_with(args.fig_name)
    mock_show.assert_called_once_with()


def test_plot_eigenvalues(mocker):
    dir_ = Path("disp_0.0")
    orb_info = mocker.Mock()
    single_point = Namespace(disp_ratio=0.0)
    ccd_init = Namespace(supercell_vbm=1.0, supercell_cbm=2.0)
    plotter = mocker.Mock()
    mock_loadfn = mocker.patch(
        "pydefect_ccd.cli.main_function.loadfn",
        side_effect=[orb_info, single_point])
    mock_plotter = mocker.patch(
        "pydefect_ccd.cli.main_function.EigenvalPlotter",
        return_value=plotter)
    args = Namespace(dirs=[dir_ / "disp_0.0"],
                     ccd_init=ccd_init,
                     y_range=[-1.0, 1.0])

    plot_eigenvalues(args)

    assert mock_loadfn.call_args_list == [
        mocker.call(dir_ / "disp_0.0" / "band_edge_orbital_infos.json"),
        mocker.call(dir_ / "disp_0.0" / "single_point.json")]
    mock_plotter.assert_called_once_with([orb_info],
                                         [0.0],
                                         1.0,
                                         2.0,
                                         y_range=[-1.0, 1.0])
    plotter.construct_plot.assert_called_once_with()
    plotter.plt.savefig.assert_called_once_with("eigenvalues.pdf")
    plotter.plt.show.assert_called_once_with()

def test_make_wswq_dirs(tmpdir, mocker):

    print(tmpdir)
    tmpdir.chdir()

    for state, c in zip(["ground", "excited"], [0, 1]):
        Path(f"{state}_original").mkdir(parents=True)
        Path(f"{state}/disp_-0.2").mkdir(parents=True)

        Path(f"{state}_original/WAVECAR").write_text("wave")
        Path(f"{state}/disp_-0.2/KPOINTS").write_text("kpoints")
        Path(f"{state}/disp_-0.2/POSCAR").write_text("poscar")
        Path(f"{state}/disp_-0.2/POTCAR").write_text("potcar")
        Path(f"{state}/disp_-0.2/WAVECAR").write_text("qqq")

        Path(f"{state}/disp_-0.2/prior_info.yaml").write_text(f"charge: {c}")

        incar = ViseIncar({"NSW": 100, "LORBIT": 11})
        incar.write_file(Path(f"{state}/disp_-0.2/INCAR"))

    Path(f"excited/disp_-0.2/wswq").mkdir()

    min_point_info1 = mocker.MagicMock()
    min_point_info2 = mocker.MagicMock()
    min_point_info1.parsed_dir = str(tmpdir / "ground_original")
    min_point_info1.charge = 0

    min_point_info2.parsed_dir = str(tmpdir / "excited_original")
    min_point_info2.charge = 1

    ccd_init = CcdInit(first_relaxed_point=min_point_info1,
                       second_relaxed_point=min_point_info2,
                       vbm=1.0, cbm=2.0,
                       supercell_vbm=1.0, supercell_cbm=2.0)

    args = Namespace(dirs=[Path(f"ground/disp_-0.2"), Path(f"excited/disp_-0.2")],
                     ccd_init=ccd_init)
    make_wswq_dirs(args)

    for state in ["ground"]:
        wswq_dir = Path(f"{state}/disp_-0.2/wswq/")
        assert Path(wswq_dir/"KPOINTS").read_text() == "kpoints"
        assert Path(wswq_dir/"POSCAR").read_text() == "poscar"
        assert Path(wswq_dir/"POTCAR").read_text() == "potcar"
        assert Path(wswq_dir/"WAVECAR.qqq").read_text() == "qqq"
        assert Path(wswq_dir/"WAVECAR").read_text() == "wave"

        actual_incar = ViseIncar.from_file(wswq_dir/"INCAR")
        expected_incar = ViseIncar({"NSW": 100,
                                    "ALGO": "None",
                                    "LWSWQ": True,
                                    "NELM": 1,
                                    "LWAVE": False})
        assert actual_incar == expected_incar

        assert Path(wswq_dir/"KPOINTS").is_symlink()
        assert Path(wswq_dir/"POSCAR").is_symlink()
        assert Path(wswq_dir/"POTCAR").is_symlink()
        assert Path(wswq_dir/"WAVECAR.qqq").is_symlink()
        assert Path(wswq_dir/"WAVECAR").is_symlink()

    assert Path("excited/disp_-0.2/wswq/KPOINTS").exists() is False


def test_make_e_p_matrix_element(mocker):
    dir_0 = Path("disp_0.0")
    dir_1 = Path("disp_0.1")
    potential_curve = Namespace(charge=0)
    base_orbital_infos = mocker.Mock()
    single_point_0 = Namespace(Q=0.0, disp_ratio=0.0)
    single_point_1 = Namespace(Q=0.1, disp_ratio=0.1)
    ep_matrix_element = mocker.Mock(index_info="b767_d766_k1_-1")
    wswq_0 = {(2, 1): {(766, 767): 3.0 + 4.0j}}
    wswq_1 = {(2, 1): {(766, 767): 30.0 + 40.0j}}
    mocker.patch("pydefect_ccd.cli.main_function._read_WSWQ",
                 side_effect=[wswq_0, wswq_1])

    def loadfn_side_effect(filename):
        data = {
            dir_0 / "single_point.json": single_point_0,
            dir_0 / "band_edge_orbital_infos.json": base_orbital_infos,
            dir_1 / "single_point.json": single_point_1}
        return data[filename]

    mocker.patch("pydefect_ccd.cli.main_function.loadfn",
                 side_effect=loadfn_side_effect)
    mock_make = mocker.patch(
        "pydefect_ccd.cli.main_function.make_e_p_matrix_element",
        return_value=ep_matrix_element)
    mocker.patch("pydefect_ccd.cli.main_function.plt.savefig")
    args = Namespace(potential_curve=potential_curve,
                     band_edge_index=767,
                     defect_band_index=766,
                     spin=Spin.down,
                     dirs=[dir_0, dir_1])

    main_make_e_p_matrix_element(args)

    mock_make.assert_called_once_with(
        charge=0,
        base_band_edge_orbital_infos=base_orbital_infos,
        band_edge_index=767,
        defect_band_index=766,
        spin=Spin.down,
        dQs=[0.0, 0.1],
        wswqs=[3.0 + 4.0j, 30.0 + 40.0j])
    ep_matrix_element.to_json_file.assert_called_once_with()


# def test_make_capture_rate(tmpdir, test_files):
#     print(tmpdir)
#     tmpdir.chdir()
#     dir_ = test_files / "NaP/Va_P1_-1__Va_P1_0/from_0_to_-1_after_make_single_point_infos"
#     args = Namespace(base_disp=0.0,
#                      potential_curve=loadfn(dir_/"potential_curve.json"),
#                      captured_carrier=Carrier.e,
#                      band_edge_index=767,
#                      defect_band_index=766,
#                      kpoint_index=1,
#                      spin=Spin.down,
#                      dirs=[dir_/"disp_0.0", dir_/"disp_0.1"],
#                      energy_diff = )
#     #    dirs=[dir_/"disp_0.0", dir_/"disp_0.1", dir_/"disp_0.2"])

    # make_e_p_matrix_element(args)
    # actual: EPMatrixElement = loadfn("e_p_matrix_element.json")
    # print(actual)
"""
TODO:
. Check if electronic SCF are converged.
. Check if displace_ratio=1 structure is the same as the counterpart.
  Consider the energy correction
"""
