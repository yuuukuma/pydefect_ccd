# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from argparse import Namespace
from pathlib import Path

from monty.serialization import loadfn
from pydefect.analyzer.unitcell import Unitcell
from pymatgen.core import Structure
from pymatgen.electronic_structure.core import Spin
from vise.input_set.incar import ViseIncar
from vise.input_set.prior_info import PriorInfo

from pydefect_ccd.ccd import SinglePointSpec, PotentialCurveSpec
from pydefect_ccd.ccd_init import CcdInit
from pydefect_ccd.cli.main_function import make_ccd_init, make_ccd, plot_ccd, \
    make_ccd_dirs, make_wswq_dirs, plot_eigenvalues, main_make_e_p_matrix_element
from pydefect_ccd.e_p_matrix_element import EPMatrixElement
from pydefect_ccd.enum import Carrier
from pydefect_ccd.relaxed_point import NearEdgeState, RelaxedPoint


def test_make_ccd_init(test_files, tmpdir):
    tmpdir.chdir()
    dir_ = test_files / "Na3AgO2"

    args = Namespace(first_dir=dir_ / "Va_O1_1",
                     second_dir=dir_ / "Va_O1_0",
                     unitcell=Unitcell.from_yaml(dir_ / "unitcell.yaml"),
                     p_state=loadfn(dir_ / "perfect_band_edge_state.json"),
                     effective_mass=loadfn(dir_ / "effective_mass.json"))
    make_ccd_init(args)
    print(loadfn("Va_O1_1⇆Va_O1_0/ccd_init.json"))


def test_make_ccd_dirs(tmpdir, ground_structure, excited_structure,
                       intermediate_structure):
    print(tmpdir)
    tmpdir.chdir()
    ccd_init = CcdInit(
        relaxed_points=[RelaxedPoint(name="test",
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
                        RelaxedPoint(name="test",
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
                        ],
        vbm=-100.0, cbm=100.0, supercell_volume=10.0,
        supercell_vbm=-100.0, supercell_cbm=100.0,
        ave_electron_mass=1.0, ave_hole_mass=1.0, ave_static_diele_const=1.0)

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
    expected = SinglePointSpec(dQ=dQ/2, disp_ratio=0.5)
    assert actual == expected

    actual = loadfn("q_0/disp_0.0/single_point_spec.json")
    expected = SinglePointSpec(dQ=0.0, disp_ratio=0.0)
    assert actual == expected

    actual = Structure.from_file("q_1/disp_1.0/POSCAR")
    assert actual == excited_structure

    actual = loadfn("q_0/potential_curve_spec.json")
    expected = PotentialCurveSpec(charge=0, correction_energy=100.0,
                                  counter_charge=1, Q_diff=dQ)
    assert actual == expected


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


def test_make_ccd(test_files, tmpdir):
    tmpdir.chdir()
    va_p1 = Path(test_files) / "NaP/Va_P1_-1__Va_P1_0"
    ground_ccd = loadfn(va_p1 / "from_-1_to_0_after_make_single_point_infos/potential_curve.json")
    excited_ccd = loadfn(va_p1 / "from_0_to_-1_after_make_single_point_infos/potential_curve.json")
    dephon_init = loadfn(va_p1 / "ccd_init.json")
    args = Namespace(ground_ccd=ground_ccd, excited_ccd=excited_ccd,
                     dephon_init=dephon_init)
    make_ccd(args)


def test_plot_ccd(ccd, tmpdir):
    args = Namespace(ccd=ccd, fig_name=tmpdir / "ccd.pdf",
                     q_range=[-1.0, 1.0],
                     quadratic_fit=True,
                     spline_fit=True)
    plot_ccd(args)


def test_plot_eigenvalues(test_files, tmpdir):
    tmpdir.chdir()
    va_p1 = Path(test_files) / "NaP/Va_P1_-1__Va_P1_0"
    dephon_init = loadfn(va_p1 / "ccd_init.json")
    dir_ = va_p1 / "from_0_to_-1_before_make_single_point_infos"
    args = Namespace(dirs=[dir_ / "disp_0.0"],
                     dephon_init=dephon_init,
                     y_range=None)
    plot_eigenvalues(args)


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

    dephon_init = CcdInit(relaxed_points=[min_point_info1, min_point_info2],
                          vbm=1.0, cbm=2.0,
                          supercell_volume=10.0,
                          supercell_vbm=1.0, supercell_cbm=2.0,
                          ave_electron_mass=1.0, ave_hole_mass=1.0,
                          ave_static_diele_const=1.0)

    args = Namespace(dirs=[Path(f"ground/disp_-0.2"), Path(f"excited/disp_-0.2")],
                     dephon_init=dephon_init)
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


def test_make_e_p_matrix_element(tmpdir, test_files):
    print(tmpdir)
    tmpdir.chdir()
    dir_ = test_files / "NaP/Va_P1_-1__Va_P1_0/from_0_to_-1_after_make_single_point_infos"
    args = Namespace(base_disp=0.0,
                     single_ccd=loadfn(dir_/"potential_curve.json"),
                     captured_carrier=Carrier.e,
                     band_edge_index=767,
                     defect_band_index=766,
                     kpoint_index=1,
                     spin=Spin.down,
                     dirs=[dir_/"disp_0.0", dir_/"disp_0.1"],
                     energy_diff=1.0)
#    dirs=[dir_/"disp_0.0", dir_/"disp_0.1", dir_/"disp_0.2"])

    main_make_e_p_matrix_element(args)
    actual: EPMatrixElement = loadfn("e_p_matrix_element_b767_d766_k1_-1.json")
    print(actual)


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