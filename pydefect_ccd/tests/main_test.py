# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from argparse import Namespace
from pathlib import Path

from pydefect.analyzer.band_edge_states import PerfectBandEdgeState
from vise.analyzer.effective_mass import EffectiveMass

from pydefect_ccd.ccd import Ccd
from pydefect_ccd.potential_curve import PotentialCurve
from pydefect_ccd.ccd_init import CcdInit
from pydefect_ccd.cli.main import parse_args_main
from pydefect_ccd.cli.main_function import make_ccd, plot_eigenvalues


def loadfn_effect(d: dict):
    def side_effect(filename):
        try:
            return d[filename]
        except KeyError:
            raise ValueError
    return side_effect


def test_main_make_ccd_init(mocker):
    mock_unitcell = mocker.patch("pydefect.cli.main.Unitcell")
    mock_p_band_edge_state = mocker.Mock(spec=PerfectBandEdgeState, autospec=True)
    mock_effective_mass = mocker.Mock(spec=EffectiveMass, autospec=True)

    side_effect = loadfn_effect(
        {"perfect_band_edge_state.json": mock_p_band_edge_state,
         "effective_mass.json": mock_effective_mass})
    mocker.patch("pydefect.cli.main.loadfn", side_effect=side_effect)
    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)
    parsed_args = parse_args_main(["mdi",
                                   "-fd", "Va_O1_0",
                                   "-sd", "Va_O1_1",
                                   "-u", "unitcell.yaml",
                                   "-pbes", "perfect_band_edge_state.json",
                                   "-em", "effective_mass.json"
                                   ])
    expected = Namespace(
        first_dir=Path("Va_O1_0"),
        second_dir=Path("Va_O1_1"),
        unitcell=mock_unitcell.from_yaml.return_value,
        p_state=mock_p_band_edge_state,
        effective_mass=mock_effective_mass,
        func=parsed_args.func)
    assert parsed_args == expected
    mock_unitcell.from_yaml.assert_called_once_with("unitcell.yaml")


def test_main_make_dirs(mocker):
    mock_ccd_init = mocker.Mock(spec=CcdInit, autospec=True)
    side_effect = loadfn_effect({"ccd_init.json": mock_ccd_init})
    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    parsed_args = parse_args_main(["mcd"])
    expected = Namespace(
        ccd_init=mock_ccd_init,
        first_to_second_div_ratios=[-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        second_to_first_div_ratios=[-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        calc_dir=Path.cwd(),
        func=parsed_args.func)
    assert parsed_args == expected

    parsed_args = parse_args_main(["mcd", "-fsr", "0.0", "-sfr", "0.1",
                                   "-d", "dirname"])
    expected = Namespace(
        ccd_init=mock_ccd_init,
        first_to_second_div_ratios=[0.0],
        second_to_first_div_ratios=[0.1],
        calc_dir=Path("dirname"),
        func=parsed_args.func)
    assert parsed_args == expected


def test_main_make_ccd_wo_args(mocker):
    mock_ccd_init = mocker.Mock(spec=CcdInit, autospec=True)
    mock_ground_single_ccd = mocker.Mock(spec=PotentialCurve, autospec=True)
    mock_excited_single_ccd = mocker.Mock(spec=PotentialCurve, autospec=True)
    side_effect = loadfn_effect(
        {"ccd_init.json": mock_ccd_init,
         "ground_single_ccd.json": mock_ground_single_ccd,
         "excited_single_ccd.json": mock_excited_single_ccd})

    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    parsed_args = parse_args_main(["mc",
                                   "-g", "ground_single_ccd.json",
                                   "-e", "excited_single_ccd.json"])
    expected = Namespace(
        ccd_init=mock_ccd_init,
        ground_ccd=mock_ground_single_ccd,
        excited_ccd=mock_excited_single_ccd,
        func=make_ccd)
    assert parsed_args == expected


def test_main_set_quadratic_fitting_q_range(mocker):
    mock_ccd = mocker.Mock(spec=Ccd, autospec=True)
    side_effect = loadfn_effect({"ccd.json": mock_ccd})
    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    parsed_args = parse_args_main(
        ["sfr", "--single_ccd_name", "ground", "--q_range", "-1.0", "1.0"])
    expected = Namespace(
        ccd=mock_ccd,
        single_ccd_name="ground",
        q_range=[-1.0, 1.0],
        func=parsed_args.func)
    assert parsed_args == expected


def test_main_plot_ccd_wo_args(mocker):
    mock_ccd = mocker.Mock(spec=Ccd, autospec=True)
    side_effect = loadfn_effect({"ccd.json": mock_ccd})
    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    parsed_args = parse_args_main(["pc"])
    expected = Namespace(
        ccd=mock_ccd,
        fig_name="ccd.pdf",
        q_range=None,
        quadratic_fit=True,
        spline_fit=True,
        func=parsed_args.func)
    assert parsed_args == expected


def test_main_plot_eigenvalues(mocker):
    mock_ccd_init = mocker.Mock(spec=CcdInit, autospec=True)
    side_effect = loadfn_effect({"ccd_init.json": mock_ccd_init})
    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    parsed_args = parse_args_main(["pe", "-d", "disp_0.0"])
    expected = Namespace(
        dirs=[Path("disp_0.0")],
        ccd_init=mock_ccd_init,
        y_range=None,
        func=plot_eigenvalues)
    assert parsed_args == expected


# def test_make_make_initial_e_p_coupling(mocker):
#     mock_dephon_init = mocker.Mock(spec=DephonInit, autospec=True)
#     mock_ccd = mocker.Mock(spec=Ccd, autospec=True)
#     mocker.patch(
#         "dephon.cli.main.loadfn",
#         side_effect=loadfn_effect({"ccd_init.json": mock_dephon_init,
#                                    "ccd.json": mock_ccd}))
#
#     parsed_args = parse_args_main(["miepc",
#                                    "-cc", "h",
#                                    "--charge_for_e_p_coupling", "1"])
#     expected = Namespace(
#         ccd_init=mock_dephon_init,
#         ccd=mock_ccd,
#         captured_carrier=Carrier.h,
#         charge_for_e_p_coupling=1,
#         func=make_initial_e_p_coupling)
#     assert parsed_args == expected

