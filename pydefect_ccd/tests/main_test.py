# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
from pathlib import Path

import pytest
from pydefect.analyzer.band_edge_states import PerfectBandEdgeState
from vise.analyzer.effective_mass import EffectiveMass

from pydefect_ccd.ccd import Ccd
from pydefect_ccd.ccd_init import CcdInit
from pydefect_ccd.cli.main import parse_args_main
from pydefect_ccd.cli.main_function import make_ccd, make_ccd_dirs, \
    make_ccd_init, make_sommerfeld_scaling, plot_ccd, plot_eigenvalues


def loadfn_effect(d: dict):
    def side_effect(filename):
        try:
            return d[filename]
        except KeyError:
            raise ValueError
    return side_effect


@pytest.mark.parametrize("command", ["make-sommerfeld-scaling", "mss"])
def test_main_make_sommerfeld_scaling(command):
    args = parse_args_main([command,
                            "-e", "0.5",
                            "-eem", "1.0",
                            "-hem", "2.0",
                            "-Ts", "10", "20"])

    assert args.epsilon0 == pytest.approx(0.5)
    assert args.electron_effective_mass == pytest.approx(1.0)
    assert args.hole_effective_mass == pytest.approx(2.0)
    assert args.temperatures == [10.0, 20.0]
    assert args.func is make_sommerfeld_scaling


@pytest.mark.parametrize("command", ["make-ccd-init", "mci"])
def test_main_make_ccd_init(mocker, command):
    mock_unitcell = mocker.patch("pydefect.cli.main.Unitcell")
    mock_p_band_edge_state = mocker.Mock(spec=PerfectBandEdgeState, autospec=True)
    mock_effective_mass = mocker.Mock(spec=EffectiveMass, autospec=True)

    side_effect = loadfn_effect(
        {"perfect_band_edge_state.json": mock_p_band_edge_state,
         "effective_mass.json": mock_effective_mass})
    mocker.patch("pydefect.cli.main.loadfn", side_effect=side_effect)
    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    args = parse_args_main([command,
                            "-fd", "Va_O1_0",
                            "-sd", "Va_O1_1",
                            "-u", "unitcell.yaml",
                            "-pbes", "perfect_band_edge_state.json",
                            "-em", "effective_mass.json"])

    assert args.first_dir == Path("Va_O1_0")
    assert args.second_dir == Path("Va_O1_1")
    assert args.unitcell == mock_unitcell.from_yaml.return_value
    assert args.p_state is mock_p_band_edge_state
    assert args.effective_mass is mock_effective_mass
    assert args.func is make_ccd_init
    mock_unitcell.from_yaml.assert_called_once_with("unitcell.yaml")


@pytest.mark.parametrize("command", ["make-ccd-dirs", "mcdir"])
def test_main_make_dirs(mocker, command):
    mock_ccd_init = mocker.Mock(spec=CcdInit, autospec=True)
    side_effect = loadfn_effect({"ccd_init.json": mock_ccd_init})
    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    args = parse_args_main([command])

    assert args.ccd_init is mock_ccd_init
    assert args.first_to_second_div_ratios == []
    assert args.second_to_first_div_ratios == []
    assert args.calc_dir == Path.cwd()
    assert args.func is make_ccd_dirs

    args = parse_args_main([command, "-fsr", "0.0", "-sfr", "0.1",
                            "-d", "dirname"])

    assert args.ccd_init is mock_ccd_init
    assert args.first_to_second_div_ratios == [0.0]
    assert args.second_to_first_div_ratios == [0.1]
    assert args.calc_dir == Path("dirname")
    assert args.func is make_ccd_dirs


@pytest.mark.parametrize("command", ["make-ccd", "mc"])
def test_main_make_ccd(mocker, command):
    mock_ccd_init = mocker.Mock(spec=CcdInit, autospec=True)
    mock_ground_single_points = mocker.Mock()
    mock_ground_spec = mocker.Mock()
    mock_ground_fitting_func = mocker.Mock()
    mock_excited_single_points = mocker.Mock()
    mock_excited_spec = mocker.Mock()
    mock_excited_fitting_func = mocker.Mock()
    side_effect = loadfn_effect(
        {"ccd_init.json": mock_ccd_init,
         "ground_single_points.json": mock_ground_single_points,
         "ground_spec.json": mock_ground_spec,
         "ground_fitting_func.json": mock_ground_fitting_func,
         "excited_single_points.json": mock_excited_single_points,
         "excited_spec.json": mock_excited_spec,
         "excited_fitting_func.json": mock_excited_fitting_func})

    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    args = parse_args_main(
        [command,
         "--ccd_init", "ccd_init.json",
         "--ground-single-points", "ground_single_points.json",
         "--ground-potential-curve-spec", "ground_spec.json",
         "--ground-fitting-func", "ground_fitting_func.json",
         "--excited-single-points", "excited_single_points.json",
         "--excited-potential-curve-spec", "excited_spec.json",
         "--excited-fitting-func", "excited_fitting_func.json"])

    assert args.ccd_init is mock_ccd_init
    assert args.ground_single_points is mock_ground_single_points
    assert args.ground_potential_curve_spec is mock_ground_spec
    assert args.ground_fitting_func is mock_ground_fitting_func
    assert args.excited_single_points is mock_excited_single_points
    assert args.excited_potential_curve_spec is mock_excited_spec
    assert args.excited_fitting_func is mock_excited_fitting_func
    assert args.func is make_ccd


@pytest.mark.parametrize("command", ["plot-ccd", "pccd"])
def test_main_plot_ccd_wo_args(mocker, command):
    mock_ccd = mocker.Mock(spec=Ccd, autospec=True)
    side_effect = loadfn_effect({"ccd.json": mock_ccd})
    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    args = parse_args_main([command])

    assert args.ccd is mock_ccd
    assert args.fig_name == "ccd.pdf"
    assert args.ground_q_range is None
    assert args.excited_q_range is None
    assert args.func is plot_ccd


@pytest.mark.parametrize("command", ["plot-eigenvalues", "peig"])
def test_main_plot_eigenvalues(mocker, command):
    mock_ccd_init = mocker.Mock(spec=CcdInit, autospec=True)
    side_effect = loadfn_effect({"ccd_init.json": mock_ccd_init})
    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    args = parse_args_main([command, "-d", "disp_0.0"])

    assert args.dirs == [Path("disp_0.0")]
    assert args.ccd_init is mock_ccd_init
    assert args.y_range is None
    assert args.func is plot_eigenvalues


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
