# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
from pathlib import Path

import pytest
from pydefect.analyzer.band_edge_states import PerfectBandEdgeState
from pymatgen.electronic_structure.core import Spin

from pydefect_ccd.ccd import Ccd
from pydefect_ccd.ccd_init import CcdInit
from pydefect_ccd.cli.main import parse_args_main
from pydefect_ccd.cli.main_function import main_make_e_p_matrix_element, \
    make_capture_rate, make_ccd, make_ccd_corrections, make_ccd_dirs, \
    make_ccd_init, make_potential_curve, make_single_points, \
    make_sommerfeld_scaling, make_total_squared_transition_moment, \
    make_wswq_dirs, plot_ccd, plot_eigenvalues
from pydefect_ccd.fitting_func import FittingFuncType


def loadfn_effect(d: dict):
    def side_effect(filename):
        try:
            return d[filename]
        except KeyError:
            raise ValueError
    return side_effect


@pytest.mark.parametrize("command", ["make-sommerfeld-scaling", "mss"])
def test_main_make_sommerfeld_scaling(mocker, command):
    mock_effective_mass = mocker.Mock()
    side_effect = loadfn_effect({"effective_mass.json": mock_effective_mass})
    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    args = parse_args_main([command,
                            "-e", "0.5",
                            "--electron-and-hole-effective-mass",
                            "1.0", "2.0",
                            "-Ts", "10", "20"])

    assert args.epsilon0 == pytest.approx(0.5)
    assert args.electron_and_hole_effective_mass == pytest.approx([1.0, 2.0])
    assert args.effective_mass_file is None
    assert args.temperatures == [10.0, 20.0]
    assert args.func is make_sommerfeld_scaling

    args = parse_args_main([command,
                            "-e", "0.5",
                            "--effective-mass-file"])

    assert args.effective_mass_file is mock_effective_mass
    assert args.electron_and_hole_effective_mass is None

    with pytest.raises(SystemExit):
        parse_args_main([command,
                         "-e", "0.5",
                         "--electron-and-hole-effective-mass",
                         "1.0", "2.0",
                         "--effective-mass-file", "effective_mass.json"])


@pytest.mark.parametrize("command", ["make-ccd-init", "mci"])
def test_main_make_ccd_init(mocker, command):
    mock_unitcell = mocker.patch("pydefect.cli.main.Unitcell")
    mock_p_band_edge_state = mocker.Mock(spec=PerfectBandEdgeState, autospec=True)

    side_effect = loadfn_effect(
        {"perfect_band_edge_state.json": mock_p_band_edge_state})
    mocker.patch("pydefect.cli.main.loadfn", side_effect=side_effect)
    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    args = parse_args_main([command,
                            "-fd", "Va_O1_0",
                            "-sd", "Va_O1_1",
                            "-u", "unitcell.yaml",
                            "-pbes", "perfect_band_edge_state.json"])

    assert args.first_dir == Path("Va_O1_0")
    assert args.second_dir == Path("Va_O1_1")
    assert args.unitcell == mock_unitcell.from_yaml.return_value
    assert args.p_state is mock_p_band_edge_state
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


@pytest.mark.parametrize("command", ["make-ccd-corrections", "mccorr"])
def test_main_make_ccd_corrections(mocker, command):
    mock_unitcell = mocker.patch("pydefect.cli.main.Unitcell")
    mock_spec = mocker.Mock()
    mock_no_disp_calc_results = mocker.Mock()
    mock_no_disp_defect_entry = mocker.Mock()
    side_effect = loadfn_effect(
        {"potential_curve_spec.json": mock_spec,
         "calc_results.json": mock_no_disp_calc_results,
         "defect_entry.json": mock_no_disp_defect_entry})
    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    args = parse_args_main(
        [command,
         "-p", "potential_curve_spec.json",
         "-d", "disp_0.0", "disp_0.1",
         "-u", "unitcell.yaml",
         "-ndcr", "calc_results.json",
         "-ndde", "defect_entry.json"])

    assert args.potential_curve_spec is mock_spec
    assert args.dirs == [Path("disp_0.0"), Path("disp_0.1")]
    assert args.unitcell == mock_unitcell.from_yaml.return_value
    assert args.no_disp_calc_results is mock_no_disp_calc_results
    assert args.no_disp_defect_entry is mock_no_disp_defect_entry
    assert args.func is make_ccd_corrections
    mock_unitcell.from_yaml.assert_called_once_with("unitcell.yaml")


@pytest.mark.parametrize("command", ["make-ccd", "mc"])
def test_main_make_ccd(mocker, command):
    mock_ccd_init = mocker.Mock(spec=CcdInit, autospec=True)
    mock_ground_pot_curve = mocker.Mock()
    mock_excited_pot_curve = mocker.Mock()
    side_effect = loadfn_effect(
        {"ccd_init.json": mock_ccd_init,
         "ground_potential_curve.json": mock_ground_pot_curve,
         "excited_potential_curve.json": mock_excited_pot_curve})

    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    args = parse_args_main(
        [command,
         "--ccd_init", "ccd_init.json",
         "--ground-pot-curve", "ground_potential_curve.json",
         "--excited-pot-curve", "excited_potential_curve.json"])

    assert args.ccd_init is mock_ccd_init
    assert args.ground_pot_curve is mock_ground_pot_curve
    assert args.excited_pot_curve is mock_excited_pot_curve
    assert args.func is make_ccd


@pytest.mark.parametrize("command", ["make-single-points", "msp"])
def test_main_make_single_points(command):
    args = parse_args_main([command, "-d", "disp_0.0", "disp_0.1"])

    assert args.dirs == [Path("disp_0.0"), Path("disp_0.1")]
    assert args.parse_ccd_correction is True
    assert args.func is make_single_points

    args = parse_args_main(
        [command, "-d", "disp_0.0", "--no-parse-ccd-correction"])

    assert args.parse_ccd_correction is False


@pytest.mark.parametrize("command", ["make-potential-curve", "mpc"])
def test_main_make_potential_curve(mocker, command):
    mock_spec = mocker.Mock()
    side_effect = loadfn_effect({"potential_curve_spec.json": mock_spec})
    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    args = parse_args_main(
        [command, "-d", "disp_0.0", "--fitting-func", "quadratic"])

    assert args.dirs == [Path("disp_0.0")]
    assert args.potential_curve_spec is mock_spec
    assert args.fitting_func is FittingFuncType.quadratic
    assert args.func is make_potential_curve


def test_main_make_potential_curve_rejects_unknown_fitting_curve(mocker):
    side_effect = loadfn_effect({"potential_curve_spec.json": mocker.Mock()})
    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    with pytest.raises(AttributeError):
        parse_args_main(["mpc", "-d", "disp_0.0", "--fitting-func", "cubic"])


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

    args = parse_args_main(
        [command,
         "--fig-name", "custom_ccd.png",
         "--ground-q-range", "-1.0", "1.0",
         "--excited-q-range", "0.0", "2.0"])

    assert args.fig_name == "custom_ccd.png"
    assert args.ground_q_range == [-1.0, 1.0]
    assert args.excited_q_range == [0.0, 2.0]


@pytest.mark.parametrize("command", ["make-total-squared-transition-moment", "mm"])
def test_main_make_total_squared_transition_moment(mocker, command):
    mock_ccd = mocker.Mock(spec=Ccd, autospec=True)
    mock_sommerfeld = mocker.Mock()
    side_effect = loadfn_effect(
        {"ccd.json": mock_ccd,
         "sommerfeld_scaling.json": mock_sommerfeld})
    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    args = parse_args_main([command, "-s", "sommerfeld_scaling.json"])

    assert args.ccd is mock_ccd
    assert args.sommerfeld is mock_sommerfeld
    assert args.func is make_total_squared_transition_moment


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

    args = parse_args_main([command, "-d", "disp_0.0", "-y", "-1.0", "2.0"])

    assert args.y_range == [-1.0, 2.0]


@pytest.mark.parametrize("command", ["make-wswq-dirs", "mwd"])
def test_main_make_wswq_dirs(mocker, command):
    mock_ccd_init = mocker.Mock(spec=CcdInit, autospec=True)
    side_effect = loadfn_effect({"ccd_init.json": mock_ccd_init})
    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    args = parse_args_main([command, "--dirs", "disp_0.0", "disp_0.1"])

    assert args.ccd_init is mock_ccd_init
    assert args.dirs == [Path("disp_0.0"), Path("disp_0.1")]
    assert args.func is make_wswq_dirs


@pytest.mark.parametrize("command", ["make-e_p_matrix_element", "mepme"])
def test_main_make_e_p_matrix_element(mocker, command):
    mock_potential_curve = mocker.Mock()
    side_effect = loadfn_effect({"potential_curve.json": mock_potential_curve})
    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    args = parse_args_main([command,
                            "--potential_curve", "potential_curve.json",
                            "--band-edge-index", "10",
                            "--defect-band-index", "11",
                            "--spin", "down",
                            "--dirs", "disp_0.0", "disp_0.1"])

    assert args.potential_curve is mock_potential_curve
    assert args.band_edge_index == 10
    assert args.defect_band_index == 11
    assert args.spin == Spin.down
    assert args.dirs == [Path("disp_0.0"), Path("disp_0.1")]
    assert args.func is main_make_e_p_matrix_element


@pytest.mark.parametrize("command", ["make-capture-rate", "mcr"])
def test_main_make_capture_rate(mocker, command):
    mock_ccd_init = mocker.Mock(spec=CcdInit, autospec=True)
    mock_ccd = mocker.Mock(spec=Ccd, autospec=True)
    mock_sommerfeld = mocker.Mock()
    mock_e_p_matrix_element = mocker.Mock()
    mock_total_moment = mocker.Mock()
    side_effect = loadfn_effect(
        {"ccd_init.json": mock_ccd_init,
         "ccd.json": mock_ccd,
         "sommerfeld_scaling.json": mock_sommerfeld,
         "e_p_matrix_element.json": mock_e_p_matrix_element,
         "total_squared_transition_moment.json": mock_total_moment})
    mocker.patch("pydefect_ccd.cli.main.loadfn", side_effect=side_effect)

    args = parse_args_main(
        [command,
         "-s", "sommerfeld_scaling.json",
         "-epme", "e_p_matrix_element.json",
         "-tm", "total_squared_transition_moment.json",
         "-d", "3"])

    assert args.ccd_init is mock_ccd_init
    assert args.ccd is mock_ccd
    assert args.sommerfeld is mock_sommerfeld
    assert args.e_p_matrix_element is mock_e_p_matrix_element
    assert args.total_moment is mock_total_moment
    assert args.degeneracy == 3
    assert args.func is make_capture_rate
