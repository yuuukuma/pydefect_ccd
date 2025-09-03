# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import argparse
import sys
import warnings
from pathlib import Path

from dephon.cli.main_function import make_ccd_init, make_ccd, \
    make_ccd_dirs, plot_ccd, plot_eigenvalues, set_quadratic_fitting_q_range, \
    make_wswq_dirs, update_single_point_infos, add_point_infos_to_single_ccd, \
    make_e_p_matrix_element, make_capture_rate, plot_capture_rate, \
    make_ccd_correction
from dephon.enum import Carrier
from dephon.version import __version__
from monty.serialization import loadfn
from pydefect.cli.main import add_sub_parser
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.inputs import UnknownPotcarWarning

warnings.simplefilter('ignore', UnknownPotcarWarning)

description = """Helper package for calculating the non-radiative carrier 
capture rates trapped by point defects."""

epilog = f"Author: Yu Kumagai Version: {__version__}"


def parse_args_main(args):

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    subparsers = parser.add_subparsers()

    unitcell_parser = add_sub_parser(argparse, name="unitcell")
    pbes_parser = add_sub_parser(argparse, name="perfect_band_edge_state")
    dirs = add_sub_parser(argparse, name="dirs")

    dephon_init = argparse.ArgumentParser(description="", add_help=False)
    dephon_init.add_argument(
            "--dephon_init", type=loadfn, default="dephon_init.json")

    ccd = argparse.ArgumentParser(description="", add_help=False)
    ccd.add_argument(
        "--ccd", type=loadfn, default="ccd.json")

    # -- make_ccd_init -----------------------------------
    parser_make_dephon_init = subparsers.add_parser(
        name="make_ccd_init",
        description="""Create a `dephon_init.json` file from two directories containing pydefect files.
 If the excited state has one more (less) charge state, n-type (p-type) is assumed.""",
        parents=[unitcell_parser, pbes_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mdi'])

    parser_make_dephon_init.add_argument(
        "-fd", "--first_dir", type=Path, required=True,
        help="First directory considered for ccd, e.g., Va_O1_0.")
    parser_make_dephon_init.add_argument(
        "-sd", "--second_dir", type=Path, required=True,
        help="Second directory considered for ccd, e.g., Va_O1_1.")
    parser_make_dephon_init.add_argument(
        "-em", "--effective_mass", type=loadfn, required=True,
        help="effective_mass.json file.")
    parser_make_dephon_init.set_defaults(func=make_ccd_init)

    # -- make_ccd_dirs -----------------------------------
    parser_add_ccd_dirs = subparsers.add_parser(
        name="make_ccd_dirs",
        description=""" Make directories to calculate configuration coordinate 
        diagrams for ground and excited states.""",
        parents=[dephon_init],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mcd'])

    parser_add_ccd_dirs.add_argument(
        "-fsr", "--first_to_second_div_ratios", type=float, nargs="+",
        default=[-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        help="Dividing ratios from first to second charge state structures.")
    parser_add_ccd_dirs.add_argument(
        "-sfr", "--second_to_first_div_ratios", type=float, nargs="+",
        default=[-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        help="Dividing ratios from second to first charge state structures.")
    parser_add_ccd_dirs.add_argument(
        "-d", "--calc_dir", type=Path, default=Path.cwd(),
        help="Directory where directories are created.")
    parser_add_ccd_dirs.set_defaults(
        func=make_ccd_dirs)

    # -- update_single_point_infos -----------------------------------
    parser_update_single_point_infos = subparsers.add_parser(
        name="update_single_point_infos",
        description="Update single_point_info.json at each directory. "
                    "Before running this command, calc_results.json, "
                    "band_edge_orbital_infos.json, and "
                    "band_edge_states.json files need to be created using "
                    "pydefect.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[dirs],
        aliases=['uspi'])

    parser_update_single_point_infos.set_defaults(
        func=update_single_point_infos)

    # -- add_point_infos_to_single_ccd -----------------------------------
    parser_add_point_infos_to_single_ccd = subparsers.add_parser(
        name="add_point_infos_to_single_ccd",
        description="""""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[dirs],
        aliases=['apsc'])

    parser_add_point_infos_to_single_ccd.set_defaults(
        func=add_point_infos_to_single_ccd)

    # -- make_ccd -----------------------------------
    parser_make_ccd = subparsers.add_parser(
        name="make_ccd",
        description="Make ccd.json file",
        parents=[dephon_init],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mc'])

    parser_make_ccd.add_argument(
        "-g", "--ground_ccd", type=loadfn,
        help="potential_curve.json file corresponding to _default_single_ccd_for_e_p_coupling ground state.")
    parser_make_ccd.add_argument(
        "-e", "--excited_ccd", type=loadfn,
        help="potential_curve.json file corresponding to an excited state.")
    parser_make_ccd.set_defaults(func=make_ccd)

    # -- make_ccd_correction -----------------------------------
    parser_make_ccd_correction = subparsers.add_parser(
        name="make_ccd_correction",
        description="Make ccd_correction.json and modify "
                    "dephon_correction.yaml",
        parents=[dirs, unitcell_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mcc'])

    parser_make_ccd_correction.add_argument(
        "-s", "--potential_curve", type=loadfn,
        help="potential_curve.json file.")
    parser_make_ccd_correction.add_argument(
        "--to_charge", type=int)
    parser_make_ccd_correction.add_argument(
        "-ndcr", "--no_disp_calc_results", required=True, type=loadfn,
        help="Path to the calc_results.json without displacement.")
    parser_make_ccd_correction.add_argument(
        "-ndde", "--no_disp_defect_entry", required=True, type=loadfn,
        help="Path to the defect_entry.json without displacement.")
    parser_make_ccd_correction.set_defaults(func=make_ccd_correction)

    # -- set_fitting_q_range -----------------------------------
    parser_set_fitting_q_range = subparsers.add_parser(
        name="set_fitting_q_range",
        description="Set the fitting range for the quadratic potential surface",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['sfr'])

    parser_set_fitting_q_range.add_argument(
        "--ccd", type=loadfn, default="ccd.json")
    parser_set_fitting_q_range.add_argument(
        "--single_ccd_name", type=str, required=True)
    parser_set_fitting_q_range.add_argument(
        "--q_range", type=float, nargs="+")
    parser_set_fitting_q_range.set_defaults(func=set_quadratic_fitting_q_range)

    # -- plot_ccd -----------------------------------
    parser_plot_ccd = subparsers.add_parser(
        name="plot_ccd",
        description="Plot cc diagram from ccd.json file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pc'])

    parser_plot_ccd.add_argument(
        "--ccd", type=loadfn, default="ccd.json")
    parser_plot_ccd.add_argument(
        "--fig_name", type=str, default="ccd.pdf")
    parser_plot_ccd.add_argument(
        "--q_range", type=float, nargs="+")
    parser_plot_ccd.add_argument(
        "--no_quadratic_fit",  dest="quadratic_fit",  action="store_false")
    parser_plot_ccd.add_argument(
        "--no_spline_fit",  dest="spline_fit",  action="store_false")

    parser_plot_ccd.set_defaults(func=plot_ccd)

    # -- plot_eigenvalues -----------------------------------
    parser_plot_eigenvalues = subparsers.add_parser(
        name="plot_eigenvalues",
        parents=[dephon_init, dirs],
        description="Plot eigenvalues as function of displacement ratio. "
                    "band_edge_orbital_infos.json is needed at each directory.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pe'])

    parser_plot_eigenvalues.add_argument(
        "-y", "--y_range", nargs=2, type=float, default=None,
        help="Energy range in y-axis")
    parser_plot_eigenvalues.set_defaults(func=plot_eigenvalues)

    # # -- make_initial_e_p_coupling -----------------------------------
    # parser_make_initial_e_p_coupling = subparsers.add_parser(
    #     name="make_initial_e_p_coupling",
    #     parents=[dephon_init, ccd],
    #     description="Make initial e_p_coupling.json file.",
    #     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    #     aliases=['miepc'])

    # parser_make_initial_e_p_coupling.add_argument(
    #     "-cc", "--captured_carrier", type=Carrier, required=True,
    #     choices=Carrier.name_list())
    # parser_make_initial_e_p_coupling.add_argument(
    #     "-d", "--disp", type=float, default=0.0)
    # parser_make_initial_e_p_coupling.add_argument(
    #     "--charge_for_e_p_coupling", type=int,
    #     help="Default is a charge with smaller absolute value.")
    # parser_make_initial_e_p_coupling.set_defaults(
    #     func=make_initial_e_p_coupling)

    # -- make_wswq_dirs -----------------------------------
    parser_make_wswq_dirs = subparsers.add_parser(
        name="make_wswq_dirs",
        description="Make directories for calculating WSWQ files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mwd'])

    parser_make_wswq_dirs.add_argument(
        "--dephon_init", type=loadfn, default="dephon_init.json")
    parser_make_wswq_dirs.add_argument(
        "--dirs", type=Path, nargs="+", default=[])

    parser_make_wswq_dirs.set_defaults(func=make_wswq_dirs)

    # -- make_e_p_matrix_element -----------------------------------
    parser_make_e_p_matrix_element = subparsers.add_parser(
        name="make_e_p_matrix_element",
        description="Make directories for calculating WSWQ files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mepme'])

    parser_make_e_p_matrix_element.add_argument(
        "--base_disp", type=float, required=True,
        help="Base displacement that must exist in potential_curve.json file.")
    parser_make_e_p_matrix_element.add_argument(
        "--potential_curve", type=loadfn, required=True,
        help="potential_curve.json filename.")
    parser_make_e_p_matrix_element.add_argument(
        "--band_edge_index", type=int)
    parser_make_e_p_matrix_element.add_argument(
        "--defect_band_index", type=int)
    parser_make_e_p_matrix_element.add_argument(
        "--kpoint_index", type=int, required=True)
    parser_make_e_p_matrix_element.add_argument(
        "--spin", type=Spin.__getattr__, required=True)
    parser_make_e_p_matrix_element.add_argument(
        "--energy_diff", type=float)
    parser_make_e_p_matrix_element.add_argument(
        "--dirs", type=Path, nargs="+", required=True)

    parser_make_e_p_matrix_element.set_defaults(func=make_e_p_matrix_element)

    # -- make_capture_rate -----------------------------------
    parser_make_capture_rate = subparsers.add_parser(
        name="make_capture_rate",
        description="Make directories for calculating WSWQ files.",
        parents=[dephon_init, ccd],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mcr'])

    parser_make_capture_rate.add_argument(
        "--captured_carrier", type=Carrier)
    parser_make_capture_rate.add_argument(
        "--e_p_matrix_elem", type=loadfn, required=True)
    parser_make_capture_rate.add_argument(
        "-t", "--temperatures", type=float, nargs="+",
        default=[t for t in range(40, 820, 20)])

    parser_make_capture_rate.set_defaults(func=make_capture_rate)

    # -- plot_capture_rate -----------------------------------
    parser_plot_capture_rate = subparsers.add_parser(
        name="plot_capture_rate",
        description="Plot capture rate",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pcr'])
    parser_plot_capture_rate.add_argument(
        "--capture_rate", type=loadfn, default="capture_rate.json",
        help="capture_rate.json filename.")

    parser_plot_capture_rate.set_defaults(func=plot_capture_rate)
    # ------------------------------------------------------------------------
    return parser.parse_args(args)


def main():
    args = parse_args_main(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()




