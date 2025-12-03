# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import argparse
import sys
import warnings
from pathlib import Path

from monty.serialization import loadfn
from pydefect.cli.main import add_sub_parser
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.inputs import UnknownPotcarWarning

from pydefect_ccd.cli.main_function import make_ccd_init, make_ccd, \
    make_ccd_dirs, plot_eigenvalues, make_wswq_dirs, \
    make_e_p_matrix_element, make_capture_rate, plot_capture_rate, \
    make_ccd_corrections, plot_ccd, make_single_points, make_potential_curve
from pydefect_ccd.version import __version__

warnings.simplefilter('ignore', UnknownPotcarWarning)

description = """Package for calculating the 1D configuration coordination 
diagram (ccd) of point defects and non-radiative carrier capture rates trapped 
by these."""

epilog = f"Author: Yu Kumagai, Version: {__version__}"


def parse_args_main(args):

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    subparsers = parser.add_subparsers()

    unitcell_parser = add_sub_parser(argparse, name="unitcell")
    pbes_parser = add_sub_parser(argparse, name="perfect_band_edge_state")
    dirs = add_sub_parser(argparse, name="dirs")

    ccd_init = argparse.ArgumentParser(description="", add_help=False)
    ccd_init.add_argument("--ccd_init", type=loadfn, default="ccd_init.json")

    ccd = argparse.ArgumentParser(description="", add_help=False)
    ccd.add_argument(
        "--ccd", type=loadfn, default="ccd.json")

    pot_curve_spec = argparse.ArgumentParser(description="", add_help=False)
    pot_curve_spec.add_argument(
        "-p", "--potential_curve_spec", type=loadfn,
        help="potential_curve.json file.", default="potential_curve_spec.json")

    # -- make_ccd_init -----------------------------------
    parser_make_ccd_init = subparsers.add_parser(
        name="make_ccd_init",
        description="""Create a `ccd_init.json` file from two directories 
containing pydefect files. If the excited state has one more (less) charge 
state, n-type (p-type) is assumed.""",
        parents=[unitcell_parser, pbes_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['ci'])

    parser_make_ccd_init.add_argument(
        "-fd", "--first_dir", type=Path, required=True,
        help="First directory considered for ccd, e.g., Va_O1_0.")
    parser_make_ccd_init.add_argument(
        "-sd", "--second_dir", type=Path, required=True,
        help="Second directory considered for ccd, e.g., Va_O1_1.")
    parser_make_ccd_init.add_argument(
        "-em", "--effective_mass", type=loadfn, default=None,
        help="effective_mass.json file.")
    parser_make_ccd_init.set_defaults(func=make_ccd_init)

    # -- make_ccd_dirs -----------------------------------
    parser_add_ccd_dirs = subparsers.add_parser(
        name="make_ccd_dirs",
        description=""" Make directories to calculate CCD between two charge 
states.""",
        parents=[ccd_init],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['cd'])

    parser_add_ccd_dirs.add_argument(
        "-fsr", "--first_to_second_div_ratios", type=float, nargs="+",
        default=[-0.2, -0.1, 0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0],
        help="Dividing ratios from first to second charge state structures.")
    parser_add_ccd_dirs.add_argument(
        "-sfr", "--second_to_first_div_ratios", type=float, nargs="+",
        default=[-0.2, -0.1, 0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0],
        help="Dividing ratios from second to first charge state structures.")
    parser_add_ccd_dirs.add_argument(
        "-d", "--calc_dir", type=Path, default=Path.cwd(),
        help="Directory where directories are created.")
    parser_add_ccd_dirs.set_defaults(func=make_ccd_dirs)

    # -- make_ccd_corrections -----------------------------------
    parser_make_ccd_correction = subparsers.add_parser(
        name="make_ccd_corrections",
        description="Make ccd_correction.json files",
        parents=[pot_curve_spec, dirs, unitcell_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mcc'])

    parser_make_ccd_correction.add_argument(
        "-ndcr", "--no_disp_calc_results", required=True, type=loadfn,
        help="Path to the calc_results.json without displacement.")
    parser_make_ccd_correction.add_argument(
        "-ndde", "--no_disp_defect_entry", required=True, type=loadfn,
        help="Path to the defect_entry.json without displacement.")
    parser_make_ccd_correction.set_defaults(func=make_ccd_corrections)

    # -- make_single_point_results -----------------------------------
    parser_make_single_point_results = subparsers.add_parser(
        name="make_single_point_results",
        description="",
        parents=[dirs],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mspr'])
    parser_make_single_point_results.set_defaults(func=make_single_points)

    # -- make_potential_curve_result -----------------------------------
    parser_make_potential_curve_result = subparsers.add_parser(
        name="make_potential_curve_result",
        description="",
        parents=[pot_curve_spec, dirs],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mpcr'])
    parser_make_potential_curve_result.set_defaults(func=make_potential_curve)

    # -- make_ccd -----------------------------------
    parser_make_ccd = subparsers.add_parser(
        name="make_ccd",
        description="",
        parents=[ccd_init],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mc'])
    parser_make_ccd.add_argument("--ground_potential_curve", type=loadfn)
    parser_make_ccd.add_argument("--excited_potential_curve", type=loadfn)

    parser_make_ccd.set_defaults(func=make_ccd)

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
        "--ground_q_range", type=float, nargs="+")
    parser_plot_ccd.add_argument(
        "--excited_q_range", type=float, nargs="+")
    parser_plot_ccd.add_argument(
        "--no_quadratic_fit",  dest="quadratic_fit",  action="store_false")
    parser_plot_ccd.add_argument(
        "--no_spline_fit",  dest="spline_fit",  action="store_false")

    parser_plot_ccd.set_defaults(func=plot_ccd)

    # -- plot_eigenvalues -----------------------------------
    parser_plot_eigenvalues = subparsers.add_parser(
        name="plot_eigenvalues",
        parents=[ccd_init, dirs],
        description="Plot eigenvalues as function of displacement ratio. "
                    "band_edge_orbital_infos.json is needed at each directory.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['pe'])

    parser_plot_eigenvalues.add_argument(
        "-y", "--y_range", nargs=2, type=float, default=None,
        help="Energy range in y-axis")
    parser_plot_eigenvalues.set_defaults(func=plot_eigenvalues)

    # -- make_wswq_dirs -----------------------------------
    parser_make_wswq_dirs = subparsers.add_parser(
        name="make_wswq_dirs",
        description="Make directories for calculating WSWQ files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mwd'])

    parser_make_wswq_dirs.add_argument(
        "--ccd_init", type=loadfn, default="ccd_init.json")
    parser_make_wswq_dirs.add_argument(
        "--dirs", type=Path, nargs="+", default=[])

    parser_make_wswq_dirs.set_defaults(func=make_wswq_dirs)

    # -- make_e_p_matrix_element -----------------------------------
    parser_make_e_p_matrix_element = subparsers.add_parser(
        name="make_e_p_matrix_element",
        description="Make directories for calculating WSWQ files.",
        parents=[ccd_init],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mepme'])

    parser_make_e_p_matrix_element.add_argument(
        "--potential_curve", type=loadfn, required=True,
        help="potential_curve.json filename.")
    parser_make_e_p_matrix_element.add_argument(
        "--band_edge_index", type=int)
    parser_make_e_p_matrix_element.add_argument(
        "--defect_band_index", type=int)
    parser_make_e_p_matrix_element.add_argument(
        "--spin", type=Spin.__getitem__, required=True)
    parser_make_e_p_matrix_element.add_argument(
        "-T", "--temperatures", type=float, nargs="+",
        help="Temperatures for calculating capture rates in K.",
        default=[t*100 for t in range(2, 9)])
    parser_make_e_p_matrix_element.add_argument(
        "--dirs", type=Path, nargs="+", required=True)

    parser_make_e_p_matrix_element.set_defaults(func=make_e_p_matrix_element)

    # -- make_capture_rate -----------------------------------
    parser_make_capture_rate = subparsers.add_parser(
        name="make_capture_rate",
        description="Make ",
        parents=[ccd_init, ccd],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=['mcr'])

    parser_make_capture_rate.add_argument(
        "--e_p_coupling", type=loadfn, required=True)
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




