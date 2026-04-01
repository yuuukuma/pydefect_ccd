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

from pydefect_ccd.cli.main_function import (
    make_ccd_init, make_ccd, make_ccd_dirs, plot_eigenvalues,
    make_wswq_dirs, main_make_e_p_matrix_element, make_capture_rate,
    make_ccd_corrections, plot_ccd, make_single_points,
    make_potential_curve, make_sommerfeld_scaling, make_total_squared_transition_moment
)
from pydefect_ccd.version import __version__

warnings.simplefilter("ignore", UnknownPotcarWarning)

description = """Package for calculating the 1D configuration coordination
diagram (ccd) of point defects and non-radiative carrier capture rates trapped
by these.
"""

epilog = f"Author: Yu Kumagai, Version: {__version__}"


def parse_args_main(args):
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show version and exit."
    )

    subparsers = parser.add_subparsers(required=True, dest="command")

    # external common parsers
    unitcell_parser = add_sub_parser(argparse, name="unitcell")
    pbes_parser = add_sub_parser(argparse, name="perfect_band_edge_state")
    dirs = add_sub_parser(argparse, name="dirs")

    # internal common parsers
    ccd_init = argparse.ArgumentParser(add_help=False)
    # NOTE: default=loadfn(...) reads file at parse-time.
    # Kept for backward-compatibility.
    ccd_init.add_argument(
        "--ccd_init", type=loadfn, default="ccd_init.json",
        help="ccd_init.json file."
    )

    ccd = argparse.ArgumentParser(add_help=False)
    ccd.add_argument(
        "--ccd", type=loadfn, default="ccd.json",
        help="ccd.json file."
    )

    pot_curve_spec = argparse.ArgumentParser(add_help=False)
    pot_curve_spec.add_argument(
        "-p", "--potential-curve-spec",
        type=loadfn,
        default="potential_curve_spec.json",
        help="potential_curve_spec.json file."
    )

    sommerfeld_scaling = argparse.ArgumentParser(add_help=False)
    sommerfeld_scaling.add_argument("--sommerfeld", "-s", type=loadfn, required=True,
                                    help="sommerfeld_scaling.json file.")

    # -- make-sommerfeld-scaling -----------------------------------
    parser_make_sommerfeld_scaling  = subparsers.add_parser(
        name="make-sommerfeld-scaling",
        description="Calculate sommerfeld scaling parameters and make its file.",
        parents=[unitcell_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=["mss"],
    )
    parser_make_sommerfeld_scaling.add_argument(
       "--electron-effective-mass", "-eem", type=float, required=True)
    parser_make_sommerfeld_scaling.add_argument(
        "--hole-effective-mass", "-hem", type=float, required=True)
    parser_make_sommerfeld_scaling.add_argument(
        "-Ts", "--temperatures", type=float, nargs="+",
        default=[T for T in range(100, 510, 10)],
        help="Temperatures for sommerfeld scaling in K.")
    parser_make_sommerfeld_scaling.set_defaults(func=make_sommerfeld_scaling)

    # -- make-ccd-init -----------------------------------
    parser_make_ccd_init = subparsers.add_parser(
        name="make-ccd-init",
        description=(
            "Create a `ccd_init.json` file from two directories containing "
            "pydefect files. If the excited state has one more (less) charge "
            "state, n-type (p-type) is assumed."
        ),
        parents=[unitcell_parser, pbes_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=["mci"],
    )
    parser_make_ccd_init.add_argument(
        "-fd", "--first-dir", type=Path, required=True,
        help="First directory considered for ccd, e.g., Va_O1_0."
    )
    parser_make_ccd_init.add_argument(
        "-sd", "--second-dir", type=Path, required=True,
        help="Second directory considered for ccd, e.g., Va_O1_1."
    )
    parser_make_ccd_init.add_argument(
        "-em", "--effective-mass", type=loadfn, default=None,
        help="effective_mass.json file."
    )
    parser_make_ccd_init.set_defaults(func=make_ccd_init)

    # -- make-ccd-dirs -----------------------------------
    parser_add_ccd_dirs = subparsers.add_parser(
        name="make-ccd-dirs",
        description="Make directories to calculate CCD between two charge states.",
        parents=[ccd_init],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=["mcdir"],
    )

    parser_add_ccd_dirs.add_argument(
        "-fsr", "--first-to-second-div-ratios", type=float, nargs="+",
        default=[],
        help="Dividing ratios from first to second charge state structures."
    )
    parser_add_ccd_dirs.add_argument(
        "-sfr", "--second-to-first-div-ratios", type=float, nargs="+",
        default=[],
        help="Dividing ratios from second to first charge state structures."
    )
    parser_add_ccd_dirs.add_argument(
        "-d", "--calc-dir", type=Path, default=Path.cwd(),
        help="Directory where directories are created. Default: current directory."
    )
    parser_add_ccd_dirs.set_defaults(func=make_ccd_dirs)

    # -- make-ccd-corrections -----------------------------
    parser_make_ccd_correction = subparsers.add_parser(
        name="make-ccd-corrections",
        description="Make ccd_correction.json files.",
        parents=[pot_curve_spec, dirs, unitcell_parser],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=["mccorr"],
    )
    parser_make_ccd_correction.add_argument(
        "-ndcr", "--no-disp-calc-results", required=True, type=loadfn,
        help="calc_results.json without displacement."
    )
    parser_make_ccd_correction.add_argument(
        "-ndde", "--no-disp-defect-entry", required=True, type=loadfn,
        help="defect_entry.json without displacement."
    )
    parser_make_ccd_correction.set_defaults(func=make_ccd_corrections)

    # -- make-single-point-results ------------------------
    parser_make_single_point_results = subparsers.add_parser(
        name="make-single-point-results",
        description="Make single point results.",
        parents=[dirs],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=["mspr"],
    )
    parser_make_single_point_results.add_argument(
        "--parse-ccd-correction",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="If parse ccd_correction.json. Use --no-parse-ccd-correction to disable."
    )
    parser_make_single_point_results.set_defaults(func=make_single_points)

    # -- make-potential-curve-result ----------------------
    parser_make_potential_curve_result = subparsers.add_parser(
        name="make-potential-curve-result",
        description="Make potential curve result.",
        parents=[pot_curve_spec, dirs],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=["mpcr"],
    )
    parser_make_potential_curve_result.set_defaults(func=make_potential_curve)

    # -- make-ccd -----------------------------------------
    parser_make_ccd = subparsers.add_parser(
        name="make-ccd",
        description="Make ccd.json from potential curves.",
        parents=[ccd_init],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=["mc"],
    )
    parser_make_ccd.add_argument("--ground-potential-curve", type=loadfn)
    parser_make_ccd.add_argument("--excited-potential-curve", type=loadfn)
    parser_make_ccd.add_argument(
        "--fixed-Q0",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Fix Q0 (default: enabled). Use --no-fixed-Q0 to disable."
    )
    parser_make_ccd.set_defaults(func=make_ccd)

    # -- plot-ccd -----------------------------------------
    parser_plot_ccd = subparsers.add_parser(
        name="plot-ccd",
        description="Plot configuration coordinate diagram from ccd.json.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=["pccd"],
    )
    parser_plot_ccd.add_argument("--ccd", type=loadfn, default="ccd.json")
    parser_plot_ccd.add_argument("--fig-name", type=str, default="ccd.pdf")
    parser_plot_ccd.add_argument("--ground-q-range", type=float, nargs="+")
    parser_plot_ccd.add_argument("--excited-q-range", type=float, nargs="+")
    parser_plot_ccd.set_defaults(func=plot_ccd)

    # -- make-total_squared_transition_moment --------------------------------
    parser_make_total_moment = subparsers.add_parser(
        name="make-total-squared-transition-moment",
        parents=[sommerfeld_scaling],
        description="Make total_squared_transition_moment.json.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=["mm"],
    )
    parser_make_total_moment.add_argument("--ccd", type=loadfn, default="ccd.json")
    parser_make_total_moment.set_defaults(func=make_total_squared_transition_moment)

    # -- plot-eigenvalues ---------------------------------
    parser_plot_eigenvalues = subparsers.add_parser(
        name="plot-eigenvalues",
        parents=[ccd_init, dirs],
        description=(
            "Plot eigenvalues as a function of displacement ratio. "
            "band_edge_orbital_infos.json is needed at each directory."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=["peig"],
    )
    parser_plot_eigenvalues.add_argument(
        "-y", "--y-range", nargs=2, type=float, default=None,
        help="Energy range in y-axis."
    )
    parser_plot_eigenvalues.set_defaults(func=plot_eigenvalues)

    # -- make-wswq-dirs -----------------------------------
    parser_make_wswq_dirs = subparsers.add_parser(
        name="make-wswq-dirs",
        parents=[ccd_init],
        description="Make directories for calculating WSWQ files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=["mwd"],
    )
    parser_make_wswq_dirs.add_argument(
        "--dirs", type=Path, nargs="+", required=True
    )
    parser_make_wswq_dirs.set_defaults(func=make_wswq_dirs)

    # -- make-e_p_matrix_element --------------------------
    parser_make_e_p_matrix_element = subparsers.add_parser(
        name="make-e_p_matrix_element",
        description="Make a file for electron-phonon matrix elements.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=["mepme"],
    )
    parser_make_e_p_matrix_element.add_argument(
        "--potential_curve", type=loadfn, required=True,
        help="potential_curve.json filename."
    )
    parser_make_e_p_matrix_element.add_argument("--band-edge-index", type=int)
    parser_make_e_p_matrix_element.add_argument("--defect-band-index", type=int)
    parser_make_e_p_matrix_element.add_argument(
        "--spin", type=Spin.__getitem__, required=True,
        help="Spin channel: up/down or 1/-1 (see pymatgen Spin)."
    )
    parser_make_e_p_matrix_element.add_argument(
        "--dirs", type=Path, nargs="+", required=True
    )
    parser_make_e_p_matrix_element.set_defaults(func=main_make_e_p_matrix_element)

    # -- make-capture-rate --------------------------------
    parser_make_capture_rate = subparsers.add_parser(
        name="make-capture-rate",
        description="Make capture_rate.json.",
        parents=[ccd_init, ccd, sommerfeld_scaling],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        aliases=["mcr"],
    )
    parser_make_capture_rate.add_argument(
        "--e_p_matrix_element", "-epme", type=loadfn, required=True,
        help="e_p_matrix_element.json file.")

    parser_make_capture_rate.add_argument(
        "--total_moment", "-tm", type=loadfn, required=True,
        help="total_squared_transition_moment.json file.")

    parser_make_capture_rate.add_argument(
        "--degeneracy", "-d", type=int, default=1,
        # help="If this is not given, the degeneracy is determined from the final site "
        #      "symmetry. "
    )
    parser_make_capture_rate.set_defaults(func=make_capture_rate)

    try:
        import argcomplete
        argcomplete.autocomplete(parser)
    except ImportError:
        pass

    return parser.parse_args(args)


def main():
    args = parse_args_main(sys.argv[1:])
    args.func(args)


if __name__ == "__main__":
    main()