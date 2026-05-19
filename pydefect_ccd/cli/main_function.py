# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
import os
from argparse import Namespace
from pathlib import Path
from typing import List, Tuple, Optional, Type, Union

import numpy as np
import yaml
from matplotlib import pyplot as plt
from monty.serialization import loadfn
from nonrad.elphon import _read_WSWQ
from nonrad.scaling import thermal_velocity
from numpy.linalg import LinAlgError
from pydefect.analyzer.band_edge_states import BandEdgeStates, \
    BandEdgeOrbitalInfos
from pydefect.analyzer.calc_results import CalcResults
from pydefect.analyzer.defect_energy import DefectEnergyInfo
from pydefect.analyzer.defect_structure_info import DefectStructureInfo
from pydefect.cli.main_functions import get_calc_results
from pydefect.cli.main_tools import parse_dirs
from pydefect.cli.vasp.make_efnv_correction import make_efnv_correction
from pydefect.corrections.site_potential_plotter import SitePotentialMplPlotter
from pymatgen.core import Structure
from pymatgen.electronic_structure.core import Spin
from vise.input_set.incar import ViseIncar
from vise.input_set.prior_info import PriorInfo
from vise.util.file_transfer import FileLink
from vise.util.logger import get_logger

from pydefect_ccd.capture_rate import \
    CaptureRate, CaptureRatePlotter
from pydefect_ccd.fitting_func import FittingFuncType, FittingFunc
from pydefect_ccd.local_enum import Carrier
from pydefect_ccd.sommerfeld_scaling import SommerfeldScaling
from pydefect_ccd.transition_moment import CalcTotalSquaredTransitionMoment, \
    PlottedTotalSquaredTransitionMoment
from pydefect_ccd.ccd import CcdPlotter, \
    NoCcdCorrection
from pydefect_ccd.potential_curve import SinglePointSpec, SinglePoint, \
    CurveTransform, PotentialCurveSpec, PotentialCurve, SinglePoints, \
    make_curve_transform
from pydefect_ccd.ccd_init import CcdInit
from pydefect_ccd.make_ccd import MakeCcd
from pydefect_ccd.make_e_p_matrix_element import make_e_p_matrix_element
from pydefect_ccd.plot_eigenvalues import EigenvalPlotter
from pydefect_ccd.relaxed_point import NearEdgeState, RelaxedPoint
from pydefect_ccd.util import spin_to_idx

logger = get_logger(__name__)


def make_sommerfeld_scaling(args: Namespace):
    epsilon0 = (args.epsilon0 if args.epsilon0 is not None
                else args.unitcell.ave_diele)
    if args.electron_and_hole_effective_mass is not None:
        electron_effective_mass, hole_effective_mass = \
            args.electron_and_hole_effective_mass
    else:
        concentration = args.effective_mass_file.concentrations[0]
        electron_effective_mass = args.effective_mass_file.average_mass(
            "n", concentration)
        hole_effective_mass = args.effective_mass_file.average_mass(
            "p", concentration)

    scaling = SommerfeldScaling(epsilon0=epsilon0,
                                electron_effective_mass=electron_effective_mass,
                                hole_effective_mass=hole_effective_mass,
                                Ts=args.temperatures)
    scaling.to_json_file()
    _plot_sommerfeld_scaling(scaling,
                             title="attractive",
                             carrier_charges=[(Carrier.e, 1),
                                              (Carrier.h, -1)],
                             filename="sommerfeld_scaling_attractive.pdf")
    _plot_sommerfeld_scaling(scaling,
                             title="repulsive",
                             carrier_charges=[(Carrier.e, -1),
                                              (Carrier.h, 1)],
                             filename="sommerfeld_scaling_repulsive.pdf")


def _plot_sommerfeld_scaling(scaling: SommerfeldScaling,
                             title: str,
                             carrier_charges: List[Tuple[Carrier, int]],
                             filename: str):
    fig, ax = plt.subplots()
    for carrier_type, defect_charge in carrier_charges:
        scaling.add_to_ax(ax=ax,
                          carrier_type=carrier_type,
                          defect_charge=defect_charge)
        ax.lines[-1].set_label(str(carrier_type))
    ax.set_title(title)
    scaling.set_label(ax)
    fig.tight_layout()
    fig.savefig(filename)
    plt.close(fig)


def _make_near_edge_states(band_edge_orbital_infos: BandEdgeOrbitalInfos,
                           spin: Spin,
                           edge_energy: float,
                           threshold: float = 0.1):
    result = []
    info_by_spin = band_edge_orbital_infos.orbital_infos[spin_to_idx(spin)]
    for k_idx, (orb_info_by_kpt, k_coords) in \
        enumerate(zip(info_by_spin, band_edge_orbital_infos.kpt_coords)):
        k_idx_from_1 = k_idx + 1
        for rel_idx, info_by_band in enumerate(orb_info_by_kpt):
            e_from_band_edge = abs(info_by_band.energy - edge_energy)
            if e_from_band_edge > threshold:
                continue
            band_idx = rel_idx + band_edge_orbital_infos.lowest_band_index + 1
            result.append(NearEdgeState(band_idx,
                                        k_coords,
                                        k_idx_from_1,
                                        info_by_band.energy,
                                        info_by_band.occupation))
    return result


def _make_relaxed_point_from_dir(_dir: Path):
    energy_info = DefectEnergyInfo.from_yaml(_dir / "defect_energy_info.yaml")
    calc_results: CalcResults = loadfn(_dir / "calc_results.json")
    defect_structure_info: DefectStructureInfo \
        = loadfn(_dir / "defect_structure_info.json")
    band_edge_states: BandEdgeStates = loadfn(_dir / "band_edge_states.json")

    band_edge_orbital_infos: BandEdgeOrbitalInfos = \
        loadfn(_dir / "band_edge_orbital_infos.json")
    localized_orbitals, valence_bands, conduction_bands = _get_band_edge_info(
        band_edge_orbital_infos, band_edge_states)
    relaxed_point = RelaxedPoint(
        name=energy_info.name,
        charge=energy_info.charge,
        structure=calc_results.structure,
        formation_energy=energy_info.defect_energy.formation_energy,
        correction_energy=energy_info.defect_energy.total_correction,
        magnetization=calc_results.magnetization,
        localized_orbitals=localized_orbitals,
        initial_site_symmetry=defect_structure_info.initial_site_sym,
        final_site_symmetry=defect_structure_info.final_site_sym,
        parsed_dir=str(_dir.absolute()),
        valence_bands=valence_bands,
        conduction_bands=conduction_bands)

    return relaxed_point


def _get_band_edge_info(band_edge_orbital_infos: BandEdgeOrbitalInfos,
                        band_edge_states: BandEdgeStates
                        ) -> Tuple[List, List, List]:
    localized_orbitals, valence_bands, conduction_bands = [], [], []
    for state, spin in zip(band_edge_states.states, [Spin.up, Spin.down]):
        los = state.localized_orbitals
        for lo in los:
            lo.band_idx += 1
        localized_orbitals.append(los)
        valence_bands.append(_make_near_edge_states(band_edge_orbital_infos,
                                                    spin,

                                                    state.vbm_info.energy))
        conduction_bands.append(_make_near_edge_states(band_edge_orbital_infos,
                                                       spin,
                                                       state.cbm_info.energy))
    return localized_orbitals, valence_bands, conduction_bands


def make_ccd_init(args: Namespace):
    min_point_1 = _make_relaxed_point_from_dir(args.first_dir)
    min_point_2 = _make_relaxed_point_from_dir(args.second_dir)

    path = Path(f"{min_point_1.full_name}_{min_point_2.charge}")
    json_file = path / "ccd_init.json"

    if json_file.exists():
        logger.info(f"{json_file} exists. Remove it first to recreate it.")
        return

    charge_diff = min_point_1.charge - min_point_2.charge
    if abs(charge_diff) != 1:
        logger.warning(f"The charge difference is {charge_diff} and not 1. Please "
                       f"ensure you understand the implications.")

    ccd_init = CcdInit(first_relaxed_point=min_point_1,
                       second_relaxed_point=min_point_2,
                       vbm=args.unitcell.vbm,
                       cbm=args.unitcell.cbm,
                       supercell_vbm=args.p_state.vbm_info.energy,
                       supercell_cbm=args.p_state.cbm_info.energy)

    if not path.exists():
        path.mkdir(parents=True)

    vise_yaml = (path / "vise.yaml")
    if not vise_yaml.exists():
        vise_yaml.write_text("""task: defect
user_incar_settings:
  NSW: 1""")

    ccd_init.to_json_file(json_file)
    logger.info(ccd_init)


def make_ccd_dirs(args: Namespace):
    os.chdir(args.calc_dir)
    s1 = args.ccd_init.first_relaxed_point.structure
    s2 = args.ccd_init.second_relaxed_point.structure
    s1_to_s2 = s1.interpolate(s2, nimages=args.first_to_second_div_ratios)
    s2_to_s1 = s2.interpolate(s1, nimages=args.second_to_first_div_ratios)
    rp1 = args.ccd_init.first_relaxed_point
    rp2 = args.ccd_init.second_relaxed_point

    Q_diff = args.ccd_init.Q

    for ratios, structures, rp_i, rp_f \
            in zip([args.first_to_second_div_ratios,
                    args.second_to_first_div_ratios],
                   [s1_to_s2, s2_to_s1],
                   [rp1, rp2],
                   [rp2, rp1]):
        dirpath = Path(f"q_{rp_i.charge}")
        dirpath.mkdir(exist_ok=True)

        spec_file = dirpath / "potential_curve_spec.json"
        if not spec_file.exists():
            spec = PotentialCurveSpec(
                rp_i.charge, rp_i.correction_energy, rp_f.charge, Q_diff)
            spec.to_json_file(spec_file)

        dQs = [args.ccd_init.Q * r for r in ratios]
        for ratio, structure, dQ in zip(ratios, structures, dQs):
            _make_ccd_dir(rp_i.charge, ratio, structure, dQ)
            if ratio == 0.0:
                # Make sure POSCAR is exactly the same as the relaxed structure.
                structure.to(filename=dirpath / "POSCAR")


def _make_ccd_dir(charge: int, ratio: float, structure: Structure, dQ: float):
    dir_ = Path(f"q_{charge}") / f"disp_{ratio}"
    try:
        dir_.mkdir(parents=True)
        logger.info(f"Directory {dir_} was created.")

        structure.to(filename=str(dir_ / "POSCAR"))
        (dir_ / "prior_info.yaml").write_text(
            yaml.dump({"charge": charge}), None)
        single_point_spec = SinglePointSpec(Q=dQ, disp_ratio=ratio)
        single_point_spec.to_json_file(dir_ / "single_point_spec.json")

    except FileExistsError:
        logger.info(f"Directory {dir_} exists, so skip it.")


def make_ccd_corrections(args):

    def _inner(dir_: Path):
        single_point_spec: SinglePointSpec \
            = loadfn(dir_ / "single_point_spec.json")

        calc_results = get_calc_results(dir_, False)
        minus_charge_diff = (args.potential_curve_spec.charge
                             - args.potential_curve_spec.counter_charge)
        effective_charge = single_point_spec.disp_ratio * minus_charge_diff

        ccd_correction = make_efnv_correction(
            effective_charge,
            calc_results,
            args.no_disp_calc_results,
            args.unitcell.effective_ionic_diele_const,
            args.no_disp_defect_entry.defect_center)

        ccd_correction.to_json_file(dir_ / "ccd_correction.json")

        plotter = SitePotentialMplPlotter.from_efnv_corr(
            title="ccd_correction", efnv_correction=ccd_correction)
        plotter.construct_plot()
        plotter.plt.savefig(fname=dir_ / "ccd_correction.pdf")
        plotter.plt.clf()

    parse_dirs(args.dirs, _inner, True, output_filename="ccd_correction.json")


def make_single_points(args: Namespace):
    def _inner(dir_: Path):
        calc_results = get_calc_results(dir_, False)
        if args.parse_ccd_correction:
            ccd_correction = loadfn(dir_ / "ccd_correction.json")
        else:
            ccd_correction = NoCcdCorrection()

        band_edge_states: BandEdgeStates = loadfn(dir_ / "band_edge_states.json")
        band_edge_orbital_infos: BandEdgeOrbitalInfos = loadfn(dir_ / "band_edge_orbital_infos.json")

        spec: SinglePointSpec = loadfn(dir_ / "single_point_spec.json")

        localized_orbitals, valence_bands, conduction_bands \
            = _get_band_edge_info(band_edge_orbital_infos, band_edge_states)

        single_point = SinglePoint(
            spec=spec,
            energy=calc_results.energy,
            ccd_correction_energy=ccd_correction.correction_energy,
            magnetization=calc_results.magnetization,
            localized_orbitals=localized_orbitals,
            valence_bands=valence_bands,
            conduction_bands=conduction_bands,
            is_shallow=band_edge_states.is_shallow)

        single_point.to_json_file(dir_ / "single_point.json")

    parse_dirs(args.dirs, _inner, verbose=True)


def _fit_curve(
        potential_curve: PotentialCurve,
        fitting_curve: Optional[Union[Type[FittingFunc], FittingFuncType]]):
    if fitting_curve:
        if isinstance(fitting_curve, FittingFuncType):
            fitting_curve = fitting_curve.cls
        potential_curve.set_fitting_curve(fitting_curve)
        return potential_curve.fitting_func
    return None


def make_potential_curve(args: Namespace):
    def _inner(dir_: Path) -> SinglePoint:
        return loadfn(dir_ / "single_point.json")

    single_points = SinglePoints(parse_dirs(args.dirs, _inner, verbose=True))
    raw_curve = PotentialCurve(args.potential_curve_spec,
                               single_points,
                               CurveTransform(0.0, False))
    curve_transform = make_curve_transform(raw_curve)
    potential_curve = PotentialCurve(args.potential_curve_spec,
                                     single_points,
                                     curve_transform)
    _fit_curve(potential_curve, args.fitting_func)
    print(potential_curve)
    potential_curve.to_json_file("potential_curve.json")
    ax = plt.gca()
    potential_curve.add_plot(ax)
    ax.set_xlabel("Q (amu$^{1/2}$ Å)")
    ax.set_ylabel("dE (eV)")
    ax.legend()
    plt.tight_layout()
    plt.savefig("potential_curve.pdf")


def make_ccd(args: Namespace):
    ccd = MakeCcd(ground_pot_curve=args.ground_pot_curve,
                  excited_pot_curve=args.excited_pot_curve,
                  vbm=args.ccd_init.vbm,
                  cbm=args.ccd_init.cbm,
                  name=args.ccd_init.name).ccd
    print(ccd)
    ccd.to_json_file()


def plot_ccd(args: Namespace):
    plotter = CcdPlotter(args.ccd,
                         plt,
                         ground_q_range=args.ground_q_range,
                         excited_q_range=args.excited_q_range)
    plotter.construct_plot()
    plt.savefig(args.fig_name)
    plt.show()


def plot_eigenvalues(args: Namespace):
    disp_ratios, orb_infos = [], []

    for d in args.dirs:
        try:
            orb_info = loadfn(Path(d) / "band_edge_orbital_infos.json")
        except FileNotFoundError:
            logger.info(f"band_edge_orbital_infos.json does not exist in {d}")
            continue
        try:
            single_point = loadfn(Path(d) / "single_point.json")
        except FileNotFoundError:
            logger.info(f"single_point_info.json does not exist in {d}")
            continue
        orb_infos.append(orb_info)
        disp_ratios.append(single_point.disp_ratio)

    vbm, cbm = args.ccd_init.supercell_vbm, args.ccd_init.supercell_cbm
    eigval_plotter = EigenvalPlotter(orb_infos, disp_ratios, vbm, cbm,
                                     y_range=args.y_range)
    eigval_plotter.construct_plot()
    eigval_plotter.plt.savefig("eigenvalues.pdf")
    eigval_plotter.plt.show()


def make_wswq_dirs(args: Namespace):
    for dir_ in args.dirs:
        _make_wswq_dir(dir_, args.ccd_init)


def _make_wswq_dir(dir_, ccd_init: CcdInit):
    wswq_dir = (dir_ / "wswq")
    if wswq_dir.exists():
        logger.info(f"Directory {wswq_dir} exists, so skip creating it.")
        return

    charge = PriorInfo.load_yaml(dir_ / "prior_info.yaml").charge
    original_dir = Path(ccd_init.relaxed_point_from_charge(charge).parsed_dir)

    wswq_dir.mkdir()
    logger.info(f"Directory {wswq_dir} was created.")

    for f_name in ["KPOINTS", "POSCAR", "POTCAR"]:
        FileLink((dir_/f_name).absolute()).transfer(wswq_dir)

    incar = ViseIncar.from_file(dir_/"INCAR")
    incar.update({"ALGO": "None", "LWSWQ": True, "NELM": 1, "LWAVE": False})
    incar.pop("LORBIT", None)
    incar.write_file(Path(wswq_dir/"INCAR"))

    os.symlink((dir_/"WAVECAR").absolute(), (wswq_dir/"WAVECAR.qqq"))

    os.symlink((original_dir/"WAVECAR").absolute(), (wswq_dir/"WAVECAR"))


def main_make_e_p_matrix_element(args: Namespace):

    charge = args.potential_curve.charge
    if charge != 0:
        raise ValueError("The e-p coupling must be evaluated with neutral defect.")

    dQs, wswqs = [], []

    band_edge_orbital_infos = None
    for d in args.dirs:
        wswq_file = d / "wswq" / "WSWQ"
        WSWQ = _read_WSWQ(wswq_file)  # keys: (spin, kpoint) and (initial, final)
        spin_kpt_index = (spin_to_idx(args.spin, True), 1)  # k-point index is 1
        band_index = tuple(sorted((args.band_edge_index, args.defect_band_index)))

        single_point: SinglePoint = loadfn(d / "single_point.json")
        dQs.append(single_point.Q)
        wswqs.append(WSWQ[spin_kpt_index][band_index])

        if single_point.disp_ratio == 0.0:
            band_edge_orbital_infos = loadfn(d / "band_edge_orbital_infos.json")

    if band_edge_orbital_infos is None:
        raise ValueError("At least one of the dirs must have disp_ratio of 0.0")

    e_p_matrix_elem = make_e_p_matrix_element(
        charge=charge,
        base_band_edge_orbital_infos=band_edge_orbital_infos,
        band_edge_index=args.band_edge_index,
        defect_band_index=args.defect_band_index,
        spin=args.spin,
        dQs=dQs,
        wswqs=wswqs)
    print(e_p_matrix_elem)
    e_p_matrix_elem.to_json_file()

    try:
        e_p_matrix_elem.plot(plt.gca())
        plt.xlabel("dQ (amu$^{1/2}$Å)")
        plt.savefig(f"e_p_matrix_element_{e_p_matrix_elem.index_info}.pdf")
        plt.show()
    except (TypeError, LinAlgError):
        logger.info("e-ph matrix element cannot be calculated.")


def make_total_squared_transition_moment(args: Namespace):
    calc = CalcTotalSquaredTransitionMoment(args.ccd.ground_curve,
                                            args.ccd.excited_curve,
                                            args.sommerfeld.Ts)
    total_moment = calc.harmonic()
    total_moment.to_json_file()
    plotter = PlottedTotalSquaredTransitionMoment(total_moment)
    plotter.plot(plt.gca())
    plt.savefig("total_squared_transition_moment.pdf")


def make_capture_rate(args: Namespace):
    ccd_init: CcdInit = args.ccd_init
    assert args.total_moment.Ts == args.sommerfeld.Ts

    # if args.degeneracy is None:
    #     f_min_info = ccd_init.relaxed_point_from_charge(args.ccd.ground_curve.charge)
    #     args.degeneracy = f_min_info.degeneracy_by_symmetry_reduction

    carrier = args.ccd.captured_carrier
    # spin_factor = 0.5 if i_min_info.is_spin_polarized else 1.0


    cap_rate = CaptureRate(args.total_moment.Ts,
                           ccd_init.volume,
                           args.sommerfeld.scaling(carrier, args.ccd.excited_curve.charge),
                           args.e_p_matrix_element.W_if_tilde,
                           args.total_moment.total_moments,
                           site_degeneracy=args.degeneracy)
    cap_rate.to_json_file()
    plotter = CaptureRatePlotter(cap_rate, plt)
    plotter.construct_plot()
    plt.savefig("capture_rate.pdf")
