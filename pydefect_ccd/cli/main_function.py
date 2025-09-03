# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import os
from argparse import Namespace
from pathlib import Path
from typing import List, Tuple

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
from pymatgen.electronic_structure.core import Spin
from vise.input_set.incar import ViseIncar
from vise.input_set.prior_info import PriorInfo
from vise.util.file_transfer import FileLink
from vise.util.logger import get_logger

from pydefect_ccd.capture_rate import calc_phonon_overlaps, CaptureRate
from pydefect_ccd.ccd_init import CcdInit
from pydefect_ccd.config_coord import SinglePointResult, CcdPlotter, \
    PotentialCurve
from pydefect_ccd.ele_phon_coupling import EPMatrixElement
from pydefect_ccd.enum import CorrectionType
from pydefect_ccd.make_config_coord import MakeCcd
from pydefect_ccd.make_e_p_matrix_element import MakeEPMatrixElement
from pydefect_ccd.plot_eigenvalues import CcdEigenvaluePlotter
from pydefect_ccd.relaxed_point import NearEdgeState, RelaxedPoint
from pydefect_ccd.util import spin_to_idx

logger = get_logger(__name__)


def _make_near_edge_states(band_edge_orbital_infos: BandEdgeOrbitalInfos,
                           spin: Spin,
                           edge_energy: float,
                           threshold: float = 0.1):
    result = []
    info_by_spin = band_edge_orbital_infos.orbital_infos[spin_to_idx(spin)]
    for k_idx, (orb_info_by_kpt, k_coords, k_weight) in \
        enumerate(zip(info_by_spin,
                      band_edge_orbital_infos.kpt_coords,
                      band_edge_orbital_infos.kpt_weights)):
        k_idx_from_1 = k_idx + 1
        for rel_idx, info_by_band in enumerate(orb_info_by_kpt):
            e_from_band_edge = abs(info_by_band.energy - edge_energy)
            if e_from_band_edge > threshold:
                continue
            band_idx = rel_idx + band_edge_orbital_infos.lowest_band_index + 1
            result.append(NearEdgeState(band_idx,
                                        k_coords,
                                        k_weight,
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
        energy=energy_info.defect_energy.formation_energy,
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

    path = Path(f"cc/{min_point_1.full_name}⇆{min_point_2.full_name}")
    json_file = path / "ccd_init.json"

    if json_file.exists():
        logger.info(f"{json_file} exists. Remove it first to recreate it.")
        return

    if abs(min_point_1.charge - min_point_2.charge) != 1:
        logger.warning("The charge difference is not 1. "
                       "Please ensure you understand the implications.")

    volume = min_point_1.structure.volume
    if args.effective_mass:
        concentration = args.effective_mass.concentrations[0]
        ave_hole_mass = args.effective_mass.average_mass("p", concentration)
        ave_electron_mass = args.effective_mass.average_mass("n", concentration)
    else:
        ave_hole_mass, ave_electron_mass = None, None

    ccd_init = CcdInit(relaxed_points=[min_point_1, min_point_2],
                          vbm=args.unitcell.vbm,
                          cbm=args.unitcell.cbm,
                          supercell_volume=volume,
                          supercell_vbm=args.p_state.vbm_info.energy,
                          supercell_cbm=args.p_state.cbm_info.energy,
                          ave_hole_mass=ave_hole_mass,
                          ave_electron_mass=ave_electron_mass,
                          ave_static_diele_const=args.unitcell.ave_ele_diele)

    if path.exists() is False:
        path.mkdir(parents=True)

    vise_yaml = (path / "vise.yaml")
    if vise_yaml.exists() is False:
        vise_yaml.write_text("""task: defect
user_incar_settings:
  NSW: 1""")

    ccd_init.to_json_file(json_file)
    logger.info(ccd_init)


def make_ccd_dirs(args: Namespace):
    os.chdir(args.calc_dir)
    d_init = args.ccd_init
    s1 = d_init.relaxed_points[0].structure
    s2 = d_init.relaxed_points[1].structure
    s1_to_s2 = s1.interpolate(s2, nimages=args.first_to_second_div_ratios)
    s2_to_s1 = s2.interpolate(s1, nimages=args.second_to_first_div_ratios)
    initial_charges = [d_init.relaxed_points[0].charge,
                       d_init.relaxed_points[1].charge]
    final_charges = [d_init.relaxed_points[1].charge,
                     d_init.relaxed_points[0].charge]
    correction_energies = [d_init.relaxed_points[0].correction_energy,
                           d_init.relaxed_points[1].correction_energy]

    for ratios, structures, initial_charge, final_charge, corr in \
            zip([args.first_to_second_div_ratios,
                 args.second_to_first_div_ratios],
                [s1_to_s2, s2_to_s1],
                initial_charges, final_charges, correction_energies):
        dQs = [args.ccd_init.dQ * r for r in ratios]
        name = f"from_{initial_charge}_to_{final_charge}"

        for ratio, structure, dQ in zip(ratios, structures, dQs):
            _make_ccd_dir(initial_charge, name, ratio, structure, dQ, corr)

        single_ccd = Path(name) / "potential_curve.json"
        if single_ccd.exists() is False:
            id_ = CcdId(name)
            PotentialCurve(id_, charge=initial_charge).to_json_file(single_ccd)


def _make_ccd_dir(charge, dirname, ratio, structure, dQ, correction):
    dir_ = Path(dirname) / f"disp_{ratio}"
    try:
        dir_.mkdir(parents=True)
        logger.info(f"Directory {dir_} was created.")

        structure.to(filename=str(dir_ / "POSCAR"))
        (dir_ / "prior_info.yaml").write_text(
            yaml.dump({"charge": charge}), None)
        single_point_info = SinglePointResult(dQ=dQ, disp_ratio=ratio)
        single_point_info.to_json_file(dir_ / "single_point_info.json")

        correction = DephonCorrection({CorrectionType.extended_FNV: correction})
        correction.to_yaml_file(dir_ / "ccd_correction.yaml")

    except FileExistsError:
        logger.info(f"Directory {dir_} exists, so skip it.")


def make_ccd_correction(args):
    file_name = "ccd_correction.json"

    def _inner(dir_: Path):
        single_point_info: SinglePointResult \
            = loadfn(dir_ / "single_point_info.json")

        if single_point_info.disp_ratio == 0.0:
            return

        calc_results = get_calc_results(dir_, False)
        minus_charge_diff = args.potential_curve.charge - args.to_charge
        effective_charge = single_point_info.disp_ratio * minus_charge_diff

        kumagai2023 = make_efnv_correction(effective_charge,
                                           calc_results,
                                           args.no_disp_calc_results,
                                           args.unitcell.effective_ionic_diele_const,
                                           args.no_disp_defect_entry.defect_center)
        kumagai2023.to_json_file(dir_ / file_name)

        title = args.potential_curve.name
        plotter = SitePotentialMplPlotter.from_efnv_corr(
            title=title, efnv_correction=kumagai2023)
        plotter.construct_plot()
        plotter.plt.savefig(fname=dir_ / "ccd_correction.pdf")
        plotter.plt.clf()

        kumagai2023 = loadfn(dir_ / file_name)

        correction = DephonCorrection.from_yaml(dir_ / "dephon_correction.yaml")
        correction.corrections[
            str(CorrectionType.kumagai2023)] = float(kumagai2023.correction_energy)
        correction.to_yaml_file(dir_ / "dephon_correction.yaml")

    parse_dirs(args.dirs, _inner, output_filename=file_name, verbose=True)


def update_single_point_infos(args: Namespace):
    def _inner(dir_: Path):
        calc_results = get_calc_results(dir_, False)
        band_edge_states: BandEdgeStates = loadfn(dir_ / "band_edge_states.json")
        band_edge_orbital_infos: BandEdgeOrbitalInfos = loadfn(dir_ / "band_edge_orbital_infos.json")

        correction = DephonCorrection.from_yaml(dir_ / "dephon_correction.yaml")

        localized_orbitals, valence_bands, conduction_bands = _get_band_edge_info(
            band_edge_orbital_infos, band_edge_states)

        sp_info: SinglePointResult = loadfn(dir_ / "single_point_info.json")
        sp_info.corrected_energy = calc_results.energy + correction.total_correction_energy
        sp_info.magnetization = calc_results.magnetization
        sp_info.localized_orbitals = localized_orbitals
        sp_info.valence_bands = valence_bands
        sp_info.conduction_bands = conduction_bands
        sp_info.is_shallow = band_edge_states.is_shallow
        sp_info.to_json_file(dir_ / "single_point_info.json")

    parse_dirs(args.dirs, _inner, verbose=True)


def add_point_infos_to_single_ccd(args: Namespace):
    def _inner(dir_: Path):
        return loadfn(dir_ / "single_point_info.json")

    single_ccd: PotentialCurve = loadfn("potential_curve.json")
    single_ccd.single_points = parse_dirs(args.dirs, _inner, verbose=True)
    single_ccd.to_json_file("potential_curve.json")


def make_ccd(args: Namespace):
    ccd = MakeCcd(args.ground_ccd, args.excited_ccd, args.dephon_init).ccd
    print(ccd)
    ccd.to_json_file()


def set_quadratic_fitting_q_range(args: Namespace):
    single_ccd: PotentialCurve = args.ccd.potential_curve(args.single_ccd_name)
    single_ccd.set_quadratic_fitting_range(args.q_range)
    print(args.ccd)
    args.ccd.to_json_file()


def plot_ccd(args: Namespace):
    plotter = CcdPlotter(args.ccd, q_range=args.q_range,
                         quadratic_fit=args.quadratic_fit,
                         spline_fit=args.spline_fit)
    plotter.construct_plot()
    plotter.plt.savefig(args.fig_name)
    plotter.plt.show()


def plot_eigenvalues(args: Namespace):
    disp_ratios, orb_infos = [], []

    for d in args.dirs:
        try:
            orb_info = loadfn(Path(d) / "band_edge_orbital_infos.json")
        except FileNotFoundError:
            logger.info(f"band_edge_orbital_infos.json does not exist in {d}")
            continue
        try:
            single_point_info = loadfn(Path(d) / "single_point_info.json")
        except FileNotFoundError:
            logger.info(f"single_point_info.json does not exist in {d}")
            continue
        orb_infos.append(orb_info)
        disp_ratios.append(single_point_info.disp_ratio)

    vbm, cbm = args.dephon_init.supercell_vbm, args.dephon_init.supercell_cbm
    eigval_plotter = CcdEigenvaluePlotter(orb_infos, disp_ratios, vbm, cbm,
                                          y_range=args.y_range)
    eigval_plotter.construct_plot()
    eigval_plotter.plt.savefig("dephon_eigenvalues.pdf")
    eigval_plotter.plt.show()


def make_wswq_dirs(args: Namespace):
    for dir_ in args.dirs:
        _make_wswq_dir(dir_, args.dephon_init)


def _make_wswq_dir(dir_, dephon_init: CcdInit):
    wswq_dir = (dir_ / "wswq")
    if wswq_dir.exists():
        logger.info(f"Directory {wswq_dir} exists, so skip creating it.")
        return

    charge = PriorInfo.load_yaml(dir_ / "prior_info.yaml").charge
    original_dir = Path(dephon_init.relaxed_point_info_from_charge(charge).parsed_dir)

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


def make_e_p_matrix_element(args: Namespace):
    dQ_wswq_pairs = []
    for d in args.dirs:
        wswq_file = d / "wswq" / "WSWQ"
        single_point_info = loadfn(d / "single_point_info.json")
        dQ_wswq_pairs.append((single_point_info.dQ, _read_WSWQ(wswq_file)))

    make_e_p_matrix_elem = MakeEPMatrixElement(
        base_disp_ratio=args.base_disp,
        single_ccd=args.single_ccd,
        band_edge_index=args.band_edge_index,
        defect_band_index=args.defect_band_index,
        kpoint_index=args.kpoint_index,
        spin=args.spin,
        dQ_wswq_pairs=dQ_wswq_pairs,
        energy_diff=args.energy_diff)
    e_p_matrix_elem = make_e_p_matrix_elem.make()
    print(e_p_matrix_elem)
    e_p_matrix_elem.to_json_file()

    try:
        e_p_matrix_elem.e_p_matrix_element(plt.gca())
        plt.xlabel("dQ (amu$^{1/2}$Å)")
        plt.savefig(f"e_p_matrix_element_{e_p_matrix_elem.index_info}.pdf")
        plt.show()
    except (TypeError, LinAlgError):
        logger.info("e-ph matrix element cannot be calculated.")


def make_capture_rate(args: Namespace):
    dephon_init: CcdInit = args.dephon_init
    e_p_matrix_elem: EPMatrixElement = args.e_p_matrix_elem
    if e_p_matrix_elem.e_p_matrix_element is None:
        raise ValueError

    f_ccd, i_ccd = args.ccd.initial_and_final_ccd_from_captured_carrier(args.captured_carrier)
    print(i_ccd)
    print(f_ccd)

    i_min_info = dephon_init.relaxed_point_info_from_charge(i_ccd.charge)
    f_min_info = dephon_init.relaxed_point_info_from_charge(f_ccd.charge)
    i_deg = i_min_info.degeneracy_by_symmetry_reduction
    f_deg = f_min_info.degeneracy_by_symmetry_reduction

    phonon_overlaps = calc_phonon_overlaps(i_ccd, f_ccd, args.temperatures)
    em = dephon_init.effective_mass(args.captured_carrier)
    velocities = thermal_velocity(np.array(args.temperatures), em)
    spin_factor = 0.5 if i_min_info.is_spin_polarized else 1.0

    cap_rate = CaptureRate(Wif=e_p_matrix_elem.e_p_matrix_element(),
                           summed_phonon_overlaps=phonon_overlaps,
                           velocities=velocities,
                           temperatures=args.temperatures,
                           site_degeneracy=f_deg / i_deg,
                           spin_selection_factor=spin_factor,
                           volume=dephon_init.volume)
    print(cap_rate)
    cap_rate.to_json_file()


def plot_capture_rate(args: Namespace):
    cap: CaptureRate = args.capture_rate
    plt.scatter(cap.temperatures, cap.capture_rate, marker='o')
    plt.gca().set_yscale("log")
    plt.savefig("capture_rate.pdf")
    plt.show()