# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass
from pathlib import Path
from typing import List

from monty.json import MSONable
from numpy import sqrt, sum
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from pymatgen.analysis.defects.ccd import get_dQ
from pymatgen.core import Structure
from tabulate import tabulate
from vise.util.logger import get_logger
from vise.util.mix_in import ToJsonFileMixIn
from vise.util.structure_symmetrizer import num_sym_op

from dephon.enum import Carrier

logger = get_logger(__name__)


def get_dR(ground: Structure, excited: Structure) -> float:
    """Summation of atomic displacement distance

    Args:
        ground (Structure): Reference structure
        excited (Structure): Target structure

    A constant offset should be removed, for example, by aligning the farthest
    atom.

    Returns:
        The Summed atomic displacement distance in float
    """
    return sqrt(sum([x.distance(y) ** 2 for x, y in zip(ground, excited)]))


@dataclass
class BandEdgeState(MSONable):
    band_index: int  # begin from 1.
    kpt_coord: List[float]
    kpt_weight: float
    kpt_index: int  # begin from 1.
    eigenvalue: float
    occupation: float

    def __str__(self):
        k_coord = " ".join([f"{x:.2f}" for x in self.kpt_coord])
        k = [f"index : {self.kpt_index}",
             f"coord: {k_coord}",
             f"weight: {self.kpt_weight}"]
        x = [f"band index: {self.band_index}",
             f"kpt info: ({', '.join(k)})",
             f"eigenvalue: {self.eigenvalue:.2f}",
             f"occupation: {self.occupation:.2f}"]
        return ", ".join(x)


def _joined_local_orbitals(localized_orbitals) -> str:
    lo_str = []
    for lo_by_spin, spin in zip(localized_orbitals, ["up", "down"]):
        for lo_by_band in lo_by_spin:
            occupation = f"{lo_by_band.occupation:.1f}"
            lo_str.append(f"{spin}-{lo_by_band.band_idx}({occupation})")
    return ", ".join(lo_str)


@dataclass
class RelaxedPointInfo(MSONable):
    """Information at the relaxed structure for a given charge state

    Attributes:
        charge: The charge state.
        structure: The atomic configuration.
        energy: Formation energy at Ef=VBM and chemical potentials
            being standard states.
        correction_energy: Correction energy estimated e.g. by eFNV method.
        magnetization: Magnetization
        localized_orbitals: List of localized orbitals at each spin channel.
            [Spin up orbitals, Spin down orbitals]
        initial_site_symmetry (str): Site symmetry before relaxing the defect
        final_site_symmetry (str): Site symmetry after relaxing the defect
        vbm: valence band maximum in the unitcell calculation
        cbm: conduction band minimum in the unitcell calculation
        parse_dir (str): Directory where the calculation results of this
            minimum point are stored. This should be an absolute path.
    """
    name: str
    charge: int
    structure: Structure
    energy: float
    correction_energy: float
    magnetization: float
    localized_orbitals: List[List[LocalizedOrbital]]
    initial_site_symmetry: str
    final_site_symmetry: str
    valence_bands: List[List[BandEdgeState]]  # by spin
    conduction_bands: List[List[BandEdgeState]]  # by spin
    parsed_dir: str # absolute dir

    @property
    def full_name(self) -> str:
        return f"{self.name}_{self.charge}"

    @property
    def corrected_energy(self) -> float:
        return self.energy + self.correction_energy

    @property
    def degeneracy_by_symmetry_reduction(self) -> float:
        initial_num_sym_op = num_sym_op[self.initial_site_symmetry]
        final_num_sym_op = num_sym_op[self.final_site_symmetry]
        return initial_num_sym_op / final_num_sym_op

    @property
    def dir_path(self) -> Path:
        return Path(self.parsed_dir)

    @property
    def is_spin_polarized(self):
        return abs(self.magnetization) > 0.95

    @property
    def relevant_band_indices(self) -> set:
        result = set()

        def add(bands):
            for i in bands:
                for j in i:
                    try:
                        result.add(j.band_index)
                    except AttributeError:
                        result.add(j.band_idx)

        add(self.valence_bands)
        add(self.localized_orbitals)
        add(self.conduction_bands)
        return result


@dataclass
class DephonInit(MSONable, ToJsonFileMixIn):
    """ Information related to configuration coordination diagram.

    Attributes:
        relaxed_points (List[RelaxedPointInfo]): List of two relaxed defects.
            the charge state difference should be 1.
        vbm (float): valence band maximum in the unitcell calculation.
        cbm (float): conduction band minimum in the unitcell calculation.
        supercell_vbm (float): vbm in the perfect supercell calculation.
        supercell_cbm (float): cbm in the perfect supercell calculation.

    """
    relaxed_points: List[RelaxedPointInfo]
    vbm: float
    cbm: float
    supercell_vbm: float
    supercell_cbm: float
    supercell_volume: float
    ave_electron_mass: float
    ave_hole_mass: float
    ave_static_diele_const: float

    def __post_init__(self):
        assert len(self.relaxed_points) == 2

    def effective_mass(self, carrier: Carrier):
        if carrier is Carrier.e:
            return self.ave_electron_mass
        return self.ave_hole_mass

    @property
    def name(self):
        return (f"{self.relaxed_points[0].full_name} "
                f"⇆ {self.relaxed_points[1].full_name}")

    @property
    def band_gap(self):
        return self.cbm - self.vbm

    def relaxed_point_info_from_charge(self, charge: int):
        for rp in self.relaxed_points:
            if rp.charge == charge:
                return rp
        else:
            raise ValueError(f"Charge {charge} does not exist.")

    @property
    def volume(self):
        assert (self.relaxed_points[0].structure.volume
                == self.relaxed_points[1].structure.volume)
        return self.relaxed_points[0].structure.volume

    @property
    def dQ(self):
        return abs(get_dQ(self.relaxed_points[0].structure,
                          self.relaxed_points[1].structure))

    @property
    def dR(self):
        return abs(get_dR(self.relaxed_points[0].structure,
                          self.relaxed_points[1].structure))

    @property
    def modal_mass(self):
        return (self.dQ / self.dR) ** 2

    def __str__(self):
        result = [f"name: {self.name}"]
        table = [["vbm", self.vbm, "supercell vbm", self.supercell_vbm],
                 ["cbm", self.cbm, "supercell cbm", self.supercell_cbm],
                 ["volume (Å^3)", self.supercell_volume],
                 ["dQ (amu^0.5 Å)", self.dQ],
                 ["dR (Å)", self.dR],
                 ["M (amu)", self.modal_mass],
                 ["electron mass (m0)", self.ave_electron_mass],
                 ["hole mass (m0)", self.ave_hole_mass],
                 ["static diele", self.ave_static_diele_const]]
        result.append(tabulate(table, tablefmt="plain", floatfmt=".3f"))

        result.append("-" * 60)

        headers = ["q", "ini symm", "final symm", "energy",
                   "correction", "corrected energy", "magnetization",
                   "localized orbitals", "ZPL"]
        table = []

        last_energy = None

        for min_info in self.relaxed_points:

            localized_state_idxs = []
            for s, spin in zip(min_info.localized_orbitals, ["up", "down"]):
                for ss in s:
                    localized_state_idxs.append(f"{spin}-{ss.band_idx}")
            table.append(
                [min_info.charge, min_info.initial_site_symmetry,
                 min_info.final_site_symmetry, min_info.energy,
                 min_info.correction_energy, min_info.corrected_energy,
                 min_info.magnetization,
                 _joined_local_orbitals(min_info.localized_orbitals)])
            if last_energy:
                table[-1].append(last_energy - min_info.corrected_energy)
            last_energy = min_info.corrected_energy

        result.append(
            tabulate(table, tablefmt="plain", headers=headers, floatfmt=".3f",
                     stralign="center"))

        result.append("-" * 60)

        for min_info in self.relaxed_points:
            result.append(f"- q={min_info.charge}")
            for vb, spin in zip(min_info.valence_bands, ["up", "down"]):
                result.append(f"-- VBM spin-{spin}")
                for v in vb:
                    result.append(str(v))
            result.append(f"")
            for cb, spin in zip(min_info.conduction_bands, ["up", "down"]):
                result.append(f"-- CBM spin-{spin}")
                for c in cb:
                    result.append(str(c))
            result.append(f"-- parse dir: {min_info.parsed_dir}")
            result.append(f"")

        return "\n".join(result)

""" 
TODO: 
plot

1. add defect position
2. consider how to handle the small difference of origin.
remove ave_static_diele_const
"""

