# -*- coding: utf-8 -*-
#  Copyright (c) 2024 Kumagai group.
from dataclasses import dataclass
from pathlib import Path
from typing import List

from monty.json import MSONable
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from pymatgen.core import Structure
from vise.util.structure_symmetrizer import num_sym_op


@dataclass
class NearEdgeState(MSONable):
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
class RelaxedPoint(MSONable):
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
        valence_bands: List of valence band edge states at each spin channel.
        conduction_bands: List of conduction band edge states at each spin channel.
        parsed_dir (str): Directory where the calculation results of this
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
    valence_bands: List[List[NearEdgeState]]  # by spin
    conduction_bands: List[List[NearEdgeState]]  # by spin
    parsed_dir: str

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
