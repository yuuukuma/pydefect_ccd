# -*- coding: utf-8 -*-
#  Copyright (c) 2024 Kumagai group.
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from monty.json import MSONable
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from pymatgen.core import Structure, IStructure
from vise.util.structure_symmetrizer import num_sym_op


@dataclass
class NearEdgeState(MSONable):
    """Information about band-edge states near the VBM or CBM.

    Note: In this code, we assume the Gamma point only.

    Attributes:
        band_index: Index of the band (begin from 1).
        kpt_coord: Coordinates of the k-point in fractional coordinates.
        kpt_index: Index of the k-point (begin from 1).
        eigenvalue: Eigenvalue of the state (in eV).
        occupation: Occupation number of the state.
    """

    band_index: int
    kpt_coord: List[float]
    kpt_index: int
    eigenvalue: float
    occupation: float

    def __str__(self):
        k_coord = " ".join([f"{x:.2f}" for x in self.kpt_coord])
        k_info = [f"index : {self.kpt_index}", f"coord: {k_coord}"]
        x = [f"band index: {self.band_index}",
             f"kpt info: ({', '.join(k_info)})",
             f"eigenvalue: {self.eigenvalue:.2f}",
             f"occupation: {self.occupation:.2f}"]
        return ", ".join(x)


def _joined_local_orbital_info(localized_orbitals: List[List[LocalizedOrbital]]
                               ) -> str:
    lo_str = []
    for lo_by_spin, spin in zip(localized_orbitals, ["up", "down"]):
        for lo_by_band in lo_by_spin:
            occupation = f"{lo_by_band.occupation:.1f}"
            lo_str.append(f"{spin}-{lo_by_band.band_idx}({occupation})")
    return ", ".join(lo_str)


@dataclass(kw_only=True)
class PointMixIn(MSONable):
    """Mix-in class that stores information for a fixed structural point.

    Indices for orbitals and bands are organized as [spin][band].
    Attributes:
        energy: Bare energy obtained from DFT calculations.
    """
    energy: float
    magnetization: float
    localized_orbitals: Optional[List[List[LocalizedOrbital]]] = field(default=None)
    valence_bands: Optional[List[List[NearEdgeState]]] = field(default=None)
    conduction_bands: Optional[List[List[NearEdgeState]]] = field(default=None)


@dataclass
class RelaxedPoint(PointMixIn):
    """Information at the relaxed structure for a given charge state

    Attributes:
        name: The defect name.
        charge: The charge state.
        correction_energy: Correction energy estimated e.g. by eFNV method.
        structure: The atomic configuration.
        energy: Formation energy at Ef=VBM and chemical potentials
            being standard states.
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
    correction_energy: float
    structure: Structure | IStructure
    initial_site_symmetry: str
    final_site_symmetry: str
    parsed_dir: str

    @property
    def full_name(self) -> str:
        return f"{self.name}_{self.charge}"

    @property
    def corrected_energy(self) -> float:
        return self.energy + self.correction_energy

    @property
    def dir_path(self) -> Path:
        return Path(self.parsed_dir)

    @property
    def is_spin_polarized(self):
        return abs(self.magnetization) > 0.95

    @property
    def related_band_indices(self) -> set:
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

    @property
    def degeneracy_by_symmetry_reduction(self) -> float:
        initial_num_sym_op = num_sym_op[self.initial_site_symmetry]
        final_num_sym_op = num_sym_op[self.final_site_symmetry]
        return initial_num_sym_op / final_num_sym_op

