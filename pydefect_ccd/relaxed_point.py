# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from monty.json import MSONable
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from pymatgen.core import Structure, IStructure
from vise.util.structure_symmetrizer import num_sym_op


@dataclass
class NearEdgeState(MSONable):
    """Band-edge state near the VBM or CBM.

    This stores the minimum information needed to identify and report a
    valence- or conduction-band state that is relevant to a defect.

    Attributes:
        band_index: Band index starting from 1.
        kpt_coord: Fractional k-point coordinate.
        kpt_index: K-point index starting from 1.
        eigenvalue: Single-particle energy in the absolute scale in eV.
        occupation: Electron occupation number.
    """
    band_index: int
    kpt_coord: List[float]
    kpt_index: int
    eigenvalue: float
    occupation: float

    def __str__(self):
        """Return a compact one-line description of the band-edge state."""
        k_coord = " ".join([f"{x:.2f}" for x in self.kpt_coord])
        k_info = [f"index : {self.kpt_index}", f"coord: {k_coord}"]
        x = [f"band index: {self.band_index}",
             f"kpt info: ({', '.join(k_info)})",
             f"eigenvalue: {self.eigenvalue:.2f}",
             f"occupation: {self.occupation:.2f}"]
        return ", ".join(x)


def _joined_local_orbital_info(localized_orbitals: List[List[LocalizedOrbital]]
                               ) -> str:
    """Return localized orbital labels grouped by spin channel."""
    lo_str = []
    for lo_by_spin, spin in zip(localized_orbitals, ["up", "down"]):
        for lo_by_band in lo_by_spin:
            occupation = f"{lo_by_band.occupation:.1f}"
            lo_str.append(f"{spin}-{lo_by_band.band_idx}({occupation})")
    return ", ".join(lo_str)


@dataclass(kw_only=True)
class OrbitalInfoMixIn(MSONable):
    """Electronic-state information attached to one fixed structure.

    The spin-resolved lists are organized as ``[spin][state]``.  The first
    spin channel is spin-up and the second is spin-down.  Consumers in this
    package generally expect these lists to be populated with empty inner lists
    when no state is present, rather than left as ``None``.

    Attributes:
        magnetization: Total magnetization of the calculation.
        localized_orbitals: Localized defect orbitals at each spin channel.
        valence_bands: Valence-band edge states at each spin channel.
        conduction_bands: Conduction-band edge states at each spin channel.
    """
    magnetization: float
    localized_orbitals: Optional[List[List[LocalizedOrbital]]] = field(default=None)
    valence_bands: Optional[List[List[NearEdgeState]]] = field(default=None)
    conduction_bands: Optional[List[List[NearEdgeState]]] = field(default=None)


@dataclass
class RelaxedPoint(OrbitalInfoMixIn):
    """Relaxed defect structure and electronic data for one charge state.

    A ``RelaxedPoint`` represents a local minimum on the configuration
    coordinate diagram.  It stores the relaxed geometry, formation energy,
    finite-size correction, symmetry labels, and band-edge information needed
    to build CCD paths between two charge states.

    The spin-resolved band and orbital lists inherited from
    ``OrbitalInfoMixIn`` are ordered as ``[spin][state]``.  The first spin
    channel is spin-up and the second is spin-down.

    Attributes:
        name: Defect name without the charge suffix.
        charge: Defect charge state.
        formation_energy: Formation energy at ``Ef = VBM`` before adding the
            correction energy.
        correction_energy: Energy correction estimated by methods such as
            eFNV.  This is added to ``formation_energy`` in ``corrected_energy``.
        structure: Relaxed atomic structure for this charge state.
        initial_site_symmetry: Site symmetry before structural relaxation.
        final_site_symmetry: Site symmetry after structural relaxation.
        parsed_dir: Directory containing the parsed calculation results.
            This should be an absolute path when possible.
        magnetization: Total magnetization of the calculation.
        localized_orbitals: List of localized orbitals at each spin channel.
        valence_bands: Valence-band edge states at each spin channel.
        conduction_bands: Conduction-band edge states at each spin channel.
    """
    name: str
    charge: int
    formation_energy: float  # formation energy
    correction_energy: float
    structure: Structure | IStructure
    initial_site_symmetry: str
    final_site_symmetry: str
    parsed_dir: str

    @property
    def full_name(self) -> str:
        """Return defect name with charge suffix."""
        return f"{self.name}_{self.charge}"

    @property
    def corrected_energy(self) -> float:
        """Return formation energy plus correction energy."""
        return self.formation_energy + self.correction_energy

    @property
    def dir_path(self) -> Path:
        """Return the parsed calculation directory as a path."""
        return Path(self.parsed_dir)

    @property
    def is_spin_polarized(self):
        """Return whether the magnetization indicates spin polarization."""
        return abs(self.magnetization) > 0.95

    @property
    def related_band_indices(self) -> set:
        """Return all band indices referenced by stored band/orbital data."""
        result = set()

        def add(bands):
            """Add band indices from one spin-resolved collection."""
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
        """Return the degeneracy from initial-to-final symmetry reduction."""
        initial_num_sym_op = num_sym_op[self.initial_site_symmetry]
        final_num_sym_op = num_sym_op[self.final_site_symmetry]
        return initial_num_sym_op / final_num_sym_op
