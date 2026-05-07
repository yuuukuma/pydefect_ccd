# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
from dataclasses import dataclass

from monty.json import MSONable
from tabulate import tabulate
from vise.util.logger import get_logger
from vise.util.mix_in import ToJsonFileMixIn

from pydefect_ccd.relaxed_point import _joined_local_orbital_info, RelaxedPoint
from pydefect_ccd.util import get_dR, get_dQ

logger = get_logger(__name__)


@dataclass
class CcdInit(MSONable, ToJsonFileMixIn):
    """Initial data for a one-dimensional configuration coordinate diagram.

    The two relaxed points define the endpoints of the CCD path.  The first
    relaxed point is used as the reference structure for derived displacement
    quantities.

    Attributes:
        first_relaxed_point: First endpoint of the CCD path.
        second_relaxed_point: Second endpoint of the CCD path.
        vbm: Valence band maximum from the unit-cell calculation.
        cbm: Conduction band minimum from the unit-cell calculation.
        supercell_vbm: VBM from the perfect-supercell calculation.
        supercell_cbm: CBM from the perfect-supercell calculation.
    """
    first_relaxed_point: RelaxedPoint
    second_relaxed_point: RelaxedPoint
    vbm: float
    cbm: float
    supercell_vbm: float
    supercell_cbm: float

    @property
    def name(self) -> str:
        """Return a readable name combining both charge states."""
        return (f"{self.first_relaxed_point.full_name} "
                f"⇆ {self.second_relaxed_point.full_name}")

    @property
    def band_gap(self) -> float:
        """Return the unit-cell band gap in eV."""
        return self.cbm - self.vbm

    def relaxed_point_from_charge(self, charge: int) -> RelaxedPoint:
        """Return the relaxed point matching the requested charge."""
        for rp in (self.first_relaxed_point, self.second_relaxed_point):
            if rp.charge == charge:
                return rp
        else:
            raise ValueError(f"Charge {charge} does not exist.")

    @property
    def volume(self) -> float:
        """ Volume of the supercell structure in Å^3. """
        assert (self.first_relaxed_point.structure.volume
                == self.second_relaxed_point.structure.volume)
        return self.first_relaxed_point.structure.volume

    @property
    def Q(self):
        """ Unit: amu^{1/2} Å. """
        return abs(get_dQ(self.first_relaxed_point.structure,
                          self.second_relaxed_point.structure))

    @property
    def R(self):
        """ Unit: Å."""
        return abs(get_dR(self.first_relaxed_point.structure,
                          self.second_relaxed_point.structure))

    @property
    def modal_mass(self):
        """ Unit: amu^{1/2}."""
        return (self.Q / self.R) ** 2

    def __str__(self):
        """Return a report of CCD inputs, relaxed points, and band edges."""
        result = [f"name: {self.name}"]
        table = [["vbm", self.vbm, "supercell vbm", self.supercell_vbm],
                 ["cbm", self.cbm, "supercell cbm", self.supercell_cbm],
                 ["volume (Å^3)", self.volume],
                 ["Q (amu^0.5 Å)", self.Q],
                 ["R (Å)", self.R],
                 ["M (amu)", self.modal_mass]]
        result.append(tabulate(table, tablefmt="plain", floatfmt=".3f"))

        result.append("-" * 60)

        headers = ["q", "ini symm.", "final symm.", "energy",
                   "correction", "corrected energy", "magnetization",
                   "localized orbitals", "ZPL"]
        table = []

        last_energy = None

        for relaxed_point in [self.first_relaxed_point, self.second_relaxed_point]:
            localized_state_idxs = []
            for s, spin in zip(relaxed_point.localized_orbitals, ["up", "down"]):
                for ss in s:
                    localized_state_idxs.append(f"{spin}-{ss.band_idx}")
            table.append(
                [relaxed_point.charge, relaxed_point.initial_site_symmetry,
                 relaxed_point.final_site_symmetry, relaxed_point.formation_energy,
                 relaxed_point.correction_energy, relaxed_point.corrected_energy,
                 relaxed_point.magnetization,
                 _joined_local_orbital_info(relaxed_point.localized_orbitals)])
            if last_energy:
                table[-1].append(last_energy - relaxed_point.corrected_energy)
            last_energy = relaxed_point.corrected_energy

        result.append(
            tabulate(table, tablefmt="plain", headers=headers, floatfmt=".3f",
                     stralign="center"))

        result.append("-" * 60)

        for relaxed_point in [self.first_relaxed_point, self.second_relaxed_point]:
            result.append(f"- q={relaxed_point.charge}")
            for vb, spin in zip(relaxed_point.valence_bands, ["up", "down"]):
                result.append(f"-- VBM spin-{spin}")
                for v in vb:
                    result.append(str(v))
            result.append(f"")
            for cb, spin in zip(relaxed_point.conduction_bands, ["up", "down"]):
                result.append(f"-- CBM spin-{spin}")
                for c in cb:
                    result.append(str(c))
            result.append(f"-- parse dir: {relaxed_point.parsed_dir}")
            result.append(f"")

        return "\n".join(result)

""" 
TODO: 
plot

1. add defect position
2. consider how to handle the small difference of origin.
remove ave_static_diele_const
"""
