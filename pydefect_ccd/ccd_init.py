# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass
from typing import List

from monty.json import MSONable
from pymatgen.analysis.defects.ccd import get_dQ
from tabulate import tabulate
from vise.util.logger import get_logger
from vise.util.mix_in import ToJsonFileMixIn

from pydefect_ccd.enum import Carrier
from pydefect_ccd.relaxed_point import _joined_local_orbitals, RelaxedPoint
from pydefect_ccd.util import get_dR

logger = get_logger(__name__)


@dataclass
class CcdInit(MSONable, ToJsonFileMixIn):
    """ Initial information related to the configuration coordination diagram.

    Attributes:
        relaxed_points (List[RelaxedPoint]): List of two relaxed defects.
            The charge state difference must be 1.
        vbm (float): valence band maximum in the unitcell calculation.
        cbm (float): conduction band minimum in the unitcell calculation.
        supercell_vbm (float): vbm in the perfect supercell calculation.
        supercell_cbm (float): cbm in the perfect supercell calculation.
        ave_static_diele_const (float): Average of the static dielectric
        ave_electron_mass (float): Average of the electron effective mass
        ave_hole_mass (float): Average of the hole effective mass
    """
    relaxed_points: List[RelaxedPoint]
    vbm: float
    cbm: float
    supercell_vbm: float
    supercell_cbm: float
    supercell_volume: float
    ave_static_diele_const: float
    ave_electron_mass: float = None
    ave_hole_mass: float = None

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

