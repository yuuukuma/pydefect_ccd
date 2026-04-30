# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass
from typing import List, Tuple

from monty.json import MSONable
from pydefect.corrections.abstract_correction import Correction
from vise.util.logger import get_logger
from vise.util.matplotlib import float_to_int_formatter
from vise.util.mix_in import ToJsonFileMixIn

from pydefect_ccd.fitting_curve import intersections
from pydefect_ccd.local_enum import Carrier
from pydefect_ccd.potential_curve import PotentialCurve

logger = get_logger(__name__)


@dataclass
class Ccd(MSONable, ToJsonFileMixIn):
    name: str
    ground_curve: PotentialCurve
    excited_curve: PotentialCurve

    @property
    def Q_diff(self) -> float:
        return self.excited_curve.Q_diff

    @property
    def dE(self) -> float:
        return self.excited_curve.lowest_energy - self.ground_curve.lowest_energy

    @property
    def captured_carrier(self) -> Carrier:
        carrier_charge = self.ground_curve.charge - self.excited_curve.charge
        return Carrier.from_carrier_charge(carrier_charge)

    def intersections(self, min_Q_mul=-2, max_Q_mul=3,
                      num_grids=2001) -> List[Tuple[float, float]]:
        ground = self.ground_curve.fitting_curve
        excited = self.excited_curve.fitting_curve
        if excited is None:
            raise ValueError("Set excited fitting curve.")
        if ground is None:
            raise ValueError("Set ground fitting curve.")

        min_Q, max_Q = self.Q_diff * min_Q_mul, self.Q_diff * max_Q_mul
        return intersections(ground, excited, [min_Q, max_Q], num_grids)

    @property
    def crossing_points(self) -> List[Tuple[float, float]]:
        return self.intersections()

    def __str__(self):
        result = [f"name: {self.name}",
                  "excited" + "-" * 50, str(self.excited_curve),
                  "ground" + "-" * 50, str(self.ground_curve),
                  "intersections" + "-" * 45]

        try:
            for q, energy in self.intersections():

                result.append(f"Q={q:.3f}, Energy={energy:.3f}")
        except (ValueError, TypeError):
            result.append("No intersections")
            pass

        return "\n".join(result)


class CcdPlotter:
    def __init__(self,
                 ccd: Ccd,
                 plt,
                 title: str = None,
                 ground_q_range: list = None,
                 excited_q_range: list = None):
        self._title = title or ""
        self._ccd = ccd
        self._ground_q_range = ground_q_range
        self._excited_q_range = excited_q_range

        self.plt = plt

    def construct_plot(self):
        self._add_ccd()
        self._set_title()
        self._set_formatter()
        self._set_labels()
        self.plt.tight_layout()

    def _add_ccd(self):
        ax = self.plt.gca()
        self._ccd.ground_curve.add_plot(ax, "red", self._ground_q_range)
        self._ccd.excited_curve.add_plot(ax, "blue", self._excited_q_range)

    def _set_labels(self):
        ax = self.plt.gca()
        ax.set_xlabel("Q (amu$^{1/2}$ Å)")
        ax.set_ylabel("Energy (eV)")
        ax.legend()

    def _set_title(self):
        self.plt.gca().set_title(self._title)

    def _set_formatter(self):
        self.plt.gca().xaxis.set_major_formatter(float_to_int_formatter)
        self.plt.gca().yaxis.set_major_formatter(float_to_int_formatter)



class NoCcdCorrection(Correction):
    @property
    def correction_energy(self) -> float:
        return 0.0


# @dataclass
# class CcdCorrection(ExtendedFnvCorrection):
#     """
#     charge: Defect charge
#     point_charge_correction: Correction energy for point charge interactions
#     defect_region_radius (float):
#         Maximum radius of a sphere touching to the lattice plane, used
#         for defining the outside region of the defect.
#     sites: Site potentials
#     defect_coords: Position of defect site in fractional coordinates
#
#     Add units of length and potential
#     """
#     charge: float
#     point_charge_correction: float
#     defect_region_radius: float
#     sites: List["PotentialSite"]
#     defect_coords: Tuple[float, float, float]
#
#     def __str__(self):
#         d = [["charge", self.charge],
#              ["pc term", self.point_charge_correction],
#              ["alignment term", self.alignment_correction],
#              ["correction energy", self.correction_energy]]
#         return tabulate(d, tablefmt='psql')
#
#     @property
#     def average_potential_diff(self):
#         return np.mean([s.diff_pot for s in self.sites
#                         if s.distance > self.defect_region_radius])
#
#     @property
#     def alignment_correction(self) -> float:
#         return - self.average_potential_diff * self.charge
#
#     @property
#     def correction_energy(self) -> float:
#         return self.point_charge_correction + self.alignment_correction
#
#     @property
#     def correction_dict(self):
#         return {"pc term": self.point_charge_correction,
#                 "alignment term": self.alignment_correction}