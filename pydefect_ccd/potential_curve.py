# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
import inspect
from copy import deepcopy
from dataclasses import dataclass
from typing import List, Tuple, Optional

import numpy as np
from monty.json import MSONable
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from pymatgen.electronic_structure.core import Spin
from scipy.optimize import curve_fit
from tabulate import tabulate
from vise.util.mix_in import ToJsonFileMixIn

# from pydefect_ccd.fitting_curve import FittingCurve
from pydefect_ccd.relaxed_point import OrbitalInfoMixIn, NearEdgeState, \
    _joined_local_orbital_info
from pydefect_ccd.util import spin_to_idx


@dataclass(frozen=True)
class SinglePointSpec(MSONable, ToJsonFileMixIn):
    """Specification of a single point calculation in ."""
    Q: float
    disp_ratio: float


@dataclass
class SinglePoint(OrbitalInfoMixIn, ToJsonFileMixIn):
    energy: float  # bare total energy
    spec: SinglePointSpec
    ccd_correction_energy: float  # CCD correction at a fixed structure
    used_for_fitting: bool
    is_shallow: bool = None
    # This needs to be here to make table_for_plot

    @property
    def Q(self) -> float: return self.spec.Q

    @property
    def disp_ratio(self) -> float: return self.spec.disp_ratio

    def near_edge_state(self, spin: Spin, band_index: int) -> NearEdgeState:
        spin_idx = spin_to_idx(spin)
        vb = self.valence_bands[spin_idx]
        cb = self.conduction_bands[spin_idx]
        lo = self.localized_orbitals[spin_idx]

        for band in vb + lo + cb:
            if getattr(band, "band_index", None) == band_index:
                return band

        raise ValueError

    def localized_orbital(self, spin: Spin, band_index: int) -> LocalizedOrbital:

        for lo in self.localized_orbitals[spin_to_idx(spin)]:
            if getattr(lo, "band_idx", None) == band_index:
                return lo
        raise ValueError

    @property
    def ccd_corrected_energy(self):
        return self.energy + self.ccd_correction_energy

    @property
    def table_headers(self):
        return ["disp ratio", "corrected energy", "is shallow?", "localized orb"]

    def table_values(self,
                     correction_energy: float = 0.0,
                     shifted_energy: float = 0.0):
        localized_orbitals = _joined_local_orbital_info(self.localized_orbitals)
        energy = self.energy + correction_energy + self.ccd_correction_energy + shifted_energy
        return [self.disp_ratio, energy, self.is_shallow, localized_orbitals]

    def __str__(self):
        return tabulate([self.table_values()], tablefmt="plain", floatfmt=".3f",
                        headers=self.table_headers)


@dataclass
class SinglePoints(MSONable):
    single_points: List[SinglePoint]
    # fitting_curve: FittingCurve

    def __iter__(self):
        return iter(self.single_points)

    def __len__(self):
        return len(self.single_points)

    @property
    def Qs(self) -> List[float]:
        return [sp.Q for sp in self]

    @property
    def corrected_energies(self) -> List[float]:
        return [sp.ccd_corrected_energy for sp in self]

    def single_point_from_disp(self, disp_ratio: float):
        for sp in self:
            if np.isclose(sp.disp_ratio, disp_ratio):
                return sp
        raise ValueError(f"No single point found for disp_ratio={disp_ratio}")

    @property
    def lowest_energy_single_point(self) -> SinglePoint:
        return min(self.single_points, key=lambda sp: sp.ccd_corrected_energy)

    def verify_Q0_has_the_lowest_energy(self):
        if not np.isclose(self.lowest_energy_single_point.Q, 0.0):
            dQ = self.lowest_energy_single_point.Q
            raise ValueError(
                f"The single point with the lowest energy has dQ={dQ:.3f}, "
                f"which is not 0. ")

    def verify_num_Q(self, f):
        sig = inspect.signature(f)
        if n_params := len(sig.parameters) > len(self):
            raise ValueError(f"The number of Q points must be >= {n_params}.")


@dataclass
class PotentialCurveSpec(MSONable, ToJsonFileMixIn):
    charge: int
    correction_energy: float  # at relaxed structure, e.g., eFNV correction
    counter_charge: int  # charge state to which the structure shifts
    Q_diff: float


@dataclass
class PotentialCurve(MSONable, ToJsonFileMixIn):
    spec: PotentialCurveSpec
    original_single_points: SinglePoints
    shifted_energy: float = 0.0
    shifted_Q: float = 0.0
    revert: bool = False

    @property
    def charge(self) -> int: return self.spec.charge

    @property
    def correction_energy(self) -> float: return self.spec.correction_energy

    @property
    def counter_charge(self) -> int: return self.spec.counter_charge

    @property
    def Q_diff(self) -> float: return self.spec.Q_diff

    def __post_init__(self):
        self.original_single_points = list(sorted(self.original_single_points, key=lambda x: x.Q))

    @property
    def lowest_energy_single_point(self) -> SinglePoint:
        return self.original_single_points.lowest_energy_single_point

    @property
    def lowest_energy(self) -> float:
        return self.lowest_energy_single_point.ccd_corrected_energy + \
                self.spec.correction_energy + self.shifted_energy

#     def dQs_and_energies(self, disp_ratio_range: Tuple[float, float] = None):
#         """
#
#         Energy for the ccd is the sum of bare DFT energy at a fixed structure,
#         FNV correction energy, shifted energy, and CCD correction energy.
#
#         Returns:
#
#         """
#         dQs, energies = [], []
#         for single_point in self.original_single_points:
#             if disp_ratio_range and not (disp_ratio_range[0]
#                                          <= single_point.disp_ratio
#                                          <= disp_ratio_range[1]):
#                 continue
#             dQs.append(single_point.dQ)
#             energy = (single_point.energy + single_point.ccd_correction_energy
#                       + self.spec.correction_energy + self.shifted_energy)
#             energies.append(energy)
#         return dQs, energies
#
#     def dQ_revert(self, fixed_Q0: bool = True) -> "PotentialCurve":
#         new_single_points = []
#         for single_point in self.original_single_points:
#             new_dis_ratio = 1.0 - single_point.disp_ratio
#             new_dQ = new_dis_ratio * self.spec.Q_diff
#             new_spec = SinglePointSpec(Q=new_dQ, disp_ratio=new_dis_ratio)
#             new_result = deepcopy(single_point)
#             new_result.spec = new_spec
#             new_single_points.append(new_result)
#
#         result = PotentialCurve(self.spec,
#                                 new_single_points,
#                                 self.shifted_energy)
#         # curve_name = self.fitting_curve.__class__.__name__
#         # new_disp_ratio_range = [1 - self.fitting_curve.disp_ratio_range[1],
#         #                         1 - self.fitting_curve.disp_ratio_range[0]]
#
# #         result.add_curve(fitting_type=FittingCurveType.from_string(curve_name),
# #                          disp_ratio_range=new_disp_ratio_range,
#
#         return result

    def add_plot(self,
                 ax,
                 color: str,
                 q_range: Optional[List[float]] = None):
        label = f"q={self.charge}"
        dQs, energies = self.dQs_and_energies(q_range)
        ax.scatter(dQs, energies, marker='o', color=color, label=label)

        if self.fitting_curve:
            q_range = q_range or [min(dQs), max(dQs)]
            self.fitting_curve.add_plot(ax, q_range, color=color)

    @property
    def table_for_plot(self):
        tables = {"charge": self.charge,
                  "lowest energy": self.lowest_energy,
                  "corr. energy": self.correction_energy,
                  "counter charge": self.counter_charge,
                  "Q diff": self.Q_diff,
                  "shifted energy": self.shifted_energy}
        return list(tables.keys()), tables.values()

    def __str__(self):
        # TODO: improve
        headers, tabulate_data = self.table_for_plot
        table = tabulate([tabulate_data], tablefmt="plain", floatfmt=".3f",
                        headers=headers)
        fc = str(self.fitting_curve) if self.fitting_curve else "fitted curve is N.A."
        d = [sp.table_values(self.correction_energy, self.shifted_energy)
             for sp in self.original_single_points]
        table_2 = tabulate(d, tablefmt="plain", floatfmt=".3f",
                           headers=self.original_single_points[0].table_headers)
        return table + "\n" + fc + "\n" + table_2



