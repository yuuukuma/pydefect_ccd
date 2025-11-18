# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from abc import ABC, abstractmethod
from copy import deepcopy
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Union

import numpy as np
from monty.json import MSONable
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from pymatgen.electronic_structure.core import Spin
from scipy import interpolate
from scipy.optimize import curve_fit
from tabulate import tabulate
from vise.util.logger import get_logger
from vise.util.matplotlib import float_to_int_formatter
from vise.util.mix_in import ToJsonFileMixIn

from pydefect_ccd.enum import Carrier
from pydefect_ccd.relaxed_point import NearEdgeState, \
    _joined_local_orbital_info, \
    PointMixIn
from pydefect_ccd.util import spin_to_idx

logger = get_logger(__name__)


@dataclass(frozen=True)
class SinglePointSpec(MSONable, ToJsonFileMixIn):
    """Specification of a single point calculation in ."""
    dQ: float
    disp_ratio: float


@dataclass
class SinglePoint(PointMixIn, ToJsonFileMixIn):
    spec: SinglePointSpec
    ccd_correction_energy: float  # CCD correction at a fixed structure
    is_shallow: bool = None
    # This needs to be here to make table_for_plot

    @property
    def dQ(self) -> float: return self.spec.dQ

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

    def localized_orbital(self,
                          spin: Spin, band_index: int) -> LocalizedOrbital:

        for lo in self.localized_orbitals[spin_to_idx(spin)]:
            if getattr(lo, "band_idx", None) == band_index:
                return lo
        raise ValueError

    @property
    def ccd_corrected_energy(self):
        return self.energy + self.ccd_correction_energy

    @property
    def table_for_plot(self):
        tables = {"corr. energy": self.ccd_corrected_energy,
                  "is shallow?": self.is_shallow,
                  "localized orb": _joined_local_orbital_info(self.localized_orbitals)}
        return list(tables.keys()), tables.values()

    def __str__(self):
        headers, tabulate_data = self.table_for_plot
        return tabulate([tabulate_data], tablefmt="plain", floatfmt=".3f",
                        headers=headers)


@dataclass
class PotentialCurveSpec(MSONable, ToJsonFileMixIn):
    charge: int
    correction_energy: float  # at relaxed structure, e.g., eFNV correction
    counter_charge: int  # charge state to which the structure shifts
    Q_diff: float


@dataclass
class Curve(ABC):
    @abstractmethod
    def __call__(self, x: Union[float, np.array]) -> Union[float, np.array]:
        pass


@dataclass
class QuadraticCurve(MSONable, Curve):
    omega: float
    Q0: float
    disp_ratio_range: Tuple[float, float]

    def __call__(self, x: Union[float, np.array]) -> Union[float, np.array]:
        return self.omega * x**2 + self.Q0

    def __str__(self):
        return (f"QuadraticCurve: omega={self.omega:.5f}, Q0={self.Q0:.5f}, "
                f"disp_ratio_range=({self.disp_ratio_range[0]:.3f}, "
                f"{self.disp_ratio_range[1]:.3f})")


@dataclass
class PotentialCurve(MSONable, ToJsonFileMixIn):
    spec: PotentialCurveSpec
    single_points: List[SinglePoint] = field(default_factory=list)
    shifted_energy: float = 0.0  # to set the zero of energy
    fitted_curve: Curve = None  # TODO: revert

    @property
    def charge(self) -> int: return self.spec.charge

    @property
    def correction_energy(self) -> float: return self.spec.correction_energy

    @property
    def counter_charge(self) -> int: return self.spec.counter_charge

    @property
    def Q_diff(self) -> float: return self.spec.Q_diff

    @property
    def dQs(self) -> List[float]:
        return [sp.dQ for sp in self.single_points]

    def __post_init__(self):
        self.single_points = list(sorted(self.single_points, key=lambda x: x.dQ))

    def single_point_from_disp(self, disp_ratio: float):
        for sp in self.single_points:
            if np.isclose(sp.disp_ratio, disp_ratio):
                return sp
        raise ValueError(f"No single point found for disp_ratio={disp_ratio}")

    def dQs_and_energies(self, disp_ratio_range: Tuple[float, float] = None):
        """

        Energy for the ccd is the sum of bare DFT energy at a fixed structure,
        FNV correction energy, shifted energy, and CCD correction energy.

        Returns:

        """
        dQs, energies = [], []
        for single_point in self.single_points:
            if disp_ratio_range and not (disp_ratio_range[0]
                                         <= single_point.disp_ratio
                                         <= disp_ratio_range[1]):
                continue
            dQs.append(single_point.dQ)
            energy = (single_point.energy + single_point.ccd_correction_energy
                      + self.spec.correction_energy + self.shifted_energy)
            energies.append(energy)
        return dQs, energies

    def add_quadratic_curve(self,
                            disp_ratio_range: Tuple[float, float] = None,
                            fixed_Q0: bool = True):
        dQs, energies = self.dQs_and_energies(disp_ratio_range)
        if len(dQs) < 3:
            raise ValueError("The number of Q points must be >= 3.")

        if fixed_Q0:
            Q0 = self.single_point_from_disp(0.0).dQ
        else:
            Q0 = None
        omega, Q0 = calc_omega_and_Q0(dQs, energies, Q0)
        self.fitted_curve = QuadraticCurve(omega, Q0, disp_ratio_range)

    def add_plot(self,
                 ax,
                 color: str,
                 q_range: Optional[List[float]] = None,
                 spline_fit: bool = True):
        dQs, energies = self.dQs_and_energies(q_range)
        ax.scatter(dQs, energies, marker='o', color=color)
        try:
            if spline_fit:
                label = f"q={self.charge}"
#                carriers = "+".join([str(c) for c in self.carriers])
#                 if carriers:
#                     label += "+" + carriers
                x, y = spline3(dQs, energies, 100, q_range)
                ax.plot(x, y, label=label, color=color)
        except TypeError as e:
            print(f"{self.charge}: {e}")
            pass

        if self.fitted_curve:
            q_max, q_min = np.max(self.dQs), np.min(self.dQs)
            q_L = q_max - q_min
            qs = np.linspace(q_min - 0.1 * q_L, q_max + 0.1 * q_L, 1000)
            ax.plot(qs, self.fitted_curve(qs))

    @property
    def table_for_plot(self):
        tables = {"charge": self.charge,
                  "corr. energy": self.correction_energy,
                  "counter charge": self.counter_charge,
                  "Q diff": self.Q_diff,
                  "shifted energy": self.shifted_energy}
        return list(tables.keys()), tables.values()

    def __str__(self):
        headers, tabulate_data = self.table_for_plot
        return tabulate([tabulate_data], tablefmt="plain", floatfmt=".3f",
                        headers=headers)


def calc_omega_and_Q0(Qs: List[float],
                      energies: List[float],
                      Q0: float) -> Tuple[float, float]:
    def f(Q, omega, Q0, dE):
        return 0.5 * omega**2 * (Q - Q0)**2 + dE

    # set bounds to restrict Q0 to the given Q0 value
    bounds = (-np.inf, np.inf) if Q0 is None else \
        ([-np.inf, Q0 - 1e-10, -np.inf], [np.inf, Q0, np.inf])
    popt, _ = curve_fit(f, Qs, energies, bounds=bounds)
    return popt[0], popt[2]


def dQ_revert(pot_curve_result: PotentialCurve) -> PotentialCurve:
    new_single_points = []
    for single_point in pot_curve_result.single_points:
        new_dis_ratio = 1.0 - single_point.disp_ratio
        new_dQ = new_dis_ratio * pot_curve_result.spec.Q_diff
        new_spec = SinglePointSpec(dQ=new_dQ, disp_ratio=new_dis_ratio)
        new_result = deepcopy(single_point)
        new_result.spec = new_spec
        new_single_points.append(new_result)

    return PotentialCurve(pot_curve_result.spec,
                          new_single_points,
                          pot_curve_result.shifted_energy)


    # def __str__(self):
    #     try:
    #         omega = f"{self.omega():.5f}"
    #     except (ValueError, RuntimeError):
    #         omega = "N.A."
    #     result = [f"name: {self.name}", f"charge: {self.charge}",
    #               f"omega: {omega}"]

#         if self.carriers:
#             carriers = " ".join([str(carrier) for carrier in self.carriers])
#             result.append(f"carriers: {carriers}")

#         table_data = [x.table_for_plot for x in self.single_points]
#         result.append(tabulate(table_data, tablefmt="plain", floatfmt=".3f",
#                                headers=_imag_headers + ["omega"]))
#         return "\n".join(result)


@dataclass
class Ccd(MSONable, ToJsonFileMixIn):
    name: str
    ground_curve: PotentialCurve
    excited_curve: PotentialCurve
    # potential_curve_results: List[PotentialCurveResult]

    @property
    def captured_carrier(self) -> Carrier:
        charge_diff = self.excited_curve.charge - self.ground_curve.charge
        carrier_charge = - charge_diff

        return Carrier.from_carrier_charge(carrier_charge)

    def __str__(self):
        result = [f"name: {self.name}"]
        result.append("-"*50)
        result.append(self.excited_curve.__str__())
        result.append(self.ground_curve.__str__())
        return "\n".join(result)


def spline3(xs, ys, num_points, xrange=None):
    """Find the B-spline representation with 3 degree of the spline.

    Args:
        xs (array_like):
        ys (array_like):
        num_points (int): Number of interpolated points including end points.

    Returns:
        Tuple of (x_spline, y_spline) as numpy arrays.
    """
    xs, ys = np.asarray(xs), np.asarray(ys)

    if num_points < 2:
        raise ValueError("num_points must be >= 2")

    # desired degree is 3 (cubic)
    max_k = max(1, min(3, len(xs) - 1))
    if len(xs) < 2:
        raise ValueError("At least two points are required.")

    # prepare spline; use s=0 to force interpolation
    tck = interpolate.splprep([xs, ys], k=max_k, s=0)[0]

    if xrange is not None:
        x_dist = float(np.max(xs) - np.min(xs))
        if x_dist == 0.0:
            raise ValueError("xs must not be all equal")
        _min = (xrange[0] - np.min(xs)) / x_dist
        _max = (xrange[1] - np.min(xs)) / x_dist
        # clamp to [0, 1]
        _min = max(0.0, min(1.0, _min))
        _max = max(0.0, min(1.0, _max))
    else:
        _min, _max = 0.0, 1.0

    u = np.linspace(_min, _max, num=num_points, endpoint=True)
    spline = interpolate.splev(u, tck)
    return np.asarray(spline[0]), np.asarray(spline[1])


class CcdPlotter:
    def __init__(self,
                 ccd: Ccd,
                 plt,
                 title: str = None,
                 set_energy_zero: bool = True,
                 quadratic_fit: bool = True,
                 spline_fit: bool = True,
                 ground_q_range: list = None,
                 excited_q_range: list = None):
        self._title = title or ""
        self._ccd = ccd
        self._set_energy_zero = set_energy_zero
        self._quadratic_fit = quadratic_fit
        self._spline_fit = spline_fit
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
        if self._q_range:
            ax.set_xlim(self._q_range[0], self._q_range[1])

        self._ccd.ground_curve.add_plot(ax, "red", self._ground_q_range,
                                        self._quadratic_fit, self._spline_fit)
        self._ccd.excited_curve.add_plot(ax, "blue", self._excited_q_range,
                                         self._quadratic_fit, self._spline_fit)

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

