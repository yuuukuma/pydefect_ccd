# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass, field
from itertools import permutations
from math import isclose
from typing import List, Optional, Dict, Tuple

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from monty.json import MSONable
from nonrad.ccd import get_omega_from_PES
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from pymatgen.electronic_structure.core import Spin
from scipy import interpolate
from tabulate import tabulate
from vise.util.logger import get_logger
from vise.util.matplotlib import float_to_int_formatter
from vise.util.mix_in import ToJsonFileMixIn

from pydefect_ccd.enum import Carrier
from pydefect_ccd.relaxed_point import NearEdgeState, _joined_local_orbitals
from pydefect_ccd.util import spin_to_idx

logger = get_logger(__name__)


@dataclass(frozen=True)
class SinglePointSpec(MSONable, ToJsonFileMixIn):
    """Specification of a single point calculation in ."""
    dQ: float
    disp_ratio: float


@dataclass
class SinglePointResult(MSONable, ToJsonFileMixIn):
    relative_energy: float  # wrt the relaxed energy.
    correction_energy: float
    magnetization: float
    localized_orbitals: List[List[LocalizedOrbital]] \
        = field(default_factory=list) # [spin][bands]
    # [spin][bands]
    valence_bands: List[List[NearEdgeState]] = field(default_factory=list)
    conduction_bands: List[List[NearEdgeState]] = field(default_factory=list)
    is_shallow: bool = None
    # This needs to be here to make table_for_plot

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
    def corrected_relative_energy(self):
        return self.relative_energy + self.correction_energy

    @property
    def table_for_plot(self):
        tables = {
                  "corr. energy": self.corrected_energy,
                  "relative energy": self.corrected_relative_energy,
                  "used for fitting?": self.used_for_fitting,
                  "is shallow?": self.is_shallow,
                  "localized orb": _joined_local_orbitals(self.localized_orbitals)}
        return list(tables.keys()), tables.values()

    def __str__(self):
        headers, tabulate_data = self.table_for_plot
        return tabulate([tabulate_data], tablefmt="plain", floatfmt=".3f",
                        headers=headers)


@dataclass
class PotentialCurve(MSONable, ToJsonFileMixIn):
    charge: int
    relaxed_energy: float  # bottom energy
    carriers: List[Carrier] = field(default_factory=list)
    single_points: Dict[SinglePointSpec, SinglePointResult] \
        = field(default_factory=dict)
    charge_to: Optional[int] = None

    def __post_init__(self):
        has_disp0 = any(isclose(spec.disp_ratio, 0.0)
                        for spec in self.single_points.keys())
        assert has_disp0, "single_points must contain disp_ratio=0.0 result"

        self.single_points = dict(sorted(self.single_points.items(),
                                         key=lambda item: item[1].spec.dQ))

    def dQs_and_energies(self, disp_ratio_range: Tuple[float, float] = None):
        dQs, energies = [], []
        for spec, result in self.single_points.items():
            if disp_ratio_range and not (disp_ratio_range[0] <= spec.disp_ratio <= disp_ratio_range[1]):
                continue
            dQs.append(spec.dQ)
            energies.append(result.corrected_relative_energy)
        return dQs, energies

    def omega(self,
              ax: Axes = None,
              plot_q_range: Optional[List[float]] = None,
              color: str = None):

        dQs, energies = self.dQs_and_energies()

        if len(dQs) < 3:
            raise ValueError("The number of Q points is not sufficient for "
                             "calculating omega.")

        q = np.linspace(plot_q_range[0], plot_q_range[1], 1000) \
            if plot_q_range else None

        return get_omega_from_PES(np.array(dQs), np.array(energies), ax=ax, q=q)

    # def dQ_reverted_single_ccd(self) -> "PotentialCurve":
    #     result = deepcopy(self)
    #     for point_info in self.single_points:
    #         if point_info.disp_ratio == 0.0:
    #             continue
    #         disp1_dQ = point_info.dQ / point_info.disp_ratio
    #         break
    #
    #     for i in result.single_points:
    #         i.dQ = (1.0 - i.disp_ratio) * disp1_dQ
    #     result.single_points.sort(key=lambda x: x.dQ)
    #     return result

    def add_plot(self,
                 ax,
                 color: str,
                 q_range: Optional[List[float]] = None,
                 quadratic_fit: bool = True,
                 spline_fit: bool = True):
        dQs, energies = self.dQs_and_energies()
        ax.scatter(dQs, energies, marker='o', color=color)
        try:
            if spline_fit:
                x, y = spline3(dQs, energies, 100, q_range)
                ax.plot(x, y, label=self.name, color=color)
        except TypeError as e:
            print(f"{self.name}: {e}")
            pass

        if quadratic_fit:
            try:
                self.omega(ax, plot_q_range=q_range)
            except (ValueError, TypeError, RuntimeError) as e:
                print(f"{self.name}: {e}")
                pass

    # def __str__(self):
    #     try:
    #         omega = f"{self.omega():.5f}"
    #     except (ValueError, RuntimeError):
    #         omega = "N.A."
    #     result = [f"name: {self.name}", f"charge: {self.charge}",
    #               f"omega: {omega}"]
    #
    #     if self.carriers:
    #         carriers = " ".join([str(carrier) for carrier in self.carriers])
    #         result.append(f"carriers: {carriers}")
    #
    #     table_data = [x.table_for_plot for x in self.single_points]
    #     result.append(tabulate(table_data, tablefmt="plain", floatfmt=".3f",
    #                            headers=_imag_headers + ["omega"]))
    #     return "\n".join(result)


def captured_carrier(initial: PotentialCurve, final: PotentialCurve):
    carrier_diff = set(initial.carriers) - set(final.carriers)
    if len(carrier_diff) != 1:
        raise CarrierDiffError
    return carrier_diff.pop()


@dataclass
class Ccd(MSONable, ToJsonFileMixIn):
    name: str
    potential_curves: List[PotentialCurve]

    # def single_ccd(self, name, charge) -> PotentialCurve:
    #     names = []
    #     for i in self.potential_curves:
    #         if i.id == single_ccd_id:
    #             return i
    #         names.append(f"{i.id_.__str__()}")
    #     raise ValueError(f"Choose state name from {'  '.join(names)}")

    def initial_and_final_ccd_from_captured_carrier(
            self, carrier: Carrier) -> (PotentialCurve, PotentialCurve):

        for i, j in permutations(self.potential_curves, 2):
            try:
                if captured_carrier(i, j) == carrier:
                    return i, j
            except CarrierDiffError:
                continue
        raise CarrierDiffError

    def __str__(self):
        result = [f"name: {self.name}"]

        for imag in self.potential_curves:
            result.append("-"*50)
            result.append(imag.__str__())
        return "\n".join(result)


def spline3(xs, ys, num_points, xrange=None):
    """Find the B-spline representation with 3 degree of the spline.

    Args:
        xs (array_like):
        ys (array_like):
        num_points (int): Number of interpolated points including end points.

    Returns:
        Tuple of
    """
    #   tck : tuple
    #         (t,c,k) _default_single_ccd_for_e_p_coupling tuple containing the vector of knots, the B-spline
    #         coefficients, and the degree of the spline.
    tck = interpolate.splprep([xs, ys], k=5)[0]

    if xrange:
        x_dist = max(xs) - min(xs)
        _min = (xrange[0] - min(xs)) / x_dist
        _max = (xrange[1] - min(xs)) / x_dist
    else:
        _min, _max = 0.0, 1.0

    u = np.linspace(_min, _max, num=num_points, endpoint=True)
    spline = interpolate.splev(u, tck)
    return spline[0], spline[1]


class CcdPlotter:
    def __init__(self, ccd: Ccd,
                 title: str = None,
                 set_energy_zero: bool = True,
                 quadratic_fit: bool = True,
                 spline_fit: bool = True,
                 q_range: list = None):
        self._title = title or ""
        self._ccd = ccd
        self._set_energy_zero = set_energy_zero
        self._quadratic_fit = quadratic_fit
        self._spline_fit = spline_fit
        self.plt = plt
        self._q_range = q_range

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
        for imag_infos, color in zip(self._ccd.potential_curves, ["red", "blue", "green"]):
            imag_infos.add_plot(ax, color, self._q_range, self._quadratic_fit, self._spline_fit)

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


class CarrierDiffError(Exception):
    pass