# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
import inspect
from copy import deepcopy
from dataclasses import dataclass
from functools import cached_property
from typing import List, Optional, Type, Tuple

import numpy as np
from monty.json import MSONable
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from pymatgen.electronic_structure.core import Spin
from scipy.optimize import curve_fit
from tabulate import tabulate
from vise.util.mix_in import ToJsonFileMixIn

from pydefect_ccd.fitting_func import FittingFunc
from pydefect_ccd.relaxed_point import OrbitalInfoMixIn, NearEdgeState, \
    _joined_local_orbital_info
from pydefect_ccd.util import spin_to_idx


@dataclass(frozen=True)
class SinglePointSpec(MSONable, ToJsonFileMixIn):
    """Position of one single-point calculation along a CCD path."""
    Q: float
    disp_ratio: float

    def flip(self, Q_diff: float) -> "SinglePointSpec":
        """Return the same point measured from the opposite endpoint."""
        flipped_Q = Q_diff - self.Q
        flipped_disp_ratio = 1.0 - self.disp_ratio
        return SinglePointSpec(Q=flipped_Q, disp_ratio=flipped_disp_ratio)


@dataclass
class SinglePoint(OrbitalInfoMixIn, ToJsonFileMixIn):
    """Single-point calculation result at a fixed displaced structure."""
    energy: float  # bare total energy
    spec: SinglePointSpec
    ccd_correction_energy: float  # CCD correction at a fixed structure
    used_for_fitting: bool
    is_shallow: bool = None
    # This needs to be here to make table_for_plot

    @property
    def Q(self) -> float:
        """Return the configuration coordinate in amu^0.5 Angstrom."""
        return self.spec.Q

    @property
    def disp_ratio(self) -> float:
        """Return the interpolation ratio along the displacement path."""
        return self.spec.disp_ratio

    def near_edge_state(self, spin: Spin, band_index: int) -> NearEdgeState:
        """Return the near-edge state matching spin and band index."""
        spin_idx = spin_to_idx(spin)
        vb = self.valence_bands[spin_idx]
        cb = self.conduction_bands[spin_idx]
        lo = self.localized_orbitals[spin_idx]

        for band in vb + lo + cb:
            if getattr(band, "band_index", None) == band_index:
                return band

        raise ValueError

    def localized_orbital(self, spin: Spin, band_index: int) -> LocalizedOrbital:
        """Return the localized orbital matching spin and band index."""

        for lo in self.localized_orbitals[spin_to_idx(spin)]:
            if getattr(lo, "band_idx", None) == band_index:
                return lo
        raise ValueError

    @property
    def ccd_corrected_energy(self):
        """Return total energy after adding the CCD correction."""
        return self.energy + self.ccd_correction_energy

    @property
    def table_headers(self):
        """Return column labels for tabulating this point."""
        return ["disp ratio", "corrected energy", "is shallow?", "localized orb"]

    def table_values(self,
                     correction_energy: float = 0.0,
                     shifted_energy: float = 0.0):
        """Return row values with optional global correction and shift."""
        localized_orbitals = _joined_local_orbital_info(self.localized_orbitals)
        energy = self.energy + correction_energy + self.ccd_correction_energy + shifted_energy
        return [self.disp_ratio, energy, self.is_shallow, localized_orbitals]

    def __str__(self):
        """Return a one-row table for this single point."""
        return tabulate([self.table_values()], tablefmt="plain", floatfmt=".3f",
                        headers=self.table_headers)


@dataclass
class Shifter(MSONable, ToJsonFileMixIn):
    """Energy and coordinate transformation applied to single points."""
    shift_energy: float
    flip: bool


@dataclass
class SinglePoints(MSONable):
    """Container for single-point calculations on one potential curve."""
    single_points: List[SinglePoint]
    # fitting_curve: FittingCurve

    def __iter__(self):
        """Iterate over contained single points."""
        return iter(self.single_points)

    def __len__(self):
        """Return the number of single points."""
        return len(self.single_points)

    @property
    def Qs(self) -> List[float]:
        """Return all Q coordinates."""
        return [sp.Q for sp in self]

    @property
    def corrected_energies(self) -> List[float]:
        """Return CCD-corrected energies for all points."""
        return [sp.ccd_corrected_energy for sp in self]

    @property
    def Qs_and_energies(self):
        """Return paired Q coordinates and CCD-corrected energies."""
        return list(zip(self.Qs, self.corrected_energies))
        # Qs, energies = [], []
        # for sp in self:
        #     # if disp_ratio_range and not (disp_ratio_range[0]
        #     #                              <= single_point.disp_ratio
        #     #                              <= disp_ratio_range[1]):
        #     #     continue
        #     Qs.append(sp.Q)
        #     energy = (sp.energy + sp.ccd_correction_energy
        #               + self.spec.correction_energy + self.shifted_energy)
        #     energies.append(energy)
        # return Qs, energies

    def single_point_from_disp(self, disp_ratio: float):
        """Return the point whose displacement ratio matches the input."""
        for sp in self:
            if np.isclose(sp.disp_ratio, disp_ratio):
                return sp
        raise ValueError(f"No single point found for disp_ratio={disp_ratio}")

    @property
    def lowest_energy_single_point(self) -> SinglePoint:
        """Return the single point with the lowest corrected energy."""
        return min(self.single_points, key=lambda sp: sp.ccd_corrected_energy)

    @property
    def lowest_energy(self) -> float:
        """Return the lowest corrected energy among the points."""
        return self.lowest_energy_single_point.ccd_corrected_energy

    # TODO: remove this
    def verify_Q0_has_the_lowest_energy(self):
        """Raise if the minimum-energy point is not at Q=0."""
        if not np.isclose(self.lowest_energy_single_point.Q, 0.0):
            dQ = self.lowest_energy_single_point.Q
            raise ValueError(
                f"The single point with the lowest energy has dQ={dQ:.3f}, "
                f"which is not 0. ")

    def verify_num_Q(self, f):
        """Raise if there are too few points for a fitting function."""
        sig = inspect.signature(f)
        n_params = len(sig.parameters) - 1
        if n_params > len(self):
            raise ValueError(f"The number of Q points must be >= {n_params}.")

    def flip(self, shifter: Shifter, Q_diff: float, total_energy_correction: float) -> "SinglePoints":
        """Return shifted points, optionally measured from the other endpoint."""
        result = []
        for sp in self:
            new_sp = deepcopy(sp)
            new_sp.energy += shifter.shift_energy + total_energy_correction
            if shifter.flip:
                new_sp.spec = sp.spec.flip(Q_diff)
            result.append(new_sp)
        return SinglePoints(result)

    # def __str__(self):
    #     sort Single points

def make_fitting_func(curve: Type[FittingFunc], single_points: SinglePoints) -> FittingFunc:
    """Fit a potential function class to a set of single points."""
    # TODO: Consider if Q0, E0 need to be fixed or not.
    vals, _ = curve_fit(curve.fitting_func,
                        single_points.Qs,
                        single_points.corrected_energies)
    # vals, _ = curve_fit(f, self.Qs, self.corrected_energies, bounds=bounds)

    kwargs = {'Q0': 0.0, 'E0': vals[0]}
    param_names = list(inspect.signature(curve.fitting_func).parameters.keys())[2:]
    for i, name in enumerate(param_names):
        kwargs[name] = vals[i + 1]
    return curve(**kwargs)


@dataclass
class PotentialCurveSpec(MSONable, ToJsonFileMixIn):
    """Metadata for one charge-state potential curve."""
    charge: int
    correction_energy: float  # at relaxed structure, e.g., eFNV correction
    counter_charge: int  # charge state to which the structure shifts
    Q_diff: float  # Q difference between two charge states

    @property
    def effective_charge(self):
        return self.charge - self.counter_charge


def make_shifter(spec: PotentialCurveSpec,
                 single_points: SinglePoints,
                 offset: float = 0.0,
                 flip: bool = False) -> Shifter:
    """Create the energy shift that places the curve minimum at the offset."""
    lowest_energy = single_points.lowest_energy + spec.correction_energy
    shift_energy = - lowest_energy + offset
    return Shifter(shift_energy, flip)


@dataclass
class PotentialCurve(MSONable, ToJsonFileMixIn):
    """Potential-energy curve for one charge state.

    ``original_single_points`` preserves the raw parsed points.  The
    ``single_points`` property applies the configured shifter and optional
    coordinate flip before fitting or plotting.

    """

    spec: PotentialCurveSpec
    original_single_points: SinglePoints # Bare energies.
    shifter: Shifter
    fitting_curve: Optional[FittingFunc] = None

    @cached_property
    def single_points(self) -> SinglePoints:
        """Return single points after applying the curve shifter."""
        if self.shifter:
            return self.original_single_points.flip(self.shifter,
                                                    self.Q_diff,
                                                    self.spec.correction_energy)
        return self.original_single_points

    @property
    def charge(self) -> int:
        """Return the charge state of this curve."""
        return self.spec.charge

    @property
    def counter_charge(self) -> int:
        """Return the charge state at the opposite endpoint."""
        return self.spec.counter_charge

    @property
    def Qs_and_energies(self) -> List[Tuple[float, float]]:
        """Return paired Q coordinates and corrected energies."""
        return self.single_points.Qs_and_energies

    @property
    def Q_diff(self) -> float:
        """Return the endpoint separation in amu^0.5 Angstrom."""
        return self.spec.Q_diff
    # TODO: consider why they need to be sorted.
    # def __post_init__(self):
    #     self.original_single_points \
    #         = list(sorted(self.original_single_points, key=lambda x: x.Q))

    @property
    def lowest_energy_single_point(self) -> SinglePoint:
        """Return the lowest-energy shifted single point."""
        return self.single_points.lowest_energy_single_point

    @property
    def lowest_energy(self) -> float:
        """Return the lowest energy after adding the curve correction."""
        return self.lowest_energy_single_point.ccd_corrected_energy + \
                self.spec.correction_energy

    def set_fitting_curve(self, curve: Type[FittingFunc]) -> None:
        """Fit and store a fitting-function instance for this curve."""
        # TODO: Consider if Q0, E0 need to be fixed or not.
        vals, _ = curve_fit(curve.fitting_func,
                            self.single_points.Qs,
                            self.single_points.corrected_energies)
        # vals, _ = curve_fit(f, self.Qs, self.corrected_energies, bounds=bounds)

        kwargs = {'Q0': 0.0, 'E0': vals[0]}
        param_names = list(inspect.signature(curve.fitting_func).parameters.keys())[2:]
        for i, name in enumerate(param_names):
            kwargs[name] = vals[i + 1]
        self.fitting_curve = curve(**kwargs)

    @property
    def table_for_plot(self):
        """Return headers and values summarizing this curve."""
        tables = {"charge": self.charge,
                  "lowest energy": self.lowest_energy,
                  "counter charge": self.counter_charge,
                  "Q diff": self.Q_diff}
        return list(tables.keys()), tables.values()

    def __str__(self):
        """Return a tabulated summary of the curve and single points."""
        headers, tabulate_data = self.table_for_plot
        fitting_curve_str = str(self.fitting_curve) if self.fitting_curve \
            else "fitted curve is N.A."
        tables = [tabulate([tabulate_data],
                           tablefmt="plain", floatfmt=".3f", headers=headers),
                  fitting_curve_str]

        if self.single_points:
            single_point_table = tabulate(
                [sp.table_values() for sp in self.single_points],
                tablefmt="plain",
                floatfmt=".3f",
                headers=headers)
        else:
            single_point_table = "single points are N.A."
        tables.append(single_point_table)

        return "\n".join(tables)

