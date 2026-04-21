# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import matplotlib.ticker as mticker

from monty.json import MSONable
from nonrad import get_C
from nonrad.constants import HBAR
from vise.util.matplotlib import float_to_int_formatter
from vise.util.mix_in import ToJsonFileMixIn

from pydefect_ccd.potential_curve import PotentialCurve
from pydefect_ccd.fitting_curve import QuadraticFittingCurve


@dataclass
class TotalSquaredTransitionMoment(MSONable, ToJsonFileMixIn):
    Ts: List[float]
    total_moments: List[float]  # amu A^2 / eV

    def add_plot(self, ax):
        ax.plot(self.Ts, self.total_moments)


class PlottedTotalSquaredTransitionMoment:
    def __init__(self,
                 total_squared_transition_moment: TotalSquaredTransitionMoment,
                 title: Optional[str] = None):
        self._transition_moment = total_squared_transition_moment
        self._title = title or ""

    def plot(self, ax: Optional[Axes] = None) -> Tuple[Figure, Axes]:
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        fmt = mticker.ScalarFormatter(useMathText=True)
        fmt.set_scientific(True)
        fmt.set_powerlimits((-2, 3))  # この範囲を外れると 1eN 表記にする

        self._add_transition_moment(ax)
        self._set_title(ax)
        self._set_formatter(ax)
        self._set_labels(ax)

        ax.yaxis.set_major_formatter(fmt)
        fig.tight_layout()
        return fig, ax

    def _add_transition_moment(self, ax: Axes) -> None:
        self._transition_moment.add_plot(ax)

    def _set_title(self, ax: Axes) -> None:
        ax.set_title(self._title)

    def _set_formatter(self, ax: Axes) -> None:
        ax.xaxis.set_major_formatter(float_to_int_formatter)
        ax.yaxis.set_major_formatter(float_to_int_formatter)

    def _set_labels(self, ax: Axes) -> None:
        ax.set_xlabel("T (K)")
        ax.set_ylabel("$Total squared moment$ (amu Å$^2$ / eV)")

        # ラベルが付いた artists があるときだけ legend 表示（空凡例回避）
        handles, labels = ax.get_legend_handles_labels()
        if labels:
            ax.legend()

@dataclass
class CalcTotalSquaredTransitionMoment:
    ground_curve: PotentialCurve
    excited_curve: PotentialCurve
    Ts: List[float]

    def harmonic(self, overlap_method: str = "HermiteGauss") -> TotalSquaredTransitionMoment:
        """Within harmonic approximation. """
        dQ = self.excited_curve.Q_diff
        dE = abs(self.excited_curve.lowest_energy - self.ground_curve.lowest_energy)

        assert isinstance(self.ground_curve.fitting_curve, QuadraticFittingCurve)
        assert isinstance(self.excited_curve.fitting_curve, QuadraticFittingCurve)

        # at Wif=1, volume=1cm^3 * HBAR / 2 / np.pi, g=1: Then, unit is in amu x Å^2 / eV.
        # dQ [amu^0.5 x Å], dE [eV], wi [eV], wf [eV], T [K]
        result = get_C(dQ=abs(dQ),
                       dE=dE,
                       wi=self.excited_curve.fitting_curve.omega_in_eV,
                       wf=self.ground_curve.fitting_curve.omega_in_eV,
                       T=np.array(self.Ts),
                       Wif=1.0,
                       volume=1e24 * HBAR / 2 / np.pi,
                       g=1,
                       overlap_method=overlap_method)
        return TotalSquaredTransitionMoment(self.Ts, result.tolist())


# def calc_summed_squared_transition_moment(
#         ground_curve: PotentialCurve,
#         excited_curve: PotentialCurve,
#         Ts: List[float],
#         overlap_method: str = "HermiteGauss") -> List[float]:
#     """Within harmonic approximation. Unit is in amu Å^2 / eV."""
#     dQ = excited_curve.Q_diff
#     dE = abs(excited_curve.lowest_energy - ground_curve.lowest_energy)

#     assert isinstance(ground_curve.fitted_curve, QuadraticCurve)
#     assert isinstance(excited_curve.fitted_curve, QuadraticCurve)

#     # at Wif=1, volume=1cm^3/ 2 * np.pi, g=1
#     result = get_C(dQ=abs(dQ),
#                    dE=dE,
#                    wi=excited_curve.fitted_curve.omega_in_eV,
#                    wf=ground_curve.fitted_curve.omega_in_eV,
#                    T=np.array(Ts),
#                    overlap_method=overlap_method,
#                    Wif=1,
#                    volume=1e8 ** 3 / (2 * np.pi),
#                    g=1)
#     return list(result)
