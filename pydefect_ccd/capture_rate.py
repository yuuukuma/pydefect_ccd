# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass
from typing import List

import numpy as np
from monty.json import MSONable
from nonrad import get_C
from vise.util.matplotlib import float_to_int_formatter
from vise.util.mix_in import ToJsonFileMixIn

from pydefect_ccd.ccd import PotentialCurve, QuadraticCurve


@dataclass
class CaptureRate(MSONable, ToJsonFileMixIn):
    Ts: List[float]
    W_if: List[float]
    summed_squared_transition_moment: List[float]  # as a function of T
    site_degeneracy: float
    velocities: List[float] = None # characteristic carrier velocity in [cm / s]
    # TODO: add spin selection factor

    @property
    def capture_rate(self) -> np.array:
        return (2 * np.pi * self.site_degeneracy
                * np.array(self.W_if) ** 2
                * np.array(self.summed_squared_transition_moment))

    # def __str__(self):
    #     header = [["Wif:", f"{self.W_if:.1e}"],
    #               ["site degeneracy:", f"{self.site_degeneracy}"]]
    #
    #     result = [tabulate(header, tablefmt="plain")]
    #
    #     table = []
    #     columns = ["T [K]", "Phonon overlap []", "C [cm3/s]", "v [cm2/s]", "c / v [cm2]"]
    #     for T, transition_moment, rate, v in zip(self.Ts,
    #                                              self.summed_squared_transition_moment,
    #                                              self.capture_rate,
    #                                              self.velocities):
    #         table.append([T, transition_moment, rate, v, rate / v])
    #
    #     result.append(
    #         tabulate(table, headers=columns, tablefmt="plain",
    #                  floatfmt=[".1f", ".1e", ".1e", ".1e", ".1e"]))
    #
    #     return "\n".join(result)


def calc_summed_squared_transition_moment(
        ground_curve: PotentialCurve,
        excited_curve: PotentialCurve,
        Ts: List[float]):
    """Within harmonic approximation"""
    dQ = excited_curve.Q_diff
    dE = abs(excited_curve.lowest_energy - ground_curve.lowest_energy)

    assert isinstance(ground_curve.fitted_curve, QuadraticCurve)
    assert isinstance(excited_curve.fitted_curve, QuadraticCurve)

    # at Wif=1, volume=1Å^3, g=1
    result = get_C(dQ=abs(dQ),
                   dE=dE,
                   wi=ground_curve.fitted_curve.omega_in_eV,
                   wf=excited_curve.fitted_curve.omega_in_eV,
                   T=np.array(Ts),
                   Wif=1, volume=1, g=1)
    print(dQ, dE, ground_curve.fitted_curve.omega_in_eV, excited_curve.fitted_curve.omega_in_eV)
    return list(result)


class CaptureRatePlotter:

    def __init__(self, capture_rate: CaptureRate, plt, title: str = None):
        self._title = title or ""
        self._capture_rate = capture_rate
        self.plt = plt

    def construct_plot(self):
        self._add_capture_rate()
        self._set_title()
        self._set_formatter()
        self._set_labels()
        self.plt.tight_layout()

    def _add_capture_rate(self):
        ax = self.plt.gca()
        ax.semilogy(self._capture_rate.Ts, self._capture_rate.capture_rate)

    def _set_labels(self):
        ax = self.plt.gca()
        ax.set_xlabel("T (K)")
        ax.set_ylabel("C$_p$ (cm$^3$/s)")
        ax.legend()

    def _set_title(self):
        self.plt.gca().set_title(self._title)

    def _set_formatter(self):
        self.plt.gca().xaxis.set_major_formatter(float_to_int_formatter)
        self.plt.gca().yaxis.set_major_formatter(float_to_int_formatter)

