# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass
from typing import List

import numpy as np
from monty.json import MSONable
from nonrad import get_C
from tabulate import tabulate
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

    @property
    def Wif(self):
        return self.W_if

    def __str__(self):
        header = [["Wif:", f"{self.W_if:.1e}"],
                  ["site degeneracy:", f"{self.site_degeneracy}"]]

        result = [tabulate(header, tablefmt="plain")]

        table = []
        columns = ["T [K]", "Phonon overlap []", "C [cm3/s]", "v [cm2/s]", "c / v [cm2]"]
        for T, transition_moment, rate, v in zip(self.Ts,
                                                 self.summed_squared_transition_moment,
                                                 self.capture_rate,
                                                 self.velocities):
            table.append([T, transition_moment, rate, v, rate / v])

        result.append(
            tabulate(table, headers=columns, tablefmt="plain",
                     floatfmt=[".1f", ".1e", ".1e", ".1e", ".1e"]))

        return "\n".join(result)


def calc_summed_squared_transition_moment(
        ground_curve: PotentialCurve,
        excited_curve: PotentialCurve,
        Ts: List[float]):
    """Within harmonic approximation"""
    dQ = excited_curve.Q_diff - ground_curve.Q_diff
    dE = excited_curve.lowest_energy - ground_curve.lowest_energy

    assert isinstance(ground_curve.fitted_curve, QuadraticCurve)
    assert isinstance(excited_curve.fitted_curve, QuadraticCurve)

    # at Wif=1, volume=1Å^3, g=1
    result = get_C(dQ=abs(dQ),
                   dE=dE,
                   wi=ground_curve.fitted_curve.omega,
                   wf=excited_curve.fitted_curve.omega,
                   T=np.array(Ts),
                   Wif=1, volume=1, g=1)
    return list(result)

