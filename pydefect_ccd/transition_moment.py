# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
from typing import List

import numpy as np
from nonrad import get_C

from pydefect_ccd.ccd import PotentialCurve, QuadraticCurve


def calc_summed_squared_transition_moment(
        ground_curve: PotentialCurve,
        excited_curve: PotentialCurve,
        Ts: List[float],
        overlap_method: str = "HermiteGauss") -> List[float]:
    """Within harmonic approximation. Unit is in amu Å^2 / eV."""
    dQ = excited_curve.Q_diff
    dE = abs(excited_curve.lowest_energy - ground_curve.lowest_energy)

    assert isinstance(ground_curve.fitted_curve, QuadraticCurve)
    assert isinstance(excited_curve.fitted_curve, QuadraticCurve)

    # at Wif=1, volume=1cm^3/ 2 * np.pi, g=1
    result = get_C(dQ=abs(dQ),
                   dE=dE,
                   wi=excited_curve.fitted_curve.omega_in_eV,
                   wf=ground_curve.fitted_curve.omega_in_eV,
                   T=np.array(Ts),
                   overlap_method=overlap_method,
                   Wif=1,
                   volume=1e8 ** 3 / (2 * np.pi),
                   g=1)
    return list(result)
