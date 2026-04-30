# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from typing import Type

from numpy.ma.testutils import assert_almost_equal
from vise.util.logger import get_logger

from pydefect_ccd.ccd import Ccd
from pydefect_ccd.fitting_curve import FittingCurve
from pydefect_ccd.local_enum import Carrier
from pydefect_ccd.potential_curve import PotentialCurve, SinglePoints, \
    PotentialCurveSpec, make_shifter

logger = get_logger(__name__)


class MakeCcd:
    """Make Ccd instance from ground and excited potential curves.

    This automatically determine the relevant band edge and carrier type.
    The total energy is shifted by setting the Fermi level at the band edge.

    The dQs of the Excited state are reverted.

    Define charge_diff = q_{excited} - q_{ground}.
    - charge_diff = +1
      - Fermi level locates at the VBM (p-type)
      - carriers are recombined via (excited + h) → ground

    - charge_diff = -1
      - Fermi level locates at the CBM (n-type)
      - carriers are recombined via  (excited + e) → ground
    """
    def __init__(self,
                 ground_single_points: SinglePoints,
                 ground_pot_curve_spec: PotentialCurveSpec,
                 ground_fitting_curve: Type[FittingCurve],
                 excited_single_points: SinglePoints,
                 excited_pot_curve_spec: PotentialCurveSpec,
                 excited_fitting_curve: Type[FittingCurve],
                 vbm: float,
                 cbm: float,
                 name: str):
        self._vbm = vbm
        self._cbm = cbm
        self._name = name

        assert ground_pot_curve_spec.counter_charge == excited_pot_curve_spec.charge
        assert excited_pot_curve_spec.counter_charge == ground_pot_curve_spec.charge
        assert_almost_equal(ground_pot_curve_spec.Q_diff, excited_pot_curve_spec.Q_diff)

        ground_shifter = make_shifter(ground_pot_curve_spec, ground_single_points)
        self._ground_curve = PotentialCurve(ground_pot_curve_spec,
                                            ground_single_points,
                                            ground_shifter)
        self._ground_curve.set_fitting_curve(ground_fitting_curve)

        excited_offset = self._shifted_energy(excited_pot_curve_spec.charge)
        excited_shifter = make_shifter(excited_pot_curve_spec,
                                       excited_single_points,
                                       excited_offset,
                                       flip=True)
        self._excited_curve = PotentialCurve(excited_pot_curve_spec,
                                             excited_single_points,
                                             excited_shifter)
        self._excited_curve.set_fitting_curve(excited_fitting_curve)

    @property
    def _charge_diff(self):
        return self._excited_curve.charge - self._ground_curve.charge

    @property
    def _carrier_in_excited_state(self):
        carrier_charge = - self._charge_diff
        return Carrier.from_carrier_charge(carrier_charge)

    # @property
    # def _band_edge_level(self) -> float:
    #     if self._charge_diff == 1:
    #         return self._cbm
    #     elif self._charge_diff == -1:
    #         return self._vbm
    #     else:
    #         raise ValueError("The charge difference must be ±1.")

    def _shifted_energy(self, charge) -> float:
        if self._carrier_in_excited_state == Carrier.e:
            band_edge = self._cbm
        elif self._carrier_in_excited_state == Carrier.h:
            band_edge = self._vbm
        else:
            raise ValueError
        return band_edge * charge

    @property
    def ccd(self) -> Ccd:
        return Ccd(name=self._name,
                   ground_curve=self._ground_curve,
                   excited_curve=self._excited_curve)

