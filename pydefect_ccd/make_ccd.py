# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from copy import deepcopy

from vise.util.logger import get_logger

from pydefect_ccd.ccd import Ccd, dQ_revert
from pydefect_ccd.potential_curve import PotentialCurve
from pydefect_ccd.local_enum import Carrier

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
                 ground_curve: PotentialCurve,
                 excited_curve: PotentialCurve,
                 vbm: float, cbm: float,
                 name: str):
        self._vbm = vbm
        self._cbm = cbm
        self._name = name

        assert ground_curve.counter_charge == excited_curve.charge
        assert excited_curve.counter_charge == ground_curve.charge
        assert ground_curve.Q_diff == excited_curve.Q_diff

        self._ground_curve = deepcopy(ground_curve)
        self._excited_curve = deepcopy(excited_curve)
        self._lowest_energy = ground_curve.lowest_energy

        self._set_shifted_energy(self._ground_curve)
        self._set_shifted_energy(self._excited_curve)

    @property
    def _charge_diff(self):
        return self._excited_curve.charge - self._ground_curve.charge

    @property
    def _carrier_in_excited_state(self):
        carrier_charge = - self._charge_diff
        return Carrier.from_carrier_charge(carrier_charge)

    @property
    def _band_edge_level(self) -> float:
        if self._charge_diff == 1:
            return self._cbm
        elif self._charge_diff == -1:
            return self._vbm
        else:
            raise ValueError("The charge difference must be ±1.")

    def _set_shifted_energy(self, curve: PotentialCurve) -> None:
        if self._carrier_in_excited_state == Carrier.e:
            band_edge = self._cbm
        elif self._carrier_in_excited_state == Carrier.h:
            band_edge = self._vbm
        else:
            raise ValueError
        curve.shifted_energy = band_edge * curve.charge - self._lowest_energy

    @property
    def ccd(self) -> Ccd:
        return Ccd(name=self._name,
                   ground_curve=self._ground_curve,
                   excited_curve=dQ_revert(self._excited_curve))

