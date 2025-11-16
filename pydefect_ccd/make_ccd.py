# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from copy import deepcopy

from vise.util.logger import get_logger

from pydefect_ccd.ccd import PotentialCurve, Ccd
from pydefect_ccd.ccd_init import CcdInit
from pydefect_ccd.enum import Carrier

logger = get_logger(__name__)


class MakeCcd:
    """
    Define charge_diff = q_{excited} - q_{ground}.
    - charge_diff = +1
      - Fermi level locates at the VBM (p-type)
      - carriers are recombined via (ground + h + e) → (excited + h) → ground
      - minority carrier is e and majority carrier is h

    - charge_diff = -1
      - Fermi level locates at the CBM (n-type)
      - carriers are recombined via (ground + h + e) → (excited + e) → ground
      - minority carrier is h and majority carrier is e
    """
    def __init__(self,
                 ground_curve: PotentialCurve,
                 excited_curve: PotentialCurve,
                 ccd_init: CcdInit):
        self.ccd_init = ccd_init
        self.ground_curve = deepcopy(ground_curve)
        self.excited_curve = deepcopy(excited_curve)

        self.ground_curve.shifted_energy \
            = ground_curve.charge * self._band_edge_level
        self.excited_curve.shifted_energy \
            = excited_curve.charge * self._band_edge_level

    @property
    def _charge_diff(self):
        return self.excited_curve.charge - self.ground_curve.charge

    @property
    def _carrier_in_excited_state(self):
        carrier_charge = - self._charge_diff
        return Carrier.from_carrier_charge(carrier_charge)

    @property
    def _band_edge_level(self) -> float:
        if self._charge_diff == 1:
            self._band_edge = "cbm"
        elif self._charge_diff == -1:
            self._band_edge = "vbm"
        else:
            raise ValueError
        return self.ccd_init.__getattribute__(str(self._band_edge))

    @property
    def _carrier_energy(self) -> float:
        if self._carrier_in_excited_state == Carrier.e:
            band_edge = self.ccd_init.cbm
        elif self._carrier_in_excited_state == Carrier.h:
            band_edge = self.ccd_init.vbm
        else:
            raise ValueError
        return - self._carrier_in_excited_state.charge * band_edge

    @property
    def ccd(self) -> Ccd:
        return Ccd(name=self.ccd_init.name,
                   ground_curve=self.ground_curve,
                   excited_curve=self.excited_curve)

    # @property
    # def ccd(self) -> Ccd:
    #     ccds = [self._ground_ccd, self._excited_ccd, self._ground_pn_ccd]
    #     return Ccd(name=self.ccd_init.name, potential_curve_results=ccds)
