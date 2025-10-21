# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from copy import deepcopy

from vise.util.logger import get_logger

from pydefect_ccd.ccd import PotentialCurveResult, Ccd, dQ_revert
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
                 ground_curve: PotentialCurveResult,
                 excited_curve: PotentialCurveResult,
                 ccd_init: CcdInit):
        self.ccd_init = ccd_init
        self.orig_ground_curve = ground_curve
        self.orig_excited_curve = excited_curve

        self.charge_diff = self.orig_excited_curve.charge - self.orig_ground_curve.charge
        if self.charge_diff == 1:
            self.band_edge = "cbm"
        elif self.charge_diff == -1:
            self.band_edge = "vbm"
        else:
            raise ValueError

        self.orig_ground_curve.shifted_energy = ground_curve.charge * self.band_edge_level
        print(excited_curve.charge * self.band_edge_level)
        self.orig_excited_curve.shifted_energy = excited_curve.charge * self.band_edge_level

    @property
    def carrier_coexisting_with_excited(self):
        carrier_charge = - self.charge_diff
        return Carrier.from_carrier_charge(carrier_charge)

    @property
    def band_edge_level(self) -> float:
        return self.ccd_init.__getattribute__(str(self.band_edge))

    def carrier_energy(self, carrier: Carrier) -> float:
        if carrier == Carrier.e:
            band_edge = self.ccd_init.cbm
        elif carrier == Carrier.h:
            band_edge = self.ccd_init.vbm
        else:
            raise ValueError
        return - carrier.charge * band_edge

    @property
    def _ground_ccd(self) -> PotentialCurveResult:
        """
        Returns: Ground state ccd with the lowest energy to be zero
        """
        return self.orig_ground_curve

    @property
    def _ground_pn_ccd(self) -> PotentialCurveResult:
        """
        Returns: Ground state ccd with an excited electron at CBM from VBM
        """
        result = deepcopy(self.orig_ground_curve)
        result.carriers = [Carrier.e, Carrier.h]
        result.shifted_energy += self.ccd_init.cbm - self.ccd_init.vbm
        return result

    @property
    def _excited_ccd(self) -> PotentialCurveResult:
        """
        Returns: Excited state ccd capturing _default_single_ccd_for_e_p_coupling minority carrier with _default_single_ccd_for_e_p_coupling majority
                 carrier. The Q values are reverted and zero is set to that of
                 ground state ccd.
        """
        try:
            result = dQ_revert(self.orig_excited_curve)
            result.carriers = [self.carrier_coexisting_with_excited]
            return result
        except ValueError:
            logger.warning("disp=1.0 is required for excited state.")
            raise

    @property
    def ccd(self) -> Ccd:
        # ccds = [self._ground_ccd, self._excited_ccd]
        return Ccd(name=self.ccd_init.name,
                   ground_curve=self._ground_ccd,
                   excited_curve=self._excited_ccd)

    # @property
    # def ccd(self) -> Ccd:
    #     ccds = [self._ground_ccd, self._excited_ccd, self._ground_pn_ccd]
    #     return Ccd(name=self.ccd_init.name, potential_curve_results=ccds)
