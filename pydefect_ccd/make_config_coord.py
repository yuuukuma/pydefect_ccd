# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from copy import deepcopy

from dephon.ccd_init import CcdInit
from dephon.config_coord import PotentialCurve, Ccd, CcdId
from dephon.enum import Carrier, BandEdge
from vise.util.logger import get_logger

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
                 ground_ccd: PotentialCurve,
                 excited_ccd: PotentialCurve,
                 ccd_init: CcdInit):
        if abs(ground_ccd.charge - excited_ccd.charge) != 1:
            raise AssertionError(
                f"The charge difference needs to be 1. Now, ground state "
                f"{ground_ccd.charge} and excited state {excited_ccd.charge}.")
        self.dephon_init = ccd_init
        self.orig_ground_ccd = ground_ccd
        self.orig_excited_ccd = excited_ccd

        self.orig_ground_ccd.shift_energy(
            ground_ccd.charge * self.band_edge_level)
        self.orig_excited_ccd.shift_energy(
            excited_ccd.charge * self.band_edge_level)

        self.ref_energy = ground_ccd.get_single_point_by_disp_ratio(0.0).corrected_energy

    @property
    def charge_diff(self) -> int:
        return self.orig_excited_ccd.charge - self.orig_ground_ccd.charge

    @property
    def carrier_coexisting_with_excited(self):
        """ This should be _default_single_ccd_for_e_p_coupling majority carrier."""
        return Carrier.from_carrier_charge(- self.charge_diff)

    @property
    def band_edge(self) -> BandEdge:
        if self.charge_diff == 1:
            return BandEdge.cbm
        else:
            return BandEdge.vbm

    @property
    def band_edge_level(self) -> float:
        return self.dephon_init.__getattribute__(str(self.band_edge))

    def carrier_energy(self, carrier: Carrier) -> float:
        if carrier == Carrier.e:
            band_edge = self.dephon_init.cbm
        elif carrier == Carrier.h:
            band_edge = self.dephon_init.vbm
        else:
            raise ValueError
        return - carrier.charge * band_edge

    @property
    def _ground_ccd(self) -> PotentialCurve:
        """
        Returns: Ground state ccd with the lowest energy to be zero
        """
        result = deepcopy(self.orig_ground_ccd)
        result.set_base_energy(self.ref_energy)
        result.id_ = CcdId("ground")
        return result

    @property
    def _ground_pn_ccd(self) -> PotentialCurve:
        """
        Returns: Ground state ccd with an excited electron at CBM from VBM
        """
        result = deepcopy(self.orig_ground_ccd)
        result.set_base_energy(self.ref_energy)
        result.id_ = CcdId("ground",
                           carriers=[Carrier.h, Carrier.e])
        result.shift_energy(self.dephon_init.cbm - self.dephon_init.vbm)
        return result

    @property
    def _excited_ccd(self) -> PotentialCurve:
        """
        Returns: Excited state ccd capturing _default_single_ccd_for_e_p_coupling minority carrier with _default_single_ccd_for_e_p_coupling majority
                 carrier. The Q values are reverted and zero is set to that of
                 ground state ccd.
        """
        try:
            print(self.orig_excited_ccd)
            result = self.orig_excited_ccd.dQ_reverted_single_ccd()
        except ValueError:
            logger.warning("disp=1.0 is required for excited state.")
            raise

        result.set_base_energy(self.ref_energy)
        result.id_ = CcdId(
            "excited", carriers=[self.carrier_coexisting_with_excited])
        return result

    @property
    def ccd(self) -> Ccd:
        ccds = [self._ground_ccd, self._excited_ccd, self._ground_pn_ccd]
        for ccd in ccds:
            ccd.set_quadratic_fitting_range()

        return Ccd(name=self.dephon_init.name, potential_curves=ccds)
