# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
from vise.util.logger import get_logger

from pydefect_ccd.ccd import Ccd
from pydefect_ccd.local_enum import Carrier
from pydefect_ccd.potential_curve import PotentialCurve

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
                 ground_pot_curve: PotentialCurve,
                 excited_pot_curve: PotentialCurve,
                 vbm: float,
                 cbm: float,
                 name: str):
        """Prepare shifted potential curves and their fitted functions."""
        self._vbm = vbm
        self._cbm = cbm
        self._ground_charge = ground_pot_curve.charge
        self._excited_charge = excited_pot_curve.charge

        offset = self._shifted_energy(excited_pot_curve.charge)

        self._ccd = Ccd(name=name,
                        ground_curve=ground_pot_curve.shifted(),
                        excited_curve=excited_pot_curve.shifted(offset,
                                                                flip=True))

    @property
    def _charge_diff(self):
        """Return excited-charge minus ground-charge."""
        return self._excited_charge - self._ground_charge

    @property
    def _carrier_in_excited_state(self):
        """Return the carrier present in the excited state."""
        carrier_charge = - self._charge_diff
        return Carrier.from_carrier_charge(carrier_charge)

    def _shifted_energy(self, charge) -> float:
        """Return the Fermi-level energy offset for a charge state."""
        if self._carrier_in_excited_state == Carrier.e:
            band_edge = self._cbm
        elif self._carrier_in_excited_state == Carrier.h:
            band_edge = self._vbm
        else:
            raise ValueError
        return band_edge * charge

    @property
    def ccd(self) -> Ccd:
        """Return the constructed CCD object."""
        return self._ccd
