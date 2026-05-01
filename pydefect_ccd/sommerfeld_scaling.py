# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
from dataclasses import dataclass, field
from typing import List, Dict, Tuple
from monty.json import MSONable

import numpy as np
from nonrad.scaling import sommerfeld_parameter
from vise.util.mix_in import ToJsonFileMixIn

from pydefect_ccd.local_enum import Carrier

Key = Tuple[str, int]  # ("e" or "h", defect_charge)


@dataclass
class SommerfeldScaling(MSONable, ToJsonFileMixIn):
    """Sommerfeld factors on a temperature grid.

    The scaling depends on the dielectric constant, carrier effective mass,
    carrier type, and defect charge.  Values are computed lazily for each
    ``(carrier, charge)`` pair and stored in ``_scaling``.  Neutral defects
    return unity scaling without calling ``nonrad``.

    Attributes:
        epsilon0: Static dielectric constant, i.e. 1 + \epsilon_ele + \epsilon_ion.
        electron_effective_mass: Electron effective mass in units of m0.
        hole_effective_mass: Hole effective mass in units of m0.
        Ts: Temperatures in K.
        _scaling: Cached arrays keyed by ``("e" or "h", defect_charge)``.
    """
    epsilon0: float
    electron_effective_mass: float
    hole_effective_mass: float
    Ts: List[float]
    _scaling: Dict[Key, np.ndarray] = field(default_factory=dict)

    def scaling(self, carrier_type: Carrier, defect_charge: int) -> np.ndarray:
        """Return scaling factors for a carrier and defect charge."""
        if defect_charge == 0:
            return np.ones_like(self.Ts, dtype=float)

        key = (str(carrier_type), defect_charge)
        if key not in self._scaling:
            self._scaling[key] = self._compute_scaling(carrier_type,
                                                       defect_charge)
        return self._scaling[key]

    def _compute_scaling(self,
                         carrier_type: Carrier,
                         defect_charge: int,
                         method: str = "Integrate") -> np.ndarray:
        """Compute a non-neutral scaling array."""
        Z = defect_charge * carrier_type.charge
        mass = (self.electron_effective_mass
                if carrier_type is Carrier.e else self.hole_effective_mass)
        Ts = np.array(self.Ts)
        return sommerfeld_parameter(Ts, Z, mass, self.epsilon0, method=method)

    def add_to_ax(self, ax, carrier_type, defect_charge, ls="--"):
        """Add one carrier-charge scaling curve to an axes."""
        y = self.scaling(carrier_type, defect_charge)
        ax.plot(self.Ts, y, linestyle=ls)

    @staticmethod
    def set_label(ax, ls="--"):
        """Set plot labels and mark 300 K with a vertical line."""
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("Sommerfeld scaling")
        ax.legend()
        ax.axvline(x=300, color='red', ls=ls, lw=0.5)
