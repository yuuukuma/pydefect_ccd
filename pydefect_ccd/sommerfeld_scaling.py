# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
from dataclasses import dataclass, field
from typing import List, Dict, Any, Tuple
from monty.json import MSONable

import numpy as np
from nonrad.scaling import sommerfeld_parameter
from vise.util.mix_in import ToJsonFileMixIn

from pydefect_ccd.local_enum import Carrier

Key = Tuple[str, int]  # ("e" or "h", defect_charge)


@dataclass
class SommerfeldScaling(MSONable, ToJsonFileMixIn):
    """Cache Sommerfeld scaling factors for carriers and defect charges."""
    epsilon0: float
    electron_effective_mass: float
    hole_effective_mass: float
    Ts: List[float]
    _scaling: Dict[Key, np.ndarray] = field(default_factory=dict)

    def scaling(self, carrier_type: Carrier, defect_charge: int) -> np.ndarray:
        """Return cached or newly computed scaling factors."""
        if defect_charge == 0:
            return np.ones_like(self.Ts)

        key = (str(carrier_type), defect_charge)
        if key not in self._scaling:
            self.get_scaling(carrier_type, defect_charge)
        return self._scaling[key]

    def get_scaling(self,
                    carrier_type: Carrier, defect_charge: int,
                    method: str = "Integrate"):
        """Compute and cache scaling factors for one carrier-charge pair."""
        Z = defect_charge * carrier_type.charge
        mass = (self.electron_effective_mass
                if carrier_type is Carrier.e else self.hole_effective_mass)
        print("mass", mass)
        Ts = np.array(self.Ts)
        self._scaling[(str(carrier_type), defect_charge)] \
            = sommerfeld_parameter(Ts, Z, mass, self.epsilon0, method=method)

    def add_to_ax(self, ax, carrier_type, defect_charge, ls="--"):
        """Plot the scaling factors on an axes."""
        y = self.scaling(carrier_type, defect_charge)
        ax.plot(self.Ts, y, linestyle=ls)

    def set_label(self, ax, ls="--"):
        """Set labels and a 300 K reference line for a scaling plot."""
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("Sommerfeld scaling")
        ax.legend()
        ax.axvline(x=300, color='red', ls=ls, lw=0.5)
