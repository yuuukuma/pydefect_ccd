# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
from typing import List

import numpy as np
from nonrad.scaling import sommerfeld_parameter

class SommerfeldScaling:

    def __init__(self,
                 dielectric_constant: float,
                 electron_effective_mass: float,
                 hole_effective_mass: float,
                 Ts: List[float] = None,
                 defect_charge_range : range = range(-2, 3)):
        self.epsilon0 = dielectric_constant
        self.e_mass = electron_effective_mass
        self.h_mass = hole_effective_mass
        self.Ts = Ts if Ts is not None else list(range(100, 1001, 10))

        self.scaling = {}
        for a, mass in [("e", self.e_mass), ("h", self.h_mass)]:
            for defect_charge in defect_charge_range:
                if defect_charge == 0:
                    self.scaling[(a, defect_charge)] = np.ones_like(Ts)
                carrier_charge = 1 if a == "h" else -1
                Z = defect_charge * carrier_charge
                self.scaling[(a, defect_charge)] \
                    = sommerfeld_parameter(np.array(self.Ts), Z, mass, self.epsilon0)

    def plot(self, ax, key, ls="--"):
        """Plot Sommerfeld scaling curves on given matplotlib Axes."""
        y = self.scaling[key]
        ax.plot(self.Ts, y, label=key, linestyle=ls)
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("Sommerfeld scaling")
        ax.legend()
        ax.axvline(x=300, color='red', ls=ls, lw=0.5)
        return ax