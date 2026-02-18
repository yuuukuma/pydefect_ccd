# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
from typing import List

import numpy as np
from nonrad.scaling import sommerfeld_parameter

from pydefect_ccd.enum import Carrier


class SommerfeldScaling:

    def __init__(self,
                 dielectric_constant: float,
                 electron_effective_mass: float,
                 hole_effective_mass: float,
                 Ts: List[float] = None,
                 max_defect_charge = 3):
        self.epsilon0 = dielectric_constant
        self.e_mass = electron_effective_mass
        self.h_mass = hole_effective_mass
        self.Ts = Ts if Ts is not None else list(range(100, 1001, 10))

        self.scaling = {}
        for a, b in [("e", self.e_mass), ("h", self.h_mass)]:
            # for c in range(-max_defect_charge, max_defect_charge + 1):
            for c in range(1, 2):
                if c == 0:
                    self.scaling[a, c] = np.ones_like(Ts)
                Z = c * (1 if a == "h" else -1)
                self.scaling[(a, c)] \
                    = sommerfeld_parameter(np.array(self.Ts), Z, self.e_mass, self.epsilon0)

    def plot(self, ax, keys):
        """Plot Sommerfeld scaling curves on given matplotlib Axes."""
        for key in keys:
            y = self.scaling[key]
            ls = "-" if key[1] >= 0 else "--"
            ax.plot(self.Ts, y, label=key, linestyle=ls)
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("Sommerfeld scaling")
        ax.legend()
        return ax