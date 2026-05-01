# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
import inspect
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Tuple, Union, List

import numpy as np
from monty.json import MSONable
from nonrad.constants import HBAR, EV2J, ANGS2M, AMU2KG
from scipy.optimize import brentq
from vise.util.enum import ExtendedEnum


@dataclass
class FittingFunc(ABC):
    """Base class for fitted potential-energy functions in CCD analysis.

    ``Q0`` is the minimum position in amu^0.5 Angstrom and ``E0`` is the
    minimum energy in eV.  Concrete classes provide both an evaluatable curve
    and a static model function used by ``scipy.optimize.curve_fit``.

    Attributes:
        Q0: Minimum position along the configuration coordinate.
        E0: Minimum energy of the fitted potential.
    """
    Q0: float
    E0: float

    @abstractmethod
    def __call__(self, x: Union[float, np.array]) -> Union[float, np.array]:
        """Evaluate the fitted function at one or more Q values."""
        pass

    @abstractmethod
    def shift(self, shift_Q, shift_energy, revert=False) -> "FittingFunc":
        """Return a copy shifted in Q and energy."""
        pass

    @staticmethod
    @abstractmethod
    def fitting_func(Q: Union[float, np.array], E: float, *params) -> Union[float, np.array]:
        """Model function for fitting."""
        pass


    @classmethod
    def n_fit_params(cls) -> int:
        """Return the number of parameters expected by ``fitting_func``."""
        sig = inspect.signature(cls.fitting_func)
        return len(sig.parameters) - 1

    def add_plot(self, ax, x_range: List[float], color):
        """Draw the fitted function on an axes over the given Q range."""
        xs = np.linspace(x_range[0], x_range[1], 1000)
        ys = self(xs)
        ax.plot(xs, ys, color=color)

    @property
    def omega_in_eV(self) -> float:
        """Return the phonon frequency converted to eV."""
        if not hasattr(self, "omega"):
            raise AttributeError(f"{self.__class__} does not have 'omega' attribute.")

        return HBAR * self.omega * np.sqrt(EV2J / (ANGS2M**2 * AMU2KG))


@dataclass
class QuadraticFittingFunc(MSONable, FittingFunc):
    """Quadratic potential ``a * (Q - Q0)^2 + E0``."""
    a: float
    # omega: float  # in amu Å^2 / eV

    def __call__(self, Q: Union[float, np.ndarray]) -> Union[float, np.array]:
        """Evaluate the quadratic potential."""
        return self.a * (Q - self.Q0)**2 + self.E0

    def shift(self, shift_Q, shift_energy, revert=False) -> "QuadraticFittingFunc":
        """Return the same curvature with shifted minimum and energy."""
        new_Q0 = self.Q0 + shift_Q
        new_dE = self.E0 + shift_energy
        return QuadraticFittingFunc(new_Q0, new_dE, a=self.a)

    # TODO: Understand why Q0 is not considered here.
    @staticmethod
    def fitting_func(Q: Union[float, np.array], E: float, a) -> Union[float, np.array]:
        """Model used for fitting when the minimum is fixed at Q=0."""
        return a*Q**2 + E

    @property
    def omega(self) -> float:
        """Return the harmonic coefficient used for frequency conversion."""
        return 0.5 * self.a

    def __str__(self):
        """Return a compact summary of the quadratic fit."""
        return (f"Quadratic Curve: omega={self.omega_in_eV:.3f} (eV), "
                f"Q0={self.Q0:.3f} (amu**0.5*Å), Emin={self.E0:.3f} (eV)")


@dataclass
class QuarticFittingFunc(MSONable, FittingFunc):
    """Quartic potential without a linear term around ``Q0``."""
    a: float  # in ??
    b: float
    c: float

    def __call__(self, Q: Union[float, np.array]) -> Union[float, np.array]:
        """Evaluate the quartic potential."""
        return (self.a * (Q - self.Q0) ** 4 + self.b * (Q - self.Q0) ** 3
                + self.c * (Q - self.Q0) ** 2 + self.E0)

    def shift(self, shift_Q, shift_energy, revert=False) -> "QuarticFittingFunc":
        """Return a shifted quartic, optionally reversing the odd term."""
        new_Q0 = self.Q0 + shift_Q
        new_dE = self.E0 + shift_energy
        new_b = -self.b if revert else self.b
        return QuarticFittingFunc(a=self.a, b=new_b, c=self.c, Q0=new_Q0, E0=new_dE)

    @staticmethod
    def fitting_func(Q: Union[float, np.array], E: float, a, b, c) -> Union[float, np.array]:
        """Model used for fitting when the minimum is fixed at Q=0."""
        return a*Q**4 + b*Q**3 + c*Q**2 + E

    @property
    def omega(self) -> float:
        """Return the quadratic coefficient used for frequency conversion."""
        return 0.5 * self.c

    def __str__(self):
        """Return a compact summary of the quartic fit."""
        return (f"QuarticCurve: {self.a}*(Q-Q0)^4 + {self.b}*(Q-Q0)^3 + "
                f"{self.c}*(Q-Q0)^2 + {self.E0} (eV), "
                f"Q0={self.Q0:.3f} (amu**0.5*Å)")


class FittingCurveType(ExtendedEnum):
    """Registry of supported fitting-function classes."""
    quadratic = ("quadratic", QuadraticFittingFunc)
    quartic = ("quartic", QuarticFittingFunc)

    # def __init__(self, _name, cls):
    #     self.name = _name
    #     self.cls = cls


def intersections(curve1: FittingFunc,
                  curve2: FittingFunc,
                  Q_range: List[float],
                  ngrids=2001) -> List[Tuple[float, float]]:
    """Return intersection points as [(Q, energy), ...] calculated numerically."""
    xs = np.linspace(Q_range[0], Q_range[1], ngrids)
    fvals = np.asarray(curve1(xs)) - np.asarray(curve2(xs))

    raw_roots = []
    for i in range(len(xs) - 1):
        y1, y2 = fvals[i], fvals[i + 1]
        if np.isclose(y1, 0.0, atol=1e-8):
            raw_roots.append(float(xs[i]))
        if y1 * y2 < 0:
            try:
                q = brentq(lambda x: float(curve1(x) - curve2(x)), xs[i], xs[i + 1])
                raw_roots.append(float(q))
            except ValueError:
                pass  # Ignore ValueError from brentq when the interval does not bracket a root

    raw_roots.sort()
    roots = []
    for r in raw_roots:
        if not any(np.isclose(r, rr, atol=1e-6) for rr in roots):
            roots.append(r)

    return [(r, float(curve1(r))) for r in roots]
