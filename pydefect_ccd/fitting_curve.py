# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
import inspect
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Tuple, Union, List

import numpy as np
from monty.json import MSONable
from nonrad.constants import HBAR, EV2J, ANGS2M, AMU2KG
from scipy.optimize import curve_fit, brentq
from vise.util.enum import ExtendedEnum

from pydefect_ccd.potential_curve import SinglePoints


@dataclass
class FittingCurve(ABC):
    """Abstract base class for fitting curves used in CCD analysis.

    Note that the class needs to be registered in FittingCurveType to be used.
    For the fitting Q0=0 is assumed and this must be included in the fitting points.

    Attributes:
    """
    Q0: float
    dE: float

    @abstractmethod
    def __call__(self, x: Union[float, np.array]) -> Union[float, np.array]:
        pass

    @abstractmethod
    def shift(self, shift_Q, shift_energy, revert=False) -> "FittingCurve":
        pass

    @staticmethod
    @abstractmethod
    def fitting_func(Q: Union[float, np.array], dE: float, *params) -> Union[float, np.array]:
        """Model function for fitting."""
        pass

    @classmethod
    def from_single_points(cls, single_points: SinglePoints):
        # Q0 is assumed to be fixed at 0.0, so only model parameters are fitted.
        # print(single_points.Qs, single_points.corrected_energies)

        single_points.verify_Q0_has_the_lowest_energy()
        single_points.verify_num_Q(cls.fitting_func)

        vals, _ = curve_fit(cls.fitting_func,
                            single_points.Qs,
                            single_points.corrected_energies)
        # vals, _ = curve_fit(f, self.Qs, self.corrected_energies, bounds=bounds)
        dE = vals[0]
        param_names = list(inspect.signature(cls.fitting_func).parameters.keys())[2:]
        kwargs = {'Q0': 0.0, 'dE': dE}
        for i, name in enumerate(param_names):
            kwargs[name] = vals[i + 1]
        return cls(**kwargs)


    @classmethod
    def n_fit_params(cls) -> int:
        sig = inspect.signature(cls.fitting_func)
        return len(sig.parameters) - 1

    def add_plot(self, ax, x_range: List[float], color):
        xs = np.linspace(x_range[0], x_range[1], 1000)
        ys = self(xs)
        ax.plot(xs, ys, color=color)

    @property
    def omega_in_eV(self) -> float:
        if not hasattr(self, "omega"):
            raise AttributeError(f"{self.__class__} does not have 'omega' attribute.")

        return HBAR * self.omega * np.sqrt(EV2J / (ANGS2M**2 * AMU2KG))


@dataclass
class QuadraticFittingCurve(MSONable, FittingCurve):
    """
    a(Q-Q0)^2 + dE = 0.5 * omega^2 * (Q-Q0)^2 + dE
    """
    a: float
    # omega: float  # in amu Å^2 / eV

    def __call__(self, Q: Union[float, np.ndarray]) -> Union[float, np.array]:
        return self.a * (Q - self.Q0)**2 + self.dE

    def shift(self, shift_Q, shift_energy, revert=False) -> "QuadraticFittingCurve":
        new_Q0 = self.Q0 + shift_Q
        new_dE = self.dE + shift_energy
        return QuadraticFittingCurve(new_Q0, new_dE, a=self.a)

    @staticmethod
    def fitting_func(Q: Union[float, np.array], dE: float, a) -> Union[float, np.array]:
        return a*Q**2 + dE

    @property
    def omega(self) -> float:
        return 0.5 * self.a

    def __str__(self):
         return (f"Quadratic Curve: omega={self.omega_in_eV:.3f} (eV), "
                 f"Q0={self.Q0:.3f} (amu**0.5*Å), Emin={self.dE:.3f} (eV)")


@dataclass
class QuarticFittingCurve(MSONable, FittingCurve):
    """ a(Q-Q0)^4 + b(Q-Q0)^3 + c(Q-Q0)^2 + d(Q-Q0) + dE """
    a: float  # in ??
    b: float
    c: float

    def __call__(self, Q: Union[float, np.array]) -> Union[float, np.array]:
        return (self.a * (Q - self.Q0) ** 4 + self.b * (Q - self.Q0) ** 3
                + self.c * (Q - self.Q0) ** 2 + self.dE)

    def shift(self, shift_Q, shift_energy, revert=False) -> "QuarticFittingCurve":
        new_Q0 = self.Q0 + shift_Q
        new_dE = self.dE + shift_energy
        new_b = -self.b if revert else self.b
        return QuarticFittingCurve(a=self.a, b=new_b, c=self.c, Q0=new_Q0, dE=new_dE)

    @staticmethod
    def fitting_func(Q: Union[float, np.array], dE: float, a, b, c) -> Union[float, np.array]:
        return a*Q**4 + b*Q**3 + c*Q**2 + dE

    @property
    def omega(self) -> float:
        return 0.5 * self.c

    def __str__(self):
        return (f"QuarticCurve: {self.a}*(Q-Q0)^4 + {self.b}*(Q-Q0)^3 + "
                f"{self.c}*(Q-Q0)^2 + {self.dE} (eV), "
                f"Q0={self.Q0:.3f} (amu**0.5*Å)")


class FittingCurveType(ExtendedEnum):
    quadratic = ("quadratic", QuadraticFittingCurve)
    quartic = ("quartic", QuarticFittingCurve)

    # def __init__(self, _name, cls):
    #     self.name = _name
    #     self.cls = cls


def intersections(curve1: FittingCurve,
                  curve2: FittingCurve,
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
