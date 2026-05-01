# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
from dataclasses import dataclass
from typing import List, Optional

import numpy as np
from monty.json import MSONable
from scipy import constants
from tabulate import tabulate
from vise.util.matplotlib import float_to_int_formatter
from vise.util.mix_in import ToJsonFileMixIn


@dataclass
class CaptureRate(MSONable, ToJsonFileMixIn):
    """Nonradiative capture rate evaluated over a temperature grid.

    The supercell volume is stored in Angstrom cubed, while the reported
    capture rate is converted to cm^3/s.
    """
    Ts: List[float]
    volume: float  # in Å^3, which is the volume of the supercell.
    sommerfeld_parameter: List[float] # as a function of T
    W_if_tilde: float  # in eV / (amu^0.5 Å)
    total_squared_transition_moment: List[float]  # as a function of T
    site_degeneracy: float
    # velocities: Optional[List[float]] = None # characteristic carrier velocity in [cm / s], which also depends on T

    # @property
    # def cross_sections(self) -> List[float]:
    #     result = list(self.capture_rate / np.array(self.velocities))
    #     return result

    @property
    def capture_rate(self) -> np.ndarray:
        """Return capture coefficients for each temperature in cm^3/s."""
        volume_cm3 = self.volume * 1e-24  # convert from Å^3 to cm^3
        hbar_eVs = constants.hbar / constants.e  # eV s
        print(2 * np.pi / hbar_eVs * volume_cm3 * self.site_degeneracy
              * self.W_if_tilde ** 2)
        return (2 * np.pi / hbar_eVs * volume_cm3 * self.site_degeneracy
                * self.W_if_tilde ** 2
                * np.array(self.sommerfeld_parameter)
                * np.array(self.total_squared_transition_moment))

    def __str__(self):
        """Return a tabulated summary of inputs and capture coefficients."""
        header = [["site degeneracy:", f"{self.site_degeneracy}"],
                  ["volume (Å^3):", f"{self.volume}"],
                  ["W_if_tilde (eV / (amu$^{0.5}$ Å)):", f"{self.W_if_tilde:.2e}"]]
        result = [tabulate(header, tablefmt="plain")]

        columns = ["T (K)", "Sommerfeld (-)", "Total moment (amu Å^2 / eV)", "C (cm^3/s)"]
        table = [[T, sommerfeld, total_moment, C] for
                 T, sommerfeld, total_moment, C,
                 in zip(self.Ts,
                        self.sommerfeld_parameter,
                        self.total_squared_transition_moment,
                        self.capture_rate)]
        fmt = [".1f", ".1e", ".1e", ".1e"]

        result.append(tabulate(table, headers=columns, tablefmt="plain", floatfmt=fmt))

        return "\n".join(result)


class CaptureRatePlotter:
    """Plot temperature-dependent capture coefficients."""

    def __init__(self, capture_rate: CaptureRate, plt, title: str = None):
        """Set capture-rate data, plotting backend, and optional title."""
        self._title = title or ""
        self._capture_rate = capture_rate
        self.plt = plt

    def construct_plot(self):
        """Build the capture-rate plot on the current axes."""
        self._add_capture_rate()
        self._set_title()
        self._set_formatter()
        self._set_labels()
        self.plt.tight_layout()

    def _add_capture_rate(self):
        """Add capture coefficients as a semilog curve."""
        ax = self.plt.gca()
        ax.semilogy(self._capture_rate.Ts, self._capture_rate.capture_rate)

    def _set_labels(self):
        """Set axis labels and legend for the plot."""
        ax = self.plt.gca()
        ax.set_xlabel("T (K)")
        ax.set_ylabel("C$_p$ (cm$^3$/s)")
        ax.legend()

    def _set_title(self):
        """Set the plot title."""
        self.plt.gca().set_title(self._title)

    def _set_formatter(self):
        """Use compact integer-style tick labels when possible."""
        self.plt.gca().xaxis.set_major_formatter(float_to_int_formatter)
        self.plt.gca().yaxis.set_major_formatter(float_to_int_formatter)
