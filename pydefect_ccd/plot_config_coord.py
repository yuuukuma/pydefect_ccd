# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.

from matplotlib import pyplot as plt
from vise.util.logger import get_logger
from vise.util.matplotlib import float_to_int_formatter

from pydefect_ccd.ccd import Ccd

logger = get_logger(__name__)


class CcdPlotter:
    """Plot a configuration coordinate diagram with both charge states."""

    def __init__(self,
                 ccd: Ccd,
                 title: str = None,
                 set_energy_zero: bool = True):
        """Set CCD data and plotting options."""
        self._title = title or ""
        self._ccd = ccd
        self._set_energy_zero = set_energy_zero
        self.plt = plt

    def construct_plot(self):
        """Draw the CCD and apply plot formatting."""
        self._add_ccd()
        self._set_title()
        self._set_formatter()
        self._set_labels()
        self.plt.tight_layout()

    def _add_ccd(self):
        """Add ground and excited curves to the current axes."""
        ax = self.plt.gca()
        self._ccd.ground_curve.add_plot(ax, "black")
        self._ccd.excited_curve.add_plot(ax, "black")

    def _set_labels(self):
        """Set axis labels and legend for the CCD plot."""
        ax = self.plt.gca()
        ax.set_xlabel("Q (amu$^{1/2}$ Å)")
        ax.set_ylabel("Energy (eV)")
        ax.legend()

    def _set_title(self):
        """Set the plot title."""
        self.plt.gca().set_title(self._title)

    def _set_formatter(self):
        """Use compact integer-style tick labels when possible."""
        self.plt.gca().xaxis.set_major_formatter(float_to_int_formatter)
        self.plt.gca().yaxis.set_major_formatter(float_to_int_formatter)
