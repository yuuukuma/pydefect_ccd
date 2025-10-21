# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.

from matplotlib import pyplot as plt
from vise.util.logger import get_logger
from vise.util.matplotlib import float_to_int_formatter

from pydefect_ccd.ccd import Ccd

logger = get_logger(__name__)


class CcdPlotter:
    def __init__(self, ccd: Ccd,
                 title: str = None,
                 set_energy_zero: bool = True):
        self._title = title or ""
        self._ccd = ccd
        self._set_energy_zero = set_energy_zero
        self.plt = plt

    def construct_plot(self):
        self._add_ccd()
        self._set_title()
        self._set_formatter()
        self._set_labels()
        self.plt.tight_layout()

    def _add_ccd(self):
        ax = self.plt.gca()
        for imag_infos in self._ccd.potential_curve_results:
            imag_infos.add_plot(ax, "black")

    def _set_labels(self):
        ax = self.plt.gca()
        ax.set_xlabel("Q (amu$^{1/2}$ Å)")
        ax.set_ylabel("Energy (eV)")
        ax.legend()

    def _set_title(self):
        self.plt.gca().set_title(self._title)

    def _set_formatter(self):
        self.plt.gca().xaxis.set_major_formatter(float_to_int_formatter)
        self.plt.gca().yaxis.set_major_formatter(float_to_int_formatter)

