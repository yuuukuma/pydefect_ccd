# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from typing import List, Optional

from matplotlib import pyplot as plt
from pydefect.analyzer.band_edge_states import BandEdgeOrbitalInfos
from vise.util.logger import get_logger
from vise.util.matplotlib import float_to_int_formatter

logger = get_logger(__name__)


class DephonEigenvaluePlotter:
    def __init__(self,
                 orb_infos: List[BandEdgeOrbitalInfos],
                 disp_ratios: List[float],
                 supercell_vbm: float = None,
                 supercell_cbm: float = None,
                 title: str = None,
                 y_range: Optional[List[float]] = None,
                 ):
        try:
            num_kpt = len(orb_infos[0].kpt_weights)
            assert num_kpt == 1
        except AssertionError:
            logger.warning(f"Currently {self.__class__} supports only single "
                           f"k-point.")
            raise

        try:
            first_num_kpt_coords = orb_infos[0].kpt_coords
            for orb_info in orb_infos:
                assert first_num_kpt_coords == orb_info.kpt_coords
        except AssertionError:
            logger.warning(f"Numbers of k-points are not the same.")

        self._title = title or ""
        self._orb_infos = orb_infos
        self._disp_ratios = disp_ratios
        self._supercell_vbm, self._supercell_cbm = supercell_vbm, supercell_cbm
        self.plt = plt
        num_spin = len(self._orb_infos[0].orbital_infos)
        self.fig, self.axs = self.plt.subplots(1, num_spin)
        if num_spin == 1:
            self.axs = [self.axs]
        self._y_range = y_range

    def construct_plot(self):
        self._add_eigenvalues()
        if self._supercell_vbm and self._supercell_cbm:
            self._add_vbm_cbm()
        self._set_title()
        self._set_formatter()
        self._set_labels()
        if self._y_range:
            print(self._y_range)
            self._set_y_range()
        self.plt.tight_layout()

    def _add_eigenvalues(self):
        num_spin = len(self._orb_infos[0].orbital_infos)

        for spin_idx in range(num_spin):
            for orb_infos, disp_ratio in zip(self._orb_infos, self._disp_ratios):
                orb_info = orb_infos.orbital_infos[spin_idx][0]
                xs = [disp_ratio] * len(orb_info)
                energies = [oi.energy for oi in orb_info]
                occupations = [oi.occupation for oi in orb_info]
                self.axs[spin_idx].scatter(xs, energies, c=occupations,
                                      cmap='viridis')
        # self.plt.colorbar(ax=axs[0])

    def _add_vbm_cbm(self):
        args = dict(linestyle="-.", linewidth=0.75)
        for ax in self.axs:
            ax.axhline(y=self._supercell_vbm, **args)
            ax.axhline(y=self._supercell_cbm, **args)

    def _set_labels(self):
        self.fig.text(0.55, 0, "Displacement ratio", ha='center')
        self.axs[0].set_ylabel(f"Energy (eV)")

    def _set_title(self):
        self.plt.gca().set_title(self._title)

    def _set_y_range(self):
        for ax in self.axs:
            ax.set_ylim(self._y_range[0], self._y_range[1])

    def _set_formatter(self):
        self.plt.gca().xaxis.set_major_formatter(float_to_int_formatter)
        self.plt.gca().yaxis.set_major_formatter(float_to_int_formatter)

