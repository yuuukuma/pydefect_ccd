# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass, field
from typing import Union, Dict, List, Optional

import numpy as np
from monty.json import MSONable
from pymatgen.electronic_structure.core import Spin
from tabulate import tabulate
from vise.util.mix_in import ToJsonFileMixIn


@dataclass
class EPMatrixElement(MSONable, ToJsonFileMixIn):
    """ Electron-phonon (E-P) matrix element between defect band and band edge

    The phonon mode is only 1D mode.

    Attributes:
        band_edge_index: The band index for band edges starting from 1.
        eigenvalue_diff: The eigenvalue difference from final to initial states
        abs_inner_prods: List of <psi_i(0) | S(0) |psi_f(Q)>
            at the given Q points.
    """
    charge: int
    band_edge_index: int
    defect_band_index: int
    spin: Union[Spin, str]
    eigenvalue_diff: float
    dQs: List[float] = field(default_factory=list)  # Å x amu^0.5
    # |\bra_{psi_i(0)} | S(0) |\ket_{psi_f(Q)}|
    abs_inner_prods: List[float] = field(default_factory=list)
    # TODO : implement fit_q_range
    fit_q_range: List[float] = None

    def __post_init__(self):
        if isinstance(self.spin, str):
            self.spin = Spin[self.spin]
        try:
            self.grad, self.const = np.polyfit(self.dQs, self.abs_inner_prods, 1)
        except:
            print(f"Cannot fit inner products vs dQ.")
            self.grad, self.const = None, None

    @property
    def W_if_tilde(self) -> float:
        """ E-P matrix element W_if tilde """
        return self.grad * self.eigenvalue_diff

    @property
    def _json_filename(self):
        return self._filename + "_" + self.index_info + ".json"

    def to_json_file(self, filename: Optional[str] = None) -> None:
        if filename is None:
            filename = self._json_filename
        super().to_json_file(filename)

    @property
    def index_info(self):
        return "_".join([f"b{self.band_edge_index}",
                         f"d{self.defect_band_index}",
                         str(self.spin)])

    def as_dict(self) -> dict:
        result = super().as_dict()
        result["spin"] = result["spin"].name
        return result

    def plot(self, ax=None) -> None:
        """ Evaluated by computing the slope of inner products"""
        ax.scatter(self.dQs, self.abs_inner_prods)
        x = np.arange(min(self.dQs), max(self.dQs), 0.01)
        y = x * self.grad + self.const
        ax.plot(x, y, alpha=0.5)

    def __str__(self):
        result = []
        table_1 = [["band edge index", self.band_edge_index],
                   ["defect band index", self.defect_band_index],
                   ["spin", self.spin.name]]
        table_2 = [["grad (amu^0.5)", self.grad],
                   ["eigenvalue difference (eV)", self.eigenvalue_diff],
                   ["W_if_tilde (eV/(Å amu^0.5))", self.W_if_tilde]]

        result.append(tabulate(table_1, tablefmt="plain"))
        result.append(tabulate(table_2, tablefmt="plain", floatfmt=".3f"))
        inner_prods = [[dQ, aip] for dQ, aip in zip(self.dQs, self.abs_inner_prods)]

        result.append(tabulate(inner_prods,
                               headers=["dQ", "inner product",
                                        "used for fitting?"],
                               tablefmt="plain", floatfmt=".3f"))

        return "\n".join(result)

    @property
    def to_W_if_tilde(self) -> "WifTilde":
        return WifTilde(W_if_tilde=self.W_if_tilde,
                        band_edge_index=self.band_edge_index,
                        charge=self.charge)


@dataclass
class WifTilde(MSONable, ToJsonFileMixIn):
    W_if_tilde: float  # eV / Angstrom / amu^0.5
    band_edge_index: int
    charge: int
    uniform_scaling_factor: float = 1.0

