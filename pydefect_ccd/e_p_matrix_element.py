# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
from dataclasses import dataclass, field
from typing import Union, Dict, List, Optional

import numpy as np
from monty.json import MSONable
from pymatgen.electronic_structure.core import Spin
from tabulate import tabulate
from vise.util.mix_in import ToJsonFileMixIn


@dataclass
class EPMatrixElement(MSONable, ToJsonFileMixIn):
    """Electron-phonon matrix element between a defect band and a band edge.

    The stored inner products are sampled along the one-dimensional
    configuration coordinate.  A linear fit against ``dQ`` gives the gradient
    used to compute ``W_if_tilde``.

    Attributes:
        charge: Charge state used for the matrix element.
        band_edge_index: Band-edge band index, starting from 1.
        defect_band_index: Defect-localized band index, starting from 1.
        spin: Spin channel for both bands.
        eigenvalue_diff: Eigenvalue difference between final and initial
            electronic states in eV.
        dQs: Configuration-coordinate displacements in amu^0.5 Angstrom.
        abs_inner_prods: Absolute overlap-like inner products at each ``dQ``.
        fit_q_range: Optional Q range for selecting fitting points.
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
        """Normalize spin input and fit the inner-product slope."""
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
        """Return the default JSON filename including band and spin indices."""
        return self._filename + "_" + self.index_info + ".json"

    def to_json_file(self, filename: Optional[str] = None) -> None:
        """Write the matrix element to JSON."""
        if filename is None:
            filename = self._json_filename
        super().to_json_file(filename)

    @property
    def index_info(self):
        """Return a compact band/spin identifier used in filenames."""
        return "_".join([f"b{self.band_edge_index}",
                         f"d{self.defect_band_index}",
                         str(self.spin)])

    def as_dict(self) -> dict:
        """Return a MSON dictionary with spin serialized by name."""
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
        """Return a tabulated summary of the fitted matrix element."""
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
        """Convert to the smaller ``WifTilde`` transfer object."""
        return WifTilde(W_if_tilde=self.W_if_tilde,
                        band_edge_index=self.band_edge_index,
                        charge=self.charge)


@dataclass
class WifTilde(MSONable, ToJsonFileMixIn):
    """Scaled electron-phonon matrix element for a band and charge state."""
    W_if_tilde: float  # eV / Angstrom / amu^0.5
    band_edge_index: int
    charge: int
    uniform_scaling_factor: float = 1.0
