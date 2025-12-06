# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass, field
from typing import Union, Dict, List, Optional

import numpy as np
from monty.json import MSONable
from nonrad.scaling import sommerfeld_parameter
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
    base_disp_ratio: float
    band_edge_index: int
    defect_band_index: int
    spin: Union[Spin, str]
    eigenvalue_diff: float
    # key dQ, val: |\bra_{psi_i(0)} | S(0) |\ket_{psi_f(Q)}|
    dQs: List[float] = field(default_factory=list)
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
                   ["W_if_tilde (eV/amu^0.5) ", self.W_if_tilde]]

        result.append(tabulate(table_1, tablefmt="plain"))
        result.append(tabulate(table_2, tablefmt="plain", floatfmt=".3f"))
        inner_prods = [[dQ, aip] for dQ, aip in zip(self.dQs, self.abs_inner_prods)]

        result.append(tabulate(inner_prods,
                               headers=["dQ", "inner product",
                                        "used for fitting?"],
                               tablefmt="plain", floatfmt=".3f"))

        return "\n".join(result)


@dataclass
class EPCoupling(MSONable, ToJsonFileMixIn):
    """ E-P coupling constants between defect band and multiple band edges

    f(T) * V * W_if^2

    To define the localized orbitals, we need to determine the charge and
    displacement (base_disp).

    Attributes:
        ave_captured_carrier_mass:
            This is needed to calculate the Sommerfeld parameter and the
            thermal velocity.
        ave_static_diele_const:
            This is needed to calculate the Sommerfeld parameter.

    """
    W_if_tilde: float
    charge: int
    T: Union[float, np.ndarray]
    volume: float
    ave_captured_carrier_mass: float = None
    ave_static_diele_const: float = None
    uniform_scaling_factor: float = 1.0

    def __str__(self):
        # todo: update
        result = []
        if self.ave_captured_carrier_mass is not None:
            mass = round(self.ave_captured_carrier_mass, 2)
        else:
            mass = "N/A"

        if self.ave_static_diele_const is not None:
            diele_const = round(self.ave_static_diele_const, 2)
        else:
            diele_const = "N/A"

        table = [["charge", self.charge],
                 ["volume", round(self.volume, 2)],
                 ["ave. carrier mass", mass],
                 ["ave. static dielectric constant", diele_const]]
        result.append(tabulate(table, tablefmt="plain"))
        header = ["T (K)", "W_if ()"]
        table = [[t, w] for t, w in zip(self.T, self.W_if)]
        result.append(tabulate(table, headers=header, tablefmt="plain"))
        return "\n".join(result)

    @property
    def f(self):
        return self.uniform_scaling_factor * self.sommerfeld_scaling_factor

    @property
    def sommerfeld_scaling_factor(self) -> Union[float, np.ndarray]:
        if self.charge == 0:
            return 1.0
        return sommerfeld_parameter(self.T,
                                    self.charge,
                                    self.ave_captured_carrier_mass,
                                    self.ave_static_diele_const)

    @property
    def W_if(self) -> List[float]:
        """ E-P coupling constant W_if """
        return [self.f * self.volume * self.W_if_tilde] * len(self.T)

