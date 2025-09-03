# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass, field
from typing import Union, Dict

import numpy as np
from monty.json import MSONable
from pymatgen.electronic_structure.core import Spin
from tabulate import tabulate
from vise.util.mix_in import ToJsonFileMixIn

from pydefect_ccd.util import spin_to_idx


@dataclass
class InnerProduct(MSONable):
    """
    Attributes:
        inner_products: \bra_{psi_i(0)} | S(0) |\ket_{psi_f(Q)}
        used_for_fitting:
    """
    abs_inner_product: float
    used_for_fitting: bool = True


@dataclass
class EPMatrixElement(MSONable, ToJsonFileMixIn):
    """ Electron-phonon (E-P) matrix element between defect band and band edge

    The phonon mode is only 1D mode.

    Attributes:
        band_edge_index: The band index for band edges starting from 1.
        eigenvalue_diff: The eigenvalue difference from final to initial states
        inner_products: List of \bra_{psi_i(0)} | S(0) |\ket_{psi_f(Q)}
            at the given Q points.
    """
    charge: int
    base_disp_ratio: float
    band_edge_index: int
    defect_band_index: int
    spin: Union[Spin, str]
    eigenvalue_diff: float
    kpt_idx: int
#    kpt_weigh: float
    # Currently, symmetry is assumed not to be changed depending on dQ.
#    kpt_coord: List[float]
    # key dQ:
    inner_products: Dict[float, InnerProduct] = field(default_factory=dict)

    def __post_init__(self):
        if isinstance(self.spin, str):
            self.spin = Spin[self.spin]
        self.inner_products = \
            {float(k): v for k, v in self.inner_products.items()}

    @property
    def _json_filename(self):
        return self._filename + "_" + self.index_info + ".json"

    @property
    def index_info(self):
        return "_".join([f"b{self.band_edge_index}",
                         f"d{self.defect_band_index}",
                         f"k{self.kpt_idx}",
                         str(self.spin)])

    def as_dict(self) -> dict:
        result = super().as_dict()
        result["spin"] = result["spin"].name
        return result

    @property
    def spin_idx(self):
        return spin_to_idx(self.spin)

    @property
    def dQs(self):
        return [dQ for dQ in self.inner_products.keys()]

    @property
    def inner_prods(self):
        return [ip.abs_inner_product for ip in self.inner_products.values()]

    @property
    def _inner_prod_vs_q(self):
        return self.dQs, self.inner_prods

    def e_p_matrix_element(self, ax=None) -> float:
        """ Evaluated by computing the slope of inner products"""
        grad, const = np.polyfit(self.dQs, self.inner_prods, 1)

        if ax:
            ax.scatter(self.dQs, self.inner_prods)

            x = np.arange(min(self.dQs), max(self.dQs), 0.01)
            y = x * grad + const
            ax.plot(x, y, alpha=0.5)

        return grad

    def __str__(self):
        result = []
        table = [["band edge index", self.band_edge_index],
                 ["defect band index", self.defect_band_index],
                 ["spin", self.spin.name],
                 ["eigenvalue difference", round(self.eigenvalue_diff, 3)],
                 ["kpoint index", self.kpt_idx],
                 ["e-p matrix element", self.e_p_matrix_element()]]

        result.append(tabulate(table, tablefmt="plain", floatfmt=".3f"))

        inner_prods = []
        for dQ, ip in self.inner_products.items():
            inner_prods.append([dQ, ip.abs_inner_product, ip.used_for_fitting])

        result.append(tabulate(inner_prods,
                               headers=["dQ", "inner product",
                                        "used for fitting?"],
                               tablefmt="plain", floatfmt=".3f"))

        return "\n".join(result)

#
# @dataclass
# class EPCoupling(MSONable, ToJsonFileMixIn):
#     """ E-P coupling constants between defect band and multiple band edges
#
#     To define the localized orbitals, we need to determine the charge and
#     displacement (base_disp).
#
#     Attributes:
#         ave_captured_carrier_mass:
#             This is needed to calculate the Sommerfeld parameter and the
#             thermal velocity.
#         ave_static_diele_const:
#             This is needed to calculate the Sommerfeld parameter.
#
#     """
#     base_disp: float
#     volume: float
#     ave_captured_carrier_mass: float
#     ave_static_diele_const: float
#     # TODO: In the future, consider multiple host and defect bands.
#     e_p_matrix_element: Optional[EPMatrixElement] = None
#
#     def reset_inner_prod(self):
#         self.e_p_matrix_element = None
#
#     def __str__(self):
#         result = []
#         mass = round(self.ave_captured_carrier_mass, 2)
#         diele_const = round(self.ave_static_diele_const, 2)
#         table = [["charge", self.charge],
#                  ["base disp", self.base_disp],
#                  ["captured carrier", self.captured_carrier],
#                  ["volume", round(self.volume, 2)],
#                  ["averaged carrier mass", mass],
#                  ["averaged static dielectric constant", diele_const]]
#         result.append(tabulate(table, tablefmt="plain", floatfmt=".3f"))
#
#         if self.e_p_matrix_element:
#             result.append("-" * 50)
#             result.append(self.e_p_matrix_element.__str__())
#
#         return "\n".join(result)
#
#     @property
#     def wif(self):
#         return self.e_p_matrix_element.e_p_matrix_element()
