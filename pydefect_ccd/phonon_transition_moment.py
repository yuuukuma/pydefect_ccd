# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from dataclasses import dataclass
from typing import Union, List, Literal

import numpy as np
from monty.json import MSONable
from nonrad import get_C
# from qmsolve import Eigenstates
from scipy.constants import constants, physical_constants
from tabulate import tabulate
from vise.util.mix_in import ToJsonFileMixIn

from pydefect_ccd.ccd import Ccd
from pydefect_ccd.fitting_curve import QuadraticFittingCurve

# @dataclass
# class PhononEigenstates(MSONable, ToJsonFileMixIn):
#     """Info about the eigenstates"""
#     energies: List[float]
#     eigenvalues: np.array
#     extent: float
#     N: int
#
#     @property
#     def num_eigenstates(self):
#         return len(self.energies)
#
#     @property
#     def array(self):
#         return self.eigenvalues
#
#     def __str__(self):
#         result = [f"Extent (Å): {self.extent}",
#                   f"Number of grids: {self.N}",
#                   tabulate(self.energies, headers=["Energy (eV)"],
#                            tablefmt="plain", floatfmt=".3f")]
#         return "\n".join(result)
#
#
# def eigen2phonon_eigen(eigen: Eigenstates) -> PhononEigenstates:
#     """
#     Units of Eigenstates:
#        energies: hartree
#        eigenstates: a.u.^(-1/2)
#
#     Args:
#         eigen:
#
#     Returns:
#
#     """
#     return PhononEigenstates(energies=[e*HARTREE2EV for e in eigen.energies],
#                              eigenvalues=eigen.array * AUL2ANGSTROM ** -0.5,
#                              extent=eigen.extent,
#                              N=eigen.N)
#
#
# @dataclass
# class PhononOverlap:
#     initial_eigenstates: Eigenstates  # excited state
#     final_eigenstates: Eigenstates  # ground state
#     dQ: float  # in Å
#     T: float  # in K
#     occ_tol: float = 1e-4
#
#     # @property
#     # def kT(self):
#     #     return (const.k / const.e) * self.T    # [(J / K) * (eV / J)] * K = eV
#
#     def phonon_overlaps_by_n(self, n):
#         return np.abs(np.conj(self.initial_eigenstates.array[n]) * self.final_eigenstates.array)
#
#     def initial_occupancies(self):
#         pass
#
#     def smearing_factor(self, dE):
#         pass
#
#     def get_R(self) -> Union[float, np.ndarray]:
#         result = 0.0
#         for n in range(self.initial_eigenstates.number):
#             result += sum(self.phonon_overlaps_by_n(n) * self.initial_occupancies())
#         return result * self.dQ / 10**10
#
#
