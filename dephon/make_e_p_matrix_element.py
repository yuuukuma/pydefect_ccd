# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from typing import Dict, Tuple, Optional, List

import numpy as np
from pymatgen.electronic_structure.core import Spin
from vise.util.logger import get_logger

from dephon.config_coord import PotentialCurve
from dephon.ele_phon_coupling import EPMatrixElement, InnerProduct
from dephon.util import spin_to_idx

logger = get_logger(__name__)


wswq_type = Dict[Optional[Tuple[int, int]], Dict[Tuple[int, int], complex]]


class MakeEPMatrixElement:
    def __init__(self,
                 base_disp_ratio: float,
                 single_ccd: PotentialCurve,
                 band_edge_index: int,
                 defect_band_index: int,
                 kpoint_index: int,
                 spin: Spin,
                 dQ_wswq_pairs: List[Tuple[float, wswq_type]],
                 energy_diff: float = None
                 ):
        self.charge = single_ccd.charge
        self.base_disp_ratio = base_disp_ratio
        self.band_edge_index = band_edge_index
        self.defect_band_index = defect_band_index

        single_point_info = single_ccd.get_single_point_by_disp_ratio(base_disp_ratio)
        if energy_diff:
            self.energy_diff = energy_diff
        else:
            band_edge_state = single_point_info.band_edge_state(spin, band_edge_index)
            defect_state = single_point_info.localized_orbital(spin, defect_band_index)
            self.energy_diff = abs(band_edge_state.eigenvalue - defect_state.ave_energy)

        self.spin = spin
        self.kpt_index = kpoint_index
        self.dQ_wswq_pairs = dQ_wswq_pairs

    def make(self):
        inner_prods = {dQ: self._add_inner_products(wswq, dQ)
                       for dQ, wswq in self.dQ_wswq_pairs}
        return EPMatrixElement(charge=self.charge,
                               base_disp_ratio=self.base_disp_ratio,
                               band_edge_index=self.band_edge_index,
                               defect_band_index=self.defect_band_index,
                               spin=self.spin,
                               eigenvalue_diff=self.energy_diff,
                               kpt_idx=self.kpt_index,
                               inner_products=inner_prods)

    def _add_inner_products(self, wswq: wswq_type, dQ: float):
        """    Returns
        -------
        dict(dict)
            a dict of dicts that takes keys (spin, kpoint) and (initial, final)
            as indices and maps it to a complex number"""
        spin_kpt_pair = (spin_to_idx(self.spin, True), self.kpt_index)
        band_indices = (self.band_edge_index, self.defect_band_index)
        braket = np.abs(wswq[spin_kpt_pair][tuple(band_indices)]) * np.sign(dQ)
        return InnerProduct(abs_inner_product=braket)



