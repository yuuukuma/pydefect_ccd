# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from typing import List

import numpy as np
from pymatgen.electronic_structure.core import Spin
from vise.util.logger import get_logger

from pydefect_ccd.ccd import SinglePoint
from pydefect_ccd.ele_phon_coupling import EPMatrixElement

logger = get_logger(__name__)


def make_ep_matrix_element(name: str,
                           base_single_point: SinglePoint,
                           band_edge_index: int,
                           defect_band_index: int,
                           spin: Spin,
                           dQs: List[float],
                           wswqs: List[complex],
                           energy_diff: float = None) -> EPMatrixElement:
    assert len(dQs) == len(wswqs)

    base_disp_ratio = base_single_point.disp_ratio

    if energy_diff:
        eigenvalue_diff = energy_diff
    else:
        band_edge_state = base_single_point.near_edge_state(spin, band_edge_index)
        defect_state = base_single_point.localized_orbital(spin, defect_band_index)
        eigenvalue_diff = abs(band_edge_state.eigenvalue - defect_state.ave_energy)

    abs_inner_prods = [float(np.abs(wswq) * np.sign(dQ))
                       for dQ, wswq in zip(dQs, wswqs)]

    return EPMatrixElement(name=name,
                           base_disp_ratio=base_disp_ratio,
                           band_edge_index=band_edge_index,
                           defect_band_index=defect_band_index,
                           spin=spin,
                           eigenvalue_diff=eigenvalue_diff,
                           dQs=dQs,
                           abs_inner_prods=abs_inner_prods)
