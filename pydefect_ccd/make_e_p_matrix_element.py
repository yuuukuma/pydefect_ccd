# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from typing import List, Optional

import numpy as np
from pydefect.analyzer.band_edge_states import BandEdgeOrbitalInfos
from pymatgen.electronic_structure.core import Spin
from vise.util.logger import get_logger

from pydefect_ccd.e_p_matrix_element import EPMatrixElement
from pydefect_ccd.util import spin_to_idx

logger = get_logger(__name__)


#TODO: eigenvalue_diff class to handle eigenvalue correction.

def make_e_p_matrix_element(charge: int,
                            base_band_edge_orbital_infos: BandEdgeOrbitalInfos,
                            band_edge_index: int,
                            defect_band_index: int,
                            spin: Spin,
                            dQs: List[float],
                            wswqs: List[complex],
                            eigenvalue_diff: Optional[float] = None) -> EPMatrixElement:
    spin_idx = spin_to_idx(spin)
    band_idx = band_edge_index - base_band_edge_orbital_infos.lowest_band_index - 1
    defect_idx = defect_band_index - base_band_edge_orbital_infos.lowest_band_index - 1
    band_orbital_info = base_band_edge_orbital_infos.orbital_infos[spin_idx][0][band_idx]
    defect_orbital_info = base_band_edge_orbital_infos.orbital_infos[spin_idx][0][defect_idx]
    eigenvalue_diff = eigenvalue_diff or defect_orbital_info.energy - band_orbital_info.energy

    abs_inner_prods = [float(np.abs(wswq) * np.sign(dQ))
                       for dQ, wswq in zip(dQs, wswqs)]

    return EPMatrixElement(charge=charge,
                           band_edge_index=band_edge_index,
                           defect_band_index=defect_band_index,
                           spin=spin,
                           eigenvalue_diff=eigenvalue_diff,
                           dQs=dQs,
                           abs_inner_prods=abs_inner_prods)