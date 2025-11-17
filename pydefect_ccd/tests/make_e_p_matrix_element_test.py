# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.

import pytest
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from pymatgen.electronic_structure.core import Spin

from pydefect_ccd.ccd import SinglePoint, SinglePointSpec
from pydefect_ccd.ele_phon_coupling import InnerProduct, EPMatrixElement
from pydefect_ccd.make_e_p_matrix_element import MakeEPMatrixElement
from pydefect_ccd.relaxed_point import NearEdgeState


@pytest.fixture
def single_point():
    cbm = NearEdgeState(band_index=104,
                        kpt_coord=[0.0]*3,
                        kpt_index=1,
                        eigenvalue=7.0,
                        occupation=0.0)
    l_orb_upper = LocalizedOrbital(
        band_idx=103, ave_energy=4.0, occupation=0.0, orbitals={})
    l_orb_lower = LocalizedOrbital(
        band_idx=102, ave_energy=2.0, occupation=1.0, orbitals={})
    vbm = NearEdgeState(band_index=101,
                        kpt_coord=[0.0] * 3,
                        kpt_index=1,
                        eigenvalue=1.0,
                        occupation=1.0)

    return SinglePoint(SinglePointSpec(dQ=1.0, disp_ratio=0.0),
                       energy=10.0,
                       magnetization=-1.0,
                       localized_orbitals=[[l_orb_lower, l_orb_upper], []],
                       valence_bands=[[vbm], []],
                       conduction_bands=[[cbm], []],
                       ccd_correction_energy=0.0)


def test_make_e_p_matrix_element(single_point):
    maker = MakeEPMatrixElement(
        name="test",
        charge=0,
        base_single_point=single_point,
        band_edge_index=101,
        defect_band_index=102,
        kpoint_index=1,
        spin=Spin.up,
        dQ_wswq_pairs=[(0.0, {(1, 1): {(101, 102): 3.0 + 4.0j}}),
                       (0.1, {(1, 1): {(101, 102): 30.0 + 40.0j}})])

    actual = maker.make()
    expected = EPMatrixElement(name="test",
                               charge=0,
                               base_disp_ratio=0.0,
                               band_edge_index=101,
                               defect_band_index=102,
                               spin=Spin.up,
                               eigenvalue_diff=1.0,
                               kpt_idx=1,
                               abs_inner_products={0.0: InnerProduct(0.0),
                                                   0.1: InnerProduct(50.0)})
    assert actual == expected


