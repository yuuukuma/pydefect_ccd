# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.

import pytest
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from pymatgen.electronic_structure.core import Spin

from pydefect_ccd.ccd import SinglePoint, SinglePointSpec
from pydefect_ccd.ele_phon_coupling import EPMatrixElement
from pydefect_ccd.make_e_p_matrix_element import make_e_p_matrix_element
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
    actual = make_e_p_matrix_element(
        name="test",
        base_single_point=single_point,
        band_edge_index=101,
        defect_band_index=102,
        spin=Spin.up,
        dQs=[0.0, 0.1],
        wswqs=[3.0 + 4.0j,  30.0 + 40.0j])

    expected = EPMatrixElement(name="test",
                               base_disp_ratio=0.0,
                               band_edge_index=101,
                               defect_band_index=102,
                               spin=Spin.up,
                               eigenvalue_diff=1.0,
                               dQs=[0.0, 0.1],
                               abs_inner_prods=[0.0, 50.0])
    assert actual == expected


