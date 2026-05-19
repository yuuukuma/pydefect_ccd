# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.

import pytest
from pydefect.analyzer.band_edge_states import BandEdgeOrbitalInfos, OrbitalInfo
from pymatgen.electronic_structure.core import Spin

from pydefect_ccd.e_p_matrix_element import EPMatrixElement
from pydefect_ccd.make_e_p_matrix_element import make_e_p_matrix_element


@pytest.fixture
def base_band_edge_orbital_infos():
    band_edge = OrbitalInfo(energy=1.0, occupation=1.0, orbitals={})
    defect = OrbitalInfo(energy=2.0, occupation=1.0, orbitals={})
    return BandEdgeOrbitalInfos(orbital_infos=[[[band_edge, defect]], [[]]],
                                kpt_coords=[[0.0, 0.0, 0.0]],
                                kpt_weights=[1.0],
                                lowest_band_index=100,
                                fermi_level=0.0)


def test_make_e_p_matrix_element(base_band_edge_orbital_infos):
    actual = make_e_p_matrix_element(
        charge=0,
        base_band_edge_orbital_infos=base_band_edge_orbital_infos,
        band_edge_index=101,
        defect_band_index=102,
        spin=Spin.up,
        dQs=[0.0, 0.1],
        wswqs=[3.0 + 4.0j,  30.0 + 40.0j])

    expected = EPMatrixElement(charge=0,
                               band_edge_index=101,
                               defect_band_index=102,
                               spin=Spin.up,
                               eigenvalue_diff=1.0,
                               dQs=[0.0, 0.1],
                               abs_inner_prods=[0.0, 50.0])
    assert actual == expected

