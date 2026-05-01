# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
import pytest
from pymatgen.electronic_structure.core import Spin
from vise.tests.helpers.assertion import assert_json_roundtrip

from pydefect_ccd.e_p_matrix_element import EPMatrixElement, WifTilde


@pytest.fixture
def e_p_matrix_elem():
    return EPMatrixElement(charge=0,
                           band_edge_index=1,
                           defect_band_index=2,
                           spin=Spin.down,
                           eigenvalue_diff=0.1,
                           dQs=[-1.0, 0.0, 1.0],
                           abs_inner_prods=[20.0, 1.0, 2.0],
                           fit_q_range=[-1.0, 1.0])



def test_e_p_coupling_to_json_file(e_p_matrix_elem, tmpdir):
    assert_json_roundtrip(e_p_matrix_elem, tmpdir)
    print(e_p_matrix_elem)

