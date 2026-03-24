# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import numpy as np
import pytest
from pymatgen.electronic_structure.core import Spin
from vise.tests.helpers.assertion import assert_json_roundtrip

from pydefect_ccd.e_p_matrix_element import EPMatrixElement, EPCoupling, WifTilde


@pytest.fixture
def e_p_matrix_elem():
    return EPMatrixElement(charge=0,
                           base_disp_ratio=0.0,
                           band_edge_index=1,
                           defect_band_index=2,
                           spin=Spin.down,
                           eigenvalue_diff=0.1,
                           dQs=[-1.0, 0.0, 1.0],
                           abs_inner_prods=[20.0, 1.0, 2.0])


@pytest.fixture
def e_p_coupling(e_p_matrix_elem):
    return EPCoupling(W_if_tilde=[WifTilde(1.0, band_edge_index=1,charge=0)],
                      T=np.array([300.0]),
                      ave_captured_carrier_mass=1.0,
                      ave_static_diele_const=2.0)


def test_e_p_coupling_to_json_file(e_p_coupling, tmpdir):
    assert_json_roundtrip(e_p_coupling, tmpdir)


def test_e_p_matrix_element(e_p_coupling):
    print(e_p_coupling)

