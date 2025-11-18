# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import numpy as np
import pytest
from matplotlib import pyplot as plt
from numpy.testing import assert_almost_equal
from pymatgen.electronic_structure.core import Spin
from vise.tests.helpers.assertion import assert_json_roundtrip

from pydefect_ccd.ele_phon_coupling import EPMatrixElement, EPCoupling


@pytest.fixture
def e_p_matrix_elem():
    return EPMatrixElement(name="ep_matrix_elem_test",
                           base_disp_ratio=0.0,
                           band_edge_index=1,
                           defect_band_index=2,
                           spin=Spin.down,
                           eigenvalue_diff=0.1,
                           abs_inner_prods={-1.0: 20.0, 0.0: 1.0, 1.0: 2.0})


@pytest.fixture
def e_p_coupling(e_p_matrix_elem):
    return EPCoupling(e_p_matrix_elements=[1.0],
                      charge=1,
                      T=np.array([300.0]),
                      volume=100.0,
                      ave_captured_carrier_mass=1.0,
                      ave_static_diele_const=2.0)


def test_e_p_coupling_to_json_file(e_p_coupling, tmpdir):
    assert_json_roundtrip(e_p_coupling, tmpdir)


def test_inner_prod_vs_q(e_p_matrix_elem):
    assert e_p_matrix_elem.dQs == (-1.0, 0.0, 1.0)
    assert e_p_matrix_elem.abs_inner_prods == (20.0, 1.0, 2.0)


def test_e_p_matrix_element_less_inner_prod(e_p_matrix_elem):
    e_p_matrix_elem.abs_inner_products = {}
    assert e_p_matrix_elem.e_p_matrix_element() is None

    e_p_matrix_elem.abs_inner_products = {0.0: 0.2}
    assert e_p_matrix_elem.e_p_matrix_element() is None


def test_e_p_matrix_element(e_p_matrix_elem):
    e_p_matrix_elem.abs_inner_products = {0.0: 0.2, 1.0: 1.2}
    ax = plt.gca()
    assert_almost_equal(e_p_matrix_elem.e_p_matrix_element(ax), 1.0)
    plt.show()


def test_e_p_matrix_element(e_p_coupling):
    print(e_p_coupling)

