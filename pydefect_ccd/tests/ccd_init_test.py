# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from math import sqrt

import pytest
from pymatgen.core import Element
from vise.tests.helpers.assertion import assert_json_roundtrip


def test_json_roundtrip(dephon_init, tmpdir):
    assert_json_roundtrip(dephon_init, tmpdir)


def test_minimum_point_info_degeneracy_by_symm_reduction(relaxed_point):
    assert relaxed_point.degeneracy_by_symmetry_reduction == 2


def test_dephon_init_dQ(dephon_init):
    expected = sqrt((0.1*10)**2*6 * Element.H.atomic_mass)
    assert dephon_init.dQ == pytest.approx(expected)


def test_dephon_init_min_point_info_from_charge(dephon_init):
    actual = dephon_init.relaxed_point_info_from_charge(charge=1)
    assert actual == dephon_init.relaxed_points[1]


def test_dephon_init_volume(dephon_init):
    assert dephon_init.volume == 1000.0


def test_ccd_string(dephon_init):
    print(dephon_init)
