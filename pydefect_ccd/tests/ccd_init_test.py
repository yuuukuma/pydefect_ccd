# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from math import sqrt

import pytest
from pymatgen.core import Element
from vise.tests.helpers.assertion import assert_json_roundtrip


def test_json_roundtrip(ccd_init, tmpdir):
    assert_json_roundtrip(ccd_init, tmpdir)


def test_dQ(ccd_init):
    expected = sqrt((0.1*10)**2*6 * Element.H.atomic_mass)
    assert ccd_init.Q == pytest.approx(expected)


def test_dR(ccd_init):
    expected = sqrt((0.1*10)**2*6)
    assert ccd_init.dR == pytest.approx(expected)


def test_modal_mass(ccd_init):
    assert ccd_init.modal_mass == pytest.approx(Element.H.atomic_mass)


def test_relaxed_point_from_charge(ccd_init):
    actual = ccd_init.relaxed_point_from_charge(charge=1)
    assert actual == ccd_init.relaxed_points[1]


def test_volume(ccd_init):
    assert ccd_init.volume == 1000.0


def test_ccd_string(ccd_init):
    print(ccd_init)
