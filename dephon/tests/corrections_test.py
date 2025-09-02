# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import pytest
from vise.tests.helpers.assertion import assert_msonable, assert_yaml_roundtrip

from dephon.corrections import DephonCorrection
from dephon.enum import CorrectionType


@pytest.fixture
def dephon_correction():
    return DephonCorrection(corrections={CorrectionType.extended_FNV: 1.0,
                                         CorrectionType.kumagai2023: 2.0})


def test_correction_energy(dephon_correction):
    assert dephon_correction.total_correction_energy == 3.0


def test_correction_msonable(dephon_correction):
    assert_msonable(dephon_correction)


def test_correction_yaml_round_trip(dephon_correction, tmpdir):
    expected_str = f"""extended_FNV: 1.0
kumagai2023: 2.0
"""
    print(dephon_correction.to_yaml())
    assert_yaml_roundtrip(dephon_correction, tmpdir, expected_str)

