# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from copy import copy
from math import sqrt

import pytest
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from pymatgen.core import Lattice, Structure, Element
from vise.tests.helpers.assertion import assert_json_roundtrip

from dephon.relaxed_point import BandEdgeState, RelaxedPoint
from dephon.util import get_dR


def test_get_dR():
    lattice = Lattice.orthorhombic(10, 20, 30)
    structure1 = Structure(lattice, ["H"], [[0.0, 0.0, 0.0]])
    structure2 = Structure(lattice, ["H"], [[0.1, 0.1, 0.1]])
    assert get_dR(structure1, structure2) == sqrt(1**2+2**2+3**2)


vb_nes = BandEdgeState(band_index=1,
                       kpt_coord=[0.0]*3,
                       kpt_weight=1.0,
                       kpt_index=1,
                       eigenvalue=1.5,
                       occupation=1.0)


def test_near_edge_state_str():
    print(vb_nes)


@pytest.fixture
def minimum_point_info(ground_structure):
    orb_info = LocalizedOrbital(band_idx=2,
                                ave_energy=2.0,
                                occupation=1.0,
                                orbitals={"O": [0.0, 1.0, 0.0]})
    cb_nes_up = BandEdgeState(band_index=2,
                              kpt_coord=[0.0]*3,
                              kpt_weight=1.0,
                              kpt_index=1,
                              eigenvalue=2.5,
                              occupation=0.0)
    cb_nes_down = copy(cb_nes_up)
    cb_nes_down.band_index = 3
    return RelaxedPoint(name="test",
                        charge=-1,
                        structure=ground_structure,
                        energy=-90.0,
                        correction_energy=1.0,
                        magnetization=1.0,
                        localized_orbitals=[[], [orb_info]],
                        valence_bands=[[vb_nes]],
                        conduction_bands=[[cb_nes_up, cb_nes_down]],
                        initial_site_symmetry="4mm",
                        final_site_symmetry="2/m",
                        parsed_dir="/path/to/min_point")


def test_relevant_band_indices(minimum_point_info):
    assert minimum_point_info.relevant_band_indices == {1, 2, 3}


def test_str(minimum_point_info):
    print(minimum_point_info)


def test_json_roundtrip(dephon_init, tmpdir):
    assert_json_roundtrip(dephon_init, tmpdir)


def test_minimum_point_info_degeneracy_by_symm_reduction(minimum_point_info):
    assert minimum_point_info.degeneracy_by_symmetry_reduction == 2


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
