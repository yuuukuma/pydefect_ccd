# -*- coding: utf-8 -*-
#  Copyright (c) 2024 Kumagai group.
from copy import copy

import pytest
from pydefect.analyzer.band_edge_states import LocalizedOrbital

from pydefect_ccd.relaxed_point import NearEdgeState, RelaxedPoint

vb_nes = NearEdgeState(band_index=1,
                       kpt_coord=[0.0]*3,
                       kpt_weight=1.0,
                       kpt_index=1,
                       eigenvalue=1.5,
                       occupation=1.0)


def test_near_edge_state_str():
    print(vb_nes)


@pytest.fixture
def relaxed_point(ground_structure):
    orb_info = LocalizedOrbital(band_idx=2,
                                ave_energy=2.0,
                                occupation=1.0,
                                orbitals={"O": [0.0, 1.0, 0.0]})
    cb_nes_up = NearEdgeState(band_index=2,
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


def test_relevant_band_indices(relaxed_point):
    assert relaxed_point.relevant_band_indices == {1, 2, 3}
