# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.

import pytest
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from pymatgen.electronic_structure.core import Spin

from pydefect_ccd.ccd import PotentialCurve, SinglePoint, CcdId
from pydefect_ccd.ele_phon_coupling import InnerProduct, EPMatrixElement
from pydefect_ccd.enum import Carrier
from pydefect_ccd.make_e_p_matrix_element import MakeEPMatrixElement
from pydefect_ccd.relaxed_point import NearEdgeState


# @pytest.fixture
# def ccd_init():
#     lattice = Lattice.cubic(2.0)
#     coords = [[0.0, 0.0, 0.0]]
#     vbm = BandEdgeState(band_index=101, kpt_coord=[0.0] * 3, kpt_weight=1.0, kpt_index=1, eigenvalue=1000.0, occupation=0.0)
#     excited = MinimumPointInfo(charge=0,
#                                structure=Structure(lattice, species=["H"], coords=coords),
#                                energy=1.0,
#                                correction_energy=1.0,
#                                magnetization=-1.0,
#                                localized_orbitals=[[]],
#                                initial_site_symmetry="1",
#                                final_site_symmetry="1",
#                                parsed_dir="/",
#                                vbm=[vbm, vbm],
#                                cbm=[])
#     ground = deepcopy(excited)
#     ground.charge = -1
#     return DephonInit(defect_name="test",
#                       min_points=[ground, excited],
#                       vbm=0.0, cbm=10.0,
#                       supercell_vbm=1.0,
#                       supercell_cbm=9.0,
#                       ave_electron_mass=0.01,
#                       ave_hole_mass=0.02,
#                       ave_static_diele_const=0.03)
#

@pytest.fixture
def single_ccd():
    cbm = NearEdgeState(band_index=104,
                        kpt_coord=[0.0]*3,
                        kpt_weight=1.0,
                        kpt_index=1,
                        eigenvalue=7.0,
                        occupation=0.0)
    l_orb_upper = LocalizedOrbital(
        band_idx=103, ave_energy=4.0, occupation=0.0, orbitals={})
    l_orb_lower = LocalizedOrbital(
        band_idx=102, ave_energy=2.0, occupation=1.0, orbitals={})
    vbm = NearEdgeState(band_index=101,
                        kpt_coord=[0.0]*3,
                        kpt_weight=1.0,
                        kpt_index=1,
                        eigenvalue=1.0,
                        occupation=1.0)

    single_point_info = \
        SinglePoint(dQ=1.0,
                    disp_ratio=0.0,
                    magnetization=-1.0,
                    localized_orbitals=[[l_orb_lower, l_orb_upper], []],
                    valence_bands=[[vbm], []],
                    conduction_bands=[[cbm], []])
    return PotentialCurve(CcdId("excited", carriers=[Carrier.e]),
                          charge=0,
                          single_points=[single_point_info])


@pytest.fixture
def make_e_p_coupling_h_capture(single_ccd):
    return MakeEPMatrixElement(base_disp_ratio=0.0,
                               potential_curve=single_ccd,
                               band_edge_index=101,
                               defect_band_index=102,
                               kpoint_index=1,
                               spin=Spin.up,
                               dQ_wswq_pairs=[(0.0, {(1, 1): {(101, 102): 3.0 + 4.0j}}),
                                              (0.1, {(1, 1): {(101, 102): 30.0 + 40.0j}})])


def test_make_e_p_matrix_element(make_e_p_coupling_h_capture):
    actual = make_e_p_coupling_h_capture.make()
    expected = EPMatrixElement(charge=0,
                               base_disp_ratio=0.0,
                               band_edge_index=101,
                               defect_band_index=102,
                               spin=Spin.up,
                               eigenvalue_diff=1.0,
                               kpt_idx=1,
                               inner_products={0.0: InnerProduct(0.0),
                                               0.1: InnerProduct(50.0)})
    assert actual == expected


