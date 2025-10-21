# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from copy import copy
from pathlib import Path

import pytest
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from pymatgen.core import Structure, Lattice
from pymatgen.electronic_structure.core import Spin

from pydefect_ccd.ccd import SinglePointResult, Ccd, PotentialCurve
from pydefect_ccd.ccd_init import CcdInit
from pydefect_ccd.ele_phon_coupling import InnerProduct, EPMatrixElement
from pydefect_ccd.enum import Carrier
from pydefect_ccd.relaxed_point import NearEdgeState, RelaxedPoint


@pytest.fixture(scope="session")
def test_files():
    return Path(__file__).parent / "test_files"


@pytest.fixture(scope="session")
def excited_structure():
    return Structure.from_str("""H He
10.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
H
6
direct 
0.1 0.0 0.0
0.9 0.0 0.0
0.0 0.1 0.0
0.0 0.9 0.0
0.0 0.0 0.1
0.0 0.0 0.9""", fmt="poscar")


@pytest.fixture(scope="session")
def ground_structure():
    return Structure.from_str("""Mg4 O4
10.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
H
6
direct
0.2 0.0 0.0
0.8 0.0 0.0
0.0 0.2 0.0
0.0 0.8 0.0
0.0 0.0 0.2
0.0 0.0 0.8""", fmt="poscar")


@pytest.fixture(scope="session")
def intermediate_structure():
    return Structure.from_str("""Mg4 O4
10.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
H
6
direct
0.15 0.0 0.0
0.85 0.0 0.0
0.0 0.15 0.0
0.0 0.85 0.0
0.0 0.0 0.15
0.0 0.0 0.85""", fmt="poscar")


vb = NearEdgeState(band_index=1,
                   kpt_coord=[0.0]*3,
                   kpt_index=1,
                   kpt_weight=1.0,
                   eigenvalue=1.0,
                   occupation=1.0)
cb = NearEdgeState(band_index=2,
                   kpt_coord=[0.0]*3,
                   kpt_index=1,
                   kpt_weight=1.0,
                   eigenvalue=3.0,
                   occupation=0.0)


@pytest.fixture
def ccd_init(ground_structure, excited_structure):
    orb_info = LocalizedOrbital(band_idx=2,
                                ave_energy=2.0,
                                occupation=1.0,
                                orbitals={"O": [0.0, 1.0, 0.0]})
    cb_w_lo = copy(cb)
    cb_w_lo.band_index = 3

    va_o1_0 = RelaxedPoint(name="Va_O1",
                           charge=0,
                           structure=ground_structure,
                           energy=11.0,
                           correction_energy=-1.0,
                           magnetization=0.0,
                           localized_orbitals=[[], []],
                           initial_site_symmetry="2mm",
                           final_site_symmetry="2",
                           parsed_dir="/path/to/Va_O1_0",
                           valence_bands=[[vb], [vb]],
                           conduction_bands=[[cb], [cb]])
    va_o1_1 = RelaxedPoint(name="Va_O1",
                           charge=1,
                           structure=excited_structure,
                           energy=12.0,
                           correction_energy=-1.0,
                           magnetization=1.0,
                           localized_orbitals=[[], [orb_info]],
                           initial_site_symmetry="2mm",
                           final_site_symmetry="2",
                           parsed_dir="/path/to/Va_O1_1",
                           valence_bands=[[vb], [vb]],
                           conduction_bands=[[cb], [cb_w_lo]])
    # transition level = -1.0 from CBM
    return CcdInit(relaxed_points=[va_o1_0, va_o1_1],
                   vbm=1.0, cbm=3.0, supercell_volume=100.0,
                   supercell_vbm=1.1, supercell_cbm=2.9,
                   ave_electron_mass=11.0, ave_hole_mass=12.0,
                   ave_static_diele_const=13.0)


@pytest.fixture
def ccd(excited_structure, ground_structure, intermediate_structure):
    return Ccd(name="test",
               potential_curve_results=[
                   PotentialCurve(name="excited", charge=0,
                                  single_points=[SinglePointResult(-1.0, -0.1, 2.1, False, used_for_fitting=True),
                                                 SinglePointResult(0.0, 0.0, 1.1, False, used_for_fitting=True, conduction_bands=[[cb]], valence_bands=[[vb]]),
                                                 SinglePointResult(1.0, 0.1, 2.2, False, used_for_fitting=True)]),
                   PotentialCurve(name="ground", charge=1,
                                  single_points=[SinglePointResult(-1.0, -0.1, 1.1, False, used_for_fitting=False),
                                                 SinglePointResult(0.0, 0.0, 0.1, False, used_for_fitting=True, conduction_bands=[[cb]], valence_bands=[[vb]]),
                                                 SinglePointResult(1.0, 0.1, 1.2, False, used_for_fitting=True)])])


@pytest.fixture(scope="session")
def sc_structure():
    lattice = Lattice.cubic(1.0)
    coords = [[0.0, 0.0, 0.0]]
    return Structure(lattice=lattice, species=["H"], coords=coords)


@pytest.fixture
def e_p_matrix_elem():
    ip_1 = InnerProduct(abs_inner_product=20.0, used_for_fitting=False)
    ip_2 = InnerProduct(abs_inner_product=1.0, used_for_fitting=True)
    ip_3 = InnerProduct(abs_inner_product=2.0, used_for_fitting=True)

    return EPMatrixElement(charge=1,
                           defect_band_index=2,
                           spin=Spin.down,
                           eigenvalue_diff=0.1,
                           kpt_idx=1,
                           inner_products={-1.0: ip_1, 0.0: ip_2, 1.0: ip_3})


@pytest.fixture
def e_p_coupling(e_p_matrix_elem):
    return EPCoupling(
        charge=1,
        base_disp=0.0,
        captured_carrier=Carrier.e,
        volume=100.0,
        ave_captured_carrier_mass=1.0,
        ave_static_diele_const=2.0,
        band_edge_index=1,
        e_p_matrix_elements=[e_p_matrix_elem])
