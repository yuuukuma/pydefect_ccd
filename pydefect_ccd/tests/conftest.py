# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from copy import copy
from pathlib import Path

import pytest
from pydefect.analyzer.band_edge_states import LocalizedOrbital
from pymatgen.core import Structure, Lattice

from pydefect_ccd.ccd_init import CcdInit
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
    # The dQ between ground and excited is 0.6 amu^1/2 A
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
                   eigenvalue=1.0,
                   occupation=1.0)
cb = NearEdgeState(band_index=2,
                   kpt_coord=[0.0]*3,
                   kpt_index=1,
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
                           formation_energy=11.0,
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
                           formation_energy=12.0,
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


@pytest.fixture(scope="session")
def sc_structure():
    lattice = Lattice.cubic(1.0)
    coords = [[0.0, 0.0, 0.0]]
    return Structure(lattice=lattice, species=["H"], coords=coords)
