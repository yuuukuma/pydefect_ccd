# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from math import sqrt
from pathlib import Path

from pymatgen.core import Lattice, Structure
from pymatgen.electronic_structure.core import Spin

from pydefect_ccd.util import spin_to_idx, idx_to_spin, reduce_wswq, get_dR


def test_spin_to_idx():
    assert spin_to_idx(Spin.up) == 0
    assert spin_to_idx(Spin.down) == 1

    assert spin_to_idx(Spin.up, count_from_1=True) == 1
    assert spin_to_idx(Spin.down, count_from_1=True) == 2


def test_idx_to_spin():
    assert idx_to_spin(0) == Spin.up
    assert idx_to_spin(1) == Spin.down


def test_reduce_wswq(tmpdir):
    tmpdir.chdir()
    text = """   spin=1, kpoint=     1
i=     1, j=     1 :     0.000000001     0.000000000
i=     2, j=     1 :     0.000000001     0.000000000
i=     3, j=     1 :     0.000000001     0.000000000
i=     4, j=     1 :     0.000000001     0.000000000
i=     1, j=     2 :     0.000000001     0.000000000
i=     2, j=     2 :     0.000000001     0.000000000
i=     3, j=     2 :     0.000000001     0.000000000
i=     4, j=     2 :     0.000000001     0.000000000
i=     1, j=     3 :     0.000000001     0.000000000
i=     2, j=     3 :     0.000000001     0.000000000
i=     3, j=     3 :     0.000000001     0.000000000
i=     4, j=     3 :     0.000000001     0.000000000
i=     1, j=     4 :     0.000000001     0.000000000
i=     2, j=     4 :     0.000000001     0.000000000
i=     3, j=     4 :     0.000000001     0.000000000
i=     4, j=     4 :     0.000000001     0.000000000
   spin=2, kpoint=     1
i=     1, j=     1 :     0.000000101     0.000000000
i=     2, j=     1 :     0.000000101     0.000000000
i=     3, j=     1 :     0.000000101     0.000000000
i=     4, j=     1 :     0.000000101     0.000000000
i=     1, j=     2 :     0.000000101     0.000000000
i=     2, j=     2 :     0.000000101     0.000000000
i=     3, j=     2 :     0.000000101     0.000000000
i=     4, j=     2 :     0.000000101     0.000000000
i=     1, j=     3 :     0.000000101     0.000000000
i=     2, j=     3 :     0.000000101     0.000000000
i=     3, j=     3 :     0.000000101     0.000000000
i=     4, j=     3 :     0.000000101     0.000000000
i=     1, j=     4 :     0.000000101     0.000000000
i=     2, j=     4 :     0.000000101     0.000000000
i=     3, j=     4 :     0.000000101     0.000000000
i=     4, j=     4 :     0.000000101     0.000000000
"""
    Path("WSWQ").write_text(text)
    reduce_wswq(filename=Path("WSWQ"), orbs=[2, 3])
    actual = Path("WSWQ").read_text()

    expected = """   spin=1, kpoint=     1
i=     2, j=     2 :     0.000000001     0.000000000
i=     3, j=     2 :     0.000000001     0.000000000
i=     2, j=     3 :     0.000000001     0.000000000
i=     3, j=     3 :     0.000000001     0.000000000
   spin=2, kpoint=     1
i=     2, j=     2 :     0.000000101     0.000000000
i=     3, j=     2 :     0.000000101     0.000000000
i=     2, j=     3 :     0.000000101     0.000000000
i=     3, j=     3 :     0.000000101     0.000000000
"""
    assert actual == expected


def test_get_dR():
    lattice = Lattice.orthorhombic(10, 20, 30)
    structure1 = Structure(lattice, ["H"], [[0.0, 0.0, 0.0]])
    structure2 = Structure(lattice, ["H"], [[0.1, 0.1, 0.1]])
    assert get_dR(structure1, structure2) == sqrt(1**2+2**2+3**2)
