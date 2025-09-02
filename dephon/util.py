# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import re
from pathlib import Path
from typing import List

from numpy import sqrt, sum
from pymatgen.core import Structure
from pymatgen.electronic_structure.core import Spin


def spin_to_idx(spin: Spin, count_from_1=False) -> int:
    result = 0 if spin == Spin.up else 1
    if count_from_1:
        result += 1
    return result


def idx_to_spin(idx: int):
    return Spin.up if idx == 0 else Spin.down


def reduce_wswq(filename: Path, orbs: List[int]) -> None:
    """

    Args:
        filename: WSWQ filename
        orbs: list of relevant band indices

    Returns:
        None
    """
    lines = []
    with filename.open(mode='r') as f:
        for line in f:
            header = re.search(r'\s*spin=(\d+), kpoint=\s*(\d+)', line)
            if header:
                lines.append(line)
            else:
                _append_ij(line, orbs, lines)

    filename.write_text("".join(lines))


def _append_ij(line, orbitals, out):
    data = re.search(r'i=\s*(\d+), '
                     r'' r'j=\s*(\d+)\s*:\s*([0-9\-.]+)\s+([0-9\-.]+)',
                     line)
    i, j = int(data.group(1)), int(data.group(2))
    if i in orbitals and j in orbitals:
        out.append(line)


def get_dR(ground: Structure, excited: Structure) -> float:
    """Summation of atomic displacement distances

    Args:
        ground (Structure): Reference structure
        excited (Structure): Target structure

    A constant offset should be removed, for example, by aligning the farthest
    atom.

    Returns:
        The Summed atomic displacement distance in float
    """
    return sqrt(sum([x.distance(y) ** 2 for x, y in zip(ground, excited)]))
