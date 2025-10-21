# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
from copy import deepcopy

from pydefect.analyzer.band_edge_states import BandEdgeOrbitalInfos, OrbitalInfo

from pydefect_ccd.plot_eigenvalues import PotCurveEigenvalPlotter


def test_eigenvalue_plotter(mocker):
    i1 = mocker.Mock(spec=BandEdgeOrbitalInfos, autospec=True)

    i1.kpt_weights = [1.0]
    i1.kpt_coords = [[0.0]*3]
    i1.orbital_infos = [
        [[OrbitalInfo(-1.1, {}, 1.0),
          OrbitalInfo(-1.0, {}, 1.0),
          OrbitalInfo(0.0, {}, 1.0),
          OrbitalInfo(1.0, {}, 0.1),
          OrbitalInfo(2.0, {}, 0.0)]],
        [[OrbitalInfo(-1.1, {}, 1.0),
          OrbitalInfo(-1.0, {}, 1.0),
          OrbitalInfo(0.0, {}, 1.0),
          OrbitalInfo(1.0, {}, 0.1),
          OrbitalInfo(2.0, {}, 0.0)]]]

    i2 = deepcopy(i1)
    i2.orbital_infos[0][0][2].energy = 0.1

    i3 = deepcopy(i1)
    i3.orbital_infos[0][0][2].energy = 0.2

    i4 = deepcopy(i1)
    i4.orbital_infos[0][0][2].energy = 0.3

    plotter = PotCurveEigenvalPlotter(orb_infos=[i1, i2, i3, i4],
                                      disp_ratios=[-0.1, 0.0, 0.4, 0.8],
                                      supercell_vbm=-1.0,
                                      supercell_cbm=1.0)
    plotter.construct_plot()
    plotter.plt.savefig("test.png")
    # plotter.plt.show()


"""
TODO
1. support spin
2. raise assertion error for inappropriate calc

"""
