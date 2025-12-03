# -*- coding: utf-8 -*-
#  Copyright (c) 2022 Kumagai group.
import numpy as np
import pytest
from nonrad import get_C


def test_summed_squared_transition_moment():
    dQ = 1.26630 # amu^(1/2) Angstrom
    dE = 0.915 # eV
    wi = 0.04582 # eV
    wf = 0.04166 # eV
    actual = get_C(dQ, dE, wi, wf, T=300, Wif=1.0, volume=1.0, g=1)
    expected = 4.41-11
    pytest.approx(actual, abs=1e-2)


# def test_extent(pot):
#     assert pot.extent == 2.0
# #
# #
# # def test_add_plot(pot):
# #     ax = plt.gca()
# #     pot.add_plot(ax, color="red", q_range=[-1.1, 2.5])
# #     plt.show()
# #
# #
# # def test_make_eigenstates(pot):
# #     actual = pot.make_eigenstates()
# #     vis = VisualizationSingleParticle1D(actual)
# #     vis.plot_eigenstate(k=0)
# #     vis.slider_plot()
# #     print(actual)
#
#
#
# """
# TODO:
# - Check if numpy msonable works for PhononEigenstates
#
# 1. Plot potential profile.
# 2.
#
#
# """
#
# #    state = pot.make_eigenstates()
#
# # def test_harmonic():
# #     def harmonic_oscillator(particle, dQ=0.0):
# #         k = 1. * eV / Å**2
# #         return 2.0 * k * particle.x**2
#
#     # actual = make_eigenstates(harmonic_oscillator, 50.0)
#     # # Visualize the Eigenstates
#     # visualization = init_visualization(actual)
#     # visualization.animate()
#     # print(actual.energies)