# import pytest
# from matplotlib import pyplot as plt
# from qmsolve import eV, Å, init_visualization, SingleParticle
# from qmsolve.visualization.single_particle_1D import \
#     VisualizationSingleParticle1D
#
# from dephon.phonon_overlap import Potential
#
#
# @pytest.fixture
# def pot():
#     return Potential(potential_energies=[0.0, 0.25, 1.0, 2.25, 4.0],
#                      dQs=[0.0, 0.5, 1.0, 1.5, 2.0])
#
#
# # def test_potential_func(pot):
# #     func = pot.potential_func()
# #     particle = SingleParticle()
# #     particle.x = 1.0
# #     assert func(particle) == 1.0
# #     particle.x = 2.0
# #     assert func(particle) == 4.0
#
#
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