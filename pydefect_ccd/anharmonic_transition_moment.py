# -*- coding: utf-8 -*-
#  Copyright (c) 2026 Kumagai group.
from typing import Callable

import numpy as np
from scipy.linalg import eigh_tridiagonal


def solve_1d_phonon_schrodinger(potential_func: Callable,
                                Q_grid,
                                n_eigs=100):
    """
    配置座標 Q 上の 1 次元 Schrödinger 方程式を有限差分で解く。

    方程式:
        [ - hbar_sq_over_2 * d^2/dQ^2 + V(Q) ] chi(Q) = E chi(Q)

    Parameters
    ----------
    potential_func : callable
        ポテンシャル関数 V(Q)。
        入力は numpy 配列 Q_grid、出力は同じ shape の配列。
    Q_grid : array-like
        等間隔の 1 次元グリッド
    n_eigs : int
        求める最低固有状態の数

    Returns
    -------
    energies : ndarray, shape (n_eigs,)
        固有値
    wavefunctions : ndarray, shape (len(Q_grid), n_eigs)
        正規化済み固有関数
    """
    Q_grid = np.asarray(Q_grid, dtype=float)

    if Q_grid.ndim != 1 or len(Q_grid) < 3:
        raise ValueError("Q_grid must be a 1D array with at least 3 points.")

    dQs = np.diff(Q_grid)
    if not np.allclose(dQs, dQs[0], rtol=1e-10, atol=1e-12):
        raise ValueError("Q_grid must be uniformly spaced.")

    dQ = dQs[0]  # spacing
    N = len(Q_grid)

    V = np.asarray(potential_func(Q_grid), dtype=float)
    t = 1.0 / (dQ ** 2)

    main_diag = 2.0 * t + V
    off_diag = -t * np.ones(N - 1)

    n_eigs = min(n_eigs, N)
    energies, wavefunctions = eigh_tridiagonal(main_diag,
                                               off_diag,
                                               select="i",
                                               select_range=(0, n_eigs - 1))

    for i in range(wavefunctions.shape[1]):
        norm = np.sqrt(np.trapezoid(np.abs(wavefunctions[:, i]) ** 2, Q_grid))
        wavefunctions[:, i] /= norm

    return energies, wavefunctions