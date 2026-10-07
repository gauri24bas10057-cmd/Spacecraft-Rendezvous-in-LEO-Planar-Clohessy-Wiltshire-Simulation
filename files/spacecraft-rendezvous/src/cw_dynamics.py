"""Planar Clohessy-Wiltshire (Hill) relative dynamics.

Frame: origin on the target; x radial (outward), y along-track (direction of
motion). State s = [x, y, x_dot, y_dot] in km and km/s.

    x_ddot =  3 n^2 x + 2 n y_dot
    y_ddot = -2 n x_dot
"""
import numpy as np


def cw_dynamics(t: float, state, n: float) -> list:
    """Right-hand side ds/dt of the planar CW equations (for solve_ivp)."""
    x, _y, x_dot, y_dot = state
    return [x_dot, y_dot, 3.0 * n**2 * x + 2.0 * n * y_dot, -2.0 * n * x_dot]


def cw_matrix_A(n: float) -> np.ndarray:
    """Constant matrix A of the linear system ds/dt = A s."""
    return np.array(
        [
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [3.0 * n**2, 0.0, 0.0, 2.0 * n],
            [0.0, 0.0, -2.0 * n, 0.0],
        ]
    )


def cw_stm(n: float, t: float) -> np.ndarray:
    """Closed-form state-transition matrix Phi(t) = expm(A t), s(t) = Phi s(0)."""
    c, s, nt = np.cos(n * t), np.sin(n * t), n * t
    return np.array(
        [
            [4.0 - 3.0 * c, 0.0, s / n, 2.0 * (1.0 - c) / n],
            [6.0 * (s - nt), 1.0, -2.0 * (1.0 - c) / n, (4.0 * s - 3.0 * nt) / n],
            [3.0 * n * s, 0.0, c, 2.0 * s],
            [6.0 * n * (c - 1.0), 0.0, -2.0 * s, 4.0 * c - 3.0],
        ]
    )
