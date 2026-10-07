"""Two-impulse CW rendezvous solver (closed form, no optimizer).

Partition Phi(tf) = [[Phi_rr, Phi_rv], [Phi_vr, Phi_vv]] (2x2 blocks).
Requiring r(tf) = 0 gives   0 = Phi_rr r0 + Phi_rv v0_plus, hence

    v0_plus = -inv(Phi_rv) @ Phi_rr @ r0          (burn 1 target velocity)
    v_f-    =  Phi_vr @ r0 + Phi_vv @ v0_plus     (arrival velocity)
    dv1 = v0_plus - v0_minus,   dv2 = -v_f-       (burn 2 nulls arrival speed)
"""
import math
from dataclasses import dataclass

import numpy as np

from .cw_dynamics import cw_stm

# det(Phi_rv) * n^2 = 8 - 8cos(nt) - 3 nt sin(nt) = 4 sin(x) (4 sin x - 3 x cos x), x = nt/2.
# Zeros: nt = 2*pi*k, and nt where tan(nt/2) = 3 nt / 4 (first one near nt = 8.9, ~131 min at 400 km).
SINGULARITY_TOL = 1e-3


@dataclass(frozen=True)
class RendezvousSolution:
    """Burn plan; all vectors are [radial, along-track] in km/s."""

    v0_plus: np.ndarray     # relative velocity right after burn 1
    dv1: np.ndarray         # burn 1 (t = 0)
    v_arrival: np.ndarray   # relative velocity just before burn 2 (t = tf)
    dv2: np.ndarray         # burn 2 (t = tf)

    @property
    def dv1_mag(self) -> float:
        return float(np.linalg.norm(self.dv1))

    @property
    def dv2_mag(self) -> float:
        return float(np.linalg.norm(self.dv2))

    @property
    def dv_total(self) -> float:
        """Sum of impulsive burn magnitudes (idealised, not a propulsion budget)."""
        return self.dv1_mag + self.dv2_mag


def solve_rendezvous(r0, v0_minus, n: float, tf: float) -> RendezvousSolution:
    """Solve for burns that bring the chaser to the origin at time tf."""
    if not (math.isfinite(tf) and tf > 0.0):
        raise ValueError(f"Transfer time must be finite and > 0 s (got {tf!r}).")
    nt = n * tf
    det_scaled = 8.0 - 8.0 * math.cos(nt) - 3.0 * nt * math.sin(nt)
    if abs(det_scaled) < SINGULARITY_TOL:
        raise ValueError(
            f"Transfer time {tf:.1f} s makes the CW position-to-velocity map singular "
            "(n*tf is a multiple of 2*pi, or tan(n*tf/2) = 3*n*tf/4): the burn is not "
            "uniquely defined. Choose a different transfer time."
        )
    r0 = np.asarray(r0, dtype=float)
    v0_minus = np.asarray(v0_minus, dtype=float)
    phi = cw_stm(n, tf)
    phi_rr, phi_rv = phi[:2, :2], phi[:2, 2:]
    phi_vr, phi_vv = phi[2:, :2], phi[2:, 2:]
    v0_plus = -np.linalg.solve(phi_rv, phi_rr @ r0)
    v_arrival = phi_vr @ r0 + phi_vv @ v0_plus
    return RendezvousSolution(
        v0_plus=v0_plus, dv1=v0_plus - v0_minus, v_arrival=v_arrival, dv2=-v_arrival
    )
