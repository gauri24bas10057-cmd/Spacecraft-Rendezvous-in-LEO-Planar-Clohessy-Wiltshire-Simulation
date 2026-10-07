"""Cross-check of the CW burn plan against full planar two-body motion.

Both spacecraft are propagated in an inertial frame under point-mass gravity
with the CW burn plan applied; the result is converted back to the target's
radial/along-track frame. This quantifies the linearisation error only.
"""
import numpy as np
from scipy.integrate import solve_ivp

from . import config
from .config import EARTH_MU_KM3_S2 as MU
from .config import Scenario
from .orbital import circular_speed_km_s, mean_motion_rad_s, orbit_radius_km
from .rendezvous import solve_rendezvous


def _two_body(_t, s):
    r1, r2 = s[0:2], s[4:6]
    a1 = -MU * r1 / np.linalg.norm(r1) ** 3
    a2 = -MU * r2 / np.linalg.norm(r2) ** 3
    return np.concatenate([s[2:4], a1, s[6:8], a2])


def nonlinear_miss(sc: Scenario) -> tuple:
    """Return (miss distance [km], relative speed [km/s]) at tf in true two-body motion."""
    n = mean_motion_rad_s(sc.altitude_km)
    r, vc = orbit_radius_km(sc.altitude_km), circular_speed_km_s(sc.altitude_km)
    sol = solve_rendezvous([sc.x0_km, sc.y0_km], [0, 0], n, sc.transfer_time_s)
    x0, y0 = sc.x0_km, sc.y0_km
    vx, vy = sol.v0_plus
    s0 = np.array([r, 0, 0, vc,                       # target
                   r + x0, y0, vx - n * y0, vc + vy + n * x0])  # chaser (inertial)
    out = solve_ivp(_two_body, (0, sc.transfer_time_s), s0, method=config.INTEGRATOR_METHOD,
                    rtol=1e-11, atol=1e-12)
    s = out.y[:, -1]
    rt, vt, rc, vcv = s[0:2], s[2:4], s[4:6], s[6:8]
    er = rt / np.linalg.norm(rt)
    et = np.array([-er[1], er[0]])
    rho, dv = rc - rt, vcv - vt
    om = (rt[0] * vt[1] - rt[1] * vt[0]) / np.dot(rt, rt)      # target angular rate
    x, y = rho @ er, rho @ et
    vxr, vyr = dv @ er + om * y, dv @ et - om * x             # rotating-frame rates
    return float(np.hypot(x, y)), float(np.hypot(vxr, vyr))
