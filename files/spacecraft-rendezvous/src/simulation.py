"""Propagate one rendezvous scenario and measure the outcome."""
import math
from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from . import config
from .config import Scenario
from .cw_dynamics import cw_dynamics, cw_stm
from .orbital import mean_motion_rad_s, orbital_period_s
from .rendezvous import RendezvousSolution, solve_rendezvous

LINEARISATION_WARN_KM = 10.0   # CW is a small-separation approximation

STATUS_OK = "RENDEZVOUS ACHIEVED"
STATUS_POS_ONLY = "POSITION MATCH ONLY (relative velocity not nulled)"
STATUS_FAIL = "RENDEZVOUS FAILED"


@dataclass
class SimulationResult:
    scenario: Scenario
    n_rad_s: float
    period_s: float
    solution: RendezvousSolution
    apply_terminal_burn: bool
    t_s: np.ndarray                 # (N,) time samples
    states: np.ndarray              # (N, 4) coast trajectory [x, y, xd, yd]
    separation_km: np.ndarray       # d(t)
    rel_speed_km_s: np.ndarray      # v_rel(t) along the coast arc
    final_separation_km: float
    arrival_speed_km_s: float       # |v| just before burn 2
    final_speed_km_s: float         # |v| after burn 2 (or arrival if no burn 2)
    status: str
    stm_vs_ode_pos_err_km: float
    stm_vs_ode_vel_err_km_s: float
    warnings: list = field(default_factory=list)


def validate_scenario(sc: Scenario, n_steps: int) -> list:
    """Raise ValueError on invalid input; return a list of soft warnings."""
    for name in ("altitude_km", "x0_km", "y0_km", "transfer_time_s"):
        if not math.isfinite(getattr(sc, name)):
            raise ValueError(f"{name} must be a finite number.")
    if sc.altitude_km <= 0:
        raise ValueError("Target altitude must be > 0 km.")
    if sc.transfer_time_s <= 0:
        raise ValueError("Transfer time must be > 0 s.")
    if n_steps < 2:
        raise ValueError("Need at least 2 output samples.")
    sep = math.hypot(sc.x0_km, sc.y0_km)
    warns = []
    if sep > LINEARISATION_WARN_KM:
        warns.append(
            f"Initial separation {sep:.1f} km exceeds {LINEARISATION_WARN_KM:.0f} km: "
            "the linearised CW model becomes inaccurate."
        )
    if sc.transfer_time_s > orbital_period_s(sc.altitude_km):
        warns.append("Transfer time exceeds one orbital period (multi-revolution).")
    return warns


def run_simulation(sc: Scenario, n_steps: int = config.DEFAULT_N_STEPS,
                   apply_terminal_burn: bool = True) -> SimulationResult:
    """Solve burns, integrate the CW ODE numerically, and report the outcome.

    The chaser starts co-orbital (zero relative velocity before burn 1).
    """
    warns = validate_scenario(sc, n_steps)
    n = mean_motion_rad_s(sc.altitude_km)
    tf = sc.transfer_time_s
    r0 = np.array([sc.x0_km, sc.y0_km])
    sol = solve_rendezvous(r0, np.zeros(2), n, tf)

    s0 = np.concatenate([r0, sol.v0_plus])
    ode = solve_ivp(cw_dynamics, (0.0, tf), s0, args=(n,),
                    t_eval=np.linspace(0.0, tf, n_steps),
                    method=config.INTEGRATOR_METHOD,
                    rtol=config.INTEGRATOR_RTOL, atol=config.INTEGRATOR_ATOL)
    if not ode.success:
        raise RuntimeError(f"Integration failed: {ode.message}")
    states = ode.y.T
    sep = np.linalg.norm(states[:, :2], axis=1)
    speed = np.linalg.norm(states[:, 2:], axis=1)

    v_arr = states[-1, 2:]
    v_final = v_arr + sol.dv2 if apply_terminal_burn else v_arr
    final_sep = float(sep[-1])
    final_speed = float(np.linalg.norm(v_final))
    pos_ok = final_sep < config.POSITION_TOL_KM
    vel_ok = final_speed < config.VELOCITY_TOL_KM_S
    status = STATUS_OK if (pos_ok and vel_ok) else (STATUS_POS_ONLY if pos_ok else STATUS_FAIL)

    err = states[-1] - cw_stm(n, tf) @ s0   # numerical vs analytic propagation
    return SimulationResult(
        scenario=sc, n_rad_s=n, period_s=orbital_period_s(sc.altitude_km),
        solution=sol, apply_terminal_burn=apply_terminal_burn,
        t_s=ode.t, states=states, separation_km=sep, rel_speed_km_s=speed,
        final_separation_km=final_sep,
        arrival_speed_km_s=float(np.linalg.norm(v_arr)),
        final_speed_km_s=final_speed, status=status,
        stm_vs_ode_pos_err_km=float(np.linalg.norm(err[:2])),
        stm_vs_ode_vel_err_km_s=float(np.linalg.norm(err[2:])), warnings=warns)
