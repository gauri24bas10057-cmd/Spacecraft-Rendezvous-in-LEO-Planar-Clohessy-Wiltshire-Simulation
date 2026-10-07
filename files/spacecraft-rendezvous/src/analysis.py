"""Parameter sweeps, trade-off curves and sensitivity metrics."""
import numpy as np
import pandas as pd

from . import config
from .config import M_PER_KM, SECONDS_PER_MINUTE, Scenario
from .cw_dynamics import cw_stm
from .nonlinear_check import nonlinear_miss
from .orbital import mean_motion_rad_s, orbital_period_s
from .simulation import run_simulation


def parameter_sweep(separations_km=config.SWEEP_INITIAL_SEPARATIONS_KM,
                    times_min=config.SWEEP_TRANSFER_TIMES_MIN,
                    altitude_km: float = config.DEFAULT_ALTITUDE_KM,
                    x0_km: float = 0.0) -> pd.DataFrame:
    """One row per (initial along-track separation, transfer time) case.

    Chaser starts at (x0, -sep) behind the target. Velocities/distances are
    reported in m/s and m for readability (internal units remain km, km/s).
    """
    n = mean_motion_rad_s(altitude_km)
    rows = []
    for sep in separations_km:
        for tmin in times_min:
            sc = Scenario(altitude_km, x0_km, -sep, tmin * SECONDS_PER_MINUTE)
            res = run_simulation(sc)
            phi_rv = cw_stm(n, sc.transfer_time_s)[:2, 2:]
            miss_nl, vrel_nl = nonlinear_miss(sc)
            s = res.solution
            rows.append({
                "separation_km": sep, "transfer_time_min": tmin,
                "dv1_m_s": s.dv1_mag * M_PER_KM, "dv2_m_s": s.dv2_mag * M_PER_KM,
                "dv_total_m_s": s.dv_total * M_PER_KM,
                "v0_plus_radial_m_s": s.v0_plus[0] * M_PER_KM,
                "v0_plus_alongtrack_m_s": s.v0_plus[1] * M_PER_KM,
                "final_separation_m": res.final_separation_km * M_PER_KM,
                "arrival_speed_m_s": res.arrival_speed_km_s * M_PER_KM,
                "final_rel_speed_m_s": res.final_speed_km_s * M_PER_KM,
                "status": res.status,
                # sigma_max(Phi_rv) [s] x 1 mm/s = worst-case miss [m] per 1 mm/s burn-1 error
                "miss_per_1mm_s_burn_error_m": np.linalg.svd(phi_rv, compute_uv=False)[0] * 1e-3,
                "nonlinear_miss_m": miss_nl * M_PER_KM,
                "nonlinear_arrival_speed_m_s": vrel_nl * M_PER_KM,
            })
    return pd.DataFrame(rows)


def dv_curve(separation_km: float, altitude_km: float = config.DEFAULT_ALTITUDE_KM,
             t_min_range=(2.0, 120.0), num: int = 400) -> pd.DataFrame:
    """Dense transfer-time vs delta-v curve (closed form; NaN at singular times)."""
    from .rendezvous import solve_rendezvous
    n = mean_motion_rad_s(altitude_km)
    rows = []
    for tmin in np.linspace(*t_min_range, num):
        try:
            s = solve_rendezvous([0.0, -separation_km], [0, 0], n, tmin * SECONDS_PER_MINUTE)
            rows.append((tmin, s.dv1_mag * M_PER_KM, s.dv2_mag * M_PER_KM, s.dv_total * M_PER_KM))
        except ValueError:
            rows.append((tmin, np.nan, np.nan, np.nan))
    return pd.DataFrame(rows, columns=["transfer_time_min", "dv1_m_s", "dv2_m_s", "dv_total_m_s"])


def to_markdown_table(df: pd.DataFrame) -> str:
    """Compact results table for README/report."""
    lines = ["| Initial separation | Transfer time | Δv1 [m/s] | Δv2 [m/s] | Δv total [m/s] "
             "| Final separation [m] | Final rel. velocity [m/s] | Status |",
             "|---|--:|--:|--:|--:|--:|--:|---|"]
    for r in df.itertuples():
        lines.append(f"| {r.separation_km:.1f} km | {r.transfer_time_min:.0f} min | "
                     f"{r.dv1_m_s:.4f} | {r.dv2_m_s:.4f} | {r.dv_total_m_s:.4f} | "
                     f"{r.final_separation_m:.2e} | {r.final_rel_speed_m_s:.2e} | {r.status} |")
    return "\n".join(lines)
