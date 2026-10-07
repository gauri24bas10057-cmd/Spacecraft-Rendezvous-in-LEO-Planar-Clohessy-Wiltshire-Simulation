"""Run the default scenario (400 km, y0 = -1 km, 30 min), print a summary, save figures."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")

from src import visualization as viz                      # noqa: E402
from src.config import M_PER_KM, Scenario                  # noqa: E402
from src.simulation import run_simulation                  # noqa: E402


def main() -> None:
    res = run_simulation(Scenario())
    s, sc = res.solution, res.scenario
    print(f"Target altitude        : {sc.altitude_km:.1f} km")
    print(f"Mean motion n          : {res.n_rad_s:.6e} rad/s  (period {res.period_s/60:.2f} min)")
    print(f"Initial chaser (x, y)  : ({sc.x0_km}, {sc.y0_km}) km, transfer {sc.transfer_time_s/60:.0f} min")
    print(f"Burn-1 velocity (x, y) : ({s.v0_plus[0]*M_PER_KM:.5f}, {s.v0_plus[1]*M_PER_KM:.5f}) m/s")
    print(f"dv1 / dv2 / total      : {s.dv1_mag*M_PER_KM:.5f} / {s.dv2_mag*M_PER_KM:.5f} / "
          f"{s.dv_total*M_PER_KM:.5f} m/s")
    print(f"Final separation       : {res.final_separation_km*M_PER_KM:.3e} m")
    print(f"Arrival speed (pre-burn 2): {res.arrival_speed_km_s*M_PER_KM:.5f} m/s")
    print(f"Final rel. velocity    : {res.final_speed_km_s*M_PER_KM:.3e} m/s")
    print(f"STM vs ODE error       : {res.stm_vs_ode_pos_err_km*M_PER_KM:.2e} m, "
          f"{res.stm_vs_ode_vel_err_km_s*M_PER_KM:.2e} m/s")
    print(f"Status                 : {res.status}")
    neg = run_simulation(sc, apply_terminal_burn=False)
    print(f"Negative control (burn 1 only): {neg.status}")

    out = ROOT / "results"
    viz.save_figure(viz.plot_trajectory(res), out / "trajectory.png")
    viz.save_figure(viz.plot_separation(res), out / "separation.png")
    viz.save_figure(viz.plot_relative_velocity(res), out / "relative_velocity.png")
    viz.make_animation(res, out / "rendezvous.gif")
    print(f"Saved figures and animation to {out}")


if __name__ == "__main__":
    main()
