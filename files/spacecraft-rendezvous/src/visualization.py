"""Plots and animation, all generated from SimulationResult / DataFrame data."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

from .config import M_PER_KM

C_TRAJ, C_TARGET, C_CHASER = "#1f4e9c", "#d62728", "#2ca02c"


def save_figure(fig, path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")


def plot_trajectory(res):
    """Hero plot: along-track (horizontal) vs radial (vertical), target at origin."""
    x, y = res.states[:, 0], res.states[:, 1]
    fig, ax = plt.subplots(figsize=(8, 4.6))
    ax.plot(y, x, color=C_TRAJ, lw=2, label="Chaser trajectory (coast arc)")
    k = len(y) // 2
    ax.annotate("", xy=(y[k + 5], x[k + 5]), xytext=(y[k], x[k]),
                arrowprops=dict(arrowstyle="->", color=C_TRAJ, lw=2))
    ax.plot(y[0], x[0], "o", color=C_CHASER, ms=10, label="Initial chaser position")
    ax.plot(0, 0, "*", color=C_TARGET, ms=18, label="Target / rendezvous point")
    ax.axhline(0, color="grey", lw=0.5), ax.axvline(0, color="grey", lw=0.5)
    ax.set_xlabel("Along-track y [km]  (direction of motion →)")
    ax.set_ylabel("Radial x [km]  (up = away from Earth)")
    ax.set_title(f"Relative trajectory in target frame, t_f = {res.scenario.transfer_time_s/60:.0f} min "
                 "(axes not equal scale)")
    ax.grid(alpha=0.3), ax.legend(loc="best")
    return fig


def plot_separation(res):
    t = res.t_s / 60
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(t, res.separation_km * M_PER_KM, color=C_TRAJ, lw=2)
    ax.plot(t[-1], res.final_separation_km * M_PER_KM, "o", color=C_TARGET,
            label=f"Final: {res.final_separation_km*M_PER_KM:.2e} m")
    ax.set_xlabel("Time [min]"), ax.set_ylabel("Separation d(t) [m]")
    ax.set_title("Separation vs time"), ax.grid(alpha=0.3), ax.legend()
    return fig


def plot_relative_velocity(res):
    """v_rel(t) including the impulsive jumps at t=0 (burn 1) and t=tf (burn 2)."""
    t = res.t_s / 60
    tt = np.concatenate([[0], t, [t[-1]]])
    vv = np.concatenate([[0], res.rel_speed_km_s, [res.final_speed_km_s]]) * M_PER_KM
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(tt, vv, color=C_TRAJ, lw=2)
    ax.plot(t[-1], res.arrival_speed_km_s * M_PER_KM, "s", color=C_CHASER,
            label=f"Arrival speed (before burn 2): {res.arrival_speed_km_s*M_PER_KM:.4f} m/s")
    ax.plot(t[-1], res.final_speed_km_s * M_PER_KM, "o", color=C_TARGET,
            label=f"Final (after burn 2): {res.final_speed_km_s*M_PER_KM:.1e} m/s")
    ax.set_xlabel("Time [min]"), ax.set_ylabel("Relative speed v_rel(t) [m/s]")
    ax.set_title("Relative speed vs time (vertical segments = impulsive burns)")
    ax.grid(alpha=0.3), ax.legend(fontsize=8)
    return fig


def plot_dv_vs_time(curves: dict, sweep_df, period_min: float):
    """Trade-off plot. curves: {separation_km: dv_curve DataFrame}."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6))
    for sep, df in curves.items():
        line, = a1.plot(df.transfer_time_min, df.dv_total_m_s, lw=2, label=f"{sep:g} km separation")
        sub = sweep_df[sweep_df.separation_km == sep]
        a1.plot(sub.transfer_time_min, sub.dv_total_m_s, "o", color=line.get_color())
    a1.axvline(period_min, color="k", ls=":", lw=1)
    a1.text(period_min, a1.get_ylim()[1] * 0.9, " 1 orbit", fontsize=8)
    a1.set_yscale("log"), a1.set_xlim(0, 120)
    a1.set_xlabel("Transfer time [min]"), a1.set_ylabel("Total Δv = Δv1 + Δv2 [m/s]")
    a1.set_title("Transfer time vs required Δv (dots = sweep cases)")
    a1.grid(alpha=0.3, which="both"), a1.legend()
    ref = sorted(curves)[len(curves) // 2]
    d = curves[ref]
    a2.plot(d.transfer_time_min, d.dv1_m_s, label="Δv1 (=Δv2 here)", lw=2)
    a2.plot(d.transfer_time_min, d.dv_total_m_s, label="Total", lw=2)
    a2.axvline(period_min, color="k", ls=":", lw=1)
    a2.set_yscale("log"), a2.set_xlim(0, 120)
    a2.set_xlabel("Transfer time [min]"), a2.set_ylabel("Δv [m/s]")
    a2.set_title(f"Burn breakdown, {ref:g} km separation"), a2.grid(alpha=0.3, which="both"), a2.legend()
    fig.tight_layout()
    return fig


def make_animation(res, path, fps: int = 20, n_frames: int = 90, hold: int = 15) -> None:
    """GIF: hold at start -> coast along the computed trajectory -> hold at rendezvous."""
    x, y, t = res.states[:, 0], res.states[:, 1], res.t_s
    idx = np.linspace(0, len(t) - 1, n_frames).astype(int)
    frames = ([("INITIAL SEPARATION", 0)] * hold + [("APPROACH", i) for i in idx]
              + [("RENDEZVOUS", len(t) - 1)] * hold)
    fig, ax = plt.subplots(figsize=(8, 4.2))
    pad_y, pad_x = 0.1 * max(np.ptp(y), 1e-3), 0.3 * max(np.ptp(x), 0.05 * np.ptp(y), 1e-3)
    ax.set_xlim(y.min() - pad_y, y.max() + pad_y), ax.set_ylim(x.min() - pad_x, x.max() + pad_x)
    ax.set_xlabel("Along-track y [km]"), ax.set_ylabel("Radial x [km]"), ax.grid(alpha=0.3)
    ax.plot(y, x, color="lightgrey", lw=1, zorder=1)
    trail, = ax.plot([], [], color=C_TRAJ, lw=2)
    ch, = ax.plot([], [], "^", color=C_CHASER, ms=13)
    ax.plot(0, 0, "*", color=C_TARGET, ms=18)
    ax.text(0, 0, "  TARGET", color=C_TARGET, va="bottom", weight="bold")
    lab = ax.text(0.01, 0.97, "", transform=ax.transAxes, va="top", family="monospace", fontsize=9)
    ctext = ax.text(y[0], x[0], " CHASER", color=C_CHASER, va="bottom", weight="bold")

    def update(f):
        phase, i = frames[f]
        trail.set_data(y[: i + 1], x[: i + 1])
        ch.set_data([y[i]], [x[i]])
        ctext.set_position((y[i], x[i]))
        if phase == "INITIAL SEPARATION":
            v = 0.0                                   # before burn 1
        elif phase == "RENDEZVOUS":
            v = res.final_speed_km_s                  # after burn 2
        else:
            v = res.rel_speed_km_s[i]
        lab.set_text(f"Phase: {phase}\nMission time: {t[i]/60:6.2f} min\n"
                     f"Distance: {res.separation_km[i]*1000:9.2f} m\nRel. speed: {v*1000:8.4f} m/s")
        return trail, ch, lab, ctext

    anim = FuncAnimation(fig, update, frames=len(frames), blit=False)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    anim.save(str(path), writer=PillowWriter(fps=fps))
    plt.close(fig)
