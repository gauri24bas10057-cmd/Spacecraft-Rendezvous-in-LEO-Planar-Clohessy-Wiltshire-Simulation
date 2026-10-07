"""Generate report/report.md from live simulation output (no hand-typed numbers).

Qualitative statements in the report are asserted against the data, so the
build fails if a claim stops being true.   Run: python scripts/build_report.py
"""
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np                                            # noqa: E402

from src import config                                        # noqa: E402
from src.analysis import dv_curve, parameter_sweep, to_markdown_table  # noqa: E402
from src.config import M_PER_KM as K, Scenario                # noqa: E402
from src.orbital import orbit_radius_km, orbital_period_s     # noqa: E402
from src.simulation import STATUS_OK, run_simulation          # noqa: E402


def main() -> None:
    df = parameter_sweep()
    res = run_simulation(Scenario())
    neg = run_simulation(Scenario(), apply_terminal_burn=False)
    s, n, T = res.solution, res.n_rad_s, res.period_s
    one = df[df.separation_km == 1.0].sort_values("transfer_time_min")
    curve = dv_curve(1.0, num=1200).dropna()
    imin = curve.dv_total_m_s.idxmin()
    analytic_T = 2 * n * 1.0 / (6 * math.pi) * K            # m/s, one-period limit (x0=0, 1 km)

    # ---- assert every qualitative claim used in the prose ----
    assert (df.status == STATUS_OK).all()
    for sep, g in df.groupby("separation_km"):
        g = g.sort_values("transfer_time_min")
        assert g.dv_total_m_s.is_monotonic_decreasing                  # dv falls with time (<=60 min)
        assert g.miss_per_1mm_s_burn_error_m.is_monotonic_increasing   # sensitivity rises
        assert g.nonlinear_miss_m.is_monotonic_increasing
    for t, g in df.groupby("transfer_time_min"):
        g = g.sort_values("separation_km")
        ratio = g.dv_total_m_s.values / g.separation_km.values
        assert np.allclose(ratio, ratio[0], rtol=1e-9)                 # dv linear in separation
        assert g.nonlinear_miss_m.is_monotonic_increasing
    assert np.allclose(df.dv1_m_s, df.dv2_m_s, rtol=1e-6)              # dv1 == dv2 for x0 = 0
    assert abs(T / 60 - curve.transfer_time_min[imin]) < 2.0           # minimum near one period
    assert abs(curve.dv_total_m_s[imin] - analytic_T) / analytic_T < 0.01
    assert (np.diff(res.separation_km) < 0).all()                      # d(t) strictly decreasing
    xmin, vmax_i = res.states[:, 0].min(), int(np.argmax(res.rel_speed_km_s))
    assert xmin < 0 and neg.final_speed_km_s > 1e-4

    pytest_out = subprocess.run([sys.executable, "-m", "pytest", "-q"], cwd=ROOT,
                                capture_output=True, text=True).stdout.strip().splitlines()[-1]
    fig = ROOT / "report" / "figures"
    for f in ("trajectory.png", "separation.png", "relative_velocity.png", "dv_vs_time.png"):
        shutil.copy(ROOT / "results" / f, fig / f)

    d10, d60 = one.dv_total_m_s.iloc[0], one.dv_total_m_s.iloc[-1]
    sens10, sens60 = one.miss_per_1mm_s_burn_error_m.iloc[0], one.miss_per_1mm_s_burn_error_m.iloc[-1]
    nl = df[(df.separation_km == 1.0) & (df.transfer_time_min == 30)].iloc[0]
    nl_min, nl_max = df.nonlinear_miss_m.min(), df.nonlinear_miss_m.max()
    sw = to_markdown_table(df)
    sens = "\n".join(
        f"| {r.transfer_time_min:.0f} min | {r.miss_per_1mm_s_burn_error_m:.2f} |" for r in one.itertuples())
    nlt = "\n".join(
        f"| {r.separation_km:.1f} km | {r.transfer_time_min:.0f} min | {r.nonlinear_miss_m:.3f} |"
        for r in df.itertuples())

    text = f"""# Design and Simulation of a Relative-Motion Spacecraft Rendezvous in Low Earth Orbit

*All numbers in this report are produced by `scripts/build_report.py` from the project's own code
(`python scripts/build_report.py` regenerates this file). Statements about trends are asserted against
the data when the report is built.*

## 1. Introduction
Rendezvous means bringing a chaser spacecraft to the position **and** velocity of a target. Close to a
target on a circular orbit, relative motion is described well by the linear Clohessy–Wiltshire (CW,
also Hill) equations. This project builds a planar CW simulator, a closed-form two-burn rendezvous
solver, a parameter study, an animation and a Streamlit application.

## 2. Problem Statement
A chaser starts {abs(res.scenario.y0_km):.1f} km behind a target in a {res.scenario.altitude_km:.0f} km
circular orbit. How do the initial relative position and the transfer time influence the velocity
correction required to rendezvous?

## 3. Mission Objective
Find, for chosen transfer times t_f, the burns that give x, y, x_dot, y_dot = 0 at t_f; quantify the
required Δv; and analyse the transfer-time / Δv trade-off from the computed results.

## 4. Rendezvous Concept
A single impulse at t = 0 has two free parameters (the two velocity components) but the goal has four
conditions (two position, two velocity). One impulse can therefore only fix the arrival *position*.
The plan used here is:

1. **Burn 1 (t = 0):** set the relative velocity so the free coast arrives at the target position at t_f.
2. **Coast** for t_f under CW dynamics.
3. **Burn 2 (t = t_f):** cancel the arrival relative velocity.

Rendezvous status is judged after burn 2. As a negative control, applying burn 1 only gives status
"{neg.status}" with a residual speed of {neg.final_speed_km_s*K:.4f} m/s, so the software does not report a
rendezvous from position agreement alone.

## 5. Reference Frame
Origin on the target. **x** is radial (away from Earth), **y** is along-track (direction of motion);
the frame rotates at the target's mean motion n. Negative y means the chaser is behind the target.
Plots show y horizontally and x vertically (the usual V-bar presentation). Axes are not equal scale in the
trajectory plot, because radial excursions are much smaller than the along-track distance.

## 6. Assumptions
1. Earth is a spherical central body (point-mass gravity).
2. The target is on a circular orbit.
3. The chaser stays close enough to the target for linearisation to hold.
4. Atmospheric drag is neglected.
5. J2 and other perturbations are neglected.
6. Spacecraft masses are negligible compared to Earth.
7. Manoeuvres are impulsive.
8. Motion is confined to the orbital plane.
9. Thrust is instantaneous (no burn duration, no finite-thrust losses).
10. Before burn 1 the chaser is co-orbital: zero relative velocity.

## 7. Theoretical Background
The CW equations were introduced for terminal rendezvous guidance (Clohessy and Wiltshire, 1960) and
follow from Hill's lunar-theory equations (Hill, 1878) by linearising two-body motion about a circular
reference orbit. Derivations appear in standard texts (see References).

## 8. Clohessy–Wiltshire Equations
In-plane form:

    x_ddot - 2 n y_dot - 3 n^2 x = 0
    y_ddot + 2 n x_dot          = 0

* x, y: radial and along-track relative position [km]; x_dot, y_dot: relative velocity [km/s].
* n = sqrt(mu / r^3): mean motion of the target [rad/s]. mu = {config.EARTH_MU_KM3_S2} km^3/s^2.

Physically, 3 n^2 x is the tidal (gravity-gradient) term that pushes a chaser above/below the reference
orbit away from it, and the 2 n terms are Coriolis effects of the rotating frame. They couple radial and
along-track motion: a chaser given a radial or along-track velocity does not move in a straight line.
In code: `cw_dynamics(t, state, n)` in `src/cw_dynamics.py`.

For this project: r = R_E + h = {config.EARTH_RADIUS_KM:.0f} + {res.scenario.altitude_km:.0f} = {orbit_radius_km(400):.0f} km,
**n = {n:.6e} rad/s**, orbital period T = 2 pi / n = **{T/60:.2f} min**.

## 9. Mathematical Formulation
**State-space form.** With s = [x, y, x_dot, y_dot], the equations become s_dot = A s:

        | 0     0    1    0   |
    A = | 0     0    0    1   |        x_ddot =  3 n^2 x + 2 n y_dot
        | 3n^2  0    0    2n  |        y_ddot = -2 n x_dot
        | 0     0   -2n   0   |

(`cw_matrix_A`). Because A is constant, the exact solution is s(t) = Phi(t) s(0) with Phi = exp(A t).
With c = cos(nt), s = sin(nt):

    x(t)     = (4-3c) x0 + (s/n) xd0 + (2/n)(1-c) yd0
    y(t)     = 6(s-nt) x0 + y0 - (2/n)(1-c) xd0 + ((4s-3nt)/n) yd0
    x_dot(t) = 3 n s x0 + c xd0 + 2 s yd0
    y_dot(t) = 6 n (c-1) x0 - 2 s xd0 + (4c-3) yd0

(`cw_stm`; verified against the matrix exponential in the tests.) Relative separation and speed:

    d(t) = sqrt(x^2 + y^2)          v_rel(t) = sqrt(x_dot^2 + y_dot^2)

## 10. Numerical Method
The coast arc is integrated numerically with `scipy.integrate.solve_ivp`, method **{config.INTEGRATOR_METHOD}**
(8th-order Dormand–Prince), rtol = {config.INTEGRATOR_RTOL:g}, atol = {config.INTEGRATOR_ATOL:g}, {config.DEFAULT_N_STEPS}
output samples by default. The analytic STM is used by the solver, and the numerical result is compared
with it as an independent check: for the default case the two final states differ by
{res.stm_vs_ode_pos_err_km*K:.2e} m and {res.stm_vs_ode_vel_err_km_s*K:.2e} m/s.

## 11. Simulation Parameters
| Quantity | Value |
|---|---|
| Earth radius R_E | {config.EARTH_RADIUS_KM:.0f} km |
| mu | {config.EARTH_MU_KM3_S2} km^3/s^2 |
| Target altitude / radius | 400 km / {orbit_radius_km(400):.0f} km |
| Default initial chaser position (x0, y0) | (0, -1) km |
| Default transfer time | 30 min |
| Sweep transfer times | {", ".join(f"{t:g}" for t in config.SWEEP_TRANSFER_TIMES_MIN)} min |
| Sweep separations (along-track, x0 = 0) | {", ".join(f"{t:g}" for t in config.SWEEP_INITIAL_SEPARATIONS_KM)} km |
| Success tolerances (after burn 2) | position < {config.POSITION_TOL_KM*K*1000:.0f} mm, speed < {config.VELOCITY_TOL_KM_S*K*1000:.0f} mm/s |

Units: km, s, km/s internally; the UI and tables show m and m/s (factor {K:.0f}).

## 12. Rendezvous Algorithm
Partition Phi(t_f) into 2x2 blocks [[Phi_rr, Phi_rv],[Phi_vr, Phi_vv]] (position/velocity rows and
columns). Requiring the final position to be zero:

    0 = Phi_rr r0 + Phi_rv v0+      =>      v0+ = - Phi_rv^-1 Phi_rr r0

* r0 = (x0, y0): initial position; v0+: relative velocity just after burn 1 (the quantity solved for).
* This is a linear 2x2 solve: no iteration, no tolerance, no black-box optimiser.
* Arrival velocity: v_f- = Phi_vr r0 + Phi_vv v0+. Burn 2 is dv2 = - v_f-.
* Burn 1 is dv1 = v0+ - v0-, with v0- = 0 (co-orbital start).
* Magnitudes: |dv1|, |dv2|, and total Δv = |dv1| + |dv2|, the sum of idealised impulsive burns.
  It is **not** a complete propulsion budget (no gravity/steering losses, margins, attitude control, etc.).

The solve is singular when det(Phi_rv) = 0. Algebraically det(Phi_rv) n^2 = 4 sin(x)(4 sin x - 3 x cos x)
with x = n t_f / 2, which vanishes at n t_f = 2 pi k and where tan(x) = 3x/4 (first such root near
n t_f = 8.9, about 131 min at 400 km). The code rejects transfer times within a small tolerance of
these with a clear error. All sweep times (<= 60 min) are far from them.

Convergence check: the solver output is verified by independent numerical integration (Section 10)
and by plugging back into Phi (unit tests, residual below 1e-12 km).

## 13. Simulation Results
Default case ({res.scenario.altitude_km:.0f} km, y0 = {res.scenario.y0_km} km, 30 min):

| Quantity | Value |
|---|---|
| Burn-1 velocity (radial, along-track) | ({s.v0_plus[0]*K:.5f}, {s.v0_plus[1]*K:.5f}) m/s |
| Δv1 / Δv2 / total | {s.dv1_mag*K:.5f} / {s.dv2_mag*K:.5f} / {s.dv_total*K:.5f} m/s |
| Arrival speed before burn 2 | {res.arrival_speed_km_s*K:.5f} m/s |
| Final separation | {res.final_separation_km*K:.2e} m |
| Final relative velocity (after burn 2) | {res.final_speed_km_s*K:.2e} m/s |
| Status | {res.status} |

Parameter sweep (all values computed by `scripts/parameter_sweep.py`; `x0 = 0`, chaser behind target):

{sw}

All {len(df)} cases meet the tolerance. Final separations of order 1e-9 to 1e-8 m reflect integrator
error, not physical accuracy (see Section 21).

![trajectory](figures/trajectory.png)

## 14. Relative Trajectory Analysis
In the default case the chaser first moves toward Earth (burn-1 radial velocity is negative), reaching a
minimum radial offset of {xmin:.4f} km ({abs(xmin)*K:.0f} m) below the target's orbit, then returns to x = 0
as it reaches y = 0. A lower orbit has a shorter period, so the chaser gains angular position on the
target and closes the along-track gap. This is the standard CW "catch-up by dipping" behaviour.

## 15. Separation Analysis
Separation d(t) decreases strictly at every computed sample from {res.separation_km[0]*K:.0f} m to
{res.final_separation_km*K:.1e} m.

![separation](figures/separation.png)

## 16. Relative Velocity Analysis
Speed is not monotonic. It jumps to {res.rel_speed_km_s[0]*K:.4f} m/s at burn 1, peaks at
{res.rel_speed_km_s[vmax_i]*K:.4f} m/s at t = {res.t_s[vmax_i]/60:.1f} min, and arrives at
{res.arrival_speed_km_s*K:.4f} m/s. Burn 2 then removes it, leaving {res.final_speed_km_s*K:.1e} m/s. In the
plot, the vertical segments are the impulsive burns.

![relative velocity](figures/relative_velocity.png)

## 17. Transfer-Time vs Δv Analysis
![dv vs time](figures/dv_vs_time.png)

* **Shorter transfers need larger Δv.** For 1 km separation, total Δv is {d10:.3f} m/s at 10 min and
  {d60:.3f} m/s at 60 min (a factor {d10/d60:.1f}). This holds for every separation in the sweep.
* **Δv scales exactly linearly with separation** for x0 = 0 (checked for each transfer time): the CW
  model is linear and the start point is on the along-track axis.
* **Δv1 = Δv2** in all computed cases (x0 = 0).
* **The decrease does not continue indefinitely.** On a dense curve for 1 km separation (2–120 min), the
  sampled minimum is {curve.dv_total_m_s[imin]:.4f} m/s at {curve.transfer_time_min[imin]:.1f} min, essentially one orbital
  period ({T/60:.2f} min). The analytic limit as t_f approaches T (x0 = 0) is 2 n |y0| / (6 pi) =
  {analytic_T:.4f} m/s, which agrees (unit-tested). Δv rises again beyond one period and diverges near the
  singularities of Section 12.
* **Sensitivity to burn errors increases with transfer time.** The worst-case position miss caused by a
  1 mm/s error in burn 1's velocity is sigma_max(Phi_rv) x 1 mm/s:

| Transfer time | Miss at t_f per 1 mm/s burn-1 error [m] |
|---|--:|
{sens}

  This is an idealised sensitivity (no navigation or control model). It grows by a factor
  {sens60/sens10:.1f} from 10 to 60 min, which is the cost of the lower Δv.
* **Practical range.** The model identifies no single optimum; it shows that Δv falls with time while
  error sensitivity and (Section 21) neglected-perturbation exposure grow. Anything beyond the sweep
  range needs more than this model can say.

## 18. Selected Rendezvous Strategy
The demonstration default is **30 min** transfer from 1 km. This is a design choice for the demo, not a
computed optimum: it gives total Δv of {one[one.transfer_time_min==30].dv_total_m_s.iloc[0]:.3f} m/s (versus {d10:.2f} m/s at 10 min and
{d60:.2f} m/s at 60 min) while keeping the coast to about a third of an orbit, where the neglected effects
are smaller and burn-error sensitivity is {one[one.transfer_time_min==30].miss_per_1mm_s_burn_error_m.iloc[0]:.1f} m per mm/s (versus {sens60:.1f} at 60 min). A different weighting
of propellant, timeline and safety would select a different time.

## 19. Interactive Demonstration
`streamlit run app/main.py` offers inputs for altitude, x0, y0, transfer time and time resolution,
shows n, burn-1 velocity, Δv1/Δv2/total, final separation, final relative velocity and status, and plots
trajectory, separation and speed. Optional animation (`results/rendezvous.gif` is the default case) goes
through INITIAL SEPARATION, APPROACH and RENDEZVOUS phases and displays mission time, distance,
relative speed and phase from the simulated states. A checkbox disables burn 2 to show a
position-only match. Invalid inputs (negative altitude, zero transfer time, singular times) show error
messages instead of results.

## 20. Validation
Automated tests (`pytest`): **{pytest_out}**.

* Mean motion and period match known values (n = {n:.4e} rad/s, T = {T/60:.2f} min at 400 km).
* `cw_dynamics` returns the correct derivative structure and equals A s.
* Closed-form STM equals the matrix exponential (atol 1e-9) and Phi(0) = I.
* Solver: final separation and speed after burn 2 are below tolerance; the position equation residual is below 1e-12.
* Negative control: burn 1 only gives "position match only".
* All {len(df)} sweep cases are processed and succeed.
* Invalid inputs raise `ValueError` with descriptive messages.
* Independent check: the burn plan is applied to full planar two-body motion (point-mass gravity,
  `src/nonlinear_check.py`). The position miss at t_f ranges from {nl_min:.3f} m to {nl_max:.2f} m
  across the sweep ({nl.nonlinear_miss_m:.3f} m for the default case). The arrival speeds agree with the CW values.

Nonlinear miss distances [m] (CW burn plan, true two-body propagation, burn 2 not included):

| Separation | Transfer time | Miss [m] |
|---|---|--:|
{nlt}

The miss grows with transfer time and with separation, as expected for a linearisation error. It is
the quantified size of the CW approximation for these cases, not an error in the CW algebra.

## 21. Limitations
The model does not capture nonlinear orbital dynamics (beyond the check above), J2, atmospheric drag,
solar radiation pressure, attitude, reaction-control dynamics, finite-duration burns, thrust-direction
errors, navigation errors, sensor noise, actuator limits, docking/contact dynamics or structural
dynamics. The reported final separations (1e-9 to 1e-8 m) are numerical residuals of an idealised model; real
performance would be dominated by the omitted effects (for example, the nonlinear check above already
shows metre-level differences). The model is also strictly planar: out-of-plane motion is not represented. The project is for preliminary relative-motion analysis and education,
not flight certification or mission planning.

## 22. Future Improvements
Out-of-plane (full 3D CW) motion; eccentric reference orbits (Tschauner–Hempel); J2 and drag
perturbations; finite-burn and thrust-limited manoeuvres; multi-impulse or optimal transfers; navigation
noise with Monte-Carlo dispersion; safety constraints such as keep-out zones; approach corridors
(V-bar/R-bar).

## 23. Conclusion
A closed-form two-burn CW solver gives a rendezvous (position and velocity) for all {len(df)} swept cases,
confirmed by independent numerical integration. Total Δv for 1 km falls from {d10:.2f} m/s (10 min) to
{d60:.2f} m/s (60 min), scales linearly with separation, and reaches a minimum of about
{curve.dv_total_m_s[imin]:.2f} m/s near one orbital period before rising again, while sensitivity to burn
errors and the linearisation error both increase with transfer time. Choosing a transfer time therefore
requires criteria (propellant, time, safety) outside this model.

## 24. References
See `references.md`. Textbook references are listed without page numbers.
"""
    (ROOT / "report" / "report.md").write_text(text)
    print("Wrote report/report.md;", pytest_out)


if __name__ == "__main__":
    main()
