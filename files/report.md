# Design and Simulation of a Relative-Motion Spacecraft Rendezvous in Low Earth Orbit

*All numbers in this report are produced by `scripts/build_report.py` from the project's own code
(`python scripts/build_report.py` regenerates this file). Statements about trends are asserted against
the data when the report is built.*

## 1. Introduction
Rendezvous means bringing a chaser spacecraft to the position **and** velocity of a target. Close to a
target on a circular orbit, relative motion is described well by the linear Clohessy–Wiltshire (CW,
also Hill) equations. This project builds a planar CW simulator, a closed-form two-burn rendezvous
solver, a parameter study, an animation and a Streamlit application.

## 2. Problem Statement
A chaser starts 1.0 km behind a target in a 400 km
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
"POSITION MATCH ONLY (relative velocity not nulled)" with a residual speed of 0.5594 m/s, so the software does not report a
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
* n = sqrt(mu / r^3): mean motion of the target [rad/s]. mu = 398600.4418 km^3/s^2.

Physically, 3 n^2 x is the tidal (gravity-gradient) term that pushes a chaser above/below the reference
orbit away from it, and the 2 n terms are Coriolis effects of the rotating frame. They couple radial and
along-track motion: a chaser given a radial or along-track velocity does not move in a straight line.
In code: `cw_dynamics(t, state, n)` in `src/cw_dynamics.py`.

For this project: r = R_E + h = 6378 + 400 = 6778 km,
**n = 1.131401e-03 rad/s**, orbital period T = 2 pi / n = **92.56 min**.

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
The coast arc is integrated numerically with `scipy.integrate.solve_ivp`, method **DOP853**
(8th-order Dormand–Prince), rtol = 1e-10, atol = 1e-12, 600
output samples by default. The analytic STM is used by the solver, and the numerical result is compared
with it as an independent check: for the default case the two final states differ by
8.51e-09 m and 1.11e-11 m/s.

## 11. Simulation Parameters
| Quantity | Value |
|---|---|
| Earth radius R_E | 6378 km |
| mu | 398600.4418 km^3/s^2 |
| Target altitude / radius | 400 km / 6778 km |
| Default initial chaser position (x0, y0) | (0, -1) km |
| Default transfer time | 30 min |
| Sweep transfer times | 10, 20, 30, 40, 60 min |
| Sweep separations (along-track, x0 = 0) | 0.5, 1, 2 km |
| Success tolerances (after burn 2) | position < 1 mm, speed < 1 mm/s |

Units: km, s, km/s internally; the UI and tables show m and m/s (factor 1000).

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
Default case (400 km, y0 = -1.0 km, 30 min):

| Quantity | Value |
|---|---|
| Burn-1 velocity (radial, along-track) | (-0.53458, 0.16481) m/s |
| Δv1 / Δv2 / total | 0.55941 / 0.55941 / 1.11882 m/s |
| Arrival speed before burn 2 | 0.55941 m/s |
| Final separation | 8.51e-09 m |
| Final relative velocity (after burn 2) | 1.11e-11 m/s |
| Status | RENDEZVOUS ACHIEVED |

Parameter sweep (all values computed by `scripts/parameter_sweep.py`; `x0 = 0`, chaser behind target):

| Initial separation | Transfer time | Δv1 [m/s] | Δv2 [m/s] | Δv total [m/s] | Final separation [m] | Final rel. velocity [m/s] | Status |
|---|--:|--:|--:|--:|--:|--:|---|
| 0.5 km | 10 min | 0.8787 | 0.8787 | 1.7573 | 7.60e-10 | 4.80e-13 | RENDEZVOUS ACHIEVED |
| 0.5 km | 20 min | 0.4510 | 0.4510 | 0.9020 | 2.71e-09 | 2.39e-12 | RENDEZVOUS ACHIEVED |
| 0.5 km | 30 min | 0.2797 | 0.2797 | 0.5594 | 4.52e-09 | 5.91e-12 | RENDEZVOUS ACHIEVED |
| 0.5 km | 40 min | 0.1825 | 0.1825 | 0.3649 | 5.43e-09 | 1.02e-11 | RENDEZVOUS ACHIEVED |
| 0.5 km | 60 min | 0.0825 | 0.0825 | 0.1650 | 7.00e-09 | 1.32e-11 | RENDEZVOUS ACHIEVED |
| 1.0 km | 10 min | 1.7573 | 1.7573 | 3.5147 | 1.36e-09 | 8.63e-13 | RENDEZVOUS ACHIEVED |
| 1.0 km | 20 min | 0.9020 | 0.9020 | 1.8041 | 5.06e-09 | 4.48e-12 | RENDEZVOUS ACHIEVED |
| 1.0 km | 30 min | 0.5594 | 0.5594 | 1.1188 | 8.51e-09 | 1.11e-11 | RENDEZVOUS ACHIEVED |
| 1.0 km | 40 min | 0.3649 | 0.3649 | 0.7298 | 1.02e-08 | 1.93e-11 | RENDEZVOUS ACHIEVED |
| 1.0 km | 60 min | 0.1650 | 0.1650 | 0.3300 | 1.33e-08 | 2.52e-11 | RENDEZVOUS ACHIEVED |
| 2.0 km | 10 min | 3.5147 | 3.5147 | 7.0294 | 2.56e-09 | 1.62e-12 | RENDEZVOUS ACHIEVED |
| 2.0 km | 20 min | 1.8041 | 1.8041 | 3.6082 | 9.75e-09 | 8.64e-12 | RENDEZVOUS ACHIEVED |
| 2.0 km | 30 min | 1.1188 | 1.1188 | 2.2376 | 1.64e-08 | 2.15e-11 | RENDEZVOUS ACHIEVED |
| 2.0 km | 40 min | 0.7298 | 0.7298 | 1.4596 | 1.97e-08 | 3.72e-11 | RENDEZVOUS ACHIEVED |
| 2.0 km | 60 min | 0.3300 | 0.3300 | 0.6600 | 2.59e-08 | 4.88e-11 | RENDEZVOUS ACHIEVED |

All 15 cases meet the tolerance. Final separations of order 1e-9 to 1e-8 m reflect integrator
error, not physical accuracy (see Section 21).

![trajectory](figures/trajectory.png)

## 14. Relative Trajectory Analysis
In the default case the chaser first moves toward Earth (burn-1 radial velocity is negative), reaching a
minimum radial offset of -0.2638 km (264 m) below the target's orbit, then returns to x = 0
as it reaches y = 0. A lower orbit has a shorter period, so the chaser gains angular position on the
target and closes the along-track gap. This is the standard CW "catch-up by dipping" behaviour.

## 15. Separation Analysis
Separation d(t) decreases strictly at every computed sample from 1000 m to
8.5e-09 m.

![separation](figures/separation.png)

## 16. Relative Velocity Analysis
Speed is not monotonic. It jumps to 0.5594 m/s at burn 1, peaks at
0.7616 m/s at t = 15.0 min, and arrives at
0.5594 m/s. Burn 2 then removes it, leaving 1.1e-11 m/s. In the
plot, the vertical segments are the impulsive burns.

![relative velocity](figures/relative_velocity.png)

## 17. Transfer-Time vs Δv Analysis
![dv vs time](figures/dv_vs_time.png)

* **Shorter transfers need larger Δv.** For 1 km separation, total Δv is 3.515 m/s at 10 min and
  0.330 m/s at 60 min (a factor 10.7). This holds for every separation in the sweep.
* **Δv scales exactly linearly with separation** for x0 = 0 (checked for each transfer time): the CW
  model is linear and the start point is on the along-track axis.
* **Δv1 = Δv2** in all computed cases (x0 = 0).
* **The decrease does not continue indefinitely.** On a dense curve for 1 km separation (2–120 min), the
  sampled minimum is 0.1199 m/s at 91.8 min, essentially one orbital
  period (92.56 min). The analytic limit as t_f approaches T (x0 = 0) is 2 n |y0| / (6 pi) =
  0.1200 m/s, which agrees (unit-tested). Δv rises again beyond one period and diverges near the
  singularities of Section 12.
* **Sensitivity to burn errors increases with transfer time.** The worst-case position miss caused by a
  1 mm/s error in burn 1's velocity is sigma_max(Phi_rv) x 1 mm/s:

| Transfer time | Miss at t_f per 1 mm/s burn-1 error [m] |
|---|--:|
| 10 min | 0.69 |
| 20 min | 1.94 |
| 30 min | 4.18 |
| 40 min | 7.37 |
| 60 min | 14.17 |

  This is an idealised sensitivity (no navigation or control model). It grows by a factor
  20.5 from 10 to 60 min, which is the cost of the lower Δv.
* **Practical range.** The model identifies no single optimum; it shows that Δv falls with time while
  error sensitivity and (Section 21) neglected-perturbation exposure grow. Anything beyond the sweep
  range needs more than this model can say.

## 18. Selected Rendezvous Strategy
The demonstration default is **30 min** transfer from 1 km. This is a design choice for the demo, not a
computed optimum: it gives total Δv of 1.119 m/s (versus 3.51 m/s at 10 min and
0.33 m/s at 60 min) while keeping the coast to about a third of an orbit, where the neglected effects
are smaller and burn-error sensitivity is 4.2 m per mm/s (versus 14.2 at 60 min). A different weighting
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
Automated tests (`pytest`): **21 passed in 0.74s**.

* Mean motion and period match known values (n = 1.1314e-03 rad/s, T = 92.56 min at 400 km).
* `cw_dynamics` returns the correct derivative structure and equals A s.
* Closed-form STM equals the matrix exponential (atol 1e-9) and Phi(0) = I.
* Solver: final separation and speed after burn 2 are below tolerance; the position equation residual is below 1e-12.
* Negative control: burn 1 only gives "position match only".
* All 15 sweep cases are processed and succeed.
* Invalid inputs raise `ValueError` with descriptive messages.
* Independent check: the burn plan is applied to full planar two-body motion (point-mass gravity,
  `src/nonlinear_check.py`). The position miss at t_f ranges from 0.007 m to 7.71 m
  across the sweep (0.389 m for the default case). The arrival speeds agree with the CW values.

Nonlinear miss distances [m] (CW burn plan, true two-body propagation, burn 2 not included):

| Separation | Transfer time | Miss [m] |
|---|---|--:|
| 0.5 km | 10 min | 0.007 |
| 0.5 km | 20 min | 0.034 |
| 0.5 km | 30 min | 0.097 |
| 0.5 km | 40 min | 0.205 |
| 0.5 km | 60 min | 0.482 |
| 1.0 km | 10 min | 0.027 |
| 1.0 km | 20 min | 0.135 |
| 1.0 km | 30 min | 0.389 |
| 1.0 km | 40 min | 0.820 |
| 1.0 km | 60 min | 1.927 |
| 2.0 km | 10 min | 0.110 |
| 2.0 km | 20 min | 0.541 |
| 2.0 km | 30 min | 1.556 |
| 2.0 km | 40 min | 3.282 |
| 2.0 km | 60 min | 7.709 |

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
A closed-form two-burn CW solver gives a rendezvous (position and velocity) for all 15 swept cases,
confirmed by independent numerical integration. Total Δv for 1 km falls from 3.51 m/s (10 min) to
0.33 m/s (60 min), scales linearly with separation, and reaches a minimum of about
0.12 m/s near one orbital period before rising again, while sensitivity to burn
errors and the linearisation error both increase with transfer time. Choosing a transfer time therefore
requires criteria (propellant, time, safety) outside this model.

## 24. References
See `references.md`. Textbook references are listed without page numbers.
