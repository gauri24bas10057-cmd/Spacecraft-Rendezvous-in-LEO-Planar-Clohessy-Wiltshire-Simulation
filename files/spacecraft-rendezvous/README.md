# Spacecraft Rendezvous in LEO — Planar Clohessy–Wiltshire Simulation

## Overview
A modular Python project that simulates a chaser spacecraft rendezvousing with a target on a circular
400 km orbit, using the planar Clohessy–Wiltshire (CW/Hill) equations. It includes a closed-form
two-burn solver, a parameter sweep, plots, an animation, a Streamlit app, tests and a technical report.
It is an educational / preliminary-analysis model, not a mission-planning tool.

## Problem statement
A chaser starts 1 km behind a target. How do initial relative position and transfer time influence the
velocity correction required to rendezvous (match position **and** velocity)?

## Objectives
- Model planar relative motion with the CW equations.
- Solve for the burns that give zero relative position and velocity at a chosen time.
- Study transfer time vs Δv for several initial separations.
- Document assumptions, validation and limitations honestly.

## Features
- Closed-form CW state-transition matrix and two-impulse solver (no black-box optimiser)
- Numerical propagation with `scipy.integrate.solve_ivp` (DOP853), cross-checked against the analytic STM
- Rendezvous judged on position **and** velocity; negative control shows position-only is rejected
- Parameter sweep (5 transfer times x 3 separations), dense Δv-vs-time curve, burn-error sensitivity
- Independent check against planar two-body motion (linearisation error)
- Animated GIF and a Streamlit simulator; 21 automated tests; report regenerated from live data

## Physics model
Frame: origin on target, x radial (outward), y along-track. State `[x, y, x_dot, y_dot]`, units km and km/s.

### CW equations
```
x_ddot =  3 n^2 x + 2 n y_dot
y_ddot = -2 n x_dot          n = sqrt(mu / r^3),  r = R_E + h
```
With Phi(t_f) partitioned into 2x2 blocks, burn 1 gives `v0+ = -inv(Phi_rv) Phi_rr r0` (arrive at the target),
and burn 2 cancels the arrival velocity. Δv is the sum of idealised impulsive burns, **not** a propulsion budget.
Full derivation: [`report/report.md`](report/report.md).

## System architecture
`config` → `orbital` → `cw_dynamics` → `rendezvous` → `simulation` → `analysis` → `visualization`;
`scripts/` and `app/` only call the top layers. `nonlinear_check.py` (an addition to the original plan)
propagates the burn plan under full two-body gravity to quantify the linearisation error.

## Project structure
```
spacecraft-rendezvous/
├── README.md  LICENSE  .gitignore  pyproject.toml  requirements.txt  requirements-dev.txt
├── src/       config  orbital  cw_dynamics  rendezvous  simulation  analysis  visualization  nonlinear_check
├── app/main.py                       Streamlit app
├── scripts/   run_simulation.py  parameter_sweep.py  build_report.py
├── tests/test_rendezvous.py
├── results/   figures, rendezvous.gif, sweep_results.csv, sweep_table.md, dv_curve_1km.csv
├── notebooks/exploration.ipynb
└── report/    report.md  references.md  figures/
```

## Installation
```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt    # or requirements-dev.txt to also get pytest
```

## Dependencies
numpy, scipy, matplotlib, pandas, streamlit (+ pytest for tests; Pillow comes with matplotlib).

## How to run
```bash
python scripts/run_simulation.py     # default case: summary + trajectory/separation/velocity plots + GIF
python scripts/parameter_sweep.py    # sweep: CSV, table, dv_vs_time.png
python scripts/build_report.py       # regenerate report/report.md from live results
pytest                               # run tests (needs requirements-dev.txt)
```

## How to launch Streamlit
```bash
streamlit run app/main.py
```

## Example simulation
Default case: 400 km, chaser at (x, y) = (0, −1) km, 30 min transfer.

| Quantity | Value |
|---|---|
| Mean motion n | 1.131401e-03 rad/s (period 92.56 min) |
| Burn-1 velocity (radial, along-track) | (-0.53458, 0.16481) m/s |
| Δv1 / Δv2 / total | 0.55941 / 0.55941 / 1.11882 m/s |
| Final separation | 8.51e-09 m |
| Final relative velocity | 1.11e-11 m/s |
| Status | RENDEZVOUS ACHIEVED |

## Results
Produced by `scripts/parameter_sweep.py` (x0 = 0, chaser behind target):

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

Key observations (details and checks in the report): total Δv falls with transfer time over the swept
range, is exactly linear in separation, has a minimum near one orbital period, and the sensitivity to
burn errors and the linearisation error both grow with transfer time. Final separations of ~1e-9 m are
numerical residuals of an idealised model.

## Visualizations
| | |
|---|---|
| ![traj](results/trajectory.png) | ![sep](results/separation.png) |
| ![vel](results/relative_velocity.png) | ![dv](results/dv_vs_time.png) |

Animation: `results/rendezvous.gif` (INITIAL SEPARATION → APPROACH → RENDEZVOUS, built from the simulated states).

## Validation
21 tests cover mean motion, CW derivative structure, STM vs matrix exponential, solver tolerances and
residuals, the position-only negative control, an analytic one-orbit Δv limit, linear Δv scaling, the
full sweep, a two-body nonlinear cross-check, and invalid-input handling.

## Assumptions
Spherical Earth; circular target orbit; small separations (linearised model); no drag; no J2; negligible
masses; impulsive, instantaneous burns; planar motion; chaser starts co-orbital (zero relative velocity).

## Limitations
No nonlinear dynamics (beyond the check), J2, drag, solar pressure, attitude, RCS dynamics, finite burns,
thrust errors, navigation errors, sensor noise, actuator limits, docking contact or structural dynamics.
Transfer times at singularities of the CW map (n·t_f = 2πk, and tan(n·t_f/2) = 3n·t_f/4) are rejected.
Not for flight certification or real mission planning.

## Future improvements
3D CW motion, eccentric orbits, J2/drag, finite-thrust burns, optimal multi-impulse transfers,
navigation noise and Monte-Carlo analysis, keep-out zones and V-bar/R-bar corridors.

## References
See [`report/references.md`](report/references.md) (Clohessy & Wiltshire 1960; Hill 1878; Vallado; Curtis; Wie; Fehse; Alfriend et al.).
Verify edition and page details before academic submission.
