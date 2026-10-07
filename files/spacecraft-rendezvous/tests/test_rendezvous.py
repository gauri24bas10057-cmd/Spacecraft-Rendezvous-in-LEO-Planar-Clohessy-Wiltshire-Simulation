"""Automated validation of the CW rendezvous project."""
import math

import numpy as np
import pytest
import scipy.linalg as sla

from src import config
from src.analysis import dv_curve, parameter_sweep
from src.config import Scenario
from src.cw_dynamics import cw_dynamics, cw_matrix_A, cw_stm
from src.nonlinear_check import nonlinear_miss
from src.orbital import mean_motion_rad_s, orbital_period_s
from src.rendezvous import solve_rendezvous
from src.simulation import STATUS_OK, STATUS_POS_ONLY, run_simulation

N = mean_motion_rad_s(400.0)


# Test 1 -- mean motion
def test_mean_motion_value_and_period():
    assert N == pytest.approx(1.1314e-3, rel=1e-3)
    assert orbital_period_s(400.0) / 60 == pytest.approx(92.56, abs=0.05)


def test_mean_motion_decreases_with_altitude():
    assert mean_motion_rad_s(800.0) < mean_motion_rad_s(400.0)


# Test 2 -- CW dynamics
def test_dynamics_structure_and_values():
    s = np.array([0.1, -0.5, 0.002, -0.001])
    d = cw_dynamics(0.0, s, N)
    assert len(d) == 4
    assert d[0] == s[2] and d[1] == s[3]                       # kinematics
    assert d[2] == pytest.approx(3 * N**2 * s[0] + 2 * N * s[3])
    assert d[3] == pytest.approx(-2 * N * s[2])
    assert np.allclose(d, cw_matrix_A(N) @ s)


def test_stm_matches_matrix_exponential_and_identity():
    assert np.allclose(cw_stm(N, 0.0), np.eye(4))
    for t in (300.0, 1800.0, 3600.0):
        assert np.allclose(cw_stm(N, t), sla.expm(cw_matrix_A(N) * t), atol=1e-9)


# Test 3 -- rendezvous solver
def test_default_rendezvous_meets_tolerances():
    res = run_simulation(Scenario())
    assert res.final_separation_km < config.POSITION_TOL_KM
    assert res.final_speed_km_s < config.VELOCITY_TOL_KM_S
    assert res.status == STATUS_OK
    assert res.stm_vs_ode_pos_err_km < 1e-6        # numerical vs analytic agree


def test_solution_satisfies_position_equation():
    sol = solve_rendezvous([0.2, -1.0], [0, 0], N, 1500.0)
    s_f = cw_stm(N, 1500.0) @ np.concatenate([[0.2, -1.0], sol.v0_plus])
    assert np.allclose(s_f[:2], 0.0, atol=1e-12)
    assert np.allclose(s_f[2:] + sol.dv2, 0.0, atol=1e-12)


def test_single_burn_does_not_report_rendezvous():
    """Negative control: position match without velocity match must not pass."""
    res = run_simulation(Scenario(), apply_terminal_burn=False)
    assert res.status == STATUS_POS_ONLY
    assert res.final_speed_km_s > 1e-4


def test_one_orbit_analytic_limit():
    """Near tf = T with x0 = 0: dv_total -> 2 n |y0| / (6 pi)."""
    T = orbital_period_s(400.0)
    sol = solve_rendezvous([0.0, -1.0], [0, 0], N, T * 0.999)
    assert sol.dv_total == pytest.approx(2 * N * 1.0 / (6 * math.pi), rel=0.01)


def test_dv_scales_linearly_with_separation():
    a = solve_rendezvous([0, -1.0], [0, 0], N, 1800).dv_total
    b = solve_rendezvous([0, -2.0], [0, 0], N, 1800).dv_total
    assert b == pytest.approx(2 * a, rel=1e-12)


def test_nonlinear_two_body_cross_check():
    miss_km, _ = nonlinear_miss(Scenario())
    assert miss_km * 1000 < 5.0        # metres: linearisation error for 1 km, 30 min


# Test 4 -- parameter sweep
def test_sweep_processes_all_cases():
    df = parameter_sweep()
    assert len(df) == 15
    assert set(df.transfer_time_min) == set(config.SWEEP_TRANSFER_TIMES_MIN)
    assert set(df.separation_km) == set(config.SWEEP_INITIAL_SEPARATIONS_KM)
    assert (df.status == STATUS_OK).all()


def test_dv_curve_has_nan_only_at_singularities():
    c = dv_curve(1.0, t_min_range=(5, 80), num=50)
    assert c.dv_total_m_s.notna().all()


# Test 5 -- invalid inputs
@pytest.mark.parametrize("sc, msg", [
    (Scenario(altitude_km=-100.0), "altitude"),
    (Scenario(altitude_km=0.0), "altitude"),
    (Scenario(transfer_time_s=0.0), "Transfer time"),
    (Scenario(transfer_time_s=-5.0), "Transfer time"),
    (Scenario(y0_km=float("nan")), "finite"),
    (Scenario(altitude_km=float("inf")), "finite"),
])
def test_invalid_scenarios_raise_clear_errors(sc, msg):
    with pytest.raises(ValueError, match=msg):
        run_simulation(sc)


def test_singular_transfer_time_raises():
    with pytest.raises(ValueError, match="singular"):
        solve_rendezvous([0, -1], [0, 0], N, 2 * math.pi / N)


def test_too_few_samples_raises():
    with pytest.raises(ValueError, match="samples"):
        run_simulation(Scenario(), n_steps=1)


def test_large_separation_warns_but_runs():
    res = run_simulation(Scenario(y0_km=-50.0))
    assert any("inaccurate" in w for w in res.warnings)
