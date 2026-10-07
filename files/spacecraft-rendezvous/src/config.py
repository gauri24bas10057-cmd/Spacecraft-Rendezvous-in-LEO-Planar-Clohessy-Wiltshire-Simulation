"""Central constants and default scenario for the CW rendezvous project.

Units (internal, used everywhere unless a name says otherwise):
    distance km | time s | velocity km/s | acceleration km/s^2
"""
from dataclasses import dataclass

# --- Physical constants -----------------------------------------------------
EARTH_RADIUS_KM: float = 6378.0
EARTH_MU_KM3_S2: float = 398600.4418

# --- Unit conversions -------------------------------------------------------
SECONDS_PER_MINUTE: float = 60.0
M_PER_KM: float = 1000.0

# --- Default scenario (brief, section 2-3) ----------------------------------
DEFAULT_ALTITUDE_KM: float = 400.0
DEFAULT_X0_KM: float = 0.0     # radial offset of chaser from target
DEFAULT_Y0_KM: float = -1.0    # along-track offset (negative = behind target)
DEFAULT_TRANSFER_TIME_MIN: float = 30.0

# --- Parameter sweep (brief, section 8) -------------------------------------
SWEEP_TRANSFER_TIMES_MIN = (10.0, 20.0, 30.0, 40.0, 60.0)
SWEEP_INITIAL_SEPARATIONS_KM = (0.5, 1.0, 2.0)

# --- Rendezvous acceptance criteria -----------------------------------------
# Checked AFTER the terminal (braking) burn; see report for rationale.
POSITION_TOL_KM: float = 1e-6        # 1 mm
VELOCITY_TOL_KM_S: float = 1e-6      # 1 mm/s

# --- Numerical integration --------------------------------------------------
INTEGRATOR_METHOD: str = "DOP853"
INTEGRATOR_RTOL: float = 1e-10
INTEGRATOR_ATOL: float = 1e-12       # km, km/s
DEFAULT_N_STEPS: int = 600           # output samples per trajectory


@dataclass(frozen=True)
class Scenario:
    """Inputs defining one rendezvous case."""
    altitude_km: float = DEFAULT_ALTITUDE_KM
    x0_km: float = DEFAULT_X0_KM
    y0_km: float = DEFAULT_Y0_KM
    transfer_time_s: float = DEFAULT_TRANSFER_TIME_MIN * SECONDS_PER_MINUTE
