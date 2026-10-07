"""Circular-orbit quantities for the target spacecraft (spherical Earth)."""
import math

from .config import EARTH_MU_KM3_S2, EARTH_RADIUS_KM


def _check_altitude(altitude_km: float) -> None:
    if not math.isfinite(altitude_km) or altitude_km <= 0.0:
        raise ValueError(
            f"Target altitude must be a finite value > 0 km (got {altitude_km!r})."
        )


def orbit_radius_km(altitude_km: float) -> float:
    """Orbit radius r = R_E + h [km]."""
    _check_altitude(altitude_km)
    return EARTH_RADIUS_KM + altitude_km


def mean_motion_rad_s(altitude_km: float) -> float:
    """Mean motion n = sqrt(mu / r^3) [rad/s]."""
    r = orbit_radius_km(altitude_km)
    return math.sqrt(EARTH_MU_KM3_S2 / r**3)


def orbital_period_s(altitude_km: float) -> float:
    """Orbital period T = 2*pi / n [s]."""
    return 2.0 * math.pi / mean_motion_rad_s(altitude_km)


def circular_speed_km_s(altitude_km: float) -> float:
    """Circular orbital speed sqrt(mu / r) [km/s] (used only by nonlinear check)."""
    return math.sqrt(EARTH_MU_KM3_S2 / orbit_radius_km(altitude_km))
