"""Streamlit interface: planar CW rendezvous simulator.  Run: streamlit run app/main.py"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st

from src import config, visualization as viz
from src.config import M_PER_KM, Scenario
from src.simulation import STATUS_OK, run_simulation

st.set_page_config(page_title="CW Rendezvous Simulator", layout="wide")
st.title("Spacecraft Rendezvous in LEO — Clohessy–Wiltshire Simulator")
st.caption("Planar, linearised relative motion about a circular orbit. Two impulsive burns. "
           "Educational/preliminary analysis only — not for mission planning.")

with st.sidebar:
    st.header("Mission parameters")
    alt = st.number_input("Target altitude [km]", value=config.DEFAULT_ALTITUDE_KM, step=10.0)
    x0 = st.number_input("Initial radial position x0 [km]", value=config.DEFAULT_X0_KM, step=0.1,
                         format="%.3f")
    y0 = st.number_input("Initial along-track position y0 [km]", value=config.DEFAULT_Y0_KM,
                         step=0.1, format="%.3f")
    tmin = st.number_input("Transfer time [min]", value=config.DEFAULT_TRANSFER_TIME_MIN, step=5.0)
    n_steps = st.slider("Time resolution (output samples)", 50, 2000, config.DEFAULT_N_STEPS, 50)
    burn2 = st.checkbox("Apply terminal braking burn (burn 2)", value=True,
                        help="Untick to see that matching position alone is not a rendezvous.")
    want_anim = st.checkbox("Generate animation (slower)", value=False)
    run = st.button("Run Simulation", type="primary")

if run:
    try:
        st.session_state["res"] = run_simulation(
            Scenario(alt, x0, y0, tmin * config.SECONDS_PER_MINUTE), n_steps, burn2)
        st.session_state["anim"] = None
        if want_anim:
            with st.spinner("Rendering animation..."):
                path = Path(tempfile.gettempdir()) / "cw_rendezvous.gif"
                viz.make_animation(st.session_state["res"], path)
                st.session_state["anim"] = path.read_bytes()
        st.session_state["err"] = None
    except (ValueError, RuntimeError) as exc:
        st.session_state["res"], st.session_state["err"] = None, str(exc)

if st.session_state.get("err"):
    st.error(st.session_state["err"])
res = st.session_state.get("res")
if res is None:
    if not st.session_state.get("err"):
        st.info("Set parameters in the sidebar and press **Run Simulation**.")
    st.stop()

for w in res.warnings:
    st.warning(w)
s = res.solution
st.subheader("Mission summary")
(st.success if res.status == STATUS_OK else st.error)(f"Status: {res.status}")
c = st.columns(4)
c[0].metric("Mean motion n", f"{res.n_rad_s:.4e} rad/s")
c[1].metric("Orbital period", f"{res.period_s/60:.2f} min")
c[2].metric("Burn-1 velocity (radial, along-track)",
            f"({s.v0_plus[0]*M_PER_KM:.4f}, {s.v0_plus[1]*M_PER_KM:.4f}) m/s")
c[3].metric("Total Δv (Δv1 + Δv2)", f"{s.dv_total*M_PER_KM:.4f} m/s")
c = st.columns(4)
c[0].metric("Δv1 (t = 0)", f"{s.dv1_mag*M_PER_KM:.4f} m/s")
c[1].metric("Δv2 (t = tf)", f"{s.dv2_mag*M_PER_KM:.4f} m/s")
c[2].metric("Final separation", f"{res.final_separation_km*M_PER_KM:.2e} m")
c[3].metric("Final relative velocity", f"{res.final_speed_km_s*M_PER_KM:.2e} m/s",
            help=f"Arrival speed before burn 2: {res.arrival_speed_km_s*M_PER_KM:.4f} m/s")
st.caption("Δv is the required velocity change under the idealised impulsive two-burn model, "
           "not a complete propulsion budget.")

left, right = st.columns([3, 2])
with left:
    st.subheader("Relative trajectory")
    fig = viz.plot_trajectory(res); st.pyplot(fig); plt.close(fig)
with right:
    st.subheader("Separation vs time")
    fig = viz.plot_separation(res); st.pyplot(fig); plt.close(fig)
    st.subheader("Relative velocity vs time")
    fig = viz.plot_relative_velocity(res); st.pyplot(fig); plt.close(fig)
if st.session_state.get("anim"):
    st.subheader("Animation")
    st.image(st.session_state["anim"])
