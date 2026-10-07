"""Sweep transfer time x initial separation; save CSV, Markdown table and trade-off plot."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")

from src import config, visualization as viz                       # noqa: E402
from src.analysis import dv_curve, parameter_sweep, to_markdown_table  # noqa: E402
from src.orbital import orbital_period_s                           # noqa: E402


def main() -> None:
    out = ROOT / "results"
    out.mkdir(exist_ok=True)

    df = parameter_sweep()

    # Save CSV using UTF-8
    df.to_csv(
        out / "sweep_results.csv",
        index=False,
        encoding="utf-8"
    )

    # Save Markdown using UTF-8
    (out / "sweep_table.md").write_text(
        to_markdown_table(df) + "\n",
        encoding="utf-8"
    )

    print(to_markdown_table(df))

    curves = {
        s: dv_curve(s)
        for s in config.SWEEP_INITIAL_SEPARATIONS_KM
    }

    ref = curves[1.0].dropna()

    i = ref.dv_total_m_s.idxmin()

    print(
        f"\nMinimum total dv on dense curve (1 km): "
        f"{ref.dv_total_m_s[i]:.4f} m/s "
        f"at {ref.transfer_time_min[i]:.1f} min"
    )

    ref.to_csv(
        out / "dv_curve_1km.csv",
        index=False,
        encoding="utf-8"
    )

    viz.save_figure(
        viz.plot_dv_vs_time(
            curves,
            df,
            orbital_period_s(config.DEFAULT_ALTITUDE_KM) / 60
        ),
        out / "dv_vs_time.png"
    )

    print(f"Saved results to {out}")


if __name__ == "__main__":
    main()
