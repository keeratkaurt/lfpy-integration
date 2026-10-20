"""Compare extracellular amplitudes across saved electrode-offset runs.

Expected input layout:
  results_offset_20/site_recordings.npz
  results_offset_50/site_recordings.npz
  results_offset_100/site_recordings.npz

The script accepts --run OFFSET=PATH repeatedly. It does not treat a single
run as a robustness comparison; at least two runs are required.

python analysis/robustness_comparison.py \
  --run 20=results_site_recordings/site_recordings.npz \
  --run 50=results_offset_50/site_recordings.npz
"""
from pathlib import Path
import argparse
import csv
import numpy as np
import matplotlib.pyplot as plt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", action="append", required=True,
                   help="OFFSET_UM=PATH_TO_SITE_RECORDINGS_NPZ")
    p.add_argument("--output-dir", type=Path,
                   default=Path("results_site_recordings/robustness"))
    args = p.parse_args()
    if len(args.run) < 2:
        p.error("Provide at least two independent offset runs")
    rows = []
    for spec in args.run:
        if "=" not in spec:
            p.error("--run must use OFFSET=PATH")
        offset, raw_path = spec.split("=", 1)
        path = Path(raw_path)
        with np.load(path, allow_pickle=False) as d:
            labels = d["electrode_region"].astype(str)
            traces_uv = d["extracellular_mV"] * 1000.0
            t = d["time_ms"]
        for label, trace in zip(labels, traces_uv):
            mask = t >= 100
            ix = np.flatnonzero(mask)[np.argmax(np.abs(trace[mask]))]
            rows.append({
                "offset_um": float(offset),
                "source": str(path),
                "region": label,
                "peak_abs_uV_100_500ms": float(abs(trace[ix])),
                "peak_time_ms": float(t[ix]),
            })
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "robustness_metrics.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    fig, ax = plt.subplots(figsize=(8, 5))
    offsets = sorted({r["offset_um"] for r in rows})
    labels = list(dict.fromkeys(r["region"] for r in rows))
    for label in labels:
        y = [next(r["peak_abs_uV_100_500ms"] for r in rows
                  if r["offset_um"] == offset and r["region"] == label)
             for offset in offsets]
        ax.plot(offsets, y, marker="o", label=label)
    ax.set(xlabel="Electrode offset (µm)",
           ylabel="Peak absolute extracellular voltage (µV)",
           title="Robustness: amplitude versus electrode offset")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output_dir / "robustness_amplitude_vs_offset.png", dpi=220)
    plt.close(fig)
    print("Saved", args.output_dir / "robustness_metrics.csv")


if __name__ == "__main__":
    main()
