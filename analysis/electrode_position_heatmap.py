"""Create an electrode-position amplitude heatmap from saved site recordings.

This is an inverse-distance visualization of the four measured probe values.
It is not a continuous physical field estimate; the measured points are shown
on top of the interpolation so the limited spatial sampling is explicit.

python analysis/electrode_position_heatmap.py
"""
from pathlib import Path
import argparse
import numpy as np
import matplotlib.pyplot as plt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path,
                   default=Path("results_site_recordings/site_recordings.npz"))
    p.add_argument("--output", type=Path,
                   default=Path("results_site_recordings/electrode_position_heatmap.png"))
    args = p.parse_args()
    with np.load(args.input, allow_pickle=False) as d:
        xyz = d["electrode_xyz_um"]
        v = d["extracellular_mV"]
        labels = d["electrode_region"].astype(str)
    peak_uv = np.max(np.abs(v), axis=1) * 1000.0
    x, y = xyz[:, 0], xyz[:, 1]
    pad = 0.15 * max(np.ptp(x), np.ptp(y), 1.0)
    gx = np.linspace(x.min() - pad, x.max() + pad, 300)
    gy = np.linspace(y.min() - pad, y.max() + pad, 300)
    xx, yy = np.meshgrid(gx, gy)
    dist = np.sqrt((xx[..., None] - x)**2 + (yy[..., None] - y)**2)
    weights = 1.0 / np.maximum(dist, 1e-9)**2
    field = np.sum(weights * peak_uv, axis=-1) / np.sum(weights, axis=-1)
    fig, ax = plt.subplots(figsize=(8, 7))
    image = ax.pcolormesh(xx, yy, field, shading="auto", cmap="viridis")
    points = ax.scatter(x, y, c=peak_uv, cmap="viridis", edgecolor="white",
                        linewidth=1.2, s=90, zorder=3)
    for px, py, label, amp in zip(x, y, labels, peak_uv):
        ax.annotate(f"{label}\n{amp:.2f} µV", (px, py), xytext=(6, 6),
                    textcoords="offset points", fontsize=8)
    fig.colorbar(image, ax=ax, label="Peak absolute voltage (µV)")
    ax.set(xlabel="x position (µm)", ylabel="y position (µm)",
           title="Electrode-position amplitude map")
    ax.set_aspect("equal", adjustable="box")
    fig.text(0.5, 0.01, "Inverse-distance visualization of four measured probes; "
             "not a continuous physical-field reconstruction.",
             ha="center", fontsize=8)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=240)
    plt.close(fig)
    print("Saved:", args.output)
    for label, amp in zip(labels, peak_uv):
        print(f"{label}: {amp:.6f} µV")


if __name__ == "__main__":
    main()
