"""One somatic stimulus, extracellular point electrodes near named BBP cell regions.

Run after compiling mechanisms in the BBP model directory:
  python lfpy_site_recordings.py MODEL_DIR --stim-file stims/5k50kInterChaoticB.csv
"""

import argparse
from pathlib import Path

import LFPy
import matplotlib.pyplot as plt
import numpy as np
from neuron import h

from lfpy_poc import build_cell


def point_to_segments(point, starts, ends, radii):
    """Approximate distance from a point to each cylindrical compartment surface."""
    direction = ends - starts
    length_sq = np.einsum("ij,ij->i", direction, direction)
    fraction = np.clip(
        np.einsum("ij,ij->i", point - starts, direction)
        / np.maximum(length_sq, 1e-12), 0, 1
    )
    nearest = starts + fraction[:, None] * direction
    return np.linalg.norm(nearest - point, axis=1) - radii


def select_probes(cell, offset_um):
    mid = np.column_stack((cell.x.mean(axis=1),
                           cell.y.mean(axis=1),
                           cell.z.mean(axis=1)))
    starts = np.column_stack((cell.x[:, 0], cell.y[:, 0], cell.z[:, 0]))
    ends = np.column_stack((cell.x[:, 1], cell.y[:, 1], cell.z[:, 1]))
    radii = np.asarray(cell.d) / 2
    targets = [
        ("Soma", "soma", 0.0),
        ("Apical dendrite", "apic", 650.0),
        ("Basal dendrite", "dend", -140.0),
        ("Axon", "axon", -550.0),
    ]
    points, labels, selected, nearest, clearances = [], [], [], [], []
    for label, group, target_y in targets:
        indices = np.asarray(cell.get_idx(group), dtype=int)
        if not len(indices):
            raise RuntimeError(f"No compartments found for {group}")
        # Prefer a section near the requested anatomical level and the soma x.
        rank = np.argsort(np.abs(mid[indices, 1] - target_y)
                          + 0.15 * np.abs(mid[indices, 0]))
        choice = None
        for idx in indices[rank[:min(80, len(indices))]]:
            for sign in (1, -1):
                for extra in (0, 5, 10, 20, 35, 50):
                    point = mid[idx].copy()
                    point[2] += sign * (offset_um + radii[idx] + extra)
                    distances = point_to_segments(point, starts, ends, radii)
                    nearest_idx = int(np.argmin(distances))
                    clearance = float(distances[nearest_idx])
                    if (clearance >= offset_um - 1
                            and (group == "soma" or nearest_idx in indices)):
                        choice = (point, idx, nearest_idx, clearance)
                        break
                if choice:
                    break
            if choice:
                break
        if choice is None:
            raise RuntimeError(f"Cannot safely place an external probe near {label}")
        point, idx, nearest_idx, clearance = choice
        points.append(point)
        labels.append(label)
        selected.append(idx)
        nearest.append(nearest_idx)
        clearances.append(clearance)
    return (np.asarray(points), labels, np.asarray(selected),
            np.asarray(nearest), np.asarray(clearances))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("model_dir", type=Path)
    parser.add_argument("--stim-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("site_recordings"))
    parser.add_argument("--dt", type=float, default=0.1, help="ms per CSV sample")
    parser.add_argument("--offset-um", type=float, default=20.0)
    parser.add_argument("--sigma", type=float, default=0.3, help="S/m")
    args = parser.parse_args()
    if args.dt <= 0 or args.offset_um <= 0:
        parser.error("dt and offset-um must be positive")
    command = np.genfromtxt(args.stim_file, dtype=float).squeeze()
    if command.ndim != 1 or not len(command) or not np.isfinite(command).all():
        parser.error("Stimulus must be one nonempty column of finite current values")
    cell = build_cell(argparse.Namespace(model_dir=args.model_dir, dt=args.dt,
                                         tstop=len(command) * args.dt,
                                         v_init=-68.0))
    soma_idx = np.asarray(cell.get_idx("soma"), dtype=int)
    assert len(soma_idx) == 1 and soma_idx[0] == 0
    points, labels, selected, nearest, clearances = select_probes(
        cell, args.offset_um)
    # Exactly one somatic current clamp. The probes below inject no current.
    clamp = LFPy.StimIntElectrode(cell, idx=int(soma_idx[0]), pptype="IClamp",
                                 amp=0.0, dur=len(command) * args.dt,
                                 delay=0.0, record_current=True)
    hoc_clamp = cell._hoc_stimlist.o(clamp.hocidx)
    vector = h.Vector().from_python(command)
    vector.play(hoc_clamp._ref_amp, args.dt)
    cell.simulate(rec_vmem=True, rec_imem=True)
    electrode = LFPy.RecExtElectrode(cell=cell, sigma=args.sigma,
                                    x=points[:, 0], y=points[:, 1],
                                    z=points[:, 2], method="linesource")
    extracellular_mV = electrode.get_transformation_matrix() @ cell.imem
    time = np.asarray(cell.tvec)
    current = np.asarray(clamp.i)
    soma_mV = np.asarray(cell.vmem[soma_idx[0]])
    if extracellular_mV.shape != (len(labels), len(time)):
        raise AssertionError("Unexpected extracellular output dimensions")
    if not np.isfinite(extracellular_mV).all():
        raise AssertionError("Extracellular recording has nonfinite samples")
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 10))
    for (xs, xe), (ys, ye) in zip(cell.x, cell.y):
        ax.plot([xs, xe], [ys, ye], color="0.55", lw=0.4)
    colors = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd"]
    for point, label, color in zip(points, labels, colors):
        ax.scatter(point[0], point[1], c=color, marker="x", s=90, zorder=5)
        ax.annotate(label, point[:2], xytext=(6, 4),
                    textcoords="offset points", color=color)
    ax.scatter(cell.somapos[0], cell.somapos[1], c="black", s=18, zorder=4)
    ax.set(xlabel="x (µm)", ylabel="y (µm)",
           title="L5_TTPC1: extracellular probe locations (x–y projection)")
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    fig.savefig(output / "site_electrode_placement.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(6, 1, figsize=(12, 12), sharex=True)
    axes[0].plot(time, current, color="black", lw=0.8)
    axes[0].set(ylabel="nA", title="One current stimulus: soma only")
    axes[1].plot(time, soma_mV, color="0.3", lw=0.9)
    axes[1].set(ylabel="mV", title="Intracellular: soma")
    for ax, trace, label, color in zip(axes[2:], extracellular_mV, labels, colors):
        ax.plot(time, trace * 1000, color=color, lw=0.9)
        ax.set(ylabel="µV", title=f"Extracellular: near {label.lower()}")
    axes[-1].set_xlabel("Time (ms)")
    for ax in axes:
        ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(output / "site_extracellular_recordings.png", dpi=200)
    plt.close(fig)

    np.savez(output / "site_recordings.npz", time_ms=time,
             somatic_stimulus_nA=current, stimulus_command_nA=command,
             soma_intracellular_mV=soma_mV,
             extracellular_mV=extracellular_mV,
             electrode_xyz_um=points, electrode_region=np.asarray(labels),
             target_compartment=selected, nearest_compartment=nearest,
             surface_clearance_um=clearances,
             dt_ms=args.dt, sigma_S_per_m=args.sigma)
    print("Cell compartments:", cell.totnsegs)
    print("Time samples:", len(time), "dt:", args.dt)
    print("Soma spike crossings:", np.sum((soma_mV[:-1] < 0) & (soma_mV[1:] >= 0)))
    for label, point, gap, trace in zip(labels, points, clearances,
                                        extracellular_mV):
        print(label, "xyz:", np.round(point, 2),
              "minimum clearance µm:", round(gap, 2),
              "peak |V| µV:", round(np.max(np.abs(trace)) * 1000, 3))
    print("Results:", output.resolve())


if __name__ == "__main__":
    main()
