"""Run an LFPy extracellular-signal proof of concept on a BBP cell model."""

import argparse
import os
import re
from pathlib import Path

import LFPy
import matplotlib.pyplot as plt
import numpy as np
from neuron import h, load_mechanisms

from extracellular import compute_extracellular


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("model_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("lfpy_results"))
    parser.add_argument("--dt", type=float, default=None,
                        help="ms; defaults to 0.1 for a CSV or 0.025 for a pulse")
    parser.add_argument("--tstop", type=float, default=500.0)
    parser.add_argument("--v-init", type=float, default=-68.0)
    parser.add_argument(
        "--stim-file",
        type=Path,
        help="CSV stimulus waveform in nA; overrides the constant-current pulse",
    )
    parser.add_argument("--stim-multiplier", type=float, default=1.0)
    parser.add_argument("--stim-dc-offset", type=float, default=0.0, help="nA")
    parser.add_argument("--stim-amp", type=float, default=0.65, help="nA")
    parser.add_argument("--stim-delay", type=float, default=100.0, help="ms")
    parser.add_argument("--stim-dur", type=float, default=300.0, help="ms")
    parser.add_argument("--sigma", type=float, default=0.3, help="S/m")
    parser.add_argument("--outside-cell", action="store_true",
                        help="Place four electrodes beyond the full x-y morphology bounds")
    parser.add_argument("--outside-margin", type=float, default=20.0,
                        help="Distance in um beyond the x-y bounds")
    return parser.parse_args()


def load_stimulus(args):
    if args.stim_file is None:
        return None

    stim_file = args.stim_file.resolve()
    if not stim_file.is_file():
        raise FileNotFoundError(f"Stimulus file not found: {stim_file}")

    values = np.genfromtxt(stim_file, dtype=np.float64)
    values = np.asarray(values).squeeze()
    if values.ndim != 1 or values.size == 0:
        raise ValueError(f"Expected one non-empty stimulus column in {stim_file}")
    if not np.isfinite(values).all():
        raise ValueError(f"Stimulus contains non-finite values: {stim_file}")

    values = values * args.stim_multiplier + args.stim_dc_offset
    args.stim_file = stim_file
    # Match DL4neurons2, which runs one dt beyond the final CSV sample.
    args.tstop = values.size * args.dt
    return values


def find_template_name(template_file):
    match = re.search(
        r"^\s*begintemplate\s+(\S+)",
        template_file.read_text(),
        flags=re.MULTILINE,
    )
    if match is None:
        raise RuntimeError(f"Could not find begintemplate in {template_file}")
    return match.group(1)


def find_morphology(model_dir):
    morphologies = sorted((model_dir / "morphology").glob("*.asc"))
    if len(morphologies) != 1:
        raise RuntimeError(
            f"Expected exactly one ASC morphology in {model_dir / 'morphology'}, "
            f"found {len(morphologies)}"
        )
    return morphologies[0]


def build_cell(args):
    model_dir = args.model_dir.resolve()
    if not load_mechanisms(str(model_dir)):
        raise RuntimeError(
            f"Compiled NEURON mechanisms were not found in {model_dir}. "
            "Run `nrnivmodl mechanisms` there first."
        )

    h.load_file("stdrun.hoc")
    h.load_file("import3d.hoc")
    h.load_file(str(model_dir / "constants.hoc"))

    template_file = model_dir / "template.hoc"
    previous_dir = Path.cwd()
    os.chdir(model_dir)
    try:
        return LFPy.TemplateCell(
            morphology=str(find_morphology(model_dir)),
            templatefile=str(template_file),
            templatename=find_template_name(template_file),
            templateargs=0,
            dt=args.dt,
            tstop=args.tstop,
            v_init=args.v_init,
            passive=False,
            delete_sections=False,
            verbose=False,
        )
    finally:
        os.chdir(previous_dir)


def main():
    args = parse_args()
    if args.dt is None:
        args.dt = 0.1 if args.stim_file is not None else 0.025
    if args.dt <= 0 or args.outside_margin <= 0:
        raise ValueError("dt and outside-margin must be positive")
    stimulus_values = load_stimulus(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cell = build_cell(args)

    if stimulus_values is None:
        stimulus = LFPy.StimIntElectrode(
            cell,
            idx=0,
            pptype="IClamp",
            amp=args.stim_amp,
            dur=args.stim_dur,
            delay=args.stim_delay,
            record_current=True,
        )
        stimulus_vector = None
    else:
        stimulus = LFPy.StimIntElectrode(
            cell,
            idx=0,
            pptype="IClamp",
            amp=0.0,
            dur=args.tstop,
            delay=0.0,
            record_current=True,
        )
        hoc_stimulus = cell._hoc_stimlist.o(stimulus.hocidx)
        stimulus_vector = h.Vector().from_python(stimulus_values)
        stimulus_vector.play(hoc_stimulus._ref_amp, args.dt)

    cell.simulate(rec_vmem=True, rec_imem=True)
    applied_stimulus = np.asarray(stimulus.i)

    soma_xyz = np.asarray(cell.somapos)
    if args.outside_cell:
        xlo, xhi = np.min(cell.x), np.max(cell.x)
        ylo, yhi = np.min(cell.y), np.max(cell.y)
        m = args.outside_margin
        electrode_x = np.array([xlo-m, xhi+m, soma_xyz[0], soma_xyz[0]])
        electrode_y = np.array([soma_xyz[1], soma_xyz[1], yhi+m, ylo-m])
        electrode_z = np.full(4, soma_xyz[2])
        labels = ["Left", "Right", "Above apical", "Below basal"]
        assert electrode_x[0] < xlo and electrode_x[1] > xhi
        assert electrode_y[2] > yhi and electrode_y[3] < ylo
    else:
        distances = np.array([20.0, 50.0, 100.0, 200.0])
        electrode_x = soma_xyz[0] + distances
        electrode_y = np.full(distances.shape, soma_xyz[1])
        electrode_z = np.full(distances.shape, soma_xyz[2])
        labels = [f"{d:.0f} µm" for d in distances]
    extracellular = compute_extracellular(
        cell, electrode_x, electrode_y, electrode_z, sigma=args.sigma
    )

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(cell.tvec, cell.vmem[0])
    ax.set(xlabel="Time (ms)", ylabel="Membrane voltage (mV)",
           title="Somatic membrane voltage")
    fig.tight_layout()
    fig.savefig(args.output_dir / "intracellular_voltage.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(cell.tvec, applied_stimulus)
    ax.set(
        xlabel="Time (ms)",
        ylabel="Injected current (nA)",
        title=(
            f"Applied stimulus: {args.stim_file.name}"
            if args.stim_file is not None
            else "Applied constant-current stimulus"
        ),
    )
    fig.tight_layout()
    fig.savefig(args.output_dir / "stimulus_current.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5))
    for trace, label in zip(extracellular, labels):
        ax.plot(cell.tvec, trace * 1_000.0, label=label)
    ax.set(xlabel="Time (ms)", ylabel="Extracellular voltage (µV)",
           title="Simulated extracellular electrode recordings")
    ax.legend(title="Electrode")
    fig.tight_layout()
    fig.savefig(args.output_dir / "extracellular_voltage.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    axes[0].plot(cell.tvec, applied_stimulus)
    axes[0].set(ylabel="Current (nA)", title="InterChaoticB LFPy proof of concept")
    axes[1].plot(cell.tvec, cell.vmem[0])
    axes[1].set(ylabel="Soma voltage (mV)")
    for trace, label in zip(extracellular, labels):
        axes[2].plot(cell.tvec, trace * 1_000.0, label=label)
    axes[2].set(xlabel="Time (ms)", ylabel="Extracellular (µV)")
    axes[2].legend(title="Electrode", ncol=4)
    fig.tight_layout()
    fig.savefig(args.output_dir / "interchaoticb_summary.png", dpi=200)
    plt.close(fig)

    if args.outside_cell:
        fig, ax = plt.subplots(figsize=(7, 10))
        for (xs, xe), (ys, ye) in zip(cell.x, cell.y):
            ax.plot([xs, xe], [ys, ye], color="0.35", linewidth=0.4)
        ax.scatter([soma_xyz[0]], [soma_xyz[1]], c="orange", s=25,
                   label="Soma", zorder=4)
        for x, y, label in zip(electrode_x, electrode_y, labels):
            ax.scatter([x], [y], c="crimson", marker="x", s=65, zorder=5)
            ax.annotate(label, (x, y), xytext=(5, 5),
                        textcoords="offset points", fontsize=8)
        ax.set(xlabel="x (µm)", ylabel="y (µm)",
               title="L5_TTPC1 morphology and outside-cell electrodes")
        ax.set_aspect("equal", adjustable="box")
        fig.tight_layout()
        fig.savefig(args.output_dir / "outside_cell_electrode_placement.png", dpi=200)
        plt.close(fig)
        fig, axes = plt.subplots(4, 1, figsize=(11, 8), sharex=True)
        for ax, trace, label in zip(axes, extracellular, labels):
            ax.plot(cell.tvec, trace * 1000)
            ax.set(ylabel="µV", title=label)
        axes[-1].set_xlabel("Time (ms)")
        fig.tight_layout()
        fig.savefig(args.output_dir / "outside_cell_voltage_traces.png", dpi=200)
        plt.close(fig)

    np.savez(
        args.output_dir / "lfpy_poc_output.npz",
        time=cell.tvec,
        stimulus_nA=applied_stimulus,
        stimulus_command_nA=(
            stimulus_values if stimulus_values is not None else applied_stimulus
        ),
        stimulus_source=str(args.stim_file) if args.stim_file is not None else "constant",
        stimulus_multiplier=args.stim_multiplier,
        stimulus_dc_offset_nA=args.stim_dc_offset,
        intracellular=cell.vmem[0],
        extracellular_mV=extracellular,
        electrode_x=electrode_x,
        electrode_y=electrode_y,
        electrode_z=electrode_z,
        electrode_label=np.asarray(labels),
        dt_ms=args.dt,
    )

    print("SUCCESS")
    print("Number of compartments:", cell.totnsegs)
    print("Transmembrane current shape:", cell.imem.shape)
    print("Extracellular voltage shape:", extracellular.shape)
    print("Stimulus current shape:", applied_stimulus.shape)
    print("Stimulus source:", args.stim_file or "constant current")
    print("Stimulus range (nA):", float(applied_stimulus.min()), "to", float(applied_stimulus.max()))
    print("Peak soma voltage (mV):", float(cell.vmem[0].max()))
    print("Peak absolute extracellular voltage (µV):", float(np.abs(extracellular).max() * 1_000.0))
    print("Output directory:", args.output_dir.resolve())
    print("Electrodes (x, y, z) um:", list(zip(electrode_x, electrode_y, electrode_z)))
    print("Peak absolute per-electrode (µV):", np.max(np.abs(extracellular), axis=1) * 1000)


if __name__ == "__main__":
    main()
