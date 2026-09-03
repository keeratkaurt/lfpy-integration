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
    parser.add_argument("--dt", type=float, default=0.025)
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

    distances = np.array([20.0, 50.0, 100.0, 200.0])
    soma_xyz = np.asarray(cell.somapos)
    electrode_x = soma_xyz[0] + distances
    electrode_y = np.full(distances.shape, soma_xyz[1])
    electrode_z = np.full(distances.shape, soma_xyz[2])
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
    for trace, distance in zip(extracellular, distances):
        ax.plot(cell.tvec, trace * 1_000.0, label=f"{distance:.0f} µm")
    ax.set(xlabel="Time (ms)", ylabel="Extracellular voltage (µV)",
           title="Simulated extracellular electrode recordings")
    ax.legend(title="Soma distance")
    fig.tight_layout()
    fig.savefig(args.output_dir / "extracellular_voltage.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    axes[0].plot(cell.tvec, applied_stimulus)
    axes[0].set(ylabel="Current (nA)", title="InterChaoticB LFPy proof of concept")
    axes[1].plot(cell.tvec, cell.vmem[0])
    axes[1].set(ylabel="Soma voltage (mV)")
    for trace, distance in zip(extracellular, distances):
        axes[2].plot(cell.tvec, trace * 1_000.0, label=f"{distance:.0f} µm")
    axes[2].set(xlabel="Time (ms)", ylabel="Extracellular (µV)")
    axes[2].legend(title="Soma distance", ncol=4)
    fig.tight_layout()
    fig.savefig(args.output_dir / "interchaoticb_summary.png", dpi=200)
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
        electrode_distance_um=distances,
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


if __name__ == "__main__":
    main()
