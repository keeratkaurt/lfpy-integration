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
    parser.add_argument("--stim-amp", type=float, default=0.65, help="nA")
    parser.add_argument("--stim-delay", type=float, default=100.0, help="ms")
    parser.add_argument("--stim-dur", type=float, default=300.0, help="ms")
    parser.add_argument("--sigma", type=float, default=0.3, help="S/m")
    return parser.parse_args()


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
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cell = build_cell(args)

    LFPy.StimIntElectrode(
        cell,
        idx=0,
        pptype="IClamp",
        amp=args.stim_amp,
        dur=args.stim_dur,
        delay=args.stim_delay,
    )
    cell.simulate(rec_vmem=True, rec_imem=True)

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

    fig, ax = plt.subplots(figsize=(10, 5))
    for trace, distance in zip(extracellular, distances):
        ax.plot(cell.tvec, trace * 1_000.0, label=f"{distance:.0f} µm")
    ax.set(xlabel="Time (ms)", ylabel="Extracellular voltage (µV)",
           title="Simulated extracellular electrode recordings")
    ax.legend(title="Soma distance")
    fig.tight_layout()
    fig.savefig(args.output_dir / "extracellular_voltage.png", dpi=200)
    plt.close(fig)

    np.savez(
        args.output_dir / "lfpy_poc_output.npz",
        time=cell.tvec,
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
    print("Peak soma voltage (mV):", float(cell.vmem[0].max()))
    print("Peak absolute extracellular voltage (µV):", float(np.abs(extracellular).max() * 1_000.0))
    print("Output directory:", args.output_dir.resolve())


if __name__ == "__main__":
    main()
