"""Audit saved site recordings without rerunning NEURON (requires NumPy only).

python verify_site_recordings.py
Writes JSON beside the input; does not modify the original arrays.
Geometry fields are checked for consistency, not independently reconstructed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path,
                        default=Path("results_site_recordings/site_recordings.npz"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    with np.load(args.input, allow_pickle=False) as saved:
        d = {key: saved[key] for key in saved.files}
    t = d["time_ms"]
    soma = d["soma_intracellular_mV"]
    voltage = d["extracellular_mV"] * 1000.0
    labels = d["electrode_region"].tolist()
    dt = float(d["dt_ms"])
    require(t.ndim == 1 and len(t) > 1, "Invalid time axis")
    require(soma.shape == t.shape, "Soma/time shape mismatch")
    require(voltage.shape == (len(labels), len(t)), "Voltage/time shape mismatch")
    for key, value in d.items():
        if np.issubdtype(value.dtype, np.number):
            require(np.isfinite(value).all(), "Nonfinite data: " + key)
    require(dt > 0 and np.all(np.diff(t) > 0), "Invalid timestep")
    require(np.allclose(np.diff(t), dt, rtol=0, atol=1e-8), "Nonuniform time")
    require(float(d["sigma_S_per_m"]) > 0, "Invalid conductivity")
    require(d["electrode_xyz_um"].shape == (len(labels), 3), "Invalid coordinates")
    require(d["surface_clearance_um"].shape == (len(labels),), "Invalid clearances")
    require(np.all(d["surface_clearance_um"] > 0), "Nonpositive saved clearance")
    for key in ("target_compartment", "nearest_compartment"):
        require(d[key].shape == (len(labels),) and
                np.issubdtype(d[key].dtype, np.integer) and np.all(d[key] >= 0),
                "Invalid saved compartment indices: " + key)
    current = d["somatic_stimulus_nA"]
    command = d["stimulus_command_nA"]
    require(current.shape == t.shape, "Current/time shape mismatch")
    require(command.ndim == 1 and len(command) == len(t) - 1,
            "Expected one fewer command sample than time samples")
    indices = np.flatnonzero((soma[:-1] < 0) & (soma[1:] >= 0)) + 1
    crossing_times = [
        float(t[i-1] - soma[i-1] / (soma[i] - soma[i-1]) * (t[i] - t[i-1]))
        for i in indices
    ]

    def metrics(trace, start, stop):
        mask = (t >= start - 1e-7) & (t <= stop + 1e-7)
        ix = np.flatnonzero(mask)
        require(len(ix) > 0, "Empty measurement window")
        peak = ix[np.argmax(np.abs(trace[ix]))]
        return dict(window_start_ms=float(start), window_end_ms=float(stop),
                    samples=int(len(ix)), peak_abs_uV=float(abs(trace[peak])),
                    peak_signed_uV=float(trace[peak]), peak_time_ms=float(t[peak]),
                    peak_to_peak_uV=float(np.ptp(trace[ix])))

    results = []
    for k, label in enumerate(labels):
        results.append(dict(
            region=label, xyz_um=d["electrode_xyz_um"][k].tolist(),
            saved_surface_clearance_um=float(d["surface_clearance_um"][k]),
            target_compartment=int(d["target_compartment"][k]),
            nearest_compartment=int(d["nearest_compartment"][k]),
            full_recording=metrics(voltage[k], t[0], t[-1]),
            displayed_interval=metrics(voltage[k], 100, 500),
            spike_windows=[
                dict(soma_crossing_ms=c, **metrics(voltage[k], c-2, c+3))
                for c in crossing_times if c-2 >= t[0] and c+3 <= t[-1]
            ]))
    report = dict(
        input=str(args.input),
        input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
        soma_zero_crossings_ms=crossing_times,
        soma_voltage_range_mV=[float(soma.min()), float(soma.max())],
        dt_ms=dt, sigma_S_per_m=float(d["sigma_S_per_m"]),
        unshifted_command_error_nA=float(np.max(np.abs(current[:-1]-command))),
        one_sample_lag_command_error_nA=float(np.max(np.abs(current[1:]-command))),
        electrodes=results,
        limitations=[
            "Spike windows are -2 to +3 ms around interpolated somatic upward "
            "zero crossings; raw window amplitudes are not baseline corrected "
            "and do not isolate spike currents.",
            "Saved clearances are positive, but morphology was not rebuilt. "
            "This does not independently verify the saved geometry.",
            "The NPZ does not record clamp count or clamp compartment; soma-only "
            "stimulation is supported by source-code inspection.",
            "No re-simulation, saved-figure pixel comparison, or noise analysis."
        ])
    output = args.output or args.input.with_name("verification_metrics.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
