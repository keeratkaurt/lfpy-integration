"""Measure baseline noise, detect somatic spikes, and summarize extracellular SNR.

python analysis/noise_spike_analysis.py
"""
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib.pyplot as plt


def window_mask(t, start, stop):
    return (t >= start) & (t <= stop)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path,
                   default=Path("results_site_recordings/site_recordings.npz"))
    p.add_argument("--output-dir", type=Path,
                   default=Path("results_site_recordings/noise_spike_analysis"))
    p.add_argument("--baseline-start-ms", type=float, default=100.0)
    p.add_argument("--baseline-stop-ms", type=float, default=280.0)
    p.add_argument("--pre-ms", type=float, default=2.0)
    p.add_argument("--post-ms", type=float, default=3.0)
    args = p.parse_args()
    with np.load(args.input, allow_pickle=False) as d:
        t = d["time_ms"]
        soma = d["soma_intracellular_mV"]
        ext_uv = d["extracellular_mV"] * 1000.0
        labels = d["electrode_region"].astype(str).tolist()
    base = window_mask(t, args.baseline_start_ms, args.baseline_stop_ms)
    if not base.any():
        raise ValueError("Baseline window is empty")
    crossing_idx = np.flatnonzero((soma[:-1] < 0) & (soma[1:] >= 0)) + 1
    crossings = [
        float(t[i-1] - soma[i-1] / (soma[i] - soma[i-1]) * (t[i] - t[i-1]))
        for i in crossing_idx
    ]
    results = []
    for label, trace in zip(labels, ext_uv):
        baseline = trace[base]
        baseline_mean = float(np.mean(baseline))
        noise_sd = float(np.std(baseline, ddof=1))
        noise_rms = float(np.sqrt(np.mean((baseline - baseline_mean)**2)))
        spike_rows = []
        for crossing in crossings:
            mask = window_mask(t, crossing - args.pre_ms, crossing + args.post_ms)
            centered = trace[mask] - baseline_mean
            peak = float(np.max(np.abs(centered)))
            spike_rows.append({
                "soma_crossing_ms": crossing,
                "peak_abs_uV_baseline_corrected": peak,
                "peak_to_peak_uV": float(np.ptp(centered)),
                "snr_peak_over_baseline_rms": peak / noise_rms if noise_rms else None,
            })
        results.append({
            "region": label,
            "baseline_window_ms": [args.baseline_start_ms, args.baseline_stop_ms],
            "baseline_mean_uV": baseline_mean,
            "baseline_sd_uV": noise_sd,
            "baseline_rms_uV": noise_rms,
            "all_recording_peak_abs_uV": float(np.max(np.abs(trace))),
            "spike_windows": spike_rows,
        })
    report = {
        "input": str(args.input),
        "baseline_definition": "100-280 ms, before the first detected somatic upward zero crossing",
        "spike_detection": "somatic membrane voltage upward zero crossing through 0 mV",
        "detected_somatic_crossings_ms": crossings,
        "electrodes": results,
        "limitations": [
            "The baseline window is a quiet pre-spike interval in this single simulation, "
            "not an empirical recording noise floor.",
            "The simulation has no added stochastic measurement noise.",
            "SNR is therefore a model-signal-to-simulation-variability ratio, not hardware SNR."
        ],
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "noise_spike_metrics.json").write_text(json.dumps(report, indent=2) + "\n")
    fig, axes = plt.subplots(len(labels), 1, figsize=(11, 2.2 * len(labels)),
                             sharex=True)
    axes = np.atleast_1d(axes)
    for ax, label, trace in zip(axes, labels, ext_uv):
        ax.plot(t, trace, lw=0.7)
        ax.axvspan(args.baseline_start_ms, args.baseline_stop_ms,
                   color="0.8", alpha=0.4, label="baseline window")
        for crossing in crossings:
            ax.axvline(crossing, color="crimson", lw=0.8, alpha=0.7)
        ax.set_ylabel("µV")
        ax.set_title(label)
        ax.grid(alpha=0.2)
    axes[-1].set_xlabel("Time (ms)")
    fig.suptitle("Baseline noise and detected somatic-spike timing")
    fig.tight_layout()
    fig.savefig(args.output_dir / "noise_spike_overview.png", dpi=220)
    plt.close(fig)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
