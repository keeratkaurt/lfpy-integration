"""Plot only extracellular signals from the validated site-recordings NPZ."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--input', type=Path,
                    default=Path('results_site_recordings/site_recordings.npz'))
parser.add_argument('--output', type=Path, default=None)
args = parser.parse_args()
destination = args.output or args.input.parent / 'extracellular_sites_only.png'
data = np.load(args.input)
time = data['time_ms']
voltage = data['extracellular_mV'] * 1000
labels = ['Perisomatic site', 'Apical dendrite site',
          'Basal dendrite site', 'Axonal site']
colors = ['#d62728', '#1f77b4', '#2ca02c', '#9467bd']
assert voltage.shape == (4, len(time))
assert np.isfinite(voltage).all()
fig, axes = plt.subplots(4, 2, figsize=(12, 9),
                         gridspec_kw={'width_ratios': [2.2, 1]})
for row, (trace, label, color) in enumerate(zip(voltage, labels, colors)):
    for col, (start, stop) in enumerate([(100, 500), (285, 305)]):
        mask = (time >= start) & (time <= stop)
        ax = axes[row, col]
        ax.plot(time[mask], trace[mask], color=color, lw=1)
        ax.set_xlim(start, stop)
        ax.grid(alpha=0.15)
        ax.spines[['top', 'right']].set_visible(False)
        ax.set_ylabel('Extracellular voltage (µV)')
        ax.set_title(label if col == 0 else 'Spike interval close-up',
                     loc='left', fontsize=11)
        if row == 3:
            ax.set_xlabel('Time (ms)')
fig.suptitle('Extracellular recordings at four sites around L5_TTPC1',
             fontsize=15, y=0.985)
fig.text(0.5, 0.947,
         'One InterChaoticB current clamp at the soma • probe colors match the location map',
         ha='center', fontsize=10)
fig.text(0.5, 0.012,
         'LFPy line-source forward model; σ = 0.3 S/m; Δt = 0.1 ms. '
         'Full recording: 0–500 ms; plot starts at 100 ms to omit initialization.',
         ha='center', fontsize=9)
fig.tight_layout(rect=[0, 0.035, 1, 0.93])
destination.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(destination, dpi=220)
plt.close(fig)
