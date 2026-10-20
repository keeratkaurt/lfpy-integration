# LFPy Extracellular Signal Integration

## Goal

Extend the existing DL4neurons2 NEURON simulation pipeline to generate simulated extracellular electrode recordings.

## Latest validation: extracellular recordings near distinct cell regions

This demonstration applies **one InterChaoticB current clamp at the soma** of
the supplied BBP L5_TTPC1 model and calculates extracellular voltage at four
passive point probes near the soma, apical dendrite, basal dendrite, and axon.
The probes do not inject current. Each probe records the field from all
transmembrane currents; its region label denotes its location, not an isolated
contribution from that region.

Install the tested environment, compile the model mechanisms as described
below, and run from the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements-lfpy-sites.txt
python lfpy_site_recordings.py /path/to/L5_TTPC1_cADpyr232_1 --stim-file stims/5k50kInterChaoticB.csv --dt 0.1 --output-dir results_site_recordings
python plot_extracellular_sites.py --input results_site_recordings/site_recordings.npz
```

Replace the model path with the actual unpacked model directory. The two new
scripts use `build_cell` from `lfpy_poc.py`; keep all three scripts together at
the repository root. NEURON 8.2.7 is required for the legacy BBP mechanisms.

The 0.1 ms interval, multiplier 1, and offset 0 match the settings in
`BBP_sbatch_def.sh`, `BBP_sbatch_submit.sh`, and `newBBP_sbatch_submit.sh` for
InterChaoticB. The 5,000 CSV samples span 500 ms. Do not infer timing from the
general `run.py` default (0.025 ms).

### Figures and numerical output

- [Probe locations](results_site_recordings/site_electrode_placement.png)
- [Extracellular voltage only, with spike-interval close-ups](results_site_recordings/extracellular_sites_only.png)
- [Supporting stimulus, soma voltage, and extracellular panels](results_site_recordings/site_extracellular_recordings.png)
- [Raw arrays and probe coordinates](results_site_recordings/site_recordings.npz)

Validated on 3,323 compartments with two somatic upward zero-crossings.
Peak absolute extracellular voltage: perisomatic 87.58 µV, apical 2.77 µV,
basal 12.39 µV, axonal 12.55 µV. Traces are stored in mV and plotted in µV.
The conductivity is 0.3 S/m and the forward model uses line sources in an
infinite homogeneous conductor. The injected-current return electrode is not
explicitly modeled, and stimulus-correlated field components are retained.

The x-y location map is a projection; probes are displaced in z and can
overlap drawn branches in that projection. Minimum clearance from modeled
compartment surfaces is approximately 20–22 µm, estimated from line segments
and compartment radii. Labels identify probe neighborhoods; the dense
perisomatic region includes nearby dendrites. The NPZ records both the target
and nearest compartment indices and exact probe coordinates.

The extracellular-only plot shows 100–500 ms to omit initialization and
285–305 ms for the close-ups. The unfiltered 0–500 ms data remain in the NPZ.
In this NEURON/LFPy combination the recorded clamp current is one sample
behind the command vector; both arrays are retained. No stimulus has been
applied to dendrites, axon, or virtual recording probes.

This completes a single-cell demonstration. HDF5 integration and NERSC
scaling remain future work.

## Existing Pipeline

The current pipeline uses NEURON to simulate detailed neuron models and generate intracellular membrane-voltage responses while varying biophysical parameters.

## LFPy Extension

The new workflow uses LFPy to calculate extracellular voltage from the same neuron simulations.

The intended flow is:

NEURON simulation  
→ transmembrane currents  
→ LFPy extracellular forward model  
→ virtual electrode recordings

## Local BBP Proof of Concept

This branch is tested with NEURON 8.2.7 and LFPy 2.3.5. The legacy BBP
probabilistic synapse mechanisms do not compile unchanged with NEURON 9.

```bash
python -m venv .venv
.venv/bin/pip install -r requirements-lfpy.txt

cd /path/to/L5_TTPC1_cADpyr232_1
CC=gcc CXX=g++ /path/to/.venv/bin/nrnivmodl mechanisms
cd /path/to/lfpy-integration

.venv/bin/python lfpy_poc.py /path/to/L5_TTPC1_cADpyr232_1 \
  --output-dir lfpy_results
```

To reproduce the test with the DL4neurons2 `InterChaoticB` stimulus:

```bash
.venv/bin/python lfpy_poc.py /path/to/L5_TTPC1_cADpyr232_1 \
  --stim-file stims/5k50kInterChaoticB.csv \
  --output-dir lfpy_results_interchaoticb
```

### Outside-cell electrode validation

```bash
.venv/bin/python lfpy_poc.py /path/to/L5_TTPC1_cADpyr232_1 \
  --stim-file stims/5k50kInterChaoticB.csv \
  --outside-cell --outside-margin 20 \
  --output-dir lfpy_results_outside
```

This run uses the BBP `L5_TTPC1_cADpyr232_1` model (morphology
`dend-C060114A2_axon-C060114A5.asc`). The placement figure uses the
anatomical x-y plane. Four point electrodes are 20 µm beyond the x-y
bounding box, at soma z: left, right, above the apical tuft, and below
the basal/axonal extent. Coordinates and labels are saved in the NPZ.

The 5,000 CSV samples are applied at **0.1 ms/sample** for this
500 ms demonstration; this is an explicit choice, not the default
of the main `run.py` pipeline (0.025 ms). Specify `--dt` to change it.
The clamp-current recorder is offset by one sample relative to the
command vector in this LFPy/NEURON combination; the saved NPZ includes both.
The extracellular plots show the recorded response, including an
initialization transient at time zero.

Validated with NEURON 8.2.7 and LFPy 2.3.5: 3,323 compartments,
two somatic zero-crossing action potentials, soma voltage -90.47 to
+39.41 mV. Peak absolute extracellular signal (µV): left 2.75,
right 2.71, above apical 1.65, below basal 2.00. The run produces
`outside_cell_electrode_placement.png` and
`outside_cell_voltage_traces.png` in addition to the common outputs.

The CSV values are interpreted as nA and played into a somatic `IClamp` at the
simulation timestep. Optional `--stim-multiplier` and `--stim-dc-offset` flags
support the scaling used by the larger data-generation pipeline; both default
to an unchanged, reproducible waveform (1.0 and 0.0 nA).

Run the Python command from a directory that does not already contain a
different `x86_64/libnrnmech.so`, because NEURON automatically loads mechanisms
from the working directory.

## Proof of Concept

The current proof-of-concept implementation:

1. Accepts an unpacked BBP model directory
2. Detects the HOC template and ASC morphology
3. Loads the compiled BBP mechanisms and creates the neuron model using LFPy
4. Applies a current stimulus
5. Records intracellular membrane voltage
6. Records transmembrane currents
7. Places virtual extracellular electrodes at multiple distances
8. Calculates extracellular voltage
9. Saves intracellular and extracellular plots
10. Saves the raw output as an NPZ file

## Files Added

### `lfpy_poc.py`

Runs the lab-specific LFPy proof of concept.

### `extracellular.py`

Contains the reusable extracellular-voltage calculation.

## Expected Outputs

When run successfully in the NERSC/lab environment, the proof of concept should generate:

- `intracellular_voltage.png`
- `extracellular_voltage.png`
- `stimulus_current.png`
- `interchaoticb_summary.png`
- `lfpy_poc_output.npz`

## Integration Plan

Once the proof of concept is validated on NERSC, the extracellular-voltage calculation will be integrated into the existing `run.py` data-generation workflow.

The existing HDF5 output contains intracellular voltage under `volts`.

The planned extension is to also save extracellular voltage under a new output such as:

`extracellular`

This will allow the existing parameter-sampling pipeline to generate both intracellular and simulated extracellular signals.

## Current Status

Validated locally with `L5_TTPC1_cADpyr232_1`:

- 3,323 compartments
- four virtual electrodes at 20, 50, 100, and 200 micrometers from the soma
- intracellular and extracellular plots generated successfully
- raw time series saved to NPZ

Pending:

- integration into the large-scale HDF5 generation pipeline
- NERSC/Perlmutter validation and scaling


## Talent Showcase analysis

The presentation work is maintained on the `talent-showcase` branch. The short
technical report is in [SHOWCASE_REPORT.md](SHOWCASE_REPORT.md). The saved-data
verification is in [verify_site_recordings.py](verify_site_recordings.py), with
its checked values in
[results_site_recordings/verification_metrics.json](results_site_recordings/verification_metrics.json).

Additional analyses are in `analysis/`:

- `electrode_position_heatmap.py` creates a labeled electrode-position amplitude map.
- `noise_spike_analysis.py` detects somatic upward zero crossings and estimates baseline variability and model SNR.
- `robustness_comparison.py` compares amplitudes across independently generated electrode-offset runs.

The heatmap and noise analysis operate on the saved NPZ without rerunning NEURON.
The robustness analysis requires at least two additional offset runs. The public-data
comparison was assessed as infeasible for the deadline because an unmatched
experimental dataset would not share the simulation's geometry, units, referencing,
and conductivity assumptions; see [PUBLIC_DATA_COMPARISON.md](PUBLIC_DATA_COMPARISON.md).

For the presentation, report whole-recording peaks as "peak absolute extracellular
voltage over 100–500 ms." The apical and basal whole-recording maxima occur around
375.9 ms, outside the two detected somatic spike windows.


The pushed 50 µm run is a placement-sensitivity comparison against the original
20 µm run. Because the current placement routine selects different nearby
compartments at different offsets, this comparison is not a controlled
distance-only experiment. The report records this limitation explicitly.
