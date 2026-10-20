# Effects of Electrode Placement on Extracellular Neural Recordings

## Question

How does the location of a virtual extracellular electrode relative to a detailed
neuron affect the amplitude and waveform of the recorded extracellular signal?

## Model and simulation

The analysis uses a detailed BBP L5_TTPC1 cortical-cell morphology in NEURON with
LFPy 2.3.5. One current waveform is applied through a single somatic IClamp.
The four recording probes are passive and are placed near the soma, apical
dendrite, basal dendrite, and axon. The probes do not inject current.

The simulation uses a 0.1 ms timestep over 500 ms, a line-source forward model,
and homogeneous conductivity of 0.3 S/m. The saved extracellular traces are in
mV; reported plot values are converted to µV by multiplying by 1,000.

## Verified checks

The saved NPZ contains 5,001 time samples from 0 to 500 ms and four extracellular
traces. All numeric values are finite. The somatic membrane voltage has two
upward zero crossings at 289.04 and 298.23 ms, with a range of -90.47 to
39.41 mV. The recorded somatic clamp current matches the command after the
documented one-sample recording offset.

The saved minimum probe clearances are approximately 20.87 µm for the soma,
21.28 µm for the apical probe, 22.16 µm for the basal probe, and 19.99 µm for
the axonal probe.

## Main results

Peak absolute extracellular voltage over the displayed 100–500 ms interval was:

| Probe region | Peak absolute voltage | Time |
|---|---:|---:|
| Perisomatic | 87.58 µV | 298.3 ms |
| Apical dendrite | 2.77 µV | 375.9 ms |
| Basal dendrite | 12.39 µV | 375.9 ms |
| Axon | 12.55 µV | 299.9 ms |

The perisomatic recording had the largest absolute signal. The apical probe had
the smallest. The basal and axonal probes were intermediate, but their peak
times and polarities differed. The apical and basal whole-recording maxima occur
outside the two detected somatic spike windows, so these values should not be
described as isolated spike amplitudes.

Within 2 ms before to 3 ms after each detected somatic crossing, the largest
raw absolute values were 87.58 µV at the soma, 2.64 µV apically, 12.04 µV
basally, and 12.55 µV axonally. These windows are not baseline-corrected.

## Contribution

My contribution was the LFPy integration and validation workflow: connecting the
NEURON cell model to an extracellular forward model, implementing passive
multi-site probe placement outside modeled compartments, saving reproducible
coordinates and raw arrays, creating the measurement audit, and preparing the
analysis scripts and documentation for electrode-location comparisons.

## Limitations

This is a single-cell deterministic simulation, not an experimental recording.
The probe labels describe neighborhoods; each extracellular trace includes
transmembrane-current contributions from the whole cell. The line-source model
assumes an infinite homogeneous conductor with conductivity 0.3 S/m, and the
stimulus return electrode is not explicitly modeled. The saved clearances are
metadata from the placement routine and should be independently checked after
rebuilding the morphology. The simulation contains no added hardware noise, so
the noise analysis is a model baseline-variability analysis rather than a
measured hardware SNR. A public-data comparison was not included because
available unmatched datasets would not provide a scientifically controlled
comparison.

## Reproduction

From the repository root, install the pinned environment:

```bash
python -m pip install -r requirements-lfpy-sites.txt
```

Run or inspect the validated site recording:

```bash
python lfpy_site_recordings.py MODEL_DIR \
  --stim-file stims/5k50kInterChaoticB.csv \
  --output-dir results_site_recordings
```

Run the saved-data audit:

```bash
python verify_site_recordings.py
```

Run the additional analyses:

```bash
python analysis/electrode_position_heatmap.py
python analysis/noise_spike_analysis.py
python analysis/robustness_comparison.py \
  --run 20=results_site_recordings/site_recordings.npz \
  --run 50=results_offset_50/site_recordings.npz
```

The robustness script requires at least two independently generated electrode
offset runs. It should not be interpreted as complete until those additional
runs exist.
