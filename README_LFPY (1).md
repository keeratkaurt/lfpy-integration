# LFPy Extracellular Signal Integration

## Goal

Extend the existing DL4neurons2 NEURON simulation pipeline to generate simulated extracellular electrode recordings.

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
