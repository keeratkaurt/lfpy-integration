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

## Proof of Concept

The current proof-of-concept implementation:

1. Loads neuron information from `cells.json`
2. Loads the existing BBP HOC templates
3. Creates the neuron model using LFPy
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
- `lfpy_poc_output.npz`

## Integration Plan

Once the proof of concept is validated on NERSC, the extracellular-voltage calculation will be integrated into the existing `run.py` data-generation workflow.

The existing HDF5 output contains intracellular voltage under `volts`.

The planned extension is to also save extracellular voltage under a new output such as:

`extracellular`

This will allow the existing parameter-sampling pipeline to generate both intracellular and simulated extracellular signals.

## Current Status

Implementation prepared.

Pending:

- NERSC account approval
- validation with the lab's BBP morphology and HOC templates
- inspection of extracellular traces
- final integration into the large-scale HDF5 generation pipeline
