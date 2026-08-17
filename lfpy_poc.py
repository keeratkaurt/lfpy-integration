import os
import json
import numpy as np
import matplotlib.pyplot as plt

from neuron import h
import LFPy


# ---------------------------------------------------------
# CONFIG
# ---------------------------------------------------------

M_TYPE = "L1_DAC"
E_TYPE = "bNAC"
CELL_I = 0

DT = 0.025
TSTOP = 500.0
V_INIT = -68.0

SIGMA = 0.3  # extracellular conductivity S/m


# Same template location currently used by DL4neurons2
TEMPLATES_DIR = "/global/cfs/cdirs/m3513/M1_Hoc_template/HocTemplate"


# ---------------------------------------------------------
# LOAD CELL INFORMATION
# ---------------------------------------------------------

with open("cells.json", "r") as f:
    cells = json.load(f)

cell_info = cells[M_TYPE][E_TYPE][CELL_I]

model_directory = cell_info["model_directory"]
template_name = cell_info["model_template"].split(":", 1)[-1]
morphology_name = cell_info["morphology"]

cell_dir = os.path.join(
    TEMPLATES_DIR,
    model_directory,
    model_directory
)

constants_file = os.path.join(cell_dir, "constants.hoc")
morphology_hoc = os.path.join(cell_dir, "morphology.hoc")
biophysics_file = os.path.join(cell_dir, "biophysics.hoc")
synapses_file = os.path.join(cell_dir, "synapses", "synapses.hoc")
template_file = os.path.join(cell_dir, "template.hoc")

morphology_file = os.path.join(
    cell_dir,
    "morphology",
    morphology_name
)


# ---------------------------------------------------------
# LOAD BBP HOC DEFINITIONS
# ---------------------------------------------------------

h.load_file("stdrun.hoc")
h.load_file("import3d.hoc")

h.load_file(constants_file)
h.load_file(morphology_hoc)
h.load_file(biophysics_file)
h.load_file(synapses_file)
h.load_file(template_file)


# ---------------------------------------------------------
# BUILD THE SAME TYPE OF CELL THROUGH LFPy
# ---------------------------------------------------------

cwd = os.getcwd()
os.chdir(cell_dir)

try:
    cell = LFPy.TemplateCell(
        morphology=morphology_file,
        templatefile=template_file,
        templatename=template_name,
        templateargs=0,        # NO_SYNAPSES in existing code
        dt=DT,
        tstop=TSTOP,
        v_init=V_INIT,
        passive=False,
        delete_sections=False,
        verbose=True
    )
finally:
    os.chdir(cwd)


# ---------------------------------------------------------
# ADD CURRENT STIMULUS
# ---------------------------------------------------------

stim = LFPy.StimIntElectrode(
    cell,
    idx=0,
    pptype="IClamp",
    amp=1.0,
    dur=300.0,
    delay=100.0
)


# ---------------------------------------------------------
# RUN SIMULATION
# ---------------------------------------------------------

cell.simulate(
    rec_vmem=True,
    rec_imem=True
)


# ---------------------------------------------------------
# VIRTUAL ELECTRODES
# ---------------------------------------------------------

electrode_x = np.array([20., 50., 100., 200.])
electrode_y = np.zeros(4)
electrode_z = np.zeros(4)

electrode = LFPy.RecExtElectrode(
    cell=cell,
    sigma=SIGMA,
    x=electrode_x,
    y=electrode_y,
    z=electrode_z,
    method="linesource"
)

M = electrode.get_transformation_matrix()

extracellular = M @ cell.imem


# ---------------------------------------------------------
# PLOT INTRACELLULAR SIGNAL
# ---------------------------------------------------------

plt.figure(figsize=(10, 4))

plt.plot(cell.tvec, cell.vmem[0])

plt.xlabel("Time (ms)")
plt.ylabel("Membrane voltage (mV)")
plt.title("Intracellular membrane voltage")

plt.tight_layout()
plt.savefig("intracellular_voltage.png", dpi=200)
plt.close()


# ---------------------------------------------------------
# PLOT EXTRACELLULAR SIGNALS
# ---------------------------------------------------------

plt.figure(figsize=(10, 5))

for i, distance in enumerate(electrode_x):
    plt.plot(
        cell.tvec,
        extracellular[i],
        label=f"{distance:.0f} µm"
    )

plt.xlabel("Time (ms)")
plt.ylabel("Extracellular voltage (mV)")
plt.title("Simulated extracellular electrode recordings")

plt.legend()
plt.tight_layout()

plt.savefig("extracellular_voltage.png", dpi=200)
plt.close()


# ---------------------------------------------------------
# SAVE RAW DATA
# ---------------------------------------------------------

np.savez(
    "lfpy_poc_output.npz",
    time=cell.tvec,
    intracellular=cell.vmem[0],
    extracellular=extracellular,
    electrode_x=electrode_x,
    electrode_y=electrode_y,
    electrode_z=electrode_z
)


print()
print("SUCCESS")
print("Number of compartments:", cell.totnsegs)
print("Transmembrane current shape:", cell.imem.shape)
print("Extracellular voltage shape:", extracellular.shape)
print()
print("Saved:")
print("  intracellular_voltage.png")
print("  extracellular_voltage.png")
print("  lfpy_poc_output.npz")
