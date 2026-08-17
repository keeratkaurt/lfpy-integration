import numpy as np
import LFPy


def compute_extracellular(
    cell,
    electrode_x,
    electrode_y,
    electrode_z,
    sigma=0.3,
    method="linesource",
):
    """
    Compute extracellular voltage from an LFPy cell simulation.
    """

    electrode = LFPy.RecExtElectrode(
        cell=cell,
        sigma=sigma,
        x=np.asarray(electrode_x),
        y=np.asarray(electrode_y),
        z=np.asarray(electrode_z),
        method=method,
    )

    transformation = electrode.get_transformation_matrix()

    return transformation @ cell.imem
