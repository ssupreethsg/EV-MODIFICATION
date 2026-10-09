"""Helpers shared by the Python tools so they use exactly the C model's data.

The OCV lookup tables are parsed from spme_v7_apex_solver.c (single source of truth),
and the C OCV proxy executable is located next to this file or in ./build.
"""
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
C_SOURCE = os.path.join(HERE, "spme_v7_apex_solver.c")

# Chen2020 (LG M50) stoichiometry windows used to tabulate the LUTs against cell SOC
ANODE_STOICH_SOC0, ANODE_STOICH_SOC100 = 0.0279, 0.9014
CATHODE_STOICH_SOC0, CATHODE_STOICH_SOC100 = 0.9287, 0.2479
NOMINAL_CAPACITY_AH = 5.0


def _load_lut(name):
    with open(C_SOURCE, encoding="utf-8") as handle:
        source = handle.read()
    match = re.search(name + r"\[101\]\s*=\s*\{([^}]*)\}", source)
    if match is None:
        raise RuntimeError(f"{name} not found in {C_SOURCE}")
    return np.array([float(v.strip().rstrip("f")) for v in match.group(1).split(",")])


LUT_SOC = np.linspace(0.0, 1.0, 101)
LUT_GRAPHITE = _load_lut("LUT_GRAPHITE_OCP_101")
LUT_NMC = _load_lut("LUT_NMC_OCP_101")


def full_cell_ocv_V(soc):
    """Full-cell OCV [V] with the same linear interpolation as the C code."""
    soc = np.clip(soc, 0.0, 1.0)
    return np.interp(soc, LUT_SOC, LUT_NMC) - np.interp(soc, LUT_SOC, LUT_GRAPHITE)


def chen2020_graphite_ocp(sto):
    return (
        1.9793 * np.exp(-39.3631 * sto)
        + 0.2482
        - 0.0909 * np.tanh(29.8538 * (sto - 0.1234))
        - 0.04478 * np.tanh(14.9159 * (sto - 0.2769))
        - 0.0205 * np.tanh(30.4444 * (sto - 0.6103))
    )


def chen2020_nmc_ocp(sto):
    return (
        -0.8090 * sto
        + 4.4875
        - 0.0428 * np.tanh(18.5138 * (sto - 0.5542))
        - 17.7326 * np.tanh(15.7890 * (sto - 0.3117))
        + 17.5842 * np.tanh(15.9308 * (sto - 0.3120))
    )


def find_ocv_proxy():
    exe = "ocv_proxy.exe" if sys.platform.startswith("win") else "ocv_proxy"
    for folder in (HERE, os.path.join(HERE, "build"), os.path.join(HERE, "build", "Release")):
        path = os.path.join(folder, exe)
        if os.path.isfile(path):
            return path
    return None
