"""Helpers shared by the Python tools so they use exactly the C model's data.

Constants and OCV lookup tables are parsed from src/spme_v7_apex_solver.{h,c}
(single source of truth). The ev_cli executable is located through the EV_CLI
environment variable or in the usual CMake build folders.

Sample convention shared by all tools (and by ev_step): current_A[i] is the current that
flowed over the interval (time_s[i-1], time_s[i]].
"""
import functools
import os
import re
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
C_HEADER = os.path.join(ROOT, "src", "spme_v7_apex_solver.h")
C_SOURCE = os.path.join(ROOT, "src", "spme_v7_apex_solver.c")
LUT_POINTS = 201


def _read(path):
    with open(path, encoding="utf-8") as handle:
        return handle.read()


def header_constant(name):
    """Numeric value of a '#define NAME <number>' in the C header."""
    match = re.search(r"^#define\s+" + name + r"\s+([-+0-9.eE]+)[fFuU]?\b", _read(C_HEADER), re.MULTILINE)
    if match is None:
        raise KeyError(f"{name} not found in {C_HEADER}")
    return float(match.group(1))


ANODE_STOICH_SOC0 = header_constant("EV_ANODE_STOICH_AT_SOC0")
ANODE_STOICH_SOC100 = header_constant("EV_ANODE_STOICH_AT_SOC100")
CATHODE_STOICH_SOC0 = header_constant("EV_CATHODE_STOICH_AT_SOC0")
CATHODE_STOICH_SOC100 = header_constant("EV_CATHODE_STOICH_AT_SOC100")
MAX_CELL_VOLTAGE_V = header_constant("EV_MAX_CELL_VOLTAGE_V")


def chen2020_graphite_ocp(sto):
    """PyBaMM graphite_LGM50_ocp_Chen2020."""
    return (
        1.9793 * np.exp(-39.3631 * sto)
        + 0.2482
        - 0.0909 * np.tanh(29.8538 * (sto - 0.1234))
        - 0.04478 * np.tanh(14.9159 * (sto - 0.2769))
        - 0.0205 * np.tanh(30.4444 * (sto - 0.6103))
    )


def chen2020_nmc_ocp(sto):
    """PyBaMM nmc_LGM50_ocp_Chen2020."""
    return (
        -0.8090 * sto
        + 4.4875
        - 0.0428 * np.tanh(18.5138 * (sto - 0.5542))
        - 17.7326 * np.tanh(15.7890 * (sto - 0.3117))
        + 17.5842 * np.tanh(15.9308 * (sto - 0.3120))
    )


def anode_stoich(soc):
    return ANODE_STOICH_SOC0 + np.asarray(soc) * (ANODE_STOICH_SOC100 - ANODE_STOICH_SOC0)


def cathode_stoich(soc):
    return CATHODE_STOICH_SOC0 + np.asarray(soc) * (CATHODE_STOICH_SOC100 - CATHODE_STOICH_SOC0)


def exact_full_cell_ocv_V(soc):
    """Chen2020 OCV without tabulation."""
    return chen2020_nmc_ocp(cathode_stoich(soc)) - chen2020_graphite_ocp(anode_stoich(soc))


def _load_lut(name):
    match = re.search(name + r"\[EV_LUT_POINTS\]\s*=\s*\{([^}]*)\}", _read(C_SOURCE))
    if match is None:
        raise RuntimeError(f"{name} not found in {C_SOURCE}")
    values = np.array([float(v.strip().rstrip("f")) for v in match.group(1).split(",") if v.strip()])
    if values.size != LUT_POINTS:
        raise RuntimeError(f"{name} has {values.size} entries, expected {LUT_POINTS}")
    return values


LUT_SOC = np.linspace(0.0, 1.0, LUT_POINTS)


@functools.lru_cache(maxsize=None)
def luts():
    """(graphite, nmc) electrode OCP tables [V] on LUT_SOC, as stored in the C source."""
    return _load_lut("LUT_GRAPHITE_OCP"), _load_lut("LUT_NMC_OCP")


def full_cell_ocv_V(soc):
    """Full-cell OCV [V] with the same linear interpolation as the C code."""
    graphite, nmc = luts()
    soc = np.clip(soc, 0.0, 1.0)
    return np.interp(soc, LUT_SOC, nmc) - np.interp(soc, LUT_SOC, graphite)


def find_ev_cli(explicit=None):
    """Path to ev_cli: the explicit argument, else $EV_CLI, else the CMake build folders.

    An explicit path or $EV_CLI that does not exist raises FileNotFoundError instead of silently
    falling back to a possibly stale build."""
    for source, path in (("--cli", explicit), ("EV_CLI", os.environ.get("EV_CLI"))):
        if path:
            if not os.path.isfile(path):
                raise FileNotFoundError(f"{source} points to {path}, which does not exist")
            return path
    exe = "ev_cli.exe" if sys.platform.startswith("win") else "ev_cli"
    for folder in ("build", os.path.join("build", "Release"), os.path.join("build", "Debug")):
        path = os.path.join(ROOT, folder, exe)
        if os.path.isfile(path):
            return path
    return None


def run_cli(cli, args, lines):
    """Runs ev_cli with one input line per entry and returns its whitespace-split output rows."""
    result = subprocess.run([cli, *args], input="\n".join(lines) + "\n", text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"ev_cli {' '.join(args)} failed ({result.returncode}): {result.stderr.strip()}")
    rows = [row.split() for row in result.stdout.splitlines() if row.strip()]
    if len(rows) != len(lines):
        raise RuntimeError(f"ev_cli returned {len(rows)} rows for {len(lines)} inputs")
    return rows
