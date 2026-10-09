"""Compares the C model's OCV + ohmic surrogate against a PyBaMM SPMe 1C charge
(Chen2020 / LG M50, 10% SOC to 4.2 V). Build ocv_proxy first (CMake or compile.bat)."""
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pybamm  # noqa: E402

import ev_model as m  # noqa: E402

INITIAL_SOC = 0.1


def main():
    proxy = m.find_ocv_proxy()
    if proxy is None:
        print("ocv_proxy executable not found: build it with CMake or compile.bat first.")
        return 1

    parameter_values = pybamm.ParameterValues("Chen2020")
    capacity_Ah = parameter_values["Nominal cell capacity [A.h]"]
    experiment = pybamm.Experiment(["Charge at 1 C until 4.2 V"])
    sim = pybamm.Simulation(pybamm.lithium_ion.SPMe(), experiment=experiment, parameter_values=parameter_values)
    sol = sim.solve(initial_soc=INITIAL_SOC)

    time_s = sol["Time [s]"].entries
    pybamm_voltage = sol["Voltage [V]"].entries
    # PyBaMM: positive current = discharge. Convert to this project's convention
    # (positive = charging).
    charge_current_A = -sol["Current [A]"].entries
    if np.mean(charge_current_A) <= 0.0:
        raise RuntimeError("Expected a charging current from the PyBaMM experiment")

    # Coulomb counting with trapezoidal integration.
    charged_Ah = np.concatenate(([0.0], np.cumsum(0.5 * (charge_current_A[1:] + charge_current_A[:-1]) * np.diff(time_s)))) / 3600.0
    soc = np.clip(INITIAL_SOC + charged_Ah / capacity_Ah, 0.0, 1.0)

    result = subprocess.run([proxy], input="\n".join(f"{s:.9f}" for s in soc) + "\n",
                            text=True, capture_output=True, check=True)
    c_ocv = np.array(result.stdout.split(), dtype=float)
    if c_ocv.size != soc.size:
        raise RuntimeError(f"ocv_proxy returned {c_ocv.size} values for {soc.size} inputs")

    # Least-squares effective resistance of the OCV + I*R surrogate. It lumps every
    # SPMe overpotential (ohmic, kinetic, solid and electrolyte concentration).
    effective_resistance_ohm = float(np.sum(charge_current_A * (pybamm_voltage - c_ocv)) / np.sum(charge_current_A ** 2))
    c_voltage = c_ocv + charge_current_A * effective_resistance_ohm
    rms_error = float(np.sqrt(np.mean((c_voltage - pybamm_voltage) ** 2)))

    print(f"Final PyBaMM voltage: {pybamm_voltage[-1]:.4f} V, final SOC: {soc[-1]:.3f}")
    print(f"Fitted effective resistance: {effective_resistance_ohm * 1000.0:.1f} mOhm")
    print(f"RMS Error against PyBaMM: {rms_error:.4f} V")

    plt.plot(time_s, pybamm_voltage, label="PyBaMM SPMe (Chen2020)")
    plt.plot(time_s, c_voltage, "--", label=f"C model OCV + I*R (R = {effective_resistance_ohm * 1000.0:.1f} mOhm)")
    plt.xlabel("Time [s]")
    plt.ylabel("Voltage [V]")
    plt.legend()
    plt.title(f"1C charge from {INITIAL_SOC:.0%} SOC (RMS error: {rms_error:.4f} V)")
    out_path = os.path.join(m.HERE, "pybamm_charge_benchmark.png")
    plt.savefig(out_path, dpi=120)
    print(f"Generated {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
