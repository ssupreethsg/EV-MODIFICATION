"""Validates the C model against PyBaMM (Chen2020 / LG M50).

Every sample convention follows ev_step(): current_A[i] is the current over (t[i-1], t[i]].

1. OCV: the C table at PyBaMM's SOC vs PyBaMM's bulk open-circuit voltage (SPMe).
2. 1-RC fit: R0, R1 and tau are fitted on a 0.5C pulse profile, then checked out of sample on a
   1C pulse profile (four cycles of up to 10 min charge, stopping at 4.2 V, + 10 min rest from
   10 % SOC, then 1C to 4.2 V).
3. EKF: the C estimator (ev_cli ekf, through ev_step) starts 0.20 SOC high on the 1C profile, with
   the fitted parameters and with R0 +/-20 % and capacity -5 % errors. From 15 min on (the end of
   the first rest) the SOC error must stay below 3 % and inside the filter's own 3-sigma band (at
   most 1 % of samples outside; the controller's SOC limit relies on that band), and no sample may
   raise a fault (a healthy cell must never trip the voltage-plausibility check).
4. Plating: the complete controller (ev_cli charge: ev_step with EKF, gradient history, voltage
   tracking and over-current supervision, a 15 A (2.9C) request and a 30 s horizon) charges a
   PyBaMM DFN with Arrhenius electrolyte and diffusion terms in a closed loop, isothermal at 2, 10
   and 25 degC, from 5 % SOC, with the DFN terminal voltage as the measurement. The DFN anode
   potential vs Li/Li+ (phi_s - phi_e) must stay >= 0 V everywhere in the electrode, no fault may
   occur, and the limiter must still charge (minimum final SOC per temperature).

Usage: python tools/run_pybamm_sim.py [--cli PATH] [--out PNG] [--skip-dfn]
"""
import argparse
import os
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pybamm  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

import ev_model as m  # noqa: E402

FARADAY = 96485.33212
INITIAL_SOC = 0.1
TEMPERATURE_C = 25.0
EKF_INITIAL_SOC_ERROR = 0.2
EKF_SETTLE_S = 900.0   # end of the first rest: the filter learns SOC mainly at low current
OCV_TOLERANCE_V = 0.001
RC_FIT_RMS_TOLERANCE_V = 0.025
EKF_SOC_TOLERANCE = 0.03
EKF_MAX_FRACTION_OUTSIDE_3SIGMA = 0.01
DFN_TEMPERATURES_C = (2.0, 10.0, 25.0)
DFN_MIN_FINAL_SOC = {2.0: 0.65, 10.0: 0.80, 25.0: 0.85}   # the limiter must keep charging, not just stay safe
DFN_START_SOC = 0.05
DFN_STEP_S = 10.0
DFN_MAX_TIME_S = 2.0 * 3600.0
DFN_REQUEST_A = 15.0
DFN_HORIZON_S = 30.0
VALIDATION_PROFILE = ["Charge at 1C for 10 minutes or until 4.2 V", "Rest for 10 minutes"] * 4 + ["Charge at 1C until 4.2 V"]
PLATING_ONSET_V = 0.0


def chen2020(temperature_C=TEMPERATURE_C, arrhenius=False):
    parameter_values = pybamm.ParameterValues("Chen2020")
    x0, x100, y100, y0 = pybamm.lithium_ion.get_min_max_stoichiometries(parameter_values)
    for name, value, header in (("x0", x0, m.ANODE_STOICH_SOC0), ("x100", x100, m.ANODE_STOICH_SOC100),
                                ("y0", y0, m.CATHODE_STOICH_SOC0), ("y100", y100, m.CATHODE_STOICH_SOC100)):
        if abs(value - header) > 1e-5:
            raise RuntimeError(f"PyBaMM stoichiometry {name} = {value:.6f} differs from the C header ({header:.6f})")
    if arrhenius:
        def factor(activation_energy, temperature):
            return pybamm.exp(activation_energy / m.header_constant("EV_GAS_J_PER_MOL_K") * (1.0 / 298.15 - 1.0 / temperature))
        kappa = parameter_values["Electrolyte conductivity [S.m-1]"]
        diffusivity = parameter_values["Electrolyte diffusivity [m2.s-1]"]
        solid_diffusivity = parameter_values["Negative particle diffusivity [m2.s-1]"]
        parameter_values.update({
            "Electrolyte conductivity [S.m-1]": lambda c, T: kappa(c, T) * factor(m.header_constant("EV_ELECTROLYTE_ACTIVATION_ENERGY_J_MOL"), T),
            "Electrolyte diffusivity [m2.s-1]": lambda c, T: diffusivity(c, T) * factor(m.header_constant("EV_ELECTROLYTE_DIFFUSIVITY_ACTIVATION_J_MOL"), T),
            "Negative particle diffusivity [m2.s-1]": lambda sto, T: solid_diffusivity * factor(m.header_constant("EV_DIFFUSION_ACTIVATION_ENERGY_J_MOL"), T),
        })
    parameter_values.update({"Ambient temperature [K]": 273.15 + temperature_C,
                             "Initial temperature [K]": 273.15 + temperature_C})
    capacity_Ah = (parameter_values["Negative electrode active material volume fraction"]
                   * parameter_values["Negative electrode thickness [m]"]
                   * parameter_values["Electrode height [m]"] * parameter_values["Electrode width [m]"]
                   * parameter_values["Maximum concentration in negative electrode [mol.m-3]"]
                   * FARADAY / 3600.0 * (x100 - x0))
    return parameter_values, capacity_Ah


def simulate_spme(steps):
    parameter_values, capacity_Ah = chen2020()
    experiment = pybamm.Experiment(steps, period="10 seconds")
    sol = pybamm.Simulation(pybamm.lithium_ion.SPMe(), experiment=experiment,
                            parameter_values=parameter_values).solve(initial_soc=INITIAL_SOC)
    time_s = sol["Time [s]"].entries
    keep = np.concatenate(([True], np.diff(time_s) > 1e-6))   # keep the ending step's boundary sample
    data = {
        "time_s": time_s[keep],
        "current_A": -sol["Current [A]"].entries[keep],   # PyBaMM: positive = discharge
        "voltage_V": sol["Voltage [V]"].entries[keep],
        "bulk_ocv_V": sol["Bulk open-circuit voltage [V]"].entries[keep],
        "soc": INITIAL_SOC - sol["Discharge capacity [A.h]"].entries[keep] / capacity_Ah,
    }
    return data, capacity_Ah


def rc_model_voltage(params, time_s, current_A, ocv_V):
    r0, r1, tau = params
    v_rc = np.zeros_like(time_s)
    for i in range(1, time_s.size):
        a = np.exp(-(time_s[i] - time_s[i - 1]) / tau)
        v_rc[i] = a * v_rc[i - 1] + r1 * (1.0 - a) * current_A[i]
    return ocv_V + v_rc + r0 * current_A


def fit_rc(data, ocv_V):
    result = least_squares(lambda p: rc_model_voltage(p, data["time_s"], data["current_A"], ocv_V) - data["voltage_V"],
                           x0=[0.02, 0.02, 100.0], bounds=([0.0, 0.0, 1.0], [0.5, 0.5, 5000.0]), x_scale=[0.01, 0.01, 100.0])
    return result.x


def c_ocv(cli, soc):
    return np.array(m.run_cli(cli, ["ocv"], [f"{s:.9f}" for s in np.clip(soc, 0.0, 1.0)]), dtype=float)[:, 0]


def run_ekf(cli, params, capacity_Ah, data):
    r0, r1, tau = params
    t, i_a, v = data["time_s"], data["current_A"], data["voltage_V"]
    lines = [f"{t[k] - t[k - 1]:.6f} {i_a[k]:.6f} {v[k]:.6f} {TEMPERATURE_C:.3f}" for k in range(1, t.size)]
    args = ["ekf", f"{capacity_Ah:.6f}", f"{r0:.6f}", f"{r1:.6f}", f"{tau:.6f}",
            f"{INITIAL_SOC + EKF_INITIAL_SOC_ERROR:.6f}", f"{TEMPERATURE_C:.3f}"]
    rows = np.array(m.run_cli(cli, args, lines), dtype=float)
    soc = np.concatenate(([INITIAL_SOC + EKF_INITIAL_SOC_ERROR], rows[:, 0]))
    sigma = np.concatenate(([0.05], rows[:, 1]))
    fault_samples = int(np.count_nonzero(rows[:, 4]))
    return soc, sigma, fault_samples


def dfn_closed_loop(cli, temperature_C, params):
    """Complete C controller charging a DFN: returns min anode potential vs Li/Li+ [V], the
    number of cycles with fault flags, and the trace (time, SOC, current, min potential)."""
    parameter_values, capacity_Ah = chen2020(temperature_C, arrhenius=True)
    parameter_values.update({"Current function [A]": "[input]"})
    parameter_values.set_initial_state(DFN_START_SOC)
    sim = pybamm.Simulation(pybamm.lithium_ion.DFN(), parameter_values=parameter_values)
    r0, r1, tau = params
    process = subprocess.Popen([cli, "charge", f"{capacity_Ah:.6f}", f"{r0:.6f}", f"{r1:.6f}", f"{tau:.6f}",
                                f"{DFN_START_SOC:.6f}", f"{temperature_C:.3f}"],
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    soc, elapsed, min_potential, faults = DFN_START_SOC, 0.0, np.inf, 0
    trace = []
    current_A = 0.0   # current over the interval that just ended (a rested start)
    try:
        sol = sim.step(DFN_STEP_S, inputs={"Current function [A]": 0.0}, save=False)
        while elapsed < DFN_MAX_TIME_S and soc < m.header_constant("EV_MAX_SAFE_SOC"):
            voltage_V = float(sol["Voltage [V]"].entries[-1])
            process.stdin.write(f"{DFN_STEP_S} {current_A:.6f} {voltage_V:.6f} {temperature_C:.3f} {DFN_REQUEST_A} {DFN_HORIZON_S}\n")
            process.stdin.flush()
            line = process.stdout.readline()
            if not line:
                raise RuntimeError("ev_cli charge stopped responding")
            fields = line.split()
            current_A, flags = float(fields[0]), int(fields[2])
            faults += int(flags != 0)
            if current_A < 0.01 and trace:
                break   # end of charge (or a fault, which is counted)
            sol = sim.step(DFN_STEP_S, inputs={"Current function [A]": -current_A}, save=False)
            potential = sol["Negative electrode surface potential difference [V]"].entries[:, -1]
            min_potential = min(min_potential, float(np.min(potential)))
            soc = DFN_START_SOC - float(sol["Discharge capacity [A.h]"].entries[-1]) / capacity_Ah
            elapsed += DFN_STEP_S
            trace.append((elapsed, soc, current_A, float(np.min(potential))))
    finally:
        process.stdin.close()
        process.wait()
    return min_potential, faults, np.array(trace)


def plot(out_path, val, rc_voltage_V, rc_rms, ekf_soc, ekf_sigma, dfn_traces):
    rows = 3 if dfn_traces else 2
    fig, axes = plt.subplots(rows, 1, figsize=(8, 3.2 * rows), constrained_layout=True)
    minutes = val["time_s"] / 60.0
    ax = axes[0]
    ax.plot(minutes, val["voltage_V"], color="#1f4e79", label="PyBaMM SPMe")
    ax.plot(minutes, rc_voltage_V, "--", color="#c55a11", label=f"C 1-RC model, fitted on 0.5C data (RMS {rc_rms * 1000:.1f} mV)")
    ax.set_ylabel("Terminal voltage [V]")
    ax.set_xlabel("Time [min]")
    ax.legend(loc="lower right", frameon=False)
    ax.grid(alpha=0.3)
    ax = axes[1]
    ax.plot(minutes, val["soc"], color="#1f4e79", label="PyBaMM SOC")
    ax.plot(minutes, ekf_soc, "--", color="#c55a11", label=f"C EKF started {EKF_INITIAL_SOC_ERROR:+.2f} SOC off")
    ax.fill_between(minutes, ekf_soc - 3 * ekf_sigma, ekf_soc + 3 * ekf_sigma, color="#c55a11", alpha=0.15, label="EKF 3-sigma")
    ax.set_ylabel("SOC [-]")
    ax.set_xlabel("Time [min]")
    ax.legend(loc="lower right", frameon=False)
    ax.grid(alpha=0.3)
    if dfn_traces:
        ax = axes[2]
        for (temperature_C, trace), color in zip(dfn_traces, ("#1f4e79", "#2e7d32", "#c55a11")):
            if trace.size:
                ax.plot(trace[:, 1], trace[:, 3] * 1000.0, color=color, label=f"{temperature_C:.0f} degC")
        ax.axhline(0.0, color="#b71c1c", linewidth=1.0, label="plating onset")
        ax.set_xlabel("SOC [-]")
        ax.set_ylabel("DFN min anode potential [mV]")
        ax.set_ylim(-5, 100)
        ax.legend(loc="upper right", frameon=False)
        ax.grid(alpha=0.3)
    fig.suptitle("EV MODIFY C model vs PyBaMM (Chen2020)")
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    fig.savefig(out_path, dpi=110)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cli", help="path to ev_cli (default: $EV_CLI or ./build)")
    parser.add_argument("--out", default=os.path.join(m.ROOT, "build", "pybamm_benchmark.png"), help="output plot")
    parser.add_argument("--skip-dfn", action="store_true", help="skip the (slow) closed-loop DFN plating check")
    args = parser.parse_args()
    try:
        cli = m.find_ev_cli(args.cli)
    except FileNotFoundError as error:
        print(f"error: {error}")
        return 2
    if cli is None:
        print("ev_cli not found: build with CMake or pass --cli / set EV_CLI.")
        return 2

    fit, _ = simulate_spme(["Charge at 0.5C for 15 minutes", "Rest for 15 minutes"] * 3 + ["Charge at 0.5C until 4.2 V"])
    params = fit_rc(fit, c_ocv(cli, fit["soc"]))
    val, capacity_Ah = simulate_spme(VALIDATION_PROFILE)
    val_ocv = c_ocv(cli, val["soc"])
    ocv_error = float(np.max(np.abs(val_ocv - val["bulk_ocv_V"])))
    rc_voltage_V = rc_model_voltage(params, val["time_s"], val["current_A"], val_ocv)
    rc_rms = float(np.sqrt(np.mean((rc_voltage_V - val["voltage_V"]) ** 2)))

    checks = [("OCV max error vs PyBaMM bulk OCV [V]", ocv_error, OCV_TOLERANCE_V),
              ("1-RC out-of-sample RMS voltage error [V]", rc_rms, RC_FIT_RMS_TOLERANCE_V)]
    settled = val["time_s"] >= EKF_SETTLE_S
    nominal = None
    cases = [("fitted", params, capacity_Ah), ("R0 +20 %", params * [1.2, 1, 1], capacity_Ah),
             ("R0 -20 %", params * [0.8, 1, 1], capacity_Ah), ("capacity -5 %", params, 0.95 * capacity_Ah)]
    for name, case_params, case_capacity in cases:
        ekf_soc, ekf_sigma, fault_samples = run_ekf(cli, case_params, case_capacity, val)
        error = np.abs(ekf_soc - val["soc"])[settled]
        outside = float(np.mean(error > 3.0 * ekf_sigma[settled]))
        checks.append((f"EKF ({name}) max SOC error after {EKF_SETTLE_S:.0f} s", float(np.max(error)), EKF_SOC_TOLERANCE))
        checks.append((f"EKF ({name}) fraction outside 3 sigma", outside, EKF_MAX_FRACTION_OUTSIDE_3SIGMA))
        checks.append((f"EKF ({name}) samples with a fault", fault_samples, 0))
        if nominal is None:
            nominal = (ekf_soc, ekf_sigma)

    dfn_traces = []
    if not args.skip_dfn:
        for temperature_C in DFN_TEMPERATURES_C:
            min_potential, faults, trace = dfn_closed_loop(cli, temperature_C, params)
            dfn_traces.append((temperature_C, trace))
            final_soc = trace[-1, 1] if trace.size else DFN_START_SOC
            if not trace.size:
                min_potential = -np.inf   # the limiter never charged: no plating evidence at all
            print(f"DFN {temperature_C:.0f} degC: charged {DFN_START_SOC:.2f} -> {final_soc:.3f} SOC in {trace.shape[0] * DFN_STEP_S / 60:.0f} min")
            checks.append((f"DFN {temperature_C:.0f} degC min anode potential (>= {PLATING_ONSET_V} V) [V]",
                           -min_potential, -PLATING_ONSET_V))
            checks.append((f"DFN {temperature_C:.0f} degC final SOC (>= {DFN_MIN_FINAL_SOC[temperature_C]})",
                           -final_soc, -DFN_MIN_FINAL_SOC[temperature_C]))
            checks.append((f"DFN {temperature_C:.0f} degC cycles with a fault", faults, 0))

    print(f"SOC window capacity {capacity_Ah:.4f} Ah; fitted 1-RC: R0 = {params[0] * 1000:.2f} mOhm, "
          f"R1 = {params[1] * 1000:.2f} mOhm, tau = {params[2]:.1f} s")
    ok = True
    for name, value, limit in checks:
        passed = value <= limit
        ok = ok and passed
        shown = -value if (name.startswith("DFN") and "fault" not in name) else value
        print(f"{'PASS' if passed else 'FAIL'}: {name} = {shown:.4f}")

    plot(args.out, val, rc_voltage_V, rc_rms, nominal[0], nominal[1], dfn_traces)
    print(f"Wrote {args.out}")
    print("BENCHMARK: " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
