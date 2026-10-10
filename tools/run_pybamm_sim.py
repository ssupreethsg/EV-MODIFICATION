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
   PyBaMM DFN in a closed loop from 5 % SOC (the DFN terminal voltage is the measurement). The DFN is
   an independent reference: its Arrhenius activation energies (electrolyte conductivity and
   diffusivity, graphite solid diffusivity, and for some cases the graphite exchange-current
   density) are NOT the controller's placeholders but a sweep: a 2x2 grid electrolyte {17, 30} kJ/mol
   x graphite D_s {35, 60} kJ/mol at 2 and 10 degC (25 degC is insensitive: one case), all with
   Chen2020's own exchange-current activation energy (35 kJ/mol), plus the harsh corner (30 / 60) with
   an exchange-current activation energy of 53.4 kJ/mol (Ecker 2015) and 60 kJ/mol (an Arrhenius
   multiplier on the negative-electrode exchange-current density relative to 298.15 K, DFN only). The
   worst case per temperature is reported. In every case the DFN anode potential vs Li/Li+
   (phi_s - phi_e) must stay >= 5 mV (the controller's own margin) at every point at the end of every 10 s step, the
   DFN terminal voltage must stay <= 4.2 V + 5 mV, and no fault may occur. The final-SOC floors (the
   limiter must still charge) are evaluated on the nominal case (17 / 35 kJ/mol). Extra cases: 2 degC
   from 75 % SOC on the harshest corner including i0 60 kJ/mol (exercises the cold CV taper); a
   lumped-thermal DFN at 25 degC with the controller fed the DFN cell temperature (max <= 45 degC);
   O'Regan2022's graphite diffusivity and exchange current (fitted together on the same cell, with
   their own temperature dependence; electrolyte 30 kJ/mol) at 2 and 25 degC and in the 0 degC
   sensor-bias case. The isothermal cases cannot show self-heating; only the lumped case does.
   Sensor errors at their limits, on the harsh corner with i0 60 kJ/mol: the DFN is at 0 degC while
   the controller reads 2 degC (EV_TEMPERATURE_SENSOR_MARGIN_C). (a) From 5 % SOC with an exact
   charger. (b) From 50 % SOC, the charger delivers 1.05 x command + 0.99 x C/50 + C/100 (just inside
   the fast over-current tolerance as measured, since the current sensor reads C/100 low): only the over-charge integral may
   fault, not before 1500 s. (c) The charger delivers 1.05 x command + 1.5 x C/50, outside the
   tolerance: EV_FAULT_OVER_CURRENT (0x800) must appear within 30 s. (d) From 5 % and from 75 % SOC,
   the charger delivers 1.05 x command + 2.45 x C/100 with the sensor C/100 low: a sustained
   excess the over-charge integral never reaches; no fault allowed. Each must stay >= 5 mV, and (a),
   (b) and (d) must charge at least 5 % SOC, so they cannot pass by never charging. In every
   case a charger obeys a stop: when the controller approves 0 A (end of charge or a fault) the
   plant delivers 0 A, and the run ends.

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
from scipy.optimize import least_squares  # noqa: E402

import ev_model as m  # noqa: E402

pybamm = None  # imported in main() after argument parsing, so that --help works without PyBaMM

FARADAY = 96485.33212
INITIAL_SOC = 0.1
TEMPERATURE_C = 25.0
EKF_INITIAL_SOC_ERROR = 0.2
EKF_SETTLE_S = 900.0  # end of the first rest: the filter learns SOC mainly at low current
OCV_TOLERANCE_V = 0.001
RC_FIT_RMS_TOLERANCE_V = 0.025
EKF_SOC_TOLERANCE = 0.03
EKF_MAX_FRACTION_OUTSIDE_3SIGMA = 0.01
DFN_ELECTROLYTE_EA_J_MOL = (
    17000.0,
    30000.0,
)  # DFN-only sweep, independent of the controller constants
DFN_SOLID_EA_J_MOL = (35000.0, 60000.0)
DFN_NOMINAL_EA_J_MOL = (
    17000.0,
    35000.0,
)  # (electrolyte, solid): the nominal case for the SOC floors
DFN_HARSH_EA_J_MOL = (30000.0, 60000.0)
DFN_I0_EA_J_MOL = (53400.0, 60000.0)  # graphite exchange-current activation energy sweep (harsh corner)
DFN_I0_NOMINAL_EA_J_MOL = 35000.0  # Chen2020's own value
DFN_COLD_CV_START_SOC = 0.75
DFN_TOLERANCE_START_SOC = 0.5
DFN_CHARGER_REL = 1.0 + m.header_constant("EV_OVER_CURRENT_TOLERANCE_REL")
DFN_TOLERANCE_C_RATE = m.header_constant("EV_OVER_CURRENT_TOLERANCE_C_RATE")  # fast over-current tolerance [1/h]
DFN_SENSOR_OFFSET_C_RATE = m.header_constant("EV_CURRENT_OFFSET_C_RATE")  # current-sensor offset requirement [1/h]
DFN_EDGE_FRACTION = 0.99  # charger excess as a fraction of the tolerance: just inside it
DFN_OUT_FRACTION = 1.5  # ... and outside it
DFN_SUSTAINED_FRACTION = 2.45  # sustained excess in units of C/100 just below what the over-charge integral can reach
FAULT_OVER_CURRENT = 0x800  # EV_FAULT_OVER_CURRENT
DFN_OVER_CURRENT_DEADLINE_S = 30.0  # outside the tolerance: the fast check (1 s)
DFN_OVER_CHARGE_EARLIEST_S = 1500.0  # just inside it: only the over-charge integral may fault, and not before this
DFN_SENSOR_BIAS_C = m.header_constant("EV_TEMPERATURE_SENSOR_MARGIN_C")  # temperature reading warmer than the cell
DFN_BIAS_MIN_CHARGE = 0.05  # SOC the sensor-bias cases must charge, so that they are not vacuous
DFN_MAX_VOLTAGE_V = 4.2 + 0.005
DFN_MAX_TEMPERATURE_C = 45.0
# The limiter must keep charging, not just stay safe. These are anti-vacuity floors, ~5 points below
# the measured result (46.8 / 65.1 / 94.9 %); every limit budgets the charger tolerance (5 % + C/50)
# and the C/100 current-sensor offset.
DFN_MIN_FINAL_SOC = {2.0: 0.42, 10.0: 0.60, 25.0: 0.90}
DFN_START_SOC = 0.05
DFN_STEP_S = 10.0
DFN_MAX_TIME_S = 2.0 * 3600.0
DFN_REQUEST_A = 15.0
DFN_HORIZON_S = 30.0
VALIDATION_PROFILE = [
    "Charge at 1C for 10 minutes or until 4.2 V",
    "Rest for 10 minutes",
] * 4 + ["Charge at 1C until 4.2 V"]
PLATING_FLOOR_V = m.header_constant("EV_MIN_SAFE_ANODE_POTENTIAL_V")  # the controller's own margin, not 0 V


def chen2020(
    temperature_C=TEMPERATURE_C,
    activation_J_mol=None,
    i0_activation_J_mol=None,
    oregan_anode=False,
    mixed_anode=False,
    measured_electrolyte=False,
):
    """Chen2020 parameters; activation_J_mol = (electrolyte, solid) adds Arrhenius terms (DFN only).
    i0_activation_J_mol multiplies the negative-electrode exchange-current density by an Arrhenius
    factor exp(Ea - 35 kJ/mol) relative to 298.15 K, i.e. replaces Chen2020's own 35 kJ/mol (DFN only).
    oregan_anode=True replaces the graphite diffusivity and exchange-current density by the pair
    O'Regan et al. (2022) measured on the same LG M50 cell, with their own temperature dependence;
    the two are fitted together, so they are swapped together (DFN only). mixed_anode=True takes
    O'Regan2022's (lower) graphite diffusivity but keeps Chen2020's (lower) exchange-current density:
    the worst pairing of the two sets, not a consistent fit (DFN only). measured_electrolyte=True
    replaces the Nyman2008 electrolyte (ideal solution, Arrhenius placeholders) by Landesfeind &
    Gasteiger's (2019) measured LiPF6 in EC:EMC 3:7 functions of concentration and temperature:
    conductivity, diffusivity, transference number and thermodynamic factor, as in PyBaMM ORegan2022
    (activation_J_mol then applies to the graphite diffusivity only)."""
    parameter_values = pybamm.ParameterValues("Chen2020")
    x0, x100, y100, y0 = pybamm.lithium_ion.get_min_max_stoichiometries(parameter_values)
    for name, value, header in (
        ("x0", x0, m.ANODE_STOICH_SOC0),
        ("x100", x100, m.ANODE_STOICH_SOC100),
        ("y0", y0, m.CATHODE_STOICH_SOC0),
        ("y100", y100, m.CATHODE_STOICH_SOC100),
    ):
        if abs(value - header) > 1e-5:
            raise RuntimeError(f"PyBaMM stoichiometry {name} = {value:.6f} differs from the C header ({header:.6f})")
    if activation_J_mol is not None:
        electrolyte_ea, solid_ea = activation_J_mol

        def factor(activation_energy, temperature):
            return pybamm.exp(
                activation_energy / m.header_constant("EV_GAS_J_PER_MOL_K") * (1.0 / 298.15 - 1.0 / temperature)
            )

        kappa = parameter_values["Electrolyte conductivity [S.m-1]"]
        diffusivity = parameter_values["Electrolyte diffusivity [m2.s-1]"]
        solid_diffusivity = parameter_values["Negative particle diffusivity [m2.s-1]"]
        parameter_values.update(
            {
                "Electrolyte conductivity [S.m-1]": lambda c, T: kappa(c, T) * factor(electrolyte_ea, T),
                "Electrolyte diffusivity [m2.s-1]": lambda c, T: diffusivity(c, T) * factor(electrolyte_ea, T),
                "Negative particle diffusivity [m2.s-1]": lambda sto, T: solid_diffusivity * factor(solid_ea, T),
            }
        )
    from pybamm.input.parameters.lithium_ion import ORegan2022

    if measured_electrolyte:
        parameter_values.update(
            {
                "Electrolyte conductivity [S.m-1]": ORegan2022.electrolyte_conductivity_EC_EMC_3_7_Landesfeind2019,
                "Electrolyte diffusivity [m2.s-1]": ORegan2022.electrolyte_diffusivity_EC_EMC_3_7_Landesfeind2019,
                "Cation transference number": ORegan2022.electrolyte_transference_number_EC_EMC_3_7_Landesfeind2019,
                "Thermodynamic factor": ORegan2022.electrolyte_TDF_EC_EMC_3_7_Landesfeind2019,
            }
        )
    if mixed_anode:
        parameter_values.update(
            {"Negative particle diffusivity [m2.s-1]": ORegan2022.graphite_LGM50_diffusivity_ORegan2022}
        )
    if oregan_anode:
        parameter_values.update(
            {
                "Negative particle diffusivity [m2.s-1]": ORegan2022.graphite_LGM50_diffusivity_ORegan2022,
                "Negative electrode exchange-current density [A.m-2]": (
                    ORegan2022.graphite_LGM50_electrolyte_exchange_current_density_ORegan2022
                ),
            }
        )
    if i0_activation_J_mol is not None:
        exchange_current = parameter_values["Negative electrode exchange-current density [A.m-2]"]
        extra_ea = i0_activation_J_mol - DFN_I0_NOMINAL_EA_J_MOL

        def exchange_current_ea(c_e, c_s_surf, c_s_max, T):
            return exchange_current(c_e, c_s_surf, c_s_max, T) * pybamm.exp(
                extra_ea / m.header_constant("EV_GAS_J_PER_MOL_K") * (1.0 / 298.15 - 1.0 / T)
            )

        parameter_values.update({"Negative electrode exchange-current density [A.m-2]": exchange_current_ea})
    parameter_values.update(
        {
            "Ambient temperature [K]": 273.15 + temperature_C,
            "Initial temperature [K]": 273.15 + temperature_C,
        }
    )
    capacity_Ah = (
        parameter_values["Negative electrode active material volume fraction"]
        * parameter_values["Negative electrode thickness [m]"]
        * parameter_values["Electrode height [m]"]
        * parameter_values["Electrode width [m]"]
        * parameter_values["Maximum concentration in negative electrode [mol.m-3]"]
        * FARADAY
        / 3600.0
        * (x100 - x0)
    )
    return parameter_values, capacity_Ah


def simulate_spme(steps):
    parameter_values, capacity_Ah = chen2020()
    experiment = pybamm.Experiment(steps, period="10 seconds")
    sol = pybamm.Simulation(
        pybamm.lithium_ion.SPMe(),
        experiment=experiment,
        parameter_values=parameter_values,
    ).solve(initial_soc=INITIAL_SOC)
    time_s = sol["Time [s]"].entries
    keep = np.concatenate(([True], np.diff(time_s) > 1e-6))  # keep the ending step's boundary sample
    data = {
        "time_s": time_s[keep],
        "current_A": -sol["Current [A]"].entries[keep],  # PyBaMM: positive = discharge
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
    result = least_squares(
        lambda p: rc_model_voltage(p, data["time_s"], data["current_A"], ocv_V) - data["voltage_V"],
        x0=[0.02, 0.02, 100.0],
        bounds=([0.0, 0.0, 1.0], [0.5, 0.5, 5000.0]),
        x_scale=[0.01, 0.01, 100.0],
    )
    return result.x


def c_ocv(cli, soc):
    return np.array(
        m.run_cli(cli, ["ocv"], [f"{s:.9f}" for s in np.clip(soc, 0.0, 1.0)]),
        dtype=float,
    )[:, 0]


def run_ekf(cli, params, capacity_Ah, data):
    r0, r1, tau = params
    t, i_a, v = data["time_s"], data["current_A"], data["voltage_V"]
    lines = [f"{t[k] - t[k - 1]:.6f} {i_a[k]:.6f} {v[k]:.6f} {TEMPERATURE_C:.3f}" for k in range(1, t.size)]
    args = [
        "ekf",
        f"{capacity_Ah:.6f}",
        f"{r0:.6f}",
        f"{r1:.6f}",
        f"{tau:.6f}",
        f"{INITIAL_SOC + EKF_INITIAL_SOC_ERROR:.6f}",
        f"{TEMPERATURE_C:.3f}",
    ]
    rows = np.array(m.run_cli(cli, args, lines), dtype=float)
    soc = np.concatenate(([INITIAL_SOC + EKF_INITIAL_SOC_ERROR], rows[:, 0]))
    sigma = np.concatenate(([0.05], rows[:, 1]))
    fault_samples = int(np.count_nonzero(rows[:, 4]))
    return soc, sigma, fault_samples


def dfn_closed_loop(
    cli,
    temperature_C,
    params,
    activation_J_mol,
    start_soc=DFN_START_SOC,
    thermal=False,
    i0_activation_J_mol=None,
    charger_offset_A=None,
    sensor_offset_A=0.0,
    temperature_bias_C=0.0,
    oregan_anode=False,
    mixed_anode=False,
    measured_electrolyte=False,
):
    """Complete C controller charging a DFN. Returns a dict: min anode potential vs Li/Li+ [V], max
    terminal voltage [V], max cell temperature [degC], cycles with fault flags, final SOC and the
    trace (time, SOC, current, min potential) and the time of the first over-current fault flag (None
    if never). With thermal=True the DFN is lumped-thermal (ambient temperature_C) and the controller
    is fed the DFN cell temperature. charger_offset_A=None: the plant receives exactly the approved
    command; otherwise it receives DFN_CHARGER_REL x command + charger_offset_A (0 A when the command
    is 0 A: the charger obeys a stop). The controller reads the current sensor_offset_A low and the
    temperature temperature_bias_C high (sensor errors in the unsafe direction)."""
    parameter_values, capacity_Ah = chen2020(
        temperature_C, activation_J_mol, i0_activation_J_mol, oregan_anode, mixed_anode, measured_electrolyte
    )
    parameter_values.update({"Current function [A]": "[input]"})
    parameter_values.set_initial_state(start_soc)
    model = pybamm.lithium_ion.DFN({"thermal": "lumped"} if thermal else {})
    sim = pybamm.Simulation(model, parameter_values=parameter_values)
    r0, r1, tau = params
    process = subprocess.Popen(
        [
            cli,
            "charge",
            f"{capacity_Ah:.6f}",
            f"{r0:.6f}",
            f"{r1:.6f}",
            f"{tau:.6f}",
            f"{start_soc:.6f}",
            f"{temperature_C + temperature_bias_C:.3f}",  # initial reading and coolant, as the BMS sees them
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True,
    )
    soc, elapsed, min_potential, faults = start_soc, 0.0, np.inf, 0
    over_current_time_s = first_fault_s = None
    fault_flags = 0
    max_voltage, max_temperature_C = -np.inf, temperature_C
    trace = []
    current_A = 0.0  # current over the interval that just ended (a rested start)
    try:
        sol = sim.step(DFN_STEP_S, inputs={"Current function [A]": 0.0}, save=False)
        while elapsed < DFN_MAX_TIME_S and soc < m.header_constant("EV_MAX_SAFE_SOC"):
            voltage_V = float(sol["Voltage [V]"].entries[-1])
            cell_C = (
                float(sol["Volume-averaged cell temperature [K]"].entries[-1]) - 273.15 if thermal else temperature_C
            ) + temperature_bias_C
            measured_A = current_A - sensor_offset_A if current_A > 0.0 else current_A
            process.stdin.write(
                f"{DFN_STEP_S} {measured_A:.6f} {voltage_V:.6f} {cell_C:.3f} {DFN_REQUEST_A} {DFN_HORIZON_S}\n"
            )
            process.stdin.flush()
            line = process.stdout.readline()
            if not line:
                raise RuntimeError("ev_cli charge stopped responding")
            fields = line.split()
            command_A, flags = float(fields[0]), int(fields[2])
            faults += int(flags != 0)
            fault_flags |= flags
            if flags and first_fault_s is None:
                first_fault_s = elapsed
            if (flags & FAULT_OVER_CURRENT) and over_current_time_s is None:
                over_current_time_s = elapsed
            if command_A < 0.01 and trace:
                break  # end of charge (or a fault, which is counted); the charger delivers 0 A
            current_A = command_A
            if charger_offset_A is not None and command_A >= 0.01:
                current_A = DFN_CHARGER_REL * command_A + charger_offset_A
            sol = sim.step(DFN_STEP_S, inputs={"Current function [A]": -current_A}, save=False)
            potential = sol["Negative electrode surface potential difference [V]"].entries  # every point and sample
            min_potential = min(min_potential, float(np.min(potential)))
            max_voltage = max(max_voltage, float(np.max(sol["Voltage [V]"].entries)))
            if thermal:
                max_temperature_C = max(
                    max_temperature_C,
                    float(np.max(sol["Volume-averaged cell temperature [K]"].entries)) - 273.15,
                )
            soc = start_soc - float(sol["Discharge capacity [A.h]"].entries[-1]) / capacity_Ah
            elapsed += DFN_STEP_S
            trace.append((elapsed, soc, current_A, float(np.min(potential))))
    finally:
        process.stdin.close()
        returncode = process.wait()
    if returncode != 0:
        raise RuntimeError(f"ev_cli charge exited with status {returncode}")
    trace = np.array(trace)
    if not trace.size:  # the limiter never charged: no plating evidence at all
        min_potential, max_voltage = -np.inf, np.inf
    return {
        "min_potential": min_potential,
        "max_voltage": max_voltage,
        "max_temperature": max_temperature_C,
        "faults": faults,
        "over_current_time_s": over_current_time_s,
        "first_fault_s": first_fault_s,
        "fault_flags": fault_flags,
        "final_soc": float(trace[-1, 1]) if trace.size else start_soc,
        "trace": trace,
    }


def plot(out_path, val, rc_voltage_V, rc_rms, ekf_soc, ekf_sigma, dfn_traces):
    rows = 3 if dfn_traces else 2
    fig, axes = plt.subplots(rows, 1, figsize=(8, 3.2 * rows), constrained_layout=True)
    minutes = val["time_s"] / 60.0
    ax = axes[0]
    ax.plot(minutes, val["voltage_V"], color="#1f4e79", label="PyBaMM SPMe")
    ax.plot(
        minutes,
        rc_voltage_V,
        "--",
        color="#c55a11",
        label=f"C 1-RC model, fitted on 0.5C data (RMS {rc_rms * 1000:.1f} mV)",
    )
    ax.set_ylabel("Terminal voltage [V]")
    ax.set_xlabel("Time [min]")
    ax.legend(loc="lower right", frameon=False)
    ax.grid(alpha=0.3)
    ax = axes[1]
    ax.plot(minutes, val["soc"], color="#1f4e79", label="PyBaMM SOC")
    ax.plot(
        minutes,
        ekf_soc,
        "--",
        color="#c55a11",
        label=f"C EKF started {EKF_INITIAL_SOC_ERROR:+.2f} SOC off",
    )
    ax.fill_between(
        minutes,
        ekf_soc - 3 * ekf_sigma,
        ekf_soc + 3 * ekf_sigma,
        color="#c55a11",
        alpha=0.15,
        label="EKF 3-sigma",
    )
    ax.set_ylabel("SOC [-]")
    ax.set_xlabel("Time [min]")
    ax.legend(loc="lower right", frameon=False)
    ax.grid(alpha=0.3)
    if dfn_traces:
        ax = axes[2]
        for (label, trace), color in zip(dfn_traces, ("#1f4e79", "#2e7d32", "#c55a11", "#6a1b9a", "#00838f")):
            if trace.size:
                ax.plot(trace[:, 1], trace[:, 3] * 1000.0, color=color, label=label)
        ax.axhline(0.0, color="#b71c1c", linewidth=1.0, label="plating onset")
        ax.set_xlabel("SOC [-]")
        ax.set_ylabel("DFN min anode potential [mV]\n(nominal Ea; worst-case in text)")
        ax.set_ylim(-5, 100)
        ax.legend(loc="upper right", frameon=False, fontsize=8)
        ax.grid(alpha=0.3)
    fig.suptitle("EV MODIFY C model vs PyBaMM (Chen2020)")
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    fig.savefig(out_path, dpi=110)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cli", help="path to ev_cli (default: $EV_CLI or ./build)")
    parser.add_argument(
        "--out",
        default=os.path.join(m.ROOT, "build", "pybamm_benchmark.png"),
        help="output plot",
    )
    parser.add_argument(
        "--skip-dfn",
        action="store_true",
        help="skip the (slow) closed-loop DFN plating check",
    )
    args = parser.parse_args()
    global pybamm
    import pybamm

    try:
        cli = m.find_ev_cli(args.cli)
    except FileNotFoundError as error:
        print(f"error: {error}")
        return 2
    if cli is None:
        print("ev_cli not found: build with CMake or pass --cli / set EV_CLI.")
        return 2

    fit, _ = simulate_spme(
        ["Charge at 0.5C for 15 minutes", "Rest for 15 minutes"] * 3 + ["Charge at 0.5C until 4.2 V"]
    )
    params = fit_rc(fit, c_ocv(cli, fit["soc"]))
    val, capacity_Ah = simulate_spme(VALIDATION_PROFILE)
    val_ocv = c_ocv(cli, val["soc"])
    ocv_error = float(np.max(np.abs(val_ocv - val["bulk_ocv_V"])))
    rc_voltage_V = rc_model_voltage(params, val["time_s"], val["current_A"], val_ocv)
    rc_rms = float(np.sqrt(np.mean((rc_voltage_V - val["voltage_V"]) ** 2)))

    checks = [
        ("OCV max error vs PyBaMM bulk OCV [V]", ocv_error, OCV_TOLERANCE_V),
        ("1-RC out-of-sample RMS voltage error [V]", rc_rms, RC_FIT_RMS_TOLERANCE_V),
    ]
    settled = val["time_s"] >= EKF_SETTLE_S
    nominal = None
    cases = [
        ("fitted", params, capacity_Ah),
        ("R0 +20 %", params * [1.2, 1, 1], capacity_Ah),
        ("R0 -20 %", params * [0.8, 1, 1], capacity_Ah),
        ("capacity -5 %", params, 0.95 * capacity_Ah),
    ]
    for name, case_params, case_capacity in cases:
        ekf_soc, ekf_sigma, fault_samples = run_ekf(cli, case_params, case_capacity, val)
        error = np.abs(ekf_soc - val["soc"])[settled]
        outside = float(np.mean(error > 3.0 * ekf_sigma[settled]))
        checks.append(
            (
                f"EKF ({name}) max SOC error after {EKF_SETTLE_S:.0f} s",
                float(np.max(error)),
                EKF_SOC_TOLERANCE,
            )
        )
        checks.append(
            (
                f"EKF ({name}) fraction outside 3 sigma",
                outside,
                EKF_MAX_FRACTION_OUTSIDE_3SIGMA,
            )
        )
        checks.append((f"EKF ({name}) samples with a fault", fault_samples, 0))
        if nominal is None:
            nominal = (ekf_soc, ekf_sigma)

    dfn_traces = []
    if not args.skip_dfn:

        def run_case(
            label,
            temperature_C,
            activation,
            expect_over_current=False,
            allow_over_charge=False,
            min_charge=DFN_BIAS_MIN_CHARGE,
            **kwargs,
        ):
            result = dfn_closed_loop(cli, temperature_C, params, activation, **kwargs)
            if min_charge > 0.0 and not expect_over_current:  # a case that never charged would pass vacuously
                charged = result["final_soc"] - kwargs.get("start_soc", DFN_START_SOC)
                checks.append((f"DFN {label} SOC charged (>= {min_charge})", -charged, -min_charge))
            over_current_s = result["over_current_time_s"]
            print(
                f"DFN case {label}: SOC {kwargs.get('start_soc', DFN_START_SOC):.2f} -> {result['final_soc']:.3f}, "
                f"min anode potential {result['min_potential'] * 1000:.1f} mV, max V {result['max_voltage']:.4f}, "
                f"max T {result['max_temperature']:.1f} degC, {result['trace'].shape[0] * DFN_STEP_S / 60:.0f} min, "
                f"faults {result['faults']}, over-current fault at "
                f"{'never' if over_current_s is None else f'{over_current_s:.0f} s'}",
                flush=True,
            )
            checks.append(
                (
                    f"DFN {label} min anode potential (>= {PLATING_FLOOR_V} V) [V]",
                    -result["min_potential"],
                    -PLATING_FLOOR_V,
                )
            )
            checks.append(
                (
                    f"DFN {label} max terminal voltage (<= {DFN_MAX_VOLTAGE_V:.3f} V) [V]",
                    result["max_voltage"],
                    DFN_MAX_VOLTAGE_V,
                )
            )
            if expect_over_current:
                checks.append(
                    (
                        f"DFN {label} time to over-current fault (<= {DFN_OVER_CURRENT_DEADLINE_S:.0f} s) [s]",
                        np.inf if over_current_s is None else over_current_s,
                        DFN_OVER_CURRENT_DEADLINE_S,
                    )
                )
            elif allow_over_charge:
                # only the over-charge integral may stop a charger inside the fast tolerance, and only late
                late = result["first_fault_s"] is None or (
                    result["fault_flags"] & ~FAULT_OVER_CURRENT == 0
                    and result["first_fault_s"] >= DFN_OVER_CHARGE_EARLIEST_S
                )
                checks.append((f"DFN {label} early or non-over-current fault (0 = none)", int(not late), 0))
            else:
                checks.append((f"DFN {label} cycles with a fault", result["faults"], 0))
            return result

        harsh = f"Ee={DFN_HARSH_EA_J_MOL[0] / 1000:.0f} Ed={DFN_HARSH_EA_J_MOL[1] / 1000:.0f}"
        for temperature_C in (2.0, 10.0, 25.0):
            # (electrolyte Ea, solid Ea, i0 Ea or None for Chen2020's own 35 kJ/mol)
            grid = (
                [(e, s_, None) for e in DFN_ELECTROLYTE_EA_J_MOL for s_ in DFN_SOLID_EA_J_MOL]
                + [(*DFN_HARSH_EA_J_MOL, i0_ea) for i0_ea in DFN_I0_EA_J_MOL]
                if temperature_C < 20.0
                else [(*DFN_NOMINAL_EA_J_MOL, None)]
            )
            worst = np.inf
            for electrolyte_ea, solid_ea, i0_ea in grid:
                label = f"{temperature_C:.0f} degC Ee={electrolyte_ea / 1000:.0f} Ed={solid_ea / 1000:.0f} kJ/mol"
                label += f" i0={(DFN_I0_NOMINAL_EA_J_MOL if i0_ea is None else i0_ea) / 1000:g} kJ/mol"
                result = run_case(label, temperature_C, (electrolyte_ea, solid_ea), i0_activation_J_mol=i0_ea)
                worst = min(worst, result["min_potential"])
                if (electrolyte_ea, solid_ea) == DFN_NOMINAL_EA_J_MOL:
                    dfn_traces.append((f"{temperature_C:.0f} degC", result["trace"]))
                    checks.append(
                        (
                            f"DFN {temperature_C:.0f} degC nominal final SOC (>= {DFN_MIN_FINAL_SOC[temperature_C]})",
                            -result["final_soc"],
                            -DFN_MIN_FINAL_SOC[temperature_C],
                        )
                    )
            print(
                f"DFN {temperature_C:.0f} degC worst-case min anode potential over {len(grid)} activation-energy "
                f"case(s): {worst * 1000:.1f} mV",
                flush=True,
            )
        i0_worst = DFN_I0_EA_J_MOL[-1]
        run_case(
            f"2 degC from SOC {DFN_COLD_CV_START_SOC:.2f} {harsh} i0={i0_worst / 1000:g} kJ/mol",
            2.0,
            DFN_HARSH_EA_J_MOL,
            start_soc=DFN_COLD_CV_START_SOC,
            i0_activation_J_mol=i0_worst,
        )
        # Sensor errors at their limits in the unsafe direction: the cell is EV_TEMPERATURE_SENSOR_MARGIN_C
        # colder than its reading and the current sensor reads C/100 low; the charger delivers 1.05 x
        # command plus an excess just inside, or outside, the fast tolerance (C/50).
        capacity = chen2020(0.0)[1]
        tolerance_A, sensor_A = DFN_TOLERANCE_C_RATE * capacity, DFN_SENSOR_OFFSET_C_RATE * capacity
        bias = f"cell {DFN_SENSOR_BIAS_C:.0f} degC below its reading"
        run_case(
            f"0 degC ({bias}) {harsh} i0={i0_worst / 1000:g} kJ/mol",
            0.0,
            DFN_HARSH_EA_J_MOL,
            i0_activation_J_mol=i0_worst,
            temperature_bias_C=DFN_SENSOR_BIAS_C,
            min_charge=DFN_BIAS_MIN_CHARGE,
        )
        for name, excess_A, expect, allow in (
            (
                f"charger excess {DFN_EDGE_FRACTION:g} x tolerance, sensor {sensor_A:.3f} A low",
                DFN_EDGE_FRACTION * tolerance_A + sensor_A,
                False,
                True,
            ),
            (f"charger excess {DFN_OUT_FRACTION:g} x tolerance (outside)", DFN_OUT_FRACTION * tolerance_A, True, False),
        ):
            run_case(
                f"0 degC ({bias}) from SOC {DFN_TOLERANCE_START_SOC:.2f} {harsh} i0={i0_worst / 1000:g} kJ/mol, {name}",
                0.0,
                DFN_HARSH_EA_J_MOL,
                expect_over_current=expect,
                allow_over_charge=allow,
                start_soc=DFN_TOLERANCE_START_SOC,
                i0_activation_J_mol=i0_worst,
                charger_offset_A=excess_A,
                sensor_offset_A=sensor_A if allow else 0.0,
                temperature_bias_C=DFN_SENSOR_BIAS_C,
                min_charge=DFN_BIAS_MIN_CHARGE,
            )
        # A sustained charger excess that never faults: the integral counts measured charge above
        # 1.05 x command + C/100 and leaks with 1 h, so it settles at the excess x 1 h and never reaches
        # its limit (C/200 x 1 h) below a measured excess of C/200. A charger at 1.05 x command +
        # 2.45 x C/100 with the sensor reading C/100 low (measured excess 0.45 x C/100) runs indefinitely. From low SOC (kinetics) and from 75 % SOC (CV taper).
        for start_soc in (DFN_START_SOC, DFN_COLD_CV_START_SOC):
            run_case(
                f"0 degC ({bias}) from SOC {start_soc:.2f} {harsh} i0={i0_worst / 1000:g} kJ/mol, "
                f"sustained charger excess {DFN_SUSTAINED_FRACTION:g} x C/100, sensor {sensor_A:.3f} A low",
                0.0,
                DFN_HARSH_EA_J_MOL,
                start_soc=start_soc,
                i0_activation_J_mol=i0_worst,
                charger_offset_A=DFN_SUSTAINED_FRACTION * sensor_A,
                sensor_offset_A=sensor_A,
                temperature_bias_C=DFN_SENSOR_BIAS_C,
                min_charge=DFN_BIAS_MIN_CHARGE,
            )
        # An independent graphite transport/kinetics pair for the same cell (O'Regan et al. 2022), with
        # its own temperature dependence; the electrolyte uses the harsh 30 kJ/mol. Cold, where its
        # diffusivity is lowest, including the sensor-bias case.
        for temperature_C, bias_C in ((2.0, 0.0), (0.0, DFN_SENSOR_BIAS_C), (25.0, 0.0)):
            run_case(
                f"{temperature_C:.0f} degC{f' ({bias})' if bias_C else ''} O'Regan2022 graphite D_s and i0, Ee=30",
                temperature_C,
                DFN_HARSH_EA_J_MOL,
                oregan_anode=True,
                temperature_bias_C=bias_C,
                min_charge=DFN_BIAS_MIN_CHARGE,
            )
        # Landesfeind & Gasteiger's measured EC:EMC 3:7 electrolyte (conductivity, diffusivity, transference
        # number and thermodynamic factor vs concentration and temperature) instead of Nyman2008 with
        # Arrhenius placeholders, combined with the worst pairing of the two graphite sets (O'Regan2022's
        # diffusivity with Chen2020's exchange current, raised to 60 kJ/mol when cold).
        for temperature_C, bias_C, start_soc, i0_ea in (
            (2.0, 0.0, DFN_START_SOC, i0_worst),
            (0.0, DFN_SENSOR_BIAS_C, DFN_START_SOC, i0_worst),
            (2.0, 0.0, DFN_COLD_CV_START_SOC, i0_worst),
            (10.0, 0.0, DFN_START_SOC, None),
            (25.0, 0.0, DFN_START_SOC, None),
        ):
            run_case(
                f"{temperature_C:.0f} degC{f' ({bias})' if bias_C else ''} from SOC {start_soc:.2f}, measured "
                f"electrolyte, O'Regan2022 D_s with Chen2020 i0" + (f" ({i0_ea / 1000:g} kJ/mol)" if i0_ea else ""),
                temperature_C,
                DFN_HARSH_EA_J_MOL if temperature_C < 20.0 else DFN_NOMINAL_EA_J_MOL,
                start_soc=start_soc,
                i0_activation_J_mol=i0_ea,
                mixed_anode=True,
                measured_electrolyte=True,
                temperature_bias_C=bias_C,
            )
        thermal = run_case("25 degC lumped-thermal", 25.0, DFN_NOMINAL_EA_J_MOL, thermal=True)
        checks.append(
            (
                f"DFN 25 degC lumped-thermal max cell temperature (<= {DFN_MAX_TEMPERATURE_C} degC) [degC]",
                thermal["max_temperature"],
                DFN_MAX_TEMPERATURE_C,
            )
        )

    print(
        f"SOC window capacity {capacity_Ah:.4f} Ah; fitted 1-RC: R0 = {params[0] * 1000:.2f} mOhm, "
        f"R1 = {params[1] * 1000:.2f} mOhm, tau = {params[2]:.1f} s"
    )
    ok = True
    for name, value, limit in checks:
        passed = value <= limit
        ok = ok and passed
        shown = (
            -value
            if (name.startswith("DFN") and ("anode potential" in name or "final SOC" in name or "SOC charged" in name))
            else value
        )
        print(f"{'PASS' if passed else 'FAIL'}: {name} = {shown:.4f}")

    plot(args.out, val, rc_voltage_V, rc_rms, nominal[0], nominal[1], dfn_traces)
    print(f"Wrote {args.out}")
    print("BENCHMARK: " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
