"""Fits the C model's lumped thermal parameters and ohmic resistance to measured data.

Usage:
    python calibrate_cell.py <test_data.csv>   # columns: time_s,current_A,temperature_C,voltage_V,initial_soc
    python calibrate_cell.py                   # synthetic self-test (parameters must be recovered)

Current is positive when charging. initial_soc is read from the first row.
"""
import sys

import numpy as np
from scipy.optimize import least_squares

import ev_model as m

# Must match spme_v7_apex_solver.h
EV_CELL_MASS_KG = 0.070
EV_CELL_COOLING_AREA_M2 = 0.0075
EV_AMBIENT_TEMPERATURE_C = 25.0

# Sensor noise used to weight the residuals (1 sigma)
TEMPERATURE_NOISE_C = 0.2
VOLTAGE_NOISE_V = 0.01

LOWER_BOUNDS = np.array([1.0, 100.0, 0.001])   # W/m2/K, J/kg/K, ohm
UPPER_BOUNDS = np.array([50.0, 2000.0, 0.100])


def simulate_cell(time_s, current_A, initial_soc, initial_temp_C, heat_transfer_W_m2_K,
                  heat_capacity_J_kg_K, resistance_ohm, capacity_Ah=m.NOMINAL_CAPACITY_AH):
    """Same coulomb counting, lumped thermal update and OCV + IR voltage as the C model."""
    thermal_conductance_W_K = heat_transfer_W_m2_K * EV_CELL_COOLING_AREA_M2
    heat_capacity_J_K = EV_CELL_MASS_KG * heat_capacity_J_kg_K

    n_steps = len(time_s)
    pred_temp = np.empty(n_steps)
    soc = np.empty(n_steps)
    pred_temp[0] = initial_temp_C
    soc[0] = initial_soc
    for i in range(1, n_steps):
        dt = time_s[i] - time_s[i - 1]
        current = current_A[i - 1]   # current applied over the interval
        soc[i] = np.clip(soc[i - 1] + current * dt / (capacity_Ah * 3600.0), 0.0, 1.0)
        steady_rise = current * current * resistance_ohm / thermal_conductance_W_K
        decay = np.exp(-thermal_conductance_W_K * dt / heat_capacity_J_K)
        pred_temp[i] = EV_AMBIENT_TEMPERATURE_C + steady_rise + (pred_temp[i - 1] - EV_AMBIENT_TEMPERATURE_C - steady_rise) * decay

    pred_voltage = m.full_cell_ocv_V(soc) + current_A * resistance_ohm
    return pred_temp, pred_voltage


def residuals(scaled, time_s, current_A, meas_temp_C, meas_voltage_V, initial_soc):
    heat_transfer, heat_capacity, resistance = scaled * (UPPER_BOUNDS - LOWER_BOUNDS) + LOWER_BOUNDS
    pred_temp, pred_voltage = simulate_cell(time_s, current_A, initial_soc, meas_temp_C[0],
                                            heat_transfer, heat_capacity, resistance)
    return np.concatenate(((pred_temp - meas_temp_C) / TEMPERATURE_NOISE_C,
                           (pred_voltage - meas_voltage_V) / VOLTAGE_NOISE_V))


def calibrate(time_s, current_A, meas_temp, meas_volt, initial_soc, initial_guess=(10.0, 800.0, 0.015)):
    guess_scaled = (np.asarray(initial_guess) - LOWER_BOUNDS) / (UPPER_BOUNDS - LOWER_BOUNDS)
    result = least_squares(residuals, guess_scaled, bounds=(0.0, 1.0),
                           args=(time_s, current_A, meas_temp, meas_volt, initial_soc), verbose=1)
    params = result.x * (UPPER_BOUNDS - LOWER_BOUNDS) + LOWER_BOUNDS
    rms_temp = np.sqrt(np.mean((result.fun[:len(time_s)] * TEMPERATURE_NOISE_C) ** 2))
    rms_volt = np.sqrt(np.mean((result.fun[len(time_s):] * VOLTAGE_NOISE_V) ** 2))
    return params, rms_temp, rms_volt


def main():
    true_params = None
    if len(sys.argv) == 2:
        data = np.genfromtxt(sys.argv[1], delimiter=",", names=True)
        time_s = data["time_s"]
        current_A = data["current_A"]
        meas_temp = data["temperature_C"]
        meas_volt = data["voltage_V"]
        initial_soc = float(data["initial_soc"][0])
    elif len(sys.argv) == 1:
        print("No CSV given: running synthetic self-test.")
        rng = np.random.default_rng(0)
        time_s = np.linspace(0.0, 3600.0, 361)        # 1 h, 10 s steps
        current_A = np.where(time_s < 2400.0, 5.0, 0.0)   # 1C charge then rest (makes R and h identifiable)
        initial_soc = 0.1
        true_params = (15.0, 900.0, 0.015)
        true_temp, true_volt = simulate_cell(time_s, current_A, initial_soc, EV_AMBIENT_TEMPERATURE_C, *true_params)
        meas_temp = true_temp + rng.normal(0.0, TEMPERATURE_NOISE_C, time_s.size)
        meas_volt = true_volt + rng.normal(0.0, VOLTAGE_NOISE_V, time_s.size)
    else:
        print(__doc__)
        return 2

    if np.any(np.diff(time_s) <= 0.0):
        print("time_s must be strictly increasing")
        return 2

    params, rms_temp, rms_volt = calibrate(time_s, current_A, meas_temp, meas_volt, initial_soc)
    print("\nCalibration Complete!")
    print(f"Heat transfer coefficient (W/m2/K): {params[0]:.3f}")
    print(f"Specific heat capacity (J/kg/K):    {params[1]:.1f}")
    print(f"Resistance (ohm):                   {params[2]:.4f}")
    print(f"RMS residual: {rms_temp:.3f} degC, {rms_volt * 1000.0:.1f} mV")

    if true_params is not None:
        rel_err = np.abs(params - np.array(true_params)) / np.array(true_params)
        ok = bool(np.all(rel_err < 0.05))
        print("CALIBRATION SELF-TEST: " + ("PASS" if ok else "FAIL") + f" (max rel. error {rel_err.max():.2%})")
        return 0 if ok else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
