"""Fits the thermal and 1-RC electrical parameters of EvCellConfig to measured cell data.

Usage:
    python tools/calibrate_cell.py DATA.csv [--capacity AH] [--ambient DEGC]
    python tools/calibrate_cell.py --self-test      # synthetic data; known parameters must be recovered

CSV columns (header row required): time_s, current_A (positive = charge; row i holds the current
over (time_s[i-1], time_s[i]], as in ev_step), temperature_C, voltage_V, initial_soc (first row). The profile should contain current steps
and rests so that R0, R1 and tau are separable, and the cell should start in thermal equilibrium.
Fitted: heat_transfer_W_m2_K, heat_capacity_J_kg_K, resistance_ohm (R0), rc_resistance_ohm (R1),
rc_time_constant_s. Exit code 0 = good fit, 1 = poor fit (residuals or parameter at a bound),
2 = bad input.
"""
import argparse
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

import ev_model as m

EV_CELL_MASS_KG = m.header_constant("EV_CELL_MASS_KG")
EV_CELL_COOLING_AREA_M2 = m.header_constant("EV_CELL_COOLING_AREA_M2")
EV_RESISTANCE_ACTIVATION_ENERGY_J_MOL = m.header_constant("EV_RESISTANCE_ACTIVATION_ENERGY_J_MOL")
GAS_CONSTANT = m.header_constant("EV_GAS_J_PER_MOL_K")

TEMPERATURE_NOISE_C = 0.2   # 1-sigma sensor noise used to weight residuals
VOLTAGE_NOISE_V = 0.005
NAMES = ["heat_transfer_W_m2_K", "heat_capacity_J_kg_K", "resistance_ohm", "rc_resistance_ohm", "rc_time_constant_s"]
LOWER = np.array([1.0, 100.0, 0.0, 0.0, 1.0])
UPPER = np.array([200.0, 5000.0, 0.5, 0.5, 5000.0])
GUESS = np.array([10.0, 900.0, 0.02, 0.02, 100.0])
REQUIRED_COLUMNS = ("time_s", "current_A", "temperature_C", "voltage_V", "initial_soc")


def resistance_scale(temperature_C):
    """R(T) / R(25 degC), as in the C model."""
    return np.exp(EV_RESISTANCE_ACTIVATION_ENERGY_J_MOL / GAS_CONSTANT
                  * (1.0 / (np.asarray(temperature_C) + 273.15) - 1.0 / 298.15))


def simulate_cell(params, time_s, current_A, initial_soc, initial_temp_C, ambient_C, capacity_Ah):
    """Same discrete equations as the C model (current_A[i] held over (t[i-1], t[i]]): coulomb
    counting, exact lumped-thermal and 1-RC updates, V = OCV(soc) + v_rc + R0(T)*I with R0 and R1
    referred to 25 degC. R0 uses the previous sample's temperature (the C EKF updates the state
    temperature after the estimator)."""
    h, cp, r0, r1, tau = params
    conductance_W_K = h * EV_CELL_COOLING_AREA_M2
    heat_capacity_J_K = EV_CELL_MASS_KG * cp
    n = time_s.size
    soc = np.empty(n)
    temp = np.empty(n)
    v_rc = np.empty(n)
    soc[0], temp[0], v_rc[0] = initial_soc, initial_temp_C, 0.0
    for i in range(1, n):
        dt = time_s[i] - time_s[i - 1]
        current = current_A[i]
        scale = resistance_scale(temp[i - 1])
        soc[i] = np.clip(soc[i - 1] + current * dt / (capacity_Ah * 3600.0), 0.0, 1.0)
        a = np.exp(-dt / tau)
        v_rc[i] = a * v_rc[i - 1] + r1 * scale * (1.0 - a) * current
        rise = current * current * (r0 + r1) * scale / conductance_W_K
        decay = np.exp(-conductance_W_K * dt / heat_capacity_J_K)
        temp[i] = ambient_C + rise + (temp[i - 1] - ambient_C - rise) * decay
    voltage = m.full_cell_ocv_V(soc) + v_rc + r0 * resistance_scale(np.r_[temp[0], temp[:-1]]) * current_A
    return temp, voltage


def integrate_cell_ode(params, time_s, current_A, initial_soc, initial_temp_C, ambient_C, capacity_Ah):
    """Independent reference for the self-test: the continuous model integrated with an adaptive
    ODE solver, interval by interval, instead of the exact discrete update above."""
    h, cp, r0, r1, tau = params
    conductance_W_K = h * EV_CELL_COOLING_AREA_M2
    heat_capacity_J_K = EV_CELL_MASS_KG * cp
    state = np.array([initial_soc, 0.0, initial_temp_C])
    soc, temp, v_rc = [initial_soc], [initial_temp_C], [0.0]
    for i in range(1, time_s.size):
        current = current_A[i]

        def rhs(_t, x):
            scale = float(resistance_scale(x[2]))
            return [current / (capacity_Ah * 3600.0),
                    (r1 * scale * current - x[1]) / tau,
                    (current * current * (r0 + r1) * scale - conductance_W_K * (x[2] - ambient_C)) / heat_capacity_J_K]

        state = solve_ivp(rhs, (time_s[i - 1], time_s[i]), state, rtol=1e-9, atol=1e-12).y[:, -1]
        soc.append(state[0])
        v_rc.append(state[1])
        temp.append(state[2])
    temp = np.array(temp)
    voltage = m.full_cell_ocv_V(np.array(soc)) + np.array(v_rc) + r0 * resistance_scale(np.r_[temp[0], temp[:-1]]) * current_A
    return temp, voltage


def calibrate(data, ambient_C, capacity_Ah):
    time_s, current_A, meas_temp, meas_volt, initial_soc = data

    def residuals(scaled):
        params = LOWER + scaled * (UPPER - LOWER)
        temp, volt = simulate_cell(params, time_s, current_A, initial_soc, meas_temp[0], ambient_C, capacity_Ah)
        return np.concatenate(((temp - meas_temp) / TEMPERATURE_NOISE_C, (volt - meas_volt) / VOLTAGE_NOISE_V))

    result = least_squares(residuals, (GUESS - LOWER) / (UPPER - LOWER), bounds=(0.0, 1.0))
    params = LOWER + result.x * (UPPER - LOWER)
    n = time_s.size
    rms_temp = float(np.sqrt(np.mean((result.fun[:n] * TEMPERATURE_NOISE_C) ** 2)))
    rms_volt = float(np.sqrt(np.mean((result.fun[n:] * VOLTAGE_NOISE_V) ** 2)))
    # 1-sigma parameter uncertainty from the Gauss-Newton covariance (weighted residuals)
    jac = result.jac * (UPPER - LOWER)[None, :] ** -1
    try:
        sigma = np.sqrt(np.diag(np.linalg.inv(jac.T @ jac)))
    except np.linalg.LinAlgError:
        sigma = np.full(params.size, np.inf)
    at_bound = (result.x < 1e-3) | (result.x > 1.0 - 1e-3)
    return params, sigma, rms_temp, rms_volt, at_bound


def load_csv(path):
    try:
        data = np.genfromtxt(path, delimiter=",", names=True, encoding="utf-8-sig", ndmin=1)
    except (OSError, ValueError) as error:
        raise ValueError(f"cannot read {path}: {error}") from error
    if data.dtype.names is None:
        raise ValueError("CSV needs a header row")
    missing = [c for c in REQUIRED_COLUMNS if c not in data.dtype.names]
    if missing:
        raise ValueError(f"missing columns: {', '.join(missing)}")
    if data.size < 10:
        raise ValueError("need at least 10 rows")
    columns = [np.asarray(data[c], dtype=float) for c in REQUIRED_COLUMNS[:4]]
    if not all(np.all(np.isfinite(c)) for c in columns):
        raise ValueError("non-finite values in the data columns")
    if np.any(np.diff(columns[0]) <= 0.0):
        raise ValueError("time_s must be strictly increasing")
    initial_soc = float(data["initial_soc"][0])
    if not 0.0 <= initial_soc <= 1.0:
        raise ValueError("initial_soc must be in [0, 1]")
    return (*columns, initial_soc)


def synthetic_data():
    rng = np.random.default_rng(0)
    time_s = np.arange(0.0, 7200.0 + 1.0, 5.0)
    current_A = np.where(((time_s - 1e-9) % 1200.0) < 600.0, 5.0, 0.0)   # 1C, 10 min on / 10 min off
    current_A[time_s > 6000.0] = 0.0
    current_A[0] = 0.0
    truth = np.array([15.0, 1000.0, 0.036, 0.015, 150.0])
    temp, volt = integrate_cell_ode(truth, time_s, current_A, 0.1, 25.0, 25.0, 5.0)
    temp = temp + rng.normal(0.0, TEMPERATURE_NOISE_C, time_s.size)
    volt = volt + rng.normal(0.0, VOLTAGE_NOISE_V, time_s.size)
    return (time_s, current_A, temp, volt, 0.1), truth


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csv", nargs="?", help="measured data")
    parser.add_argument("--self-test", action="store_true", help="fit synthetic data with known parameters")
    parser.add_argument("--capacity", type=float, default=5.0, help="capacity between the OCV limits [Ah]")
    parser.add_argument("--ambient", type=float, help="ambient/coolant temperature [degC] (default: first sample)")
    args = parser.parse_args()
    if args.self_test == (args.csv is not None):
        parser.error("give either a CSV file or --self-test")
    if not (np.isfinite(args.capacity) and args.capacity > 0.0):
        parser.error("--capacity must be finite and positive")
    if args.ambient is not None and not np.isfinite(args.ambient):
        parser.error("--ambient must be finite")

    truth = None
    if args.self_test:
        data, truth = synthetic_data()
    else:
        try:
            data = load_csv(args.csv)
        except ValueError as error:
            print(f"error: {error}")
            return 2
    ambient_C = data[2][0] if args.ambient is None else args.ambient

    params, sigma, rms_temp, rms_volt, at_bound = calibrate(data, ambient_C, args.capacity)
    print("Calibration result (value +/- 1 sigma):")
    for name, value, s, bound in zip(NAMES, params, sigma, at_bound):
        print(f"  {name:22s} {value:12.5g} +/- {s:.2g}{'   AT BOUND' if bound else ''}")
    print(f"RMS residual: {rms_temp:.3f} degC, {rms_volt * 1000.0:.2f} mV")

    good_fit = (rms_temp <= 2.0 * TEMPERATURE_NOISE_C) and (rms_volt <= 2.0 * VOLTAGE_NOISE_V) and not np.any(at_bound)
    if not good_fit:
        print("WARNING: poor fit (residuals above 2x sensor noise or a parameter at its bound)")
    if truth is not None:
        rel_err = np.abs(params - truth) / truth
        recovered = bool(np.all(rel_err < 0.05))
        print(f"CALIBRATION SELF-TEST: {'PASS' if recovered and good_fit else 'FAIL'} (max rel. error {rel_err.max():.2%})")
        return 0 if recovered and good_fit else 1
    return 0 if good_fit else 1


if __name__ == "__main__":
    sys.exit(main())
