import numpy as np
from scipy.optimize import least_squares
import csv
import sys
import os

# Physical constants from the C model
EV_CELL_MASS_KG = 0.070
EV_CELL_COOLING_AREA_M2 = 0.0075
EV_AMBIENT_TEMPERATURE_C = 25.0
EV_FARADAY_C_PER_MOL = 96485.33212
EV_GAS_J_PER_MOL_K = 8.31446262
EV_REFERENCE_TEMPERATURE_K = 298.15
EV_KINETIC_ACTIVATION_ENERGY_J_MOL = 30000.0
EV_EXCHANGE_CURRENT_REF_A_M2 = 1.5
EV_PARTICLE_RADIUS_M = 5.8e-6
EV_GEOMETRIC_ELECTRODE_AREA_M2 = 0.100
EV_ANODE_THICKNESS_M = 75.0e-6

def ev_clamp(val, min_val, max_val):
    return max(min_val, min(max_val, val))

def ev_graphite_ocp_V(theta):
    bounded_theta = ev_clamp(theta, 0.001, 0.999)
    return 0.0630 + 0.8000 * np.exp(-40.0 * bounded_theta) \
        + 0.0650 * np.exp(-0.20 * (bounded_theta - 0.35)) \
        - 0.0100 * np.log10(bounded_theta / (1.0 - bounded_theta))

def ev_nmc_ocp_V(soc):
    x = ev_clamp(soc, 0.001, 0.999)
    return 3.43 + 0.67 * x - 0.45 * x * x

def ev_full_cell_ocv_V(soc):
    return ev_nmc_ocp_V(soc) - ev_graphite_ocp_V(soc)

def simulate_cell(time_s, current_A, initial_soc, heat_transfer_W_m2_K, heat_capacity_J_kg_K, resistance_ohm):
    """
    Simulates the cell temperature and voltage given a current profile.
    Returns predicted temperature and voltage arrays.
    """
    dt = time_s[1] - time_s[0]
    n_steps = len(time_s)
    
    pred_temp = np.zeros(n_steps)
    pred_voltage = np.zeros(n_steps)
    
    current_temp_C = EV_AMBIENT_TEMPERATURE_C
    current_soc = initial_soc
    
    thermal_conductance_W_K = heat_transfer_W_m2_K * EV_CELL_COOLING_AREA_M2
    heat_capacity_J_K = EV_CELL_MASS_KG * heat_capacity_J_kg_K
    
    for i in range(n_steps):
        I = current_A[i]
        
        # 1. Update SOC
        current_soc += (I * dt) / (5.0 * 3600.0) # Nominal capacity 5Ah
        current_soc = ev_clamp(current_soc, 0.0, 1.0)
        
        # 2. Update Temperature
        joule_heat_W = I * I * resistance_ohm
        temp_rise_steady = joule_heat_W / thermal_conductance_W_K
        decay = np.exp(-(thermal_conductance_W_K * dt) / heat_capacity_J_K)
        current_temp_C = EV_AMBIENT_TEMPERATURE_C + temp_rise_steady + (current_temp_C - EV_AMBIENT_TEMPERATURE_C - temp_rise_steady) * decay
        
        pred_temp[i] = current_temp_C
        
        # 3. Predict Voltage (Simplified for calibration: OCV + Ohmic)
        # Note: A full DFN surrogate model includes polarization, but 
        # OCV + Ohmic dominates the macro prediction.
        pred_voltage[i] = ev_full_cell_ocv_V(current_soc) + I * resistance_ohm
        
    return pred_temp, pred_voltage

def objective_function(params, time_s, current_A, meas_temp_C, meas_voltage_V, initial_soc):
    heat_transfer = params[0]
    heat_capacity = params[1]
    resistance_ohm = params[2]
    
    pred_temp, pred_voltage = simulate_cell(
        time_s, current_A, initial_soc, 
        heat_transfer, heat_capacity, resistance_ohm
    )
    
    # We want to minimize the error in both temperature and voltage
    temp_error = (pred_temp - meas_temp_C) / np.std(meas_temp_C)
    volt_error = (pred_voltage - meas_voltage_V) / np.std(meas_voltage_V)
    
    return np.concatenate((temp_error, volt_error))

def main():
    if len(sys.argv) != 2:
        print("Usage: python calibrate_cell.py <test_data.csv>")
        print("Using synthetic data for self-test...")
        
        # Generate synthetic data
        time_s = np.linspace(0, 3600, 360) # 1 hour, 10s steps
        current_A = np.full_like(time_s, 5.0) # 1C charge
        
        true_ht = 15.0
        true_hc = 900.0
        initial_soc = 0.1
        res_ohm = 0.015
        
        true_temp, true_volt = simulate_cell(time_s, current_A, initial_soc, true_ht, true_hc, res_ohm)
        # Add noise
        meas_temp = true_temp + np.random.normal(0, 0.2, len(true_temp))
        meas_volt = true_volt + np.random.normal(0, 0.01, len(true_volt))
    else:
        # Load CSV
        filename = sys.argv[1]
        data = np.genfromtxt(filename, delimiter=',', names=True)
        time_s = data['time_s']
        current_A = data['current_A']
        meas_temp = data['temperature_C']
        meas_volt = data['voltage_V']
        initial_soc = data['initial_soc'][0]
        res_ohm = data['resistance_ohm'][0]

    lower_bounds = np.array([1.0, 100.0, 0.001])
    upper_bounds = np.array([50.0, 2000.0, 0.100])
    
    def scale_params(p):
        return (p - lower_bounds) / (upper_bounds - lower_bounds)
        
    def unscale_params(p_scaled):
        return p_scaled * (upper_bounds - lower_bounds) + lower_bounds

    initial_guess_scaled = scale_params(np.array([10.0, 800.0, 0.015]))
    
    def scaled_objective(p_scaled, *args):
        p_real = unscale_params(p_scaled)
        return objective_function(p_real, *args)
        
    bounds_scaled = ([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
    
    print("Starting optimization (with scaled parameters)...")
    res_scaled = least_squares(
        scaled_objective, 
        initial_guess_scaled, 
        bounds=bounds_scaled,
        args=(time_s, current_A, meas_temp, meas_volt, initial_soc),
        verbose=1
    )
    
    optimal_params = unscale_params(res_scaled.x)
    
    print("\nCalibration Complete!")
    print(f"Optimal Heat Transfer (W/m2/K): {optimal_params[0]:.3f}")
    print(f"Optimal Heat Capacity (J/kg/K): {optimal_params[1]:.1f}")
    print(f"Optimal Resistance (ohm): {optimal_params[2]:.4f}")
    print(f"Cost (MSE): {res_scaled.cost:.4f}")

if __name__ == '__main__':
    main()
