import pybamm
import numpy as np
import matplotlib.pyplot as plt
import subprocess

def main():
    model = pybamm.lithium_ion.SPMe()
    parameter_values = pybamm.ParameterValues("Chen2020")
    
    experiment = pybamm.Experiment(
        ["Charge at 1 C until 4.2 V"]
    )
    
    sim = pybamm.Simulation(model, experiment=experiment, parameter_values=parameter_values)
    sol = sim.solve(initial_soc=0.1)
    
    pybamm_time = sol["Time [s]"].entries
    pybamm_voltage = sol["Terminal voltage [V]"].entries
    pybamm_current = sol["Current [A]"].entries
    
    c_sim_voltage = []
    soc = 0.1
    
    # We pipe all SOCs to ocv_proxy to run the C physics model efficiently
    # Batch processing to avoid calling subprocess repeatedly
    soc_list = []
    for i in range(len(pybamm_time)):
        dt = pybamm_time[i] - pybamm_time[i-1] if i > 0 else 0
        I = pybamm_current[i]
        soc += (I * dt) / (5.0 * 3600.0)
        soc = max(0.0, min(1.0, soc))
        soc_list.append(soc)
        
    input_str = "\n".join(f"{s}" for s in soc_list) + "\n"
    result = subprocess.run(["ocv_proxy.exe"], input=input_str, text=True, capture_output=True)
    c_ocvs = [float(v) for v in result.stdout.strip().split("\n")]
    
    for i in range(len(pybamm_time)):
        I = pybamm_current[i]
        v = c_ocvs[i] + I * 0.015
        c_sim_voltage.append(v)
        
    c_sim_voltage = np.array(c_sim_voltage)
    rms_error = np.sqrt(np.mean((c_sim_voltage - pybamm_voltage)**2))
    
    plt.plot(pybamm_time, pybamm_voltage, label="PyBaMM SPMe")
    plt.plot(pybamm_time, c_sim_voltage, label="C Model Surrogate")
    plt.xlabel("Time [s]")
    plt.ylabel("Voltage [V]")
    plt.legend()
    plt.title(f"Voltage Comparison (RMS Error: {rms_error:.4f} V)")
    plt.savefig("pybamm_charge_benchmark.png")
    print(f"RMS Error against PyBaMM: {rms_error:.4f} V")
    print("Generated pybamm_charge_benchmark.png")
        
if __name__ == "__main__":
    main()
