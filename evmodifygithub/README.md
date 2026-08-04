# EV MODIFY

**Copyright (c) 2026. All Rights Reserved.**
**License**: PROPRIETARY (see `LICENSE` file for details). You may view the code, but you must ask for explicit permission to use or modify it.

---

EV MODIFY is a BMS firmware simulator for conservative, physics-informed charge-current limiting in lithium-ion battery modules.

It features a numerically stable Joseph Form Extended Kalman Filter (EKF), robust calendric aging degradation tracking, and entirely bounds-checked floating-point safety without reliance on NaNs or heavy transcendental math during real-time loops.

## What is included

- `spme_v7_apex_solver.c` — deterministic C charge-limit demonstrator with
  persistent SOC/temperature states, input validation, a pack-current filter,
  per-cell prediction, and bisection-based current limiting.
- `run_pybamm_sim.py` — correctly signed PyBaMM **charge** comparison from
  10% SOC to 4.2 V. It acts as a reference for voltage/temperature.
- `calibrate_cell.py` — calibration script to fit model parameters (heat transfer, heat capacity, resistance) to empirical data.
- `diagrams/` — editable SVG architecture, model-scope, and validation-roadmap
  figures used by the executive report.

## Reproducible verification

Install the pinned Python dependencies, then run:

```powershell
python .\run_pybamm_sim.py
```

The PyBaMM script confirms the run ends at the charge-voltage limit and writes `pybamm_charge_benchmark.png` beside itself. RMS error between the PyBaMM output and the C-model surrogate is output to the console and displayed on the plot.

## Model boundary

The reduced model has explicit parameter units, state persistence, and input
checks, but its electrode parameters and safety threshold require calibration
for the specific cell, ageing state, sensors, cooling system, and operating
envelope. Its use on hardware requires at minimum:

1. Cell characterization across SOC, temperature, current, and ageing.
2. Comparison with a validated DFN/SPMe or measured reference data.
3. Hardware-in-the-loop, fault-injection, numerical-range, and timing tests on
   the target microcontroller.
4. Requirements, traceability, verification, and an ISO 26262 safety case if
   the controller will affect a vehicle.

## Commercially credible next step

The credible proposition is a calibrated
software feature that uses existing current/temperature/voltage signals to
reduce fast-charge degradation and warranty cost. The first commercial proof
point should be a measured reduction in plating/degradation compared with a
baseline charger at equal charging time and safety constraints.
