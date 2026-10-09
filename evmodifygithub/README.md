# EV MODIFY

**Copyright (c) 2026. All Rights Reserved.**
**License**: PROPRIETARY (see `LICENSE` file for details). You may view the code, but you must ask for explicit permission to use or modify it.

---

EV MODIFY is a BMS firmware simulator for conservative, physics-informed charge-current limiting in lithium-ion battery modules.

Each control cycle it updates a per-cell SOC estimate (Extended Kalman Filter with Joseph-form covariance update, OCV lookup tables from PyBaMM Chen2020) together with calendar/cycling ageing, then searches for the largest charge current that keeps every series cell inside three limits over a 30 s horizon: anode potential >= 5 mV vs Li (plating margin), SOC <= 98 %, temperature <= 45 °C. All inputs are checked for finiteness and plausibility, and an invalid or out-of-range cell fails closed (FAULT, 0 A). The controller recovers automatically once the measurements are valid again.

Sign convention: positive current = charging.

## What is included

- `spme_v7_apex_solver.c/.h`: the C11 controller. It holds persistent per-cell state, validates inputs, runs the EKF and ageing update, and evaluates a reduced electrochemical predictor per cell (lumped thermal model, solid-diffusion surface-SOC correction, Butler-Volmer kinetics, electrolyte ohmic and concentration losses). A coarse upward sweep plus 24-step bisection finds the safe current. Electrode and electrolyte parameters follow the Chen2020 (LG M50) set; activation energies and ageing rates are calibration placeholders.
- `test_solver.c`: demo simulator with a simple plant model, plus `--self-test`. The self-test covers closed-loop limiting, physical trends (current/temperature), full-cell blocking, over-temperature fault and recovery, and invalid inputs (NaN, NULL, bad config).
- `ocv_proxy.c`: stdin/stdout bridge to the C OCV model, used by the Python tools.
- `ev_model.py`: shared Python helpers. It parses the OCV tables directly from the C source, so there is a single source of truth.
- `test_ocv.py`: checks the C OCV tables against the PyBaMM Chen2020 OCP functions, and the built proxy against Python.
- `run_pybamm_sim.py`: PyBaMM SPMe 1C **charge** from 10 % SOC to 4.2 V (PyBaMM current sign converted to charge-positive). It compares the voltage against the C OCV + I·R surrogate with a least-squares effective resistance.
- `calibrate_cell.py`: fits the heat-transfer coefficient, specific heat capacity and resistance of the C model to measured data (or a synthetic self-test that must recover known parameters).
- `diagrams/`: editable SVG architecture, model-scope and validation-roadmap figures.

## Build and verify

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure     # C self-test
./build/ev_modify_solver                       # demo run (CSV to stdout)
```

On Windows with MinGW you can run `compile.bat` instead.

```bash
pip install -r requirements.txt
python test_ocv.py          # OCV tables vs Chen2020 and C proxy
python calibrate_cell.py    # synthetic calibration self-test
python run_pybamm_sim.py    # writes pybamm_charge_benchmark.png beside the script
```

The PyBaMM benchmark compares terminal voltage only. The I·R surrogate lumps every SPMe overpotential into one fitted resistance (about 51 mΩ for this run, about 35 mV RMS error). It is a sanity check of the OCV/SOC model, not a validation of the plating predictor.

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
