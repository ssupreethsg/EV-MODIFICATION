# EV MODIFY

**Copyright (c) 2026 Supreeth S G. All rights reserved.**
**License:** proprietary, see [`LICENSE`](LICENSE). The code may be viewed. Building, running, modifying or any other use requires the author's prior written permission. The build instructions below are for permitted users.

EV MODIFY is a C11 battery-management function that limits the charge current of a series lithium-ion module. The limit is conservative and based on physics. Every control cycle it:

1. validates all inputs and the per-cell configuration;
2. updates each cell's state:
   - a 2-state Extended Kalman Filter tracks SOC and the 1-RC polarisation voltage, with a Joseph-form update, a current-dependent measurement noise and a 5σ innovation gate;
   - a voltage-tracking check that catches a sensor frozen at an in-range value, with or without noise. Over a window, the change of each reading is compared with the change of an open-loop model (coulomb counting, 1-RC and R0, started from the EKF state). Allowances cover 30 % R0, 100 % R1 (applied to the current steps in the window and to recent steps whose RC transients are still decaying) and 20 % OCV-change errors, and the EKF's 3σ SOC band. Above 10 mV of mismatch the EKF stops using the reading. Above 20 mV for 3 s (and at least 2 cycles) the cell faults, and the fault stays latched until `ev_clear_latched_faults()`, so the charger stop it causes cannot clear it. Short glitches (for example 0.3 s at a 0.1 s period) do not latch. A reading that stays suspect, never validating, while 10 % of the capacity is charged also latches, so pauses caused by other faults cannot hide it;
   - an SOC-divergence backstop: since the reading was last shown to track the model, the EKF may fall below the coulomb count only by its 3σ uncertainty, 2 % SOC and a 20 % capacity error on the charge counted. A reading that drifts, steps or freezes pulls the EKF further and latches `EV_FAULT_SOC_DIVERGENCE`;
   - optionally, `ev_step_with_pack_voltage()` compares the sum of the cell readings with an independent string-voltage measurement (see Sensor limits below);
   - temperature, capacity fade (SOH) and resistance growth are updated alongside it (open-loop ageing models, not estimated from data);
3. detects faults, with hysteresis and a 5-cycle hold. This includes a check that the measured current follows the previous approved current, allowing the charger 2 s to ramp down;
4. finds the largest charge current up to the request that keeps **every** cell inside all constraints for max(horizon, 10 s maximum control period + the 2 s charger ramp-down). The search is a 20-step upward sweep followed by 24 bisection steps.

| Constraint (per cell) | Limit | Flag |
|---|---|---|
| Graphite potential vs Li/Li+ at the end of the horizon (plating margin) | ≥ 5 mV | `EV_CONSTRAINT_ANODE_POTENTIAL` |
| Terminal voltage, bounded over the whole horizon (CV taper) | ≤ 4.2 V | `EV_CONSTRAINT_CELL_VOLTAGE` |
| Upper SOC bound: the larger of EKF SOC + 3σ and the coulomb count since the reading was last validated (less its 3σ and a 5 % capacity error on the charge counted), plus the charge over the horizon. | ≤ 98 % | `EV_CONSTRAINT_SOC` |
| Predicted cell temperature + 2 °C sensor margin | ≤ 45 °C | `EV_CONSTRAINT_TEMPERATURE_HIGH` |
| Cell temperature − 2 °C sensor margin while charging (3 °C release hysteresis) | ≥ 0 °C | `EV_CONSTRAINT_TEMPERATURE_LOW` |
| Electrode-averaged electrolyte depletion (keeps c_e ≥ 10 % at the collector) | < 60 % | `EV_CONSTRAINT_ELECTROLYTE_DEPLETION` |

The decision is **ALLOW** (full request), **LIMIT** (a constraint binds; the binding cell and constraint bits are reported) or **FAULT** (0 A, with `EV_FAULT_*` bits). FAULT is returned for:

| Fault | Condition |
|---|---|
| `EV_FAULT_INPUT` | a NULL pointer, or a pack-level input that is non-finite or out of range: current beyond ±1000 A, dt outside [0.001, 10] s, horizon outside (0, 600] s, a non-finite request, or (with `ev_step_with_pack_voltage()`) a pack tolerance or path resistance out of range, or a pack reading that went missing while current flowed (until it is valid at rest). Finite requests above 1000 A are clamped instead. |
| `EV_FAULT_CONFIG` | an `EvCellConfig` field is out of range. |
| `EV_FAULT_VOLTAGE_SENSOR`, `EV_FAULT_TEMPERATURE_SENSOR` | a reading is outside its plausibility window. An intermittent voltage reading (3 invalid samples within about 20 s) is latched until `ev_clear_latched_faults()`. |
| `EV_FAULT_OVER_VOLTAGE` | measured cell voltage > 4.25 V; releases at 4.20 V. |
| `EV_FAULT_VOLTAGE_PLAUSIBILITY` | 10 consecutive EKF innovations beyond 5σ, or a reading that does not follow the model by more than 20 mV beyond the allowances for 3 s (a frozen or stuck sensor). The tracking fault is latched until `ev_clear_latched_faults()`. |
| `EV_FAULT_SOC_DIVERGENCE` | the EKF SOC fell below the coulomb count by more than the backstop allows (a drifting, offset or frozen reading). Latched until `ev_clear_latched_faults()`. |
| `EV_FAULT_PACK_VOLTAGE` | only with `ev_step_with_pack_voltage()`: the sum of the cell readings differs from the pack reading by more than the given tolerance for 3 s. Latched until `ev_clear_latched_faults()`. |
| `EV_FAULT_OVER_CURRENT` | measured charge current above the allowance for 1 s: the previous approved current, plus a linear 2 s ramp-down after a lower command, plus 5 % + C/50 (charger not following the limit); or charge beyond 1.05 × the expected current + C/100, integrated with a 1 h leak, above 0.5 % of capacity (a sustained excess inside the fast tolerance, such as a charger floor); or, at a 0 A command after the ramp-down, a reading above 4.21 V, or an EKF SOC or coulomb-count lower bound above 99 % (a charger that does not stop). |
| `EV_FAULT_OVER_TEMPERATURE`, `EV_FAULT_UNDER_TEMPERATURE` | cell above 45 °C or below -30 °C; releases 3 °C inside the limit. |
| `EV_FAULT_STATE` | a corrupted estimator state was reset (SOC σ = 0.5): a numeric field out of range or a flag byte other than 0 or 1 (read as set). A state with another `layout_version` faults until `ev_cell_init()`. |
| `EV_FAULT_END_OF_LIFE` | SOH has reached 50 %. |
| `EV_FAULT_HOLD` | not a fault itself: a fault has cleared and the 5-cycle hold is running. |

A prediction that is not finite sets `EV_CONSTRAINT_INVALID_INPUT` and is never treated as safe.

Every fault that is not latched, including pack-level input faults, is held until 5 consecutive good cycles have passed. Latched faults need `ev_clear_latched_faults()`, a service action after the sensor is repaired; it keeps SOH and resistance growth and resets the SOC σ to 0.5, so the EKF re-learns the SOC from the repaired reading. Persist the whole `EvCellState` across power cycles and call `ev_cell_init()` only for a new cell, so latches and ageing survive a reboot. `ev_step()` checks a restored state's `layout_version`, its numeric ranges and its flag bytes, but a flag flipped between valid values cannot be detected, so protect the stored copy with a CRC. When the current or dt measurement is invalid, or a cell's configuration is invalid, the SOC mean is kept and its variance is widened, so uncounted charge can never make the limit optimistic. Positive current means charging.

## Layout

```
src/ev_modify.{h,c}         controller library (public API documented in the header)
tests/test_solver.c         C self-test (208 checks) and demo
tests/test_ocv.py           OCV table / stoichiometry checks
tests/cli/                  ev_cli exit-code and output tests (run by CTest through cmake -P)
tests/consumer/             find_package / add_subdirectory consumer used by CI
tests/cm4/                  link check for the Cortex-M4F CI job
tools/ev_cli.c              stdin/stdout bridge to the C model (ocv, limit, ekf and charge modes)
tools/ev_model.py           shared Python helpers; parse constants and tables from the C sources
tools/generate_ocv_lut.py   regenerates (or --check) the OCV tables in the C source
tools/calibrate_cell.py     fits thermal + 1-RC parameters to measured data
tools/run_pybamm_sim.py     validation against PyBaMM SPMe and DFN (Chen2020)
docs/                       design notes, diagrams and the benchmark plot
.github/workflows/ci.yml    GCC, Clang, ASan/UBSan, MSVC, Cortex-M4F, cppcheck, coverage, ruff, PyBaMM
```

## Build and test

Requires CMake ≥ 3.21 and a C11 compiler. The Python checks need Python 3 (validated on 3.13) with `numpy` and `scipy`.

```bash
pip install -r requirements.txt
cmake -S . -B build -DEV_WARNINGS_AS_ERRORS=ON
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
./build/test_solver            # demo: 4 cells, 30 A request, CSV to stdout (build/Release/test_solver.exe with MSVC)
```

CTest runs the C self-test, the demo, 23 `ev_cli` exit-code and output tests (22 where `/dev/full` is missing), `tests/test_ocv.py`, the OCV-table freshness check and the calibration self-test. If numpy or scipy is missing, the Python checks are not registered and CMake prints a warning. Options:

- `-DEV_SANITIZE=ON`: AddressSanitizer and UndefinedBehaviorSanitizer (GCC/Clang). For testing only, so it turns the install rules off.
- `-DEV_WARNINGS_AS_ERRORS=ON`: `-Werror` / `/WX`. Off by default; CI turns it on.
- `-DEV_BUILD_TESTS=OFF`: library only. This is the default when the project is used through `add_subdirectory`.
- `-DEV_PYTHON_TESTS=OFF`: skip the Python checks.
- `-DEV_INSTALL=OFF`: no install rules. This is the default under `add_subdirectory`.

The library refuses to compile with `-ffast-math`, `-ffinite-math-only` or `/fp:fast`, because its NaN checks are safety relevant. It is also cross-compiled in CI for a Cortex-M4F (about 19 KB of code and constants, excluding the soft-double and libm runtime; CI fails above 22 000 bytes and links a small program against newlib). The EKF and SOC accumulation use double precision on purpose, which is software-emulated on single-precision FPUs.

To use it from another CMake project, either embed it:

```cmake
add_subdirectory(EV-MODIFICATION)
target_link_libraries(app PRIVATE ev_modify::ev_modify)
```

or install it (`cmake --install build --prefix <dir>`) and use `find_package(ev_modify 2.0 REQUIRED)`.

## Using the API

```c
#include "ev_modify.h"

EvCellConfig config[EV_NUM_CELLS];   /* capacity, R0, R1, tau (25 degC), thermal, coolant, maker's max charge current */
EvCellState  state[EV_NUM_CELLS];
for (unsigned i = 0; i < EV_NUM_CELLS; ++i) {
    config[i] = (EvCellConfig){5.0f, 0.041f, 0.0135f, 190.0f, 10.0f, 1100.0f, 25.0f, 3.5f};
}
if (!nvm_restore(state)) {           /* your CRC-checked copy; ev_step() also checks state[i].layout_version */
    for (unsigned i = 0; i < EV_NUM_CELLS; ++i) {
        ev_cell_init(&state[i], soc_from_rested_ocv[i], temperature_C[i]);   /* new cell only */
    }
}
/* every control cycle (1 ms <= dt <= 10 s); current = average over the cycle that just ended: */
EvPackResult r = ev_step(config, state, measured_current_A, requested_current_A,
                         cell_voltage_V, cell_temperature_C, dt_s, 30.0f);
set_charger_current(r.approved_charge_current_A);
nvm_save(state);                     /* periodically and at shutdown */
```

With an independent measurement of the string voltage (a second AFE channel or the charger's output voltage, sampled in the same cycle), call `ev_step_with_pack_voltage(..., pack_voltage_V, tolerance_V, path_resistance_ohm)` instead: the tolerance in (0, 0.5] V from both measurement accuracies, and the resistance in [0, 0.1] Ω of the busbars between the cell taps and the pack sense point. The resistance is trusted; a value set too high hides a cell reading low by up to the current times the excess.

`ev_cell_max_charge_current()` gives the same limit for a single cell. Every struct field's unit, range and role is documented in `src/ev_modify.h`. The library is single-threaded, allocates no memory and keeps no global state.

## Validation

`tools/run_pybamm_sim.py` (needs `requirements-benchmark.txt`) checks the C code against PyBaMM's Chen2020 LG M50 cell. The CI job runs it on every push to `main`, on pull requests and on manual runs.

| Check | Result | Limit |
|---|---|---|
| C OCV table vs PyBaMM bulk OCV (SPMe, 1C pulse profile) | 0.2 mV max | 1 mV |
| 1-RC model fitted on 0.5C pulses (R0 41.3 mΩ, R1 13.5 mΩ, τ 193 s), checked on 1C pulses | 17.4 mV RMS | 25 mV |
| C EKF started 0.20 SOC off; fitted / R0 +20 % / R0 −20 % / capacity −5 % | 1.1 / 1.3 / 1.7 / 1.6 % max SOC error from 15 min on (end of the first rest) | 3 % |
| Same four cases: samples outside the EKF's own 3σ band | 0 % | 1 % |
| Same four cases: samples with any fault flag (a healthy cell must never trip the plausibility checks) | 0 | 0 |
| Closed loop: the **complete controller** (`ev_step`: EKF on the DFN terminal voltage, gradient memory, voltage tracking, SOC-divergence and over-current checks) charging a PyBaMM **DFN** (15 A = 2.9C request with no maker cap, to stress the electrochemical limits; 30 s horizon, from 5 % SOC, 2 h). At 2 and 10 °C the DFN's own activation energies are swept independently of the controller's: electrolyte 17 / 30 kJ/mol × graphite diffusion 35 / 60 kJ/mol at exchange-current 35 kJ/mol, plus the harshest corner at exchange-current 53.4 and 60 kJ/mol. Worst-case minimum anode potential vs Li/Li+, over every point of the electrode at the end of every 10 s step | +38.0 / +41.6 / +38.1 mV at 2 / 10 / 25 °C | ≥ 5 mV |
| Same sweep: maximum DFN terminal voltage | 3.807 / 3.945 / 4.172 V | ≤ 4.205 V |
| Cold constant-voltage taper: 2 °C from 75 % SOC, harshest corner | +41.0 mV, 4.113 V, 81.5 % SOC | ≥ 5 mV, ≤ 4.205 V |
| Sensor errors at their limits: cell at 0 °C while its sensor reads 2 °C, harshest corner, from 5 % SOC | +28.9 mV, 41.4 % SOC, no fault | ≥ 5 mV, ≥ 5 % SOC charged |
| Same sensors, from 50 % SOC, charger delivering 1.05 × command + 0.99 × C/50 + C/100 while the current sensor reads C/100 low (just inside the fast tolerance as measured) | +28.8 mV; only the over-charge integral may fault, not before 1500 s | ≥ 5 mV, no other fault |
| Same, charger at 1.05 × command + 1.5 × C/50 (outside the tolerance) | over-current fault at the first 10 s step, +31.0 mV | fault ≤ 30 s, ≥ 5 mV |
| Same sensors, charger at 1.05 × command + 2.45 × C/100 (a sustained excess the over-charge integral never reaches), from 5 % and from 75 % SOC | +17.3 / +16.8 mV, no fault | ≥ 5 mV, 0 faults |
| O'Regan2022's graphite diffusivity and exchange current (a second parameterisation of the same cell, with their own temperature dependence; electrolyte 30 kJ/mol) at 2 °C, at 0 °C with the 2 °C sensor error, and at 25 °C | +98.8 / +95.1 / +58.0 mV | ≥ 5 mV |
| Landesfeind & Gasteiger's measured electrolyte (conductivity, diffusivity, transference number, thermodynamic factor) with the worst graphite pairing (O'Regan2022 diffusivity, Chen2020 exchange current at 60 kJ/mol when cold): 2 °C, 0 °C with the 2 °C sensor error, 2 °C from 75 % SOC, 10 °C, 25 °C | +33.7 / +24.7 / +38.5 / +37.5 / +6.6 mV | ≥ 5 mV |
| Lumped-thermal DFN at 25 °C, controller fed the DFN cell temperature | +38.1 mV, peak 34.3 °C | ≥ 5 mV, ≤ 45 °C |
| Nominal cases: SOC reached in 2 h, so a limiter that never charges cannot pass | 44.0 / 63.1 / 94.7 % | ≥ 42 / 60 / 90 % |
| Every other closed-loop case: control cycles with any fault flag | 0 | 0 |

![PyBaMM benchmark](docs/pybamm_charge_benchmark.png)

The DFN uses PyBaMM's ideal-solution electrolyte (thermodynamic factor 1), so the controller's measured-electrolyte diffusion-potential factor shows up as extra margin here. The DFN loop measures only what a BMS measures (current, terminal voltage, temperature) and starts from the true SOC, as after a rest. The C self-test (208 checks, ≥ 98 % line coverage of the controller) uses the controller's own model as its plant, including parameter-mismatch and fault-injection cases. PyBaMM is the independent reference: its electrode and electrolyte parameters are the published Chen2020 set, and its temperature dependences (electrolyte, solid diffusion, exchange current) are swept rather than copied from the controller, and its plating floor is the controller's own 5 mV margin rather than 0 V. Neither replaces experimental plating validation.

## Model and parameter provenance

- **Electrochemistry (anode).** PyBaMM `Chen2020` (LG M50 21700). The LG M50 negative electrode is graphite with a small SiOx fraction; Chen2020, and this model, treat it as graphite only:
  - stoichiometry windows from `get_min_max_stoichiometries`: OCV 2.5 V at SOC 0 and 4.2 V at SOC 1, with a SOC-window capacity of 5.153 Ah;
  - particle radius, thickness, active fraction, diffusivity and maximum concentration;
  - the exchange-current coefficient and its 35 kJ/mol activation energy.

  The electrolyte values (conductivity, diffusivity, t+) are Nyman2008 at 1 mol/L, with 17 kJ/mol activation energies as in PyBaMM `OKane2022` (Ecker2015). This placeholder is not a bound by construction; the DFN check also runs the electrolyte at 30 kJ/mol and with the measured functions. Against Landesfeind & Gasteiger's (2019) measured EC:EMC 3:7 electrolyte the controller's conductivity, diffusivity and depletion rate are conservative down to 0 °C (docs/design-notes.md). The concentration overpotential uses (1 − t+)(1 + d ln f/d ln c) = 1.85, the upper bound of measured LiPF6 in EC:EMC 3:7 at ≤ 1 mol/L and −10 to 45 °C (Landesfeind & Gasteiger 2019, via PyBaMM `ORegan2022`). Chen2020's ideal-solution value of 0.74 understates the term by 2.1–2.5×. The OCV tables are generated from the Chen2020 OCP functions on a 201-point grid. Interpolation error is ≤ 2.4 mV overall and ≤ 0.2 mV above 10 % SOC.
- **Placeholders without a source.**
  - The activation energy of graphite diffusion. Chen2020 gives none; reported values range from about 17 (O'Regan2022) to 43 kJ/mol (Mohtat2020). The controller uses 60 kJ/mol below 25 °C and 17 kJ/mol above, the smallest diffusivity on either side.
  - The graphite diffusivity itself: published LG M50 values differ by up to 100× (Chen2020 3.3e-14 m²/s; O'Regan2022 1e-16 to 2e-14 m²/s, stoichiometry-dependent). The controller uses half of Chen2020's value, which keeps the 5 mV margin in the DFN for the worst pairing of the two sets (see Validation).
  - The graphite exchange-current activation energy. Reported values range from 35 kJ/mol (Chen2020) to about 53 kJ/mol (Ecker2015); the controller uses 60 kJ/mol below 25 °C and 35 kJ/mol above, the smallest exchange current on either side. This and the diffusion bound are why cold charging is slow (≈0.19C at 5 °C).
  - The 17 kJ/mol temperature scaling of R0/R1. These are fitted at 25 °C; characterise the cell at several temperatures to replace it.
  - The ageing model: linear rates (≈2 %/year calendar, ≈20 % after 1000 full equivalent cycles at 25 °C, 30 kJ/mol) and resistance growth of 2.5 × the relative capacity loss.
  - The ageing mode. The anode stoichiometry is bounded by (x0 + SOC · window · SOH) / (1 − f) and the interfacial current density is raised by 1 / (1 − f), where f is the larger of 1 − SOH and `EvCellConfig.anode_material_loss`. That bound holds for any mix of lithium loss and anode-material loss at any lithiation (checked on 90 combinations, docs/design-notes.md). Anode material lost while the anode still had spare capacity does not show in SOH, so `anode_material_loss` must come from a degradation-mode diagnostic (dV/dQ or ICA); 0 is right only for a fresh cell. The OCV table remains the fresh cell's.
- **Plating predictor.** A reduced, conservative surrogate:
  - a steady-state polynomial surface-concentration correction;
  - inverted Butler–Volmer kinetics;
  - the full-thickness electrolyte ohmic drop;
  - steady parabolic electrolyte depletion with the exact profile-averaged concentration overpotential, since the depletion time constant is ≈33 s;
  - an SEI film drop equal to the full resistance growth, scaled with temperature like the electrolyte (17 kJ/mol). SEI resistance is usually more temperature-sensitive, so this term is optimistic in the cold for an aged cell;
  - interfacial current density raised by 1/SOH;
  - concentration-gradient memory: after a step-down the solid surface excess and the electrolyte depletion take tens of seconds to minutes to relax, so their overpotentials are evaluated at the larger of the candidate current and a per-cell history current that decays with the slowest diffusion mode (≈52 s for the particles and ≈33 s for the electrolyte at 25 °C, slower in the cold). The history follows the measured current capped at the previous approval, so an in-tolerance current-sensor gain error cannot make the limit oscillate. The depletion limit itself is checked at the candidate current;
  - evaluation at the 3σ upper SOC, and at the colder of the present temperature and the zero-current prediction, minus the 2 °C sensor margin. No credit is taken for Joule self-heating, so every constraint is monotone in current and the approved current never falls when the request rises.

  With the example parameters (5 Ah, 30 s horizon, SOC known exactly, no recent charge current) it approves ≈0.58C at 25 °C and 50 % SOC, ≈0.32C at 90 % SOC (≈0.30C with a 1 % SOC σ), and ≈0.19C at 5 °C and 50 % SOC (`ev_cli limit`). Each approval is checked at the highest current that can flow without a fault: 5 % more plus the C/50 fast over-current tolerance, with the current sensor reading low by its C/100 offset. A charger anywhere inside the tolerance therefore never takes the cell past the limits; a sustained excess above 5 % + C/100 is also caught by the over-charge integral. `max_charge_current_A` caps the delivered-current bound at the cell maker's limit (`EV_CONSTRAINT_CURRENT_LIMIT`); update it from the maker's temperature derating when the coolant temperature changes.
- **Terminal voltage.** The bound over the horizon is the larger of the model (OCV at the 3σ upper SOC, the RC voltage's worst point and R0·I) and the latest reading plus the R0 step and the predicted OCV and RC rise, so a model sitting below the reading cannot hide the peak.
- **Thermal model.** Lumped. Mass 69 g, cooling area 53.1 cm² (PyBaMM's 21700 value). Heat source I²(R0 + R1) at the cell temperature.

## Sensor limits

With one voltage channel per cell, some sensor faults look the same as model error, so `ev_step()` can only bound them. The self-test checks each bound, and prints the measured values for the frozen-taper, offset, drift and charger-tolerance cases (`test_solver --self-test`):

- an offset of 50 mV or more is detected at once; a 30 mV offset is inside the 30 % R0 allowance and is not detected by the tracking threshold, but the reading becomes suspect, the growing SOC uncertainty ends the charge early and the suspect-charge latch raises `EV_FAULT_VOLTAGE_PLAUSIBILITY` (peak true SOC about 55 % in the self-test);
- a reading frozen during the constant-voltage taper, where the current changes slowly, can stay inside the tracking tolerance. It cannot restart the charge, and the coulomb-count floor on the upper SOC bound stops it: the worst case in the self-test peaks at about 101 % true SOC (OCV ≈ 4.22 V, inside the LG M50's 4.20 ± 0.05 V charge-voltage specification);
- a slow drift (−20 µV/s) is caught only by the SOC-divergence backstop, at about 105 % true SOC (OCV ≈ 4.29 V; the test bounds it below 106 %). With one voltage channel per cell this is a known gap: use the pack-voltage cross-check below or a redundant AFE measurement;
- a current-sensor gain error inside the 5 % over-current tolerance is not detected; the limits are computed at the current the charger may deliver with that error, so it does not erode the plating margin;
- the current sensor must be zero-calibrated to within C/100 (0.05 A for the example cell). This is a **safety requirement** in the under-reading direction: the limits budget exactly that offset, so a sensor that reads low by more hides current the limits do not cover. An over-reading offset fails safe (over-current or plausibility faults; up to C/100 is tolerated). A production BMS should cross-check the charger's own current measurement.

The tracking and divergence checks also allow for the fitted OCV table being imperfect: in the self-test a healthy cell whose OCV differs from the table by a sinusoidal error (10 mV over 0.33 SOC, or 5 mV over 0.2 SOC) charges without a latched fault across the parameter-mismatch corners. Larger OCV-shape errors (for example 10 mV over 0.2 SOC with 15 % less capacity, 30 % more R0 and twice R1) can latch a false plausibility fault, so the OCV table must be fitted to that accuracy for the cell in use.

`ev_step_with_pack_voltage()` closes the first three (the frozen-taper case peaks at about 97 %): in the self-test it detects the 30 mV offset at once and the drift at about 69 % SOC, with no false fault across the parameter-mismatch corners or with a 20 mΩ path between the cell taps and the pack sense point (pass that resistance in; the figures hold for a tolerance of 30 mV or less). An unavailable pack reading (for example a charger output that is not live yet) is skipped at rest; if it goes missing while current flows, the charge stays blocked until a valid reading is seen at rest. A production BMS should provide that redundant measurement.

## Limits of the evidence

This is research software, not a qualified product. Use on a vehicle requires at minimum:

1. cell-specific calibration of the OCV, 1-RC (at several temperatures), thermal and ageing parameters (`tools/calibrate_cell.py`);
2. experimental plating validation (for example three-electrode cells or post-mortem analysis);
3. HIL, fault-injection and timing tests on the target MCU;
4. requirements traceability and an ISO 26262 safety case.

The diagrams in `docs/diagrams/` summarise the architecture, the model scope and the validation roadmap.
