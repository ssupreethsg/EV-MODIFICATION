# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/).

## [2.0.0] - 2026-10-10

Reworked from the earlier firmware prototype into a validated charge limiter. Review rounds 1–5.

### Changed (breaking)
- The library files are now `src/ev_modify.h` and `src/ev_modify.c` (formerly `spme_v7_apex_solver.*`, a name left over from an earlier project). Include `ev_modify.h`.
- `ev_cli ekf` prints a fifth column with the fault flags.

### Added
- Two-state EKF (SOC, RC voltage) with Joseph-form update and innovation gate; per-cell SOH and resistance-growth tracking.
- Voltage-tracking check that catches a sensor frozen at an in-range value, with or without noise; allowances for parameter errors and for RC transients of recent current steps, so a healthy PyBaMM cell never trips it.
- Over-current fault: the measured current must follow the previous approval, with a 2 s ramp-down allowance counted in seconds.
- Concentration-gradient memory that decays with the slowest solid and electrolyte diffusion modes (temperature dependent), driven by the measured current capped at the previous approval.
- Terminal-voltage floor from the latest reading includes the OCV and RC rise over the horizon.
- Latched sensor faults: the voltage-tracking fault (confirmed over 3 s, so short glitches do not latch) stays set until `ev_cell_init()`, so the charger stop it causes cannot clear it.
- `EV_FAULT_SOC_DIVERGENCE`: backstop that faults when the EKF falls below the coulomb count kept since the reading last validated; the same count also floors the upper SOC bound.
- `ev_step_with_pack_voltage()` and `EV_FAULT_PACK_VOLTAGE`: optional cross-check of the cell readings against an independent string-voltage measurement, which catches small offsets and slow drift that one channel per cell cannot.
- Over-current also integrates charge beyond the expected current (1 h leak, 2 % of capacity), so a charger floor below the absolute tolerance is caught.
- `ev_cli charge`: the complete controller in a closed loop, used by the PyBaMM DFN benchmark.
- CI: GCC, Clang, ASan/UBSan, MSVC, Cortex-M4F cross-compile, cppcheck, coverage ≥ 95 %, ruff, PyBaMM validation. Dependabot for actions and pip.

### Fixed
- The depletion limit is checked at the candidate current; evaluating it at the previous current blocked charging even at 0 A and made the approval alternate with a 1 % current-sensor gain error.
- No credit for Joule self-heating in the electrochemistry and 0 °C checks, so every constraint is monotone in current.
- Concentration overpotential uses the measured LiPF6 diffusion-potential factor bound (1.85, Landesfeind & Gasteiger 2019) instead of the ideal-solution 0.74.
- Thermal prediction exact at equilibrium.
- Tools: the PyBaMM validation profile stops at 4.2 V; calibration uses the measured temperature for R0 and R1, as the controller does; `ev_cli limit` range-checks the request and horizon; sanitizer builds no longer generate install rules.

### Documented
- Sensor limits of a single voltage channel per cell, with the bound the self-test measures for each.
- Ageing assumption and its limits (SOH cannot bound loss of partly lithiated anode material).
- Assumptions behind the quoted C-rates.
