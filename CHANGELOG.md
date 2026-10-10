# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/).

## [Unreleased]

## [2.0.0] - 2026-10-10

Reworked from the earlier firmware prototype into a validated charge limiter.

### Changed
- **Breaking:** the library files are now `src/ev_modify.h` and `src/ev_modify.c` (formerly `spme_v7_apex_solver.*`, a name left over from an earlier project). Include `ev_modify.h`.
- **Breaking:** `ev_cli ekf` prints a fifth column with the fault flags; `ev_cli` exits 3 for dt outside [0.001, 10] s and 4 when writing stdout fails.
- **Breaking:** latched faults are cleared with `ev_clear_latched_faults()`; `ev_cell_init()` is for a new cell only (persist `EvCellState` across power cycles).
- Every limit is computed at the highest current the charger may deliver indefinitely without a fault (command × 1.05 + C/100), and the over-charge integral faults on a sustained excess above it.
- **Breaking:** `EvCellConfig` has a new `max_charge_current_A` field (the cell maker's limit, `EV_CONSTRAINT_CURRENT_LIMIT`); a configuration without it is invalid.
- The graphite exchange-current activation energy is a worst-case bound (60 kJ/mol below 25 °C, 35 kJ/mol above), so cold charging is slower and the plating margin holds over the published range.
- The PyBaMM reference also sweeps the exchange-current activation energy (35 / 53.4 / 60 kJ/mol) and runs a charger at and outside its tolerance.
- The PyBaMM DFN reference uses activation energies independent of the controller (a 2 × 2 sweep at 2 and 10 °C), and adds a lumped-thermal run, a cold constant-voltage taper run and a terminal-voltage check.
- Documentation: sensor limits of a single voltage channel per cell with the bound the self-test measures for each, the ageing assumption and its limits, and the assumptions behind the quoted C-rates.

### Added
- Two-state EKF (SOC, RC voltage) with Joseph-form update and innovation gate; per-cell SOH and resistance-growth models.
- Voltage-tracking check that catches a sensor frozen at an in-range value, with or without noise. It is confirmed over 3 s, so short glitches do not latch, and it stays latched, so the charger stop it causes cannot clear it.
- `EV_FAULT_SOC_DIVERGENCE`: backstop that faults when the EKF falls below the coulomb count kept since the reading last validated; the same count also floors the upper SOC bound.
- `ev_step_with_pack_voltage()` and `EV_FAULT_PACK_VOLTAGE`: optional cross-check of the cell readings against an independent string-voltage measurement, with a path-resistance term, which catches small offsets and slow drift that one channel per cell cannot.
- Over-current fault: the measured current must follow the previous approval, with a 2 s ramp-down allowance; a leaky integral of the charge beyond the delivered-current bound catches a sustained excess inside the absolute tolerance.
- Concentration-gradient memory that decays with the slowest solid and electrolyte diffusion modes (temperature dependent).
- Terminal-voltage floor from the latest reading includes the OCV and RC rise over the horizon.
- `ev_cli charge`: the complete controller in a closed loop, used by the PyBaMM DFN benchmark.
- CI: GCC, Clang, ASan/UBSan, MSVC, Cortex-M4F cross-compile, cppcheck, coverage ≥ 95 %, ruff lint and format, CLI exit-code tests, `find_package` and `add_subdirectory` consumer builds, PyBaMM validation. Dependabot for actions and pip.

### Fixed
- The depletion limit is checked at the candidate current; evaluating it at the previous current blocked charging even at 0 A and made the approval alternate with a 1 % current-sensor gain error.
- A surface stoichiometry beyond the OCP table is a plating violation, and the history depletion is no longer capped at the limit, so neither can under-predict the overpotential.
- An invalid voltage sample no longer restarts the tracking window or clears a suspicion; an intermittent voltage reading is a latched sensor fault, and while a reading is suspect the SOC sigma grows with the charge counted.
- A pack reading that goes missing during a charge keeps the charge blocked until it is valid at rest, instead of letting it duty-cycle.
- No credit for Joule self-heating in the electrochemistry and 0 °C checks, so every constraint is monotone in current.
- Concentration overpotential uses the measured LiPF6 diffusion-potential factor bound (1.85, Landesfeind & Gasteiger 2019) instead of the ideal-solution 0.74.
- Thermal prediction exact at equilibrium.
- Tools: the PyBaMM validation profile stops at 4.2 V; calibration uses the measured temperature for R0 and R1, as the controller does; `ev_cli` range-checks every value before converting it and rejects NUL bytes; sanitizer builds no longer generate install rules.

[Unreleased]: https://github.com/ssupreethsg/EV-MODIFICATION/compare/v2.0.0...HEAD
[2.0.0]: https://github.com/ssupreethsg/EV-MODIFICATION/compare/fb3d4ff...v2.0.0
