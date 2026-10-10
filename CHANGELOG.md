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
- Every limit is computed at the highest current that can flow without a fault: command × 1.05 plus the C/50 fast over-current tolerance, with the current sensor reading low by its C/100 offset. The fast tolerance (formerly 5 % + 0.5 A) scales with capacity, so the limits and the fault logic use one allowance; a sustained excess above command × 1.05 + C/100 is also caught by the over-charge integral.
- **Breaking:** `EvCellState` starts with `layout_version` (`EV_STATE_LAYOUT_VERSION` 4); a restored state with another version faults until `ev_cell_init()`. Flag bytes other than 0 or 1 are a state fault and read as set.
- **Breaking:** `EvCellConfig` has a new `max_charge_current_A` field (the cell maker's limit, `EV_CONSTRAINT_CURRENT_LIMIT`); a configuration without it is invalid.
- The graphite exchange-current and solid-diffusion activation energies are worst-case bounds (60 kJ/mol below 25 °C; 35 and 17 kJ/mol above), so cold charging is slower and the plating margin holds over the published range. Resistances use the same 2 °C cold sensor margin as the electrochemistry.
- The PyBaMM reference also sweeps the exchange-current activation energy (35 / 53.4 / 60 kJ/mol), runs a 0 °C cell whose sensor reads 2 °C warm with a charger just inside (current sensor reading C/100 low), outside, and at the highest current that never faults, and runs O'Regan2022's graphite diffusivity and exchange current for the same cell, cold and warm. Its plating floor is the controller's 5 mV margin, sampled at every point and time step.
- The PyBaMM DFN reference uses activation energies independent of the controller (a 2 × 2 sweep at 2 and 10 °C), and adds a lumped-thermal run, a cold constant-voltage taper run and a terminal-voltage check.
- Documentation: sensor limits of a single voltage channel per cell with the bound the self-test measures for each, the ageing assumption and its limits, and the assumptions behind the quoted C-rates.

### Added
- Two-state EKF (SOC, RC voltage) with Joseph-form update and innovation gate; per-cell SOH and resistance-growth models.
- Voltage-tracking check that catches a sensor frozen at an in-range value, with or without noise. It is confirmed over 3 s, so short glitches do not latch, and it stays latched, so the charger stop it causes cannot clear it.
- `EV_FAULT_SOC_DIVERGENCE`: backstop that faults when the EKF falls below the coulomb count kept since the reading last validated; the same count also floors the upper SOC bound.
- `ev_step_with_pack_voltage()` and `EV_FAULT_PACK_VOLTAGE`: optional cross-check of the cell readings against an independent string-voltage measurement, with a path-resistance term, which catches small offsets and slow drift that one channel per cell cannot.
- Over-current fault: the measured current must follow the previous approval, with a 2 s ramp-down allowance; a leaky integral of the charge beyond command × 1.05 + C/100 catches a sustained excess inside the fast tolerance.
- `EV_MODIFY_VERSION_MAJOR/MINOR/PATCH` in the header, checked against the CMake project version.
- Concentration-gradient memory that decays with the slowest solid and electrolyte diffusion modes (temperature dependent).
- Terminal-voltage floor from the latest reading includes the OCV and RC rise over the horizon.
- `ev_cli charge`: the complete controller in a closed loop, used by the PyBaMM DFN benchmark.
- CI: GCC, Clang, ASan/UBSan, MSVC, Cortex-M4F cross-compile with a size budget and a newlib link check, cppcheck, coverage ≥ 95 %, ruff lint and format (pinned in `requirements-dev.txt`), CLI exit-code and output tests, `find_package` and `add_subdirectory` consumer builds, PyBaMM validation. Dependabot for actions and pip.

### Fixed
- The depletion limit is checked at the candidate current; evaluating it at the previous current blocked charging even at 0 A and made the approval alternate with a 1 % current-sensor gain error.
- A surface stoichiometry beyond the OCP table is a plating violation, and the history depletion is no longer capped at the limit, so neither can under-predict the overpotential.
- An invalid voltage sample no longer restarts the tracking window or clears a suspicion; an intermittent voltage reading is a latched sensor fault, and while a reading is suspect the SOC sigma grows with the charge counted.
- A pack reading that goes missing during a charge keeps the charge blocked until it is valid at rest, instead of letting it duty-cycle; plain `ev_step()` reports the block too.
- A reading that stays suspect while 10 % of the capacity is charged latches the plausibility fault, so charge pauses caused by invalid samples cannot keep resetting the confirmation timer.
- The gradient memory is driven by the measured current plus the C/100 sensor offset, as the limits are, and every check spans at least dt plus the 2 s charger ramp-down.
- `ev_clear_latched_faults()` documents everything it clears (intermittent-sensor latch, pack block, over-charge integral).
- An estimator reset repairs every validated field (pack-mismatch timer, intermittency score, hold counter), so one corrupted field no longer faults forever.
- No credit for Joule self-heating in the electrochemistry and 0 °C checks, so every constraint is monotone in current.
- Concentration overpotential uses the measured LiPF6 diffusion-potential factor bound (1.85, Landesfeind & Gasteiger 2019) instead of the ideal-solution 0.74.
- Thermal prediction exact at equilibrium.
- Tools: the PyBaMM validation profile stops at 4.2 V; calibration uses the measured temperature for R0 and R1, as the controller does; `ev_cli` range-checks every value before converting it and rejects NUL bytes; sanitizer builds no longer generate install rules.

[Unreleased]: https://github.com/ssupreethsg/EV-MODIFICATION/compare/v2.0.0...HEAD
[2.0.0]: https://github.com/ssupreethsg/EV-MODIFICATION/compare/fb3d4ff...v2.0.0
