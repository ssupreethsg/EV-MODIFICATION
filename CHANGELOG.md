# Changelog

## 2026-10-10 — review rounds 3 and 4

- Voltage-tracking check replaces the frozen-sensor dead band: catches a sensor frozen at an in-range value with or without noise, keeps the EKF from following it and latches the fault until the reading tracks again.
- Over-current fault: the measured current must follow the previous approved current, with a 2 s ramp-down allowance counted in seconds, not cycles.
- No credit for Joule self-heating in the electrochemistry and 0 °C checks, so every constraint is monotone in current and the approved current never falls when the request rises.
- Concentration overpotential uses the measured LiPF6 diffusion-potential factor bound (1.85, Landesfeind & Gasteiger 2019) instead of the ideal-solution 0.74.
- Concentration gradients remember the previous charge current after a step-down; thermal prediction is exact at equilibrium.
- Ageing assumption documented with its limits (SOH cannot bound loss of partly lithiated anode material).
- Tools: the DFN plating check fails when the limiter does not charge, `ev_cli` range-checks its input, calibration uses the C model's temperature convention, `generate_ocv_lut.py` rejects unknown flags; CMake keeps sanitizer flags out of installed targets and adds `EV_INSTALL`; CI lints the Python tools.

## 2026-10-10 — 8cf4af1

- Reworked into a validated charge limiter with CI (review rounds 1–2).
