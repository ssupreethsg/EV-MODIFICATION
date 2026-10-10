# Design notes: v2.0.0 review

This note records the open issues found in the v2.0.0 review. For each one it gives the cause, the existing work it draws on, the fix chosen and the evidence for it. The numbers come from `tests/test_solver.c --self-test` and `tools/run_pybamm_sim.py`.

## 1. Voltage reading frozen near full charge (single channel)

**Cause.** In the constant-voltage taper the cell voltage hardly moves. A frozen reading therefore matches the model within the tracking tolerance, and the reading never validates again. Before this release, the SOC limit used only the *lower* bound of the coulomb count as a floor. In the self-test the true SOC peaked at 101.2 %. A slow drift (−20 µV/s) re-validated window after window and reached 105 % (OCV ≈ 4.29 V) before the divergence backstop fired.

**Existing practice.** Production BMS software, foxBMS for example (BSD-3), cross-checks every cell reading against a redundant measurement: a second AFE ADC, or the string voltage from the current sensor. AFEs built for ISO 26262 provide that redundant ADC path. This library's `ev_step_with_pack_voltage()` is the same pattern.

**Tried and not adopted.** We tried also holding the 98 % SOC limit against the coulomb count's *upper* bound once the reading stops validating. In the self-test that stopped the frozen taper at 94 % and the drift at 98 %. Against the DFN, however, a healthy reading seldom re-validates against the 1-RC model. The anchor uncertainty then ended a 25 °C charge at 78 % SOC instead of 95 %. With one channel, a frozen or drifting reading near full charge cannot be told apart from a healthy one without that cost.

**Decision.** The single-channel bounds stay as documented: about 101 % true SOC (OCV 4.22 V, inside the LG M50's 4.20 ± 0.05 V charge specification) for a reading frozen in the taper, and about 105 % for a slow drift. Safety-rated use needs the redundant measurement: `ev_step_with_pack_voltage()` (frozen taper 97 %, drift caught at 68 %) or a redundant AFE ADC.

## 2. Graphite transport and kinetics parameters

**Cause.** Two published parameter sets for the same LG M50 cell disagree. Chen2020 gives a constant D_s of 3.3e-14 m²/s. O'Regan2022 gives a D_s that depends on stoichiometry and drops to about 3e-16 m²/s, together with a larger exchange current. Each set is self-consistent. Their worst pairing, O'Regan's D_s with Chen's i0 and the measured electrolyte, left only +3.8 mV in the DFN at 25 °C, below the 5 mV margin.

**Fix.** The surrogate uses half of Chen2020's diffusivity (`EV_SOLID_DIFFUSIVITY_UNCERTAINTY`).

**Evidence.**
- The worst pairing improves to +6.5 mV at 25 °C, and the cold cases keep more than +20 mV.
- The nominal final SOC moves by less than 1 %.
- The DFN check now runs the worst pairing at 0, 2, 10 and 25 °C.

## 3. Electrolyte temperature dependence

**Cause.** The controller scales the Nyman2008 electrolyte transport with a 17 kJ/mol Arrhenius placeholder, which is not a proven bound.

**Existing data.** Landesfeind & Gasteiger (2019) measured LiPF6 in EC:EMC 3:7 from −10 to 50 °C, and PyBaMM's `ORegan2022` set ships those functions. At 1 mol/L and 0–2 °C the measurements give:

| Quantity | Measured (Landesfeind 2019) | Controller |
|---|---|---|
| Conductivity κ | 0.55–0.58 S/m | 0.51–0.53 S/m |
| Diffusivity D_e | 1.5–1.6e-10 m²/s | 0.95–0.99e-10 m²/s |
| (1−t+)·TDF | ≤ 1.82 | 1.85 |
| Depletion rate per unit current | | ≈ 7 % above the measured value |

The controller is the more conservative side on every row, so the measured electrolyte is milder than the controller's assumption.

**Fix.** The DFN check now also runs the measured electrolyte (all four functions, non-ideal), not only the ideal-solution Nyman2008 set.

## 4. Ageing mode (loss of anode material)

**Cause.** Capacity fade alone cannot reveal anode material lost while the anode still has spare capacity (LAM_NE with N/P excess). Such a loss pushes the graphite to a higher stoichiometry at the same SOC. In the old mapping, 20 % LAM_NE delithiated at low SOC moved the stoichiometry up by as much as 0.23 with SOH still 1.0.

**Fix.**
- A new config field `anode_material_loss` takes LAM_NE as measured by a degradation-mode diagnostic (dV/dQ or ICA, following Birkl 2017 and Dubarry).
- The bulk stoichiometry is bounded by (x0 + SOC · window · SOH) / (1 − f), with f = max(`anode_material_loss`, 1 − SOH). The interfacial current density is scaled by 1 / (1 − f).

**Evidence.**
- The bound was checked against 90 combinations: LAM_NE 5–20 %, lost at SOC 0–1, lithiated, delithiated or half, with 0 or 5 % LLI. It is never exceeded, and is tight in the worst case.
- With f = 1 − SOH (no diagnostic) it reduces to the previous mapping plus x0(1/SOH − 1).

## 5. Charger trickle at a 0 A command

**Cause.** At a 0 A command the over-current tolerances accept a measured C/100 (a true C/50 with the sensor offset). No limit budgets that current, so a charger that never stops could overcharge a full cell without any fault.

**Fix.** At a 0 A command, once the ramp-down has ended, a measured voltage above 4.21 V, or an EKF SOC or coulomb-count lower bound above 99 %, raises `EV_FAULT_OVER_CURRENT`.

**Evidence.** A C/50 trickle into a 90 % cell faults at a true SOC of 98.9 %, and a full cell at rest does not fault.

## 6. Smaller corrections

- **Voltage-tracking false fault.** A steady charge in the 15 %-low-capacity corner raised a false fault in 2 of 128 runs of a wider sweep. The cause was the window model relaxing an EKF polarisation offset that the plant does not have.
  - The part of that offset relaxed so far is now allowed for, and the OCV spread also scales the window by the capacity error.
  - The sweep now gives 0 of 128, and the detection tests are unchanged.
- **Measured voltage floor after a step down.** It now uses 0.7 × R0, the smallest plausible value.
- **Temperature sensor margin.** The 2 °C margin now applies to the 45 °C charge limit as well as the cold limit.
- **Prediction horizon.** The floor is EV_MAX_CONTROL_DT_S + 2 s, so a late next cycle stays covered.
- **Suspect reading.** An estimator reset keeps a valid suspect flag.

## Remaining known gaps

These need hardware or cell data rather than code:

- a redundant voltage measurement (pack-voltage cross-check or redundant AFE ADC) for safety-rated use; one channel per cell bounds a frozen or drifting reading only to about 101 % and 105 % true SOC;
- experimental plating validation, for example three-electrode cells;
- cell-specific OCV, 1-RC, thermal and ageing calibration;
- a measured LAM_NE value for aged cells;
- HIL and timing tests on the target MCU;
- an ISO 26262 safety case.
