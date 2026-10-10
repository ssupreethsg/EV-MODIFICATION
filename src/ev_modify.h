/*
 * EV MODIFY - conservative charge-current limiter for a series lithium-ion module.
 *
 * Copyright (c) 2026 Supreeth S G. All rights reserved. See LICENSE.
 *
 * Usage (one call per control cycle):
 *
 *     EvCellConfig config[EV_NUM_CELLS];   // filled once by the integrator
 *     EvCellState  state[EV_NUM_CELLS];
 *     if (!nvm_restore(state)) {            // integrator's CRC-checked copy (layout_version is checked too)
 *         for (i = 0; i < EV_NUM_CELLS; ++i) { ev_cell_init(&state[i], initial_soc[i], initial_temp_C[i]); }
 *     }                                      // new cell only: ev_cell_init() clears latches and ageing
 *     ...
 *     EvPackResult r = ev_step(config, state, pack_current_A, requested_A,
 *                              cell_voltages_V, cell_temperatures_C, dt_s, horizon_s);
 *     command_charger(r.approved_charge_current_A);   // 0 A on EV_DECISION_FAULT
 *     nvm_save(state);                                // periodically and at shutdown
 *
 * Sign convention: positive current = charging. Units are SI unless the name says otherwise.
 * The library is single-threaded and keeps no global state; all state lives in EvCellState.
 */
#ifndef EV_MODIFY_H
#define EV_MODIFY_H

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

#define EV_NUM_CELLS 4U
#define EV_MODIFY_VERSION_MAJOR 2   /* library version; CMake checks it against project(VERSION) */
#define EV_MODIFY_VERSION_MINOR 0
#define EV_MODIFY_VERSION_PATCH 0
#define EV_STATE_LAYOUT_VERSION 4U   /* EvCellState.layout_version; changes whenever EvCellState changes */

/* ---- Physical constants ---------------------------------------------------------------- */
#define EV_FARADAY_C_PER_MOL 96485.33212f
#define EV_GAS_J_PER_MOL_K 8.31446262f
#define EV_KELVIN_OFFSET 273.15f
#define EV_REFERENCE_TEMPERATURE_K 298.15f

/* ---- Safety limits ---------------------------------------------------------------------- */
#define EV_MAX_CELL_VOLTAGE_V 4.20f          /* LG M50 charge voltage */
#define EV_MAX_SAFE_SOC 0.980f               /* applied to the 3-sigma upper SOC bound */
#define EV_MIN_SAFE_ANODE_POTENTIAL_V 0.005f /* plating margin vs Li/Li+ */
#define EV_MIN_CHARGE_TEMPERATURE_C 0.0f     /* LG M50 datasheet: no charging below 0 degC */
#define EV_MAX_CELL_TEMPERATURE_C 45.0f      /* charge limit and over-temperature fault */
#define EV_MIN_CELL_TEMPERATURE_C -30.0f     /* under-temperature fault */
#define EV_TEMPERATURE_HYSTERESIS_C 3.0f     /* a faulted cell must come this far back inside */
#define EV_FAULT_HOLD_CYCLES 5U              /* consecutive good cycles needed to clear a fault */
#define EV_TEMPERATURE_SENSOR_MARGIN_C 2.0f  /* cold bias for the plating check and the 0 degC limit */
#define EV_OVER_VOLTAGE_MARGIN_V 0.05f       /* measured V above 4.2 V + margin => over-voltage fault */
#define EV_VOLTAGE_HYSTERESIS_V 0.05f        /* a faulted cell must come back to 4.2 V */
#define EV_SOC_SIGMA_MULTIPLIER 3.0          /* SOC constraints use soc + k * sigma */
#define EV_INNOVATION_GATE_SIGMA 5.0         /* EKF rejects voltage innovations beyond this */
#define EV_INNOVATION_FAULT_CYCLES 10U       /* consecutive rejections => voltage plausibility fault */
/* Voltage tracking (frozen or stuck sensor, open wire): over a window the change of the reading
 * must match the change of an open-loop model (coulomb counting + 1-RC + R0, started from the EKF
 * state at the window start), after allowing for 30 % R0, 100 % R1 and 20 % OCV-change errors and
 * for the OCV slope anywhere within the EKF's 3-sigma SOC band. Noise on a stuck
 * reading does not help it pass: the test compares changes, not equality. */
#define EV_TRACK_SUSPECT_V 0.010f            /* mismatch above this: the EKF ignores the reading until the window validates */
#define EV_TRACK_FAULT_V 0.020f              /* mismatch above this ... */
#define EV_INTERMITTENT_SAMPLES 2.5f         /* decaying count of invalid voltage samples (3 within about 20 s) ... */
#define EV_INTERMITTENT_WINDOW_S 60.0f       /* ... decaying with this time constant => EV_FAULT_VOLTAGE_SENSOR, latched */
#define EV_TRACK_FAULT_CONFIRM_S 3.0f        /* ... sustained this long (seconds, at least 2 cycles) => plausibility fault, latched until ev_clear_latched_faults() */
#define EV_TRACK_SUSPECT_LATCH_SOC 0.10       /* charge counted while suspect, never validating, beyond this fraction of capacity => plausibility fault, latched */
#define EV_TRACK_VALIDATE_V 0.030f           /* window validated (and restarted) once the reading has moved this much beyond the tolerance, as modelled */
#define EV_TRACK_SPIKE_V 0.050f              /* consecutive readings further apart than this are not used to (re)start a window */
#define EV_TRACK_R1_UNCERTAINTY 1.0f         /* relative R1 error allowed for; with any time constant the plant RC voltage lies in [0, (1 + u) R1 dI] */
#define EV_TRACK_TAU_MIN_FACTOR 0.3f        /* the plant RC time constant may be this fraction of the configured one */
#define EV_TRACK_SOC_BIAS 0.05             /* SOC error beyond the EKF's 3 sigma that resistance errors can cause (R0 +30 %, R1 x2) */
#define EV_TRACK_CAPACITY_UNCERTAINTY 0.2f   /* relative error of the open-loop OCV change (capacity, OCV slope) */
#define EV_TRACK_REST_C_RATE 0.02f           /* below this |current| / capacity [1/h] the cell counts as at rest */
/* Measured charge current must follow the approved command (charger ignoring the limit). Every
 * limit is computed at the highest current that can flow without a fault: the command plus this
 * tolerance plus the current-sensor offset, so no excess the fault logic accepts is unbudgeted. */
#define EV_OVER_CURRENT_TOLERANCE_C_RATE 0.02f /* allowed excess over the approved current, as a C-rate ... */
#define EV_OVER_CURRENT_TOLERANCE_REL 0.05f  /* ... plus this fraction of it */
#define EV_OVER_CURRENT_RAMP_S 2.0f          /* after a lower command the charger may ramp down linearly over this time */
#define EV_OVER_CURRENT_CONFIRM_S 1.0f       /* excess sustained this long (in seconds, not cycles) => fault */
#define EV_OVER_CHARGE_LIMIT_FRACTION 0.005f /* charge beyond expected x 1.05 + offset (a sustained excess inside the tolerance), leaky over ... */
#define EV_CURRENT_OFFSET_C_RATE 0.01f       /* current-sensor offset after zero calibration, as a C-rate (requirement) */
#define EV_OVER_CHARGE_LEAK_S 3600.0f        /* ... this time constant, above this fraction of capacity => over-current fault */
#define EV_SOC_DIVERGENCE_LIMIT 0.02         /* EKF below the coulomb count since the sensor last worked, beyond 3 sigma => fault */
#define EV_COUNT_CAPACITY_UNCERTAINTY 0.05   /* capacity error of the coulomb count (SOH model) for its lower SOC bound */
#define EV_SOC_DIVERGENCE_ARM_SIGMA 0.03     /* the divergence check is armed once the EKF SOC sigma is below this */
#define EV_CURRENT_TOLERANCE_A 1.0e-3f

/* ---- Input plausibility (outside => FAULT) ----------------------------------------------- */
#define EV_MAX_PLAUSIBLE_CURRENT_A 1000.0f   /* |measured current|; requests are clamped to it */
#define EV_MIN_PLAUSIBLE_VOLTAGE_V 2.0f
#define EV_MAX_PLAUSIBLE_VOLTAGE_V 4.5f
#define EV_MIN_PLAUSIBLE_TEMPERATURE_C -40.0f
#define EV_MAX_PLAUSIBLE_TEMPERATURE_C 100.0f
#define EV_MIN_CONTROL_DT_S 0.001f
#define EV_MAX_CONTROL_DT_S 10.0f
#define EV_MAX_PREDICTION_HORIZON_S 600.0f
#define EV_MAX_PACK_VOLTAGE_TOLERANCE_V 0.5f /* upper limit of the pack-voltage cross-check tolerance */
#define EV_MAX_PACK_PATH_RESISTANCE_OHM 0.1f /* upper limit of the resistance between cell taps and pack sense */

/* ---- Electrochemical parameters: PyBaMM Chen2020 (LG M50 21700) -------------------------- */
/* Stoichiometry windows from pybamm.lithium_ion.get_min_max_stoichiometries (2.5-4.2 V,
 * lithium conserving). SOC is linear in stoichiometry between these limits. */
#define EV_ANODE_STOICH_AT_SOC0 0.026346f
#define EV_ANODE_STOICH_AT_SOC100 0.910618f
#define EV_CATHODE_STOICH_AT_SOC0 0.853975f
#define EV_CATHODE_STOICH_AT_SOC100 0.263845f
#define EV_REFERENCE_CAPACITY_AH 5.1532f     /* SOC-window capacity of the electrode area below */
#define EV_GEOMETRIC_ELECTRODE_AREA_M2 0.1027f
#define EV_ANODE_THICKNESS_M 85.2e-6f
#define EV_PARTICLE_RADIUS_M 5.86e-6f
#define EV_MAX_CONCENTRATION_MOL_M3 33133.0f
#define EV_SOLID_DIFFUSIVITY_REF_M2_S 3.3e-14f
#define EV_ANODE_ACTIVE_FRACTION 0.75f
/* i0 = k * sqrt(theta * (1 - theta) * c_e / c_e0), k = m_ref * sqrt(c_e0) * c_s_max (Chen2020) */
#define EV_EXCHANGE_CURRENT_REF_A_M2 0.679f
/* Graphite exchange-current activation energy: reported values span 35 (Chen2020) to about
 * 53 kJ/mol (Ecker2015) and more. The bound on i0 uses the high end below 25 degC and the low
 * end above it, so it is the smallest i0 the range allows on either side. */
#define EV_KINETIC_ACTIVATION_ENERGY_COLD_J_MOL 60000.0f
#define EV_KINETIC_ACTIVATION_ENERGY_WARM_J_MOL 35000.0f
/* Graphite solid-diffusivity activation energy: Chen2020 gives none; reported values span about
 * 17 (O'Regan2022) to 43 kJ/mol (Mohtat2020), and the DFN check sweeps up to 60. The same rule as
 * for i0 gives the smallest diffusivity (largest surface rise, slowest relaxation) on either side. */
#define EV_DIFFUSION_ACTIVATION_ENERGY_COLD_J_MOL 60000.0f
#define EV_DIFFUSION_ACTIVATION_ENERGY_WARM_J_MOL 17000.0f

/* Electrolyte: LiPF6 in EC:EMC, Nyman2008 values at 1 mol/L */
#define EV_ELECTROLYTE_CONCENTRATION_MOL_M3 1000.0f
#define EV_ELECTROLYTE_POROSITY 0.25f
#define EV_BRUGGEMAN_EXPONENT 1.5f
#define EV_TRANSFERENCE_NUMBER 0.2594f
/* Diffusion-potential factor (1 - t+) * (1 + dln f/dln c) of the concentration overpotential.
 * Chen2020/Nyman2008 assume an ideal solution (thermodynamic factor 1, product 0.74), but
 * measured LiPF6 in EC:EMC 3:7 (Landesfeind & Gasteiger 2019, J. Electrochem. Soc. 166 A3079)
 * gives 1.58-1.82 at 1 mol/L between -10 and 45 degC and less at lower concentration. The
 * upper bound is used, so the term is never smaller than the measured electrolyte's. */
#define EV_DIFFUSION_POTENTIAL_FACTOR 1.85f
#define EV_ELECTROLYTE_CONDUCTIVITY_S_M 0.9487f
#define EV_ELECTROLYTE_DIFFUSIVITY_M2_S 1.769e-10f
/* Nyman2008 has no temperature dependence; 17 kJ/mol as in PyBaMM OKane2022 (Ecker2015) */
#define EV_ELECTROLYTE_ACTIVATION_ENERGY_J_MOL 17000.0f
#define EV_ELECTROLYTE_DIFFUSIVITY_ACTIVATION_J_MOL 17000.0f
/* R0 and R1 are fitted at 25 degC and scaled like the electrolyte conductivity (placeholder until
 * the cell is characterised at several temperatures) */
#define EV_RESISTANCE_ACTIVATION_ENERGY_J_MOL 17000.0f

/* Lumped thermal model: 21700 can, 21 mm x 70 mm (PyBaMM "Cell cooling surface area") */
#define EV_CELL_MASS_KG 0.069f
#define EV_CELL_COOLING_AREA_M2 0.00531f

/* Ageing: capacity fade raises the anode interfacial current density by 1/SOH (as if it were all
 * loss of anode active material) while the fresh SOC -> stoichiometry mapping is kept. That is
 * exact only for anode material lost fully lithiated, and conservative for loss of lithium
 * inventory (the usual dominant mode: anode potential +2 to +11 mV at 5-20 % LLI, Chen2020).
 * Anode material lost partly lithiated raises the anode stoichiometry instead: 10 % lost at
 * theta = 0.7 / 0.5 / 0.3 (SOH 0.93 / 0.95 / 0.97) puts the anode 11 / 20 / 38 mV below this
 * model at 1C, and SOH alone cannot bound it. Use with aged cells requires the degradation
 * modes from dV/dQ or three-electrode diagnostics.
 * Linear placeholder rates of typical NMC/graphite magnitude at 25 degC (about 2 %
 * capacity per year calendar, 20 % after 1000 full equivalent cycles), Arrhenius scaled.
 * They must be replaced by cell-specific ageing data before use. */
#define EV_CALENDAR_LOSS_PER_S 6.3e-10
#define EV_CYCLE_LOSS_PER_FEC 2.0e-4
#define EV_AGING_ACTIVATION_ENERGY_J_MOL 30000.0f
#define EV_RESISTANCE_GROWTH_PER_LOSS 2.5    /* R0 rises 2.5 x relative capacity loss */
#define EV_MIN_SOH 0.5

/* ---- Types ------------------------------------------------------------------------------- */

/* Per-cell configuration, written by the integrator and validated on every ev_step(). */
typedef struct {
    float nominal_capacity_Ah;     /* [0.1, 1000] capacity between the 2.5 V and 4.2 V OCV limits */
    float resistance_ohm;          /* [0, 1]  ohmic resistance R0 at 25 degC, beginning of life */
    float rc_resistance_ohm;       /* [0, 1]  polarisation resistance R1 at 25 degC (1-RC model) */
    float rc_time_constant_s;      /* [1, 1e4] R1*C1 time constant */
    float heat_transfer_W_m2_K;    /* [1, 1000] convective coefficient to the coolant */
    float heat_capacity_J_kg_K;    /* [100, 5000] specific heat */
    float coolant_temperature_C;   /* [-40, 100] ambient/coolant temperature */
    float max_charge_current_A;    /* (0, 1000] cell maker's maximum charge current (LG M50: 0.7C); set it from the
                                    * maker's temperature derating when the coolant temperature changes */
} EvCellConfig;

/* Per-cell persistent state. Initialise with ev_cell_init(); afterwards only the library's functions (ev_step(), ev_step_with_pack_voltage(), ev_clear_latched_faults()) write it. */
typedef struct {
    uint32_t layout_version;       /* EV_STATE_LAYOUT_VERSION, set by ev_cell_init(); any other value => EV_FAULT_STATE until ev_cell_init() */
    double soc;                    /* [0, 1] EKF estimate */
    double covariance[3];          /* EKF covariance of [soc, v_rc]: P00, P01, P11 */
    double soh;                    /* [EV_MIN_SOH, 1] capacity retention, never increases */
    double throughput_Ah;          /* total |charge| processed */
    float rc_voltage_V;            /* EKF estimate of the polarisation voltage */
    float temperature_C;           /* latest plausible measurement (model value if the sensor fails) */
    float resistance_growth_ohm;   /* [0, 1] ageing-related R0 increase, treated as anode SEI film */
    float last_voltage_V;          /* latest plausible measured terminal voltage (0 = none yet) */
    float last_current_A;          /* current of the latest cycle */
    float solid_history_current_A; /* decaying charge-current memory for the solid surface gradient */
    float electrolyte_history_current_A; /* the same for the electrolyte concentration gradient */
    uint16_t innovation_fail_cycles; /* consecutive EKF innovations outside the gate */
    float track_start_V;           /* reading at the start of the voltage-tracking window (0 = none) */
    float track_start_model_V;     /* open-loop model voltage at the start of the window */
    float track_start_ocv_V;       /* OCV at the start of the window */
    double track_start_soc;        /* EKF SOC at the start of the window */
    float track_rc_voltage_V;      /* open-loop polarisation voltage since the start of the window */
    double track_soc;              /* open-loop (coulomb-counted) SOC since the start of the window */
    float track_start_current_A;   /* current at the start of the window */
    float track_step_filtered_A;   /* |current - start current| filtered with the fastest plausible RC time constant */
    float track_prior_step_A;      /* recent steps before the window started, decaying with the RC time constant */
    float recent_current_step_A;   /* current steps decaying with the RC time constant (RC transients in progress) */
    uint16_t track_fail_cycles;    /* consecutive cycles with a mismatch above EV_TRACK_FAULT_V */
    float track_fail_time_s;       /* time since the first of those cycles (latched at EV_TRACK_FAULT_CONFIRM_S) */
    float rest_time_s;             /* continuous time at rest (|I| < EV_TRACK_REST_C_RATE * C) */
    double reference_soc;          /* coulomb count since the voltage sensor was last shown to track the model */
    double reference_sigma;        /* EKF SOC sigma when that count started */
    double reference_throughput;   /* |charge| / capacity counted since then */
    bool soc_divergence;           /* EKF and coulomb count disagreed: latched until ev_clear_latched_faults() */
    bool sensor_fault_latched;     /* voltage-tracking plausibility fault: latched until ev_clear_latched_faults() */
    float excess_charge_Ah;        /* leaky integral of charge beyond the expected current */
    uint16_t pack_mismatch_cycles; /* consecutive cycles with the cell sum outside the pack-voltage tolerance */
    float pack_mismatch_time_s;    /* time since the first of those cycles */
    bool pack_voltage_fault_latched; /* pack-voltage cross-check failed: latched until ev_clear_latched_faults() */
    bool pack_unavailable_block;   /* pack reading missing while current flowed: charging blocked until it is valid at rest */
    float invalid_voltage_score;   /* invalid voltage samples, decaying with EV_INTERMITTENT_WINDOW_S */
    bool voltage_sensor_latched;   /* intermittent voltage reading: latched until ev_clear_latched_faults() */
    bool track_suspect;            /* mismatch above EV_TRACK_SUSPECT_V: EKF ignores the reading until it validates */
    double suspect_charge;         /* |charge| / capacity counted while suspect (latched at EV_TRACK_SUSPECT_LATCH_SOC) */
    float commanded_current_A;     /* approved current of the previous ev_step() (0 after a FAULT) */
    float ramp_from_A;             /* level the charger may still be ramping down from */
    float ramp_age_s;              /* time since that ramp started */
    float over_current_time_s;     /* time the measured current has exceeded the allowance */
    uint16_t fault_hold_cycles;    /* > 0 while a fault is being held */
    uint16_t fault_flags;          /* EV_FAULT_* bits seen on the latest cycle */
    bool cold_charge_block;        /* below 0 degC (with margin); released 3 degC above */
} EvCellState;

/* Constraint bits (EvCellPrediction.violated, EvPackResult.limiting_constraints) */
#define EV_CONSTRAINT_SOC 0x01U
#define EV_CONSTRAINT_TEMPERATURE_HIGH 0x02U
#define EV_CONSTRAINT_TEMPERATURE_LOW 0x04U
#define EV_CONSTRAINT_ANODE_POTENTIAL 0x08U
#define EV_CONSTRAINT_CELL_VOLTAGE 0x10U
#define EV_CONSTRAINT_ELECTROLYTE_DEPLETION 0x20U
#define EV_CONSTRAINT_INVALID_INPUT 0x40U
#define EV_CONSTRAINT_CURRENT_LIMIT 0x80U /* delivered-current bound above max_charge_current_A */

/* Fault bits (EvCellState.fault_flags, EvPackResult.fault_flags) */
#define EV_FAULT_INPUT 0x01U          /* pack-level input: pointer, current, dt, horizon, request */
#define EV_FAULT_CONFIG 0x02U         /* EvCellConfig out of range */
#define EV_FAULT_VOLTAGE_SENSOR 0x04U
#define EV_FAULT_TEMPERATURE_SENSOR 0x08U
#define EV_FAULT_OVER_TEMPERATURE 0x10U
#define EV_FAULT_UNDER_TEMPERATURE 0x20U
#define EV_FAULT_STATE 0x40U          /* non-finite estimator state (reset applied) */
#define EV_FAULT_HOLD 0x80U           /* condition cleared, waiting EV_FAULT_HOLD_CYCLES */
#define EV_FAULT_OVER_VOLTAGE 0x100U  /* measured cell voltage above 4.2 V + margin */
#define EV_FAULT_VOLTAGE_PLAUSIBILITY 0x200U /* voltage inconsistent with the model (e.g. stuck sensor) */
#define EV_FAULT_END_OF_LIFE 0x400U   /* SOH reached EV_MIN_SOH: model no longer valid */
#define EV_FAULT_OVER_CURRENT 0x800U  /* measured charge current above the approved command */
#define EV_FAULT_SOC_DIVERGENCE 0x1000U /* EKF SOC left the coulomb count (a drifting, offset or frozen voltage reading) */
#define EV_FAULT_PACK_VOLTAGE 0x2000U /* sum of cell readings disagrees with the independent pack-voltage reading */

typedef uint8_t EvDecision;
#define EV_DECISION_ALLOW 0U          /* full request approved */
#define EV_DECISION_LIMIT 1U          /* approved < requested because a constraint binds */
#define EV_DECISION_FAULT 2U          /* approved = 0 because an input, config or cell is faulted */

#define EV_NO_CELL 0xFFU

typedef struct {
    EvDecision decision;
    float requested_charge_current_A;  /* request after clamping to [0, EV_MAX_PLAUSIBLE_CURRENT_A] */
    float approved_charge_current_A;
    /* Binding cell: on LIMIT the cell whose constraint stops a higher current, on ALLOW the cell
     * with the lowest anode potential, on FAULT the first faulted cell (EV_NO_CELL for pack-level
     * input faults). The predictions below are for that cell at the approved current. */
    uint8_t limiting_cell_index;
    uint8_t limiting_constraints;      /* EV_CONSTRAINT_* bits that bind (LIMIT only) */
    uint16_t fault_flags;              /* EV_FAULT_* bits, union over cells (FAULT only) */
    float limiting_anode_potential_V;
    float limiting_terminal_voltage_V;
    float limiting_temperature_C;
    float limiting_soc;                /* 3-sigma upper bound at the end of the horizon */
} EvPackResult;

typedef struct {
    float anode_potential_V;           /* graphite potential vs Li/Li+ at the end of the horizon */
    float terminal_voltage_V;          /* upper bound over the whole horizon */
    float predicted_temperature_C;
    float predicted_soc;               /* 3-sigma upper bound */
    uint8_t violated;                  /* EV_CONSTRAINT_* bits */
    bool safe;                         /* violated == 0 */
} EvCellPrediction;

/* Sets a fresh state for a NEW cell: given SOC (clamped to [0, 1]) with sigma 0.05, zero
 * polarisation, soh = 1 and no faults. Call it at commissioning only: across power cycles,
 * save the whole EvCellState to non-volatile memory with a CRC and restore it, so SOH, resistance
 * growth and latched faults survive a reboot. ev_step() checks a restored state's layout_version,
 * numeric ranges and flag bytes; a flag flipped between valid values needs the CRC to be caught.
 * Pass a SOC from a rested OCV reading; if it is unknown pass NaN
 * (sigma 0.5: charging stays blocked until the EKF has converged). */
void ev_cell_init(EvCellState *state, double initial_soc, float initial_temperature_C);

/* Service action after the faulty sensor has been repaired: clears the latched faults
 * (EV_FAULT_VOLTAGE_PLAUSIBILITY tracking latch, EV_FAULT_SOC_DIVERGENCE, EV_FAULT_PACK_VOLTAGE)
 * and keeps SOH, resistance growth and throughput. The SOC sigma is set to 0.5, so the EKF
 * re-learns the SOC from the repaired reading (at rest within a cycle) before charge is allowed
 * on it. Does nothing to an invalid state. */
void ev_clear_latched_faults(EvCellState *state);

/* True when every EvCellConfig field is finite and inside its documented range. */
bool ev_config_is_valid(const EvCellConfig *config);

/*
 * One control cycle:
 *  1. validate inputs; measured current and dt describe the interval that just elapsed, during
 *     which the previous call's approved current was in force (more than that => over-current);
 *  2. update every cell estimate (EKF, temperature, ageing) with those measurements;
 *  3. apply fault detection with hysteresis and hold;
 *  4. return the largest charge current <= request that keeps every cell inside all
 *     constraints for max(prediction_horizon_s, control_dt_s).
 * State is updated even on FAULT so that estimates keep tracking; when the current or dt
 * measurement itself is invalid the SOC mean is kept and its variance is inflated instead.
 */
EvPackResult ev_step(const EvCellConfig config[EV_NUM_CELLS], EvCellState state[EV_NUM_CELLS],
                     float measured_pack_current_A, float requested_charge_current_A,
                     const float measured_voltages_V[EV_NUM_CELLS], const float measured_temperatures_C[EV_NUM_CELLS],
                     float control_dt_s, float prediction_horizon_s);

/*
 * ev_step() plus a cross-check of the cell readings against an independent measurement of the
 * series-string voltage taken in the same cycle (a second AFE channel or the charger's output
 * voltage). When |sum of cell readings - measured_pack_voltage_V| exceeds pack_voltage_tolerance_V
 * for EV_TRACK_FAULT_CONFIRM_S (and at least 2 cycles), every cell latches EV_FAULT_PACK_VOLTAGE
 * until ev_clear_latched_faults() or ev_cell_init(); plain ev_step() keeps reporting the latch.
 * The expected pack reading is sum(cells) + current * pack_path_resistance_ohm (busbars and
 * connectors between the cell taps and the pack sense point, in [0, EV_MAX_PACK_PATH_RESISTANCE_OHM]).
 * A pack reading outside [EV_NUM_CELLS * EV_MIN_PLAUSIBLE_VOLTAGE_V, EV_NUM_CELLS *
 * EV_MAX_PLAUSIBLE_VOLTAGE_V] is "unavailable" (for example a charger output that is not live
 * yet): it never latches the pack-voltage fault, but once it is missing while |current| >= C/50
 * the charge is blocked (EV_FAULT_INPUT) until a valid pack reading is seen at rest, so the
 * check cannot be lost during a charge. The path resistance is trusted: a value set too high
 * hides a cell reading low in proportion to the current, up to current * the excess.
 * This is the check that catches a single cell channel that drifts slowly,
 * freezes during the constant-voltage taper or carries an offset smaller than the model-tracking
 * allowance; ev_step() alone cannot separate those from model error (see README, Limits).
 * Set the tolerance from both measurement accuracies and the path-resistance uncertainty times
 * the maximum current, in (0, EV_MAX_PACK_VOLTAGE_TOLERANCE_V]; the README's detection figures
 * hold for 30 mV or less. An out-of-range tolerance or path resistance is EV_FAULT_INPUT.
 */
EvPackResult ev_step_with_pack_voltage(const EvCellConfig config[EV_NUM_CELLS], EvCellState state[EV_NUM_CELLS],
                                       float measured_pack_current_A, float requested_charge_current_A,
                                       const float measured_voltages_V[EV_NUM_CELLS], const float measured_temperatures_C[EV_NUM_CELLS],
                                       float control_dt_s, float prediction_horizon_s,
                                       float measured_pack_voltage_V, float pack_voltage_tolerance_V,
                                       float pack_path_resistance_ohm);

/* Predicts one cell at a constant charge current. Returns false (and safe = false,
 * violated = EV_CONSTRAINT_INVALID_INPUT) for invalid arguments. */
bool ev_predict_cell(const EvCellConfig *config, const EvCellState *state, float charge_current_A,
                     float horizon_s, EvCellPrediction *prediction);

/* Largest charge current in [0, requested_A] that keeps one cell inside all constraints for
 * horizon_s (same search as ev_step). Writes the binding EV_CONSTRAINT_* bits (0 when the full
 * request is safe) and the prediction at the returned current. Returns 0 for invalid arguments. */
float ev_cell_max_charge_current(const EvCellConfig *config, const EvCellState *state, float requested_A,
                                 float horizon_s, uint8_t *binding_constraints, EvCellPrediction *at_limit);

/* Full-cell OCV [V] and dOCV/dSOC [V] at a SOC (clamped to [0, 1]). Non-finite SOC gives NaN.
 * Either output pointer may be NULL. */
void ev_full_cell_ocv_V_and_dv_dsoc(float soc, float *ocv, float *dv_dsoc);

#ifdef __cplusplus
}
#endif

#endif /* EV_MODIFY_H */
