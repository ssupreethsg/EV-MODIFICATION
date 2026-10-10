/*
 * EV MODIFY - conservative charge-current limiter for a series lithium-ion module.
 *
 * Copyright (c) 2026 Supreeth S G. All rights reserved. See LICENSE.
 *
 * Usage (one call per control cycle):
 *
 *     EvCellConfig config[EV_NUM_CELLS];   // filled once by the integrator
 *     EvCellState  state[EV_NUM_CELLS];
 *     for (i = 0; i < EV_NUM_CELLS; ++i) { ev_cell_init(&state[i], initial_soc[i], initial_temp_C[i]); }
 *     ...
 *     EvPackResult r = ev_step(config, state, pack_current_A, requested_A,
 *                              cell_voltages_V, cell_temperatures_C, dt_s, horizon_s);
 *     command_charger(r.approved_charge_current_A);   // 0 A on EV_DECISION_FAULT
 *
 * Sign convention: positive current = charging. Units are SI unless the name says otherwise.
 * The library is single-threaded and keeps no global state; all state lives in EvCellState.
 */
#ifndef SPME_V7_APEX_SOLVER_H
#define SPME_V7_APEX_SOLVER_H

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

#define EV_NUM_CELLS 4U

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
#define EV_CURRENT_TOLERANCE_A 1.0e-3f

/* ---- Input plausibility (outside => FAULT) ----------------------------------------------- */
#define EV_MAX_PLAUSIBLE_CURRENT_A 1000.0f   /* |measured current|; requests are clamped to it */
#define EV_MIN_PLAUSIBLE_VOLTAGE_V 2.0f
#define EV_MAX_PLAUSIBLE_VOLTAGE_V 4.5f
#define EV_MIN_PLAUSIBLE_TEMPERATURE_C -40.0f
#define EV_MAX_PLAUSIBLE_TEMPERATURE_C 100.0f
#define EV_MAX_CONTROL_DT_S 10.0f
#define EV_MAX_PREDICTION_HORIZON_S 600.0f

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
#define EV_KINETIC_ACTIVATION_ENERGY_J_MOL 35000.0f      /* Chen2020 */
#define EV_DIFFUSION_ACTIVATION_ENERGY_J_MOL 35000.0f    /* placeholder (Chen2020 gives none) */

/* Electrolyte: LiPF6 in EC:EMC, Nyman2008 values at 1 mol/L */
#define EV_ELECTROLYTE_CONCENTRATION_MOL_M3 1000.0f
#define EV_ELECTROLYTE_POROSITY 0.25f
#define EV_BRUGGEMAN_EXPONENT 1.5f
#define EV_TRANSFERENCE_NUMBER 0.2594f
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

/* Ageing: linear placeholder rates of typical NMC/graphite magnitude at 25 degC (about 2 %
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
} EvCellConfig;

/* Per-cell persistent state. Initialise with ev_cell_init(); afterwards only ev_step() writes it. */
typedef struct {
    double soc;                    /* [0, 1] EKF estimate */
    double covariance[3];          /* EKF covariance of [soc, v_rc]: P00, P01, P11 */
    double soh;                    /* [EV_MIN_SOH, 1] capacity retention, never increases */
    double throughput_Ah;          /* total |charge| processed */
    float rc_voltage_V;            /* EKF estimate of the polarisation voltage */
    float temperature_C;           /* latest plausible measurement (model value if the sensor fails) */
    float resistance_growth_ohm;   /* [0, 1] ageing-related R0 increase, treated as anode SEI film */
    float last_voltage_V;          /* latest plausible measured terminal voltage (0 = none yet) */
    float last_current_A;          /* current of the latest cycle */
    uint16_t innovation_fail_cycles; /* consecutive EKF innovations outside the gate */
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

/* Sets a fresh state: given SOC (clamped to [0, 1]) with sigma 0.05, zero polarisation,
 * soh = 1 and no faults. Pass a SOC from a rested OCV reading; if it is unknown pass NaN
 * (sigma 0.5: charging stays blocked until the EKF has converged). */
void ev_cell_init(EvCellState *state, double initial_soc, float initial_temperature_C);

/* True when every EvCellConfig field is finite and inside its documented range. */
bool ev_config_is_valid(const EvCellConfig *config);

/*
 * One control cycle:
 *  1. validate inputs; measured current and dt describe the interval that just elapsed;
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

#endif /* SPME_V7_APEX_SOLVER_H */
