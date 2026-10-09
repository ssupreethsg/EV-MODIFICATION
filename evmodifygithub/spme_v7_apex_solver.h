#ifndef SPME_V7_APEX_SOLVER_H
#define SPME_V7_APEX_SOLVER_H

#include <stdint.h>
#include <stdbool.h>

/* Sign convention: positive current = charging. */

#define EV_NUM_CELLS 4U

/* Physical constants */
#define EV_FARADAY_C_PER_MOL 96485.33212f
#define EV_GAS_J_PER_MOL_K 8.31446262f
#define EV_KELVIN_OFFSET 273.15f
#define EV_REFERENCE_TEMPERATURE_K 298.15f

/* Operating envelope and safety limits */
#define EV_AMBIENT_TEMPERATURE_C 25.0f
#define EV_MAX_CELL_TEMPERATURE_C 45.0f
#define EV_MIN_CELL_TEMPERATURE_C -30.0f
#define EV_TEMPERATURE_SENSOR_MARGIN_C 2.0f   /* cold bias applied to the plating check */
#define EV_MIN_SAFE_ANODE_POTENTIAL_V 0.005f
#define EV_MAX_SAFE_SOC 0.980f
#define EV_CURRENT_TOLERANCE_A 1.0e-3f

/* Sensor plausibility windows (outside => sensor fault) */
#define EV_MAX_PLAUSIBLE_CURRENT_A 1000.0f
#define EV_MIN_PLAUSIBLE_VOLTAGE_V 2.0f
#define EV_MAX_PLAUSIBLE_VOLTAGE_V 4.5f
#define EV_MIN_PLAUSIBLE_TEMPERATURE_C -40.0f
#define EV_MAX_PLAUSIBLE_TEMPERATURE_C 100.0f

/* Negative electrode (graphite), Chen2020 / LG M50 parameter set */
#define EV_ANODE_STOICH_AT_SOC0 0.0279f
#define EV_ANODE_STOICH_AT_SOC100 0.9014f
#define EV_ANODE_THICKNESS_M 85.2e-6f
#define EV_PARTICLE_RADIUS_M 5.86e-6f
#define EV_MAX_CONCENTRATION_MOL_M3 33133.0f
#define EV_SOLID_DIFFUSIVITY_REF_M2_S 3.3e-14f
#define EV_DIFFUSION_ACTIVATION_ENERGY_J_MOL 35000.0f   /* calibration placeholder */
#define EV_GEOMETRIC_ELECTRODE_AREA_M2 0.1027f
/* i0 = k * sqrt(theta * (1 - theta)), k = m_ref * sqrt(c_e) * c_s_max (Chen2020) */
#define EV_EXCHANGE_CURRENT_REF_A_M2 0.679f
#define EV_KINETIC_ACTIVATION_ENERGY_J_MOL 35000.0f

/* Electrolyte (LiPF6 EC:EMC, Nyman2008 at 1 mol/L) */
#define EV_ELECTROLYTE_CONCENTRATION_MOL_M3 1000.0f
#define EV_ELECTROLYTE_POROSITY 0.25f
#define EV_BRUGGEMAN_EXPONENT 1.5f
#define EV_TRANSFERENCE_NUMBER 0.2594f
#define EV_ELECTROLYTE_CONDUCTIVITY_S_M 0.9487f
#define EV_ELECTROLYTE_ACTIVATION_ENERGY_J_MOL 15000.0f       /* calibration placeholder */
#define EV_ELECTROLYTE_DIFFUSIVITY_REF_M2_S 1.769e-10f
#define EV_ELECTROLYTE_DIFFUSIVITY_ACTIVATION_J_MOL 10000.0f  /* calibration placeholder */

/* Lumped thermal model */
#define EV_CELL_MASS_KG 0.070f
#define EV_CELL_COOLING_AREA_M2 0.0075f

/* Ageing (calibration placeholders, Arrhenius-scaled to 1 at 25 degC) */
#define EV_SEI_GROWTH_RATE_PER_AH 1.5e-8f
#define EV_CALENDRIC_AGING_RATE_PER_S 2.0e-10f
#define EV_AGING_ACTIVATION_ENERGY_J_MOL 30000.0f
#define EV_SEI_RESISTANCE_FACTOR_OHM 0.015f   /* resistance rise per unit capacity loss */
#define EV_MIN_SOH 0.5f

typedef uint8_t EvDecision;
#define EV_DECISION_ALLOW 0U
#define EV_DECISION_LIMIT 1U
#define EV_DECISION_FAULT 2U

typedef struct {
    double soc;                   /* 0..1, EKF estimate */
    double soc_variance;          /* EKF covariance */
    float temperature_C;          /* latest plausible cell temperature */
    float soh;                    /* 0..1 capacity-retention factor */
    float throughput_Ah;          /* total charge throughput */
    float nominal_capacity_Ah;
    float resistance_ohm;
    bool sensor_valid;            /* refreshed by ev_step() from measurement plausibility */
    float heat_transfer_W_m2_K;
    float heat_capacity_J_kg_K;
    float active_solid_fraction;
} EvCellState;

typedef struct {
    EvDecision decision;
    float requested_charge_current_A;
    float approved_charge_current_A;
    float limiting_anode_potential_V;
    float limiting_temperature_C;
    float limiting_soc;
    uint8_t limiting_cell_index;
} EvPackResult;

#define EV_PREDICTION_OK 0
#define EV_PREDICTION_INVALID_INPUT 1
#define EV_PREDICTION_ELECTROLYTE_DEPLETION 2

typedef struct {
    float anode_potential_V;
    float predicted_temperature_C;
    float predicted_soc;
    bool safe;
    int32_t error_code;
} EvCellPrediction;

/*
 * One control cycle: update every cell estimate with the measurements taken over the
 * last control_dt_s, then find the largest charge current <= request that keeps every
 * cell inside its limits over prediction_horizon_s.
 */
EvPackResult ev_step(EvCellState *cells, float measured_pack_current_A, float requested_charge_current_A,
                     const float measured_voltages_V[EV_NUM_CELLS], const float measured_temperatures_C[EV_NUM_CELLS],
                     float control_dt_s, float prediction_horizon_s);

bool ev_predict_cell(const EvCellState *cell, float charge_current_A, float horizon_s, EvCellPrediction *prediction);
void ev_full_cell_ocv_V_and_dv_dsoc(float soc, float *ocv, float *dv_dsoc);

#endif /* SPME_V7_APEX_SOLVER_H */
