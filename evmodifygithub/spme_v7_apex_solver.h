#ifndef SPME_V7_APEX_SOLVER_H
#define SPME_V7_APEX_SOLVER_H

#include <stdint.h>
#include <stdbool.h>

#define EV_NUM_CELLS 4U
#define EV_FARADAY_C_PER_MOL 96485.33212f
#define EV_GAS_J_PER_MOL_K 8.31446262f
#define EV_REFERENCE_TEMPERATURE_K 298.15f
#define EV_AMBIENT_TEMPERATURE_C 25.0f
#define EV_MAX_CELL_TEMPERATURE_C 45.0f
#define EV_MIN_CELL_TEMPERATURE_C -30.0f
#define EV_MIN_SAFE_ANODE_POTENTIAL_V 0.005f
#define EV_MAX_SAFE_SOC 0.980f\n#define EV_CURRENT_TOLERANCE_A 1.0e-3f
#define EV_ELECTROLYTE_ACTIVATION_ENERGY_J_MOL 15000.0f
#define EV_TRANSFERENCE_NUMBER 0.36f
#define EV_ELECTROLYTE_DIFFUSIVITY_REF_M2_S 2.0e-10f
#define EV_ELECTROLYTE_DIFFUSIVITY_ACTIVATION_J_MOL 10000.0f
#define EV_CALENDRIC_AGING_RATE 2.0e-10f
#define EV_SOLID_DIFFUSIVITY_REF_M2_S 3.9e-14f
#define EV_DIFFUSION_ACTIVATION_ENERGY_J_MOL 35000.0f
#define EV_KINETIC_ACTIVATION_ENERGY_J_MOL 30000.0f
#define EV_MAX_CONCENTRATION_MOL_M3 33133.0f
#define EV_ANODE_THICKNESS_M 75.0e-6f

#define EV_ELECTROLYTE_POROSITY 0.25f
#define EV_ELECTROLYTE_CONDUCTIVITY_S_M 1.10f
#define EV_PARTICLE_RADIUS_M 5.8e-6f
#define EV_GEOMETRIC_ELECTRODE_AREA_M2 0.100f
#define EV_EXCHANGE_CURRENT_REF_A_M2 1.5f
#define EV_CELL_MASS_KG 0.070f
#define EV_CELL_COOLING_AREA_M2 0.0075f
#define EV_SEI_GROWTH_RATE_CONSTANT 1.5e-8f
#define EV_SEI_RESISTANCE_FACTOR 0.015f
typedef uint8_t EvDecision;
#define EV_DECISION_ALLOW 0U
#define EV_DECISION_LIMIT 1U
#define EV_DECISION_FAULT 2U
typedef struct {
    double soc;
    double soc_variance;
    float temperature_C;          /* persisted lumped temperature state */
    float soh;                    /* 0..1 capacity-retention factor (Degradation Tracking) */
    float throughput_Ah;          /* Total charge throughput for SEI aging */
    float nominal_capacity_Ah;
    float resistance_ohm;
    bool sensor_valid;
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
typedef struct {
    float anode_potential_V;
    float predicted_temperature_C;
    float predicted_soc;
    bool safe;
    int32_t error_code;
} EvCellPrediction;
EvPackResult ev_step(EvCellState *cells, float measured_pack_current_A, float requested_charge_current_A, const float measured_voltages_V[EV_NUM_CELLS], const float measured_temperatures_C[EV_NUM_CELLS], float control_dt_s, float prediction_horizon_s);

bool ev_predict_cell(const EvCellState *cell, float charge_current_A, float horizon_s, EvCellPrediction *prediction);
void ev_full_cell_ocv_V_and_dv_dsoc(float soc, float *ocv, float *dv_dsoc);

#endif // SPME_V7_APEX_SOLVER_H
