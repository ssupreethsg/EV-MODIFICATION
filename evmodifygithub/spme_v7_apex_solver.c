/*
 * Copyright (c) 2026. All Rights Reserved.
 * PROPRIETARY AND CONFIDENTIAL
 */

#include <stdint.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>

#include "spme_v7_apex_solver.h"

#ifndef STATIC
#define STATIC static
#endif

const float LUT_GRAPHITE_OCP_101[101] = {
    1.06374f, 0.87128f, 0.73440f, 0.63668f, 0.56629f, 0.51459f, 0.47510f, 0.44281f, 0.41375f, 0.38512f, 0.35565f, 0.32600f, 0.29828f, 0.27474f, 0.25649f, 0.24327f, 0.23404f, 0.22760f, 0.22294f, 0.21927f, 0.21605f, 0.21291f, 0.20960f, 0.20593f, 0.20181f, 0.19718f, 0.19208f, 0.18658f, 0.18082f, 0.17497f, 0.16922f, 0.16376f, 0.15873f, 0.15423f, 0.15030f, 0.14695f, 0.14415f, 0.14185f, 0.13998f, 0.13848f, 0.13729f, 0.13635f, 0.13560f, 0.13502f, 0.13457f, 0.13422f, 0.13395f, 0.13374f, 0.13357f, 0.13344f, 0.13334f, 0.13326f, 0.13320f, 0.13314f, 0.13309f, 0.13303f, 0.13295f, 0.13284f, 0.13266f, 0.13237f, 0.13190f, 0.13113f, 0.12988f, 0.12794f, 0.12506f, 0.12110f, 0.11616f, 0.11075f, 0.10558f, 0.10125f, 0.09800f, 0.09576f, 0.09430f, 0.09339f, 0.09284f, 0.09250f, 0.09231f, 0.09219f, 0.09212f, 0.09208f, 0.09205f, 0.09204f, 0.09203f, 0.09203f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f
};
const float LUT_NMC_OCP_101[101] = {
    3.54498f, 3.55049f, 3.55600f, 3.56150f, 3.56701f, 3.57252f, 3.57803f, 3.58354f, 3.58904f, 3.59455f, 3.60006f, 3.60557f, 3.61108f, 3.61658f, 3.62209f, 3.62760f, 3.63311f, 3.63862f, 3.64413f, 3.64964f, 3.65515f, 3.66066f, 3.66617f, 3.67169f, 3.67720f, 3.68272f, 3.68824f, 3.69377f, 3.69930f, 3.70483f, 3.71038f, 3.71593f, 3.72150f, 3.72708f, 3.73269f, 3.73832f, 3.74399f, 3.74970f, 3.75547f, 3.76132f, 3.76725f, 3.77330f, 3.77950f, 3.78588f, 3.79248f, 3.79935f, 3.80654f, 3.81412f, 3.82213f, 3.83064f, 3.83967f, 3.84923f, 3.85931f, 3.86983f, 3.88069f, 3.89175f, 3.90285f, 3.91383f, 3.92455f, 3.93491f, 3.94485f, 3.95437f, 3.96347f, 3.97222f, 3.98067f, 3.98892f, 3.99704f, 4.00511f, 4.01322f, 4.02144f, 4.02983f, 4.03847f, 4.04740f, 4.05668f, 4.06633f, 4.07637f, 4.08678f, 4.09753f, 4.10852f, 4.11962f, 4.13064f, 4.14133f, 4.15141f, 4.16056f, 4.16850f, 4.17500f, 4.18001f, 4.18365f, 4.18633f, 4.18870f, 4.19159f, 4.19591f, 4.20247f, 4.21187f, 4.22436f, 4.23986f, 4.25796f, 4.27805f, 4.29944f, 4.32140f, 4.34332f
};


STATIC float ev_clamp(float value, float lower, float upper) {
    if (value < lower) return lower;
    if (value > upper) return upper;
    return value;
}

STATIC double ev_clamp_d(double value, double lower, double upper) {
    if (value < lower) return lower;
    if (value > upper) return upper;
    return value;
}

STATIC bool ev_is_valid_cell(const EvCellState *cell) {
    return (cell != NULL) && (cell->soc >= 0.0 && cell->soc <= 1.0)
        && (cell->soh > 0.0f) && (cell->soh <= 1.0f)
        && (cell->nominal_capacity_Ah > 0.0f)
        && (cell->resistance_ohm >= 0.0f) && (cell->soc_variance >= 0.0)
        && (cell->temperature_C >= EV_MIN_CELL_TEMPERATURE_C)
        && (cell->temperature_C <= EV_MAX_CELL_TEMPERATURE_C)
        && (cell->sensor_valid == true);
}

void ev_full_cell_ocv_V_and_dv_dsoc(float soc, float *ocv, float *dv_dsoc) {
    float bounded_soc = ev_clamp(soc, 0.0f, 1.0f);
    float f_idx = bounded_soc * 100.0f;
    uint8_t idx = (uint8_t)f_idx;
    if (idx >= 100U) { idx = 99U; }
    float t = f_idx - (float)idx;
    
    float graphite_val = LUT_GRAPHITE_OCP_101[idx] + t * (LUT_GRAPHITE_OCP_101[idx+1] - LUT_GRAPHITE_OCP_101[idx]);
    float nmc_val = LUT_NMC_OCP_101[idx] + t * (LUT_NMC_OCP_101[idx+1] - LUT_NMC_OCP_101[idx]);
    
    float graphite_slope = (LUT_GRAPHITE_OCP_101[idx+1] - LUT_GRAPHITE_OCP_101[idx]) * 100.0f;
    float nmc_slope = (LUT_NMC_OCP_101[idx+1] - LUT_NMC_OCP_101[idx]) * 100.0f;
    
    if (ocv != NULL) { *ocv = nmc_val - graphite_val; }
    if (dv_dsoc != NULL) { *dv_dsoc = nmc_slope - graphite_slope; }
}

STATIC float ev_full_cell_ocv_V(float soc) {
    float ocv = 0.0f;
    ev_full_cell_ocv_V_and_dv_dsoc(soc, &ocv, NULL);
    return ocv;
}

STATIC float ev_graphite_ocp_V(float soc) {
    float bounded_soc = ev_clamp(soc, 0.0f, 1.0f);
    float f_idx = bounded_soc * 100.0f;
    uint8_t idx = (uint8_t)f_idx;
    if (idx >= 100U) { idx = 99U; }
    float t = f_idx - (float)idx;
    return LUT_GRAPHITE_OCP_101[idx] + t * (LUT_GRAPHITE_OCP_101[idx+1] - LUT_GRAPHITE_OCP_101[idx]);
}

STATIC float ev_diffusivity_m2_s(float temperature_K) {
    float inv_T = 1.0f / ev_clamp(temperature_K, 243.15f, 333.15f);
    float exponent = (-EV_DIFFUSION_ACTIVATION_ENERGY_J_MOL / EV_GAS_J_PER_MOL_K) * (inv_T - (1.0f / EV_REFERENCE_TEMPERATURE_K));
    return EV_SOLID_DIFFUSIVITY_REF_M2_S * expf(ev_clamp(exponent, -20.0f, 20.0f));
}

STATIC float ev_predict_temperature_C(const EvCellState *cell, float charge_current_A, float dt_s) {
    float thermal_conductance_W_K = cell->heat_transfer_W_m2_K * EV_CELL_COOLING_AREA_M2;
    float heat_capacity_J_K = EV_CELL_MASS_KG * cell->heat_capacity_J_kg_K;
    float joule_heat_W = charge_current_A * charge_current_A * cell->resistance_ohm;
    float temperature_rise_steady_C = joule_heat_W / thermal_conductance_W_K;
    float exponent = -(thermal_conductance_W_K * dt_s) / heat_capacity_J_K;
    float decay = expf(ev_clamp(exponent, -20.0f, 0.0f));

    return EV_AMBIENT_TEMPERATURE_C + temperature_rise_steady_C
        + (cell->temperature_C - EV_AMBIENT_TEMPERATURE_C - temperature_rise_steady_C) * decay;
}

bool ev_predict_cell(const EvCellState *cell, float charge_current_A, float horizon_s, EvCellPrediction *prediction) {
    prediction->anode_potential_V = 0.0f;
    prediction->predicted_temperature_C = 0.0f;
    prediction->predicted_soc = 0.0f;
    prediction->safe = false;
    prediction->error_code = 1;

    if (!ev_is_valid_cell(cell) || (charge_current_A < 0.0f) || (horizon_s <= 0.0f)) {
        return false;
    }
    
    float effective_capacity_Ah = cell->nominal_capacity_Ah * cell->soh;

    prediction->predicted_soc = (float)(cell->soc + ((double)charge_current_A * (double)horizon_s) / (effective_capacity_Ah * 3600.0));
    prediction->predicted_temperature_C = ev_predict_temperature_C(cell, charge_current_A, horizon_s);
    
    float temperature_K = (prediction->predicted_temperature_C - 2.0f) + 273.15f;
    if (temperature_K < 243.15f || temperature_K > 333.15f) {
        return false;
    }

    float diffusivity_m2_s = ev_diffusivity_m2_s(temperature_K);
    float geometric_current_density_A_m2 = charge_current_A / EV_GEOMETRIC_ELECTRODE_AREA_M2;
    float specific_surface_area_m_inv = (3.0f * cell->active_solid_fraction) / EV_PARTICLE_RADIUS_M;
    float molar_flux_mol_m2_s = -geometric_current_density_A_m2 / (specific_surface_area_m_inv * EV_ANODE_THICKNESS_M * EV_FARADAY_C_PER_MOL);
    float surface_soc = prediction->predicted_soc - (molar_flux_mol_m2_s * EV_PARTICLE_RADIUS_M) / (5.0f * diffusivity_m2_s * EV_MAX_CONCENTRATION_MOL_M3);
    surface_soc = ev_clamp(surface_soc, 0.001f, 0.999f);
    
    /* Strictly bounded math */
    float exchange_current_A_m2 = EV_EXCHANGE_CURRENT_REF_A_M2 * sqrtf(ev_clamp(surface_soc * (1.0f - surface_soc), 0.0f, 1.0f));
    float inv_T = 1.0f / temperature_K;
    float kin_exponent = (-EV_KINETIC_ACTIVATION_ENERGY_J_MOL / EV_GAS_J_PER_MOL_K) * (inv_T - (1.0f / EV_REFERENCE_TEMPERATURE_K));
    exchange_current_A_m2 *= expf(ev_clamp(kin_exponent, -20.0f, 20.0f));
    
    float effective_conductivity_S_m = EV_ELECTROLYTE_CONDUCTIVITY_S_M * powf(EV_ELECTROLYTE_POROSITY, 1.5f);
    float elec_exponent = (-EV_ELECTROLYTE_ACTIVATION_ENERGY_J_MOL / EV_GAS_J_PER_MOL_K) * (inv_T - (1.0f / EV_REFERENCE_TEMPERATURE_K));
    effective_conductivity_S_m *= expf(ev_clamp(elec_exponent, -20.0f, 20.0f));
    
    float diff_exponent = (-EV_ELECTROLYTE_DIFFUSIVITY_ACTIVATION_J_MOL / EV_GAS_J_PER_MOL_K) * (inv_T - (1.0f / EV_REFERENCE_TEMPERATURE_K));
    float electrolyte_diffusivity = EV_ELECTROLYTE_DIFFUSIVITY_REF_M2_S * expf(ev_clamp(diff_exponent, -20.0f, 20.0f));
    float concentration_gradient_factor = (1.0f - EV_TRANSFERENCE_NUMBER) / (EV_FARADAY_C_PER_MOL * electrolyte_diffusivity);
    float argument = (geometric_current_density_A_m2 * EV_ANODE_THICKNESS_M * concentration_gradient_factor) / 1000.0f;
    
    /* Safely bound logarithm to avoid PDE divergence or C0 discontinuities without tanhf duct-tape */
    if (argument <= -0.99f) {
        prediction->error_code = 2; /* Lithium depletion limit */
        return false;
    }
    float ce_polarization_V = (2.0f * EV_GAS_J_PER_MOL_K * temperature_K / EV_FARADAY_C_PER_MOL) * log1pf(argument);
    
    float asinh_arg = (molar_flux_mol_m2_s * EV_FARADAY_C_PER_MOL) / (2.0f * exchange_current_A_m2);
    float eta_kinetic_V = (2.0f * EV_GAS_J_PER_MOL_K * temperature_K / EV_FARADAY_C_PER_MOL) * asinhf(ev_clamp(asinh_arg, -100.0f, 100.0f));
    float electrolyte_ohmic_drop_V = -(geometric_current_density_A_m2 * EV_ANODE_THICKNESS_M) / (2.0f * effective_conductivity_S_m);
    
    prediction->anode_potential_V = ev_graphite_ocp_V(surface_soc) + eta_kinetic_V + electrolyte_ohmic_drop_V - ce_polarization_V;
    prediction->safe = (prediction->predicted_soc <= EV_MAX_SAFE_SOC)
        && (prediction->predicted_temperature_C <= EV_MAX_CELL_TEMPERATURE_C)
        && (prediction->anode_potential_V >= EV_MIN_SAFE_ANODE_POTENTIAL_V);
    prediction->error_code = 0;
    return true;
}

STATIC bool ev_current_is_safe(const EvCellState *cells, float current_A, float horizon_s, uint8_t *limiting_index, EvCellPrediction *limiting_prediction) {
    uint8_t index;
    EvCellPrediction worst = {FLT_MAX, -FLT_MAX, -FLT_MAX, true, 0};
    uint8_t worst_index = 0U;

    for (index = 0U; index < EV_NUM_CELLS; ++index) {
        EvCellPrediction prediction;
        if (!ev_predict_cell(&cells[index], current_A, horizon_s, &prediction) || !prediction.safe) {
            if (limiting_index != NULL) { *limiting_index = index; }
            if (limiting_prediction != NULL) { *limiting_prediction = prediction; }
            return false;
        }
        if (prediction.anode_potential_V < worst.anode_potential_V) {
            worst = prediction;
            worst_index = index;
        }
    }
    if (limiting_index != NULL) { *limiting_index = worst_index; }
    if (limiting_prediction != NULL) { *limiting_prediction = worst; }
    return true;
}

STATIC float ev_find_safe_charge_current(const EvCellState *cells, float requested_current_A, float horizon_s, uint8_t *limiting_index, EvCellPrediction *limiting_prediction) {
    /* Replace ticking time bomb bisection with a deterministic linear sweep upwards. 
       Starts from 0 and increases until it hits a limit. This is mathematically safe 
       even if the temperature response is non-monotonic. */
    float safe_current_A = 0.0f;
    float candidate_current_A = 0.0f;
    float step_A = requested_current_A / 20.0f;
    if (step_A < 0.1f) step_A = 0.1f;
    
    if (!ev_current_is_safe(cells, 0.0f, horizon_s, limiting_index, limiting_prediction)) {
        return 0.0f;
    }
    
    while (candidate_current_A <= requested_current_A) {
        if (ev_current_is_safe(cells, candidate_current_A, horizon_s, limiting_index, limiting_prediction)) {
            safe_current_A = candidate_current_A;
            candidate_current_A += step_A;
        } else {
            break;
        }
    }
    
    (void)ev_current_is_safe(cells, safe_current_A, horizon_s, limiting_index, limiting_prediction);
    return safe_current_A;
}

STATIC void ev_ekf_and_sei_update(EvCellState *cell, float measured_current_A, float measured_voltage_V, float measured_temperature_C, float dt_s) {
    float Ah_processed = (fabsf(measured_current_A) * dt_s) / 3600.0f;
    cell->throughput_Ah += Ah_processed;
    
    float temp_C = ev_clamp(cell->temperature_C, EV_MIN_CELL_TEMPERATURE_C, EV_MAX_CELL_TEMPERATURE_C);
    float aging_exponent = temp_C / 25.0f;
    float cycling_loss = EV_SEI_GROWTH_RATE_CONSTANT * Ah_processed * expf(ev_clamp(aging_exponent, -5.0f, 5.0f));
    float calendric_loss = EV_CALENDRIC_AGING_RATE * dt_s * expf(ev_clamp(aging_exponent, -5.0f, 5.0f));
    float incremental_capacity_loss = cycling_loss + calendric_loss;
    cell->soh = ev_clamp(cell->soh - incremental_capacity_loss, 0.5f, 1.0f);
    cell->resistance_ohm += EV_SEI_RESISTANCE_FACTOR * incremental_capacity_loss;

    const double Q_process_noise = 1.0e-7; 
    const double R_meas_noise = 1.0e-2; 
    
    const double effective_capacity_Ah = (double)(cell->nominal_capacity_Ah * cell->soh);
    double soc_pred = cell->soc + ((double)measured_current_A * (double)dt_s) / (effective_capacity_Ah * 3600.0);
    soc_pred = ev_clamp_d(soc_pred, 0.0, 1.0);
    
    double P_pred = cell->soc_variance + Q_process_noise;
    
    float ocv, dv_dsoc;
    ev_full_cell_ocv_V_and_dv_dsoc((float)soc_pred, &ocv, &dv_dsoc);
    
    /* Analytical Jacobian: No floating point truncation errors */
    double H = (double)dv_dsoc;
    double V_pred = (double)ocv + (double)(measured_current_A * cell->resistance_ohm);
    
    double S = H * P_pred * H + R_meas_noise;
    double K = (P_pred * H) / S;
    double innovation = (double)measured_voltage_V - V_pred;
    
    /* S is strictly >= R_meas_noise (1e-2) by construction, no need for arbitrary 0.0001 check */
    cell->soc = ev_clamp_d(soc_pred + K * innovation, 0.0, 1.0);
    double I_KH = 1.0 - K * H;
    cell->soc_variance = I_KH * P_pred * I_KH + K * R_meas_noise * K;
    
    if (measured_temperature_C > -100.0f && measured_temperature_C < 200.0f) {
        cell->temperature_C = ev_clamp(measured_temperature_C, EV_MIN_CELL_TEMPERATURE_C, EV_MAX_CELL_TEMPERATURE_C);
    } else {
        cell->temperature_C = ev_clamp(ev_predict_temperature_C(cell, measured_current_A, dt_s), EV_MIN_CELL_TEMPERATURE_C, EV_MAX_CELL_TEMPERATURE_C);
    }
}

EvPackResult ev_step(EvCellState *cells, float measured_pack_current_A, float requested_charge_current_A, const float measured_voltages_V[EV_NUM_CELLS], const float measured_temperatures_C[EV_NUM_CELLS], float control_dt_s, float prediction_horizon_s) {
    EvPackResult result = {EV_DECISION_FAULT, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0U};
    EvCellPrediction limiting_prediction = {0.0f, 0.0f, 0.0f, false, 0};
    uint8_t limiting_index = 0U;
    uint8_t index;

    if ((cells == NULL) || (measured_voltages_V == NULL) || (control_dt_s <= 0.0f) || (prediction_horizon_s <= 0.0f)) {
        return result;
    }
    for (index = 0U; index < EV_NUM_CELLS; ++index) {
        if (!ev_is_valid_cell(&cells[index])) {
            return result;
        }
    }

    requested_charge_current_A = requested_charge_current_A > 0.0f ? requested_charge_current_A : 0.0f;
    result.requested_charge_current_A = requested_charge_current_A;
    result.approved_charge_current_A = ev_find_safe_charge_current(cells, requested_charge_current_A, prediction_horizon_s, &limiting_index, &limiting_prediction);
    result.limiting_cell_index = limiting_index;
    result.limiting_anode_potential_V = limiting_prediction.anode_potential_V;
    result.limiting_temperature_C = limiting_prediction.predicted_temperature_C;
    result.limiting_soc = limiting_prediction.predicted_soc;

    if (result.approved_charge_current_A < 0.0f) {
        return result;
    }
    result.decision = ((result.approved_charge_current_A + EV_CURRENT_TOLERANCE_A) < requested_charge_current_A) ? EV_DECISION_LIMIT : EV_DECISION_ALLOW;
        
    for (index = 0U; index < EV_NUM_CELLS; ++index) {
        ev_ekf_and_sei_update(&cells[index], measured_pack_current_A, measured_voltages_V[index], measured_temperatures_C[index], control_dt_s);
    }
    return result;
}
