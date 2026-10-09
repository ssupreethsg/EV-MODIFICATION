/*
 * Copyright (c) 2026. All Rights Reserved.
 * PROPRIETARY AND CONFIDENTIAL
 */

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "spme_v7_apex_solver.h"

#define EV_LUT_SEGMENTS 100U
#define EV_COARSE_SEARCH_STEPS 20U
#define EV_BISECTION_STEPS 24U
#define EV_SOC_PROCESS_NOISE_PER_S 1.0e-7
#define EV_VOLTAGE_MEASUREMENT_NOISE_V2 1.0e-2

/* Electrode OCPs from PyBaMM Chen2020, tabulated against cell SOC (0..1, 101 points)
 * via the stoichiometry windows. Checked by test_ocv.py. */
static const float LUT_GRAPHITE_OCP_101[101] = {
    1.06374f, 0.87128f, 0.73440f, 0.63668f, 0.56629f, 0.51459f, 0.47510f, 0.44281f, 0.41375f, 0.38512f, 0.35565f, 0.32600f, 0.29828f, 0.27474f, 0.25649f, 0.24327f, 0.23404f, 0.22760f, 0.22294f, 0.21927f, 0.21605f, 0.21291f, 0.20960f, 0.20593f, 0.20181f, 0.19718f, 0.19208f, 0.18658f, 0.18082f, 0.17497f, 0.16922f, 0.16376f, 0.15873f, 0.15423f, 0.15030f, 0.14695f, 0.14415f, 0.14185f, 0.13998f, 0.13848f, 0.13729f, 0.13635f, 0.13560f, 0.13502f, 0.13457f, 0.13422f, 0.13395f, 0.13374f, 0.13357f, 0.13344f, 0.13334f, 0.13326f, 0.13320f, 0.13314f, 0.13309f, 0.13303f, 0.13295f, 0.13284f, 0.13266f, 0.13237f, 0.13190f, 0.13113f, 0.12988f, 0.12794f, 0.12506f, 0.12110f, 0.11616f, 0.11075f, 0.10558f, 0.10125f, 0.09800f, 0.09576f, 0.09430f, 0.09339f, 0.09284f, 0.09250f, 0.09231f, 0.09219f, 0.09212f, 0.09208f, 0.09205f, 0.09204f, 0.09203f, 0.09203f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f, 0.09202f
};
static const float LUT_NMC_OCP_101[101] = {
    3.54498f, 3.55049f, 3.55600f, 3.56150f, 3.56701f, 3.57252f, 3.57803f, 3.58354f, 3.58904f, 3.59455f, 3.60006f, 3.60557f, 3.61108f, 3.61658f, 3.62209f, 3.62760f, 3.63311f, 3.63862f, 3.64413f, 3.64964f, 3.65515f, 3.66066f, 3.66617f, 3.67169f, 3.67720f, 3.68272f, 3.68824f, 3.69377f, 3.69930f, 3.70483f, 3.71038f, 3.71593f, 3.72150f, 3.72708f, 3.73269f, 3.73832f, 3.74399f, 3.74970f, 3.75547f, 3.76132f, 3.76725f, 3.77330f, 3.77950f, 3.78588f, 3.79248f, 3.79935f, 3.80654f, 3.81412f, 3.82213f, 3.83064f, 3.83967f, 3.84923f, 3.85931f, 3.86983f, 3.88069f, 3.89175f, 3.90285f, 3.91383f, 3.92455f, 3.93491f, 3.94485f, 3.95437f, 3.96347f, 3.97222f, 3.98067f, 3.98892f, 3.99704f, 4.00511f, 4.01322f, 4.02144f, 4.02983f, 4.03847f, 4.04740f, 4.05668f, 4.06633f, 4.07637f, 4.08678f, 4.09753f, 4.10852f, 4.11962f, 4.13064f, 4.14133f, 4.15141f, 4.16056f, 4.16850f, 4.17500f, 4.18001f, 4.18365f, 4.18633f, 4.18870f, 4.19159f, 4.19591f, 4.20247f, 4.21187f, 4.22436f, 4.23986f, 4.25796f, 4.27805f, 4.29944f, 4.32140f, 4.34332f
};
static float ev_clamp(float value, float lower, float upper) {
    if (value < lower) { return lower; }
    if (value > upper) { return upper; }
    return value;
}

static double ev_clamp_d(double value, double lower, double upper) {
    if (value < lower) { return lower; }
    if (value > upper) { return upper; }
    return value;
}

static bool ev_is_finite(float value) {
    return (value == value) && (value <= FLT_MAX) && (value >= -FLT_MAX);
}

/* exp(-Ea/R * (1/T - 1/Tref)), exponent bounded */
static float ev_arrhenius(float activation_energy_J_mol, float temperature_K) {
    float exponent = (-activation_energy_J_mol / EV_GAS_J_PER_MOL_K) * ((1.0f / temperature_K) - (1.0f / EV_REFERENCE_TEMPERATURE_K));
    return expf(ev_clamp(exponent, -20.0f, 20.0f));
}

/* Parameters that must hold before any state update divides by them. */
static bool ev_is_valid_config(const EvCellState *cell) {
    return (cell != NULL)
        && (cell->nominal_capacity_Ah > 0.0f)
        && (cell->soh > 0.0f) && (cell->soh <= 1.0f)
        && (cell->resistance_ohm >= 0.0f)
        && (cell->heat_transfer_W_m2_K > 0.0f)
        && (cell->heat_capacity_J_kg_K > 0.0f)
        && (cell->active_solid_fraction > 0.0f) && (cell->active_solid_fraction < 1.0f)
        && (cell->soc_variance >= 0.0)
        && (cell->soc >= 0.0) && (cell->soc <= 1.0);
}

static bool ev_is_valid_cell(const EvCellState *cell) {
    return ev_is_valid_config(cell)
        && (cell->temperature_C >= EV_MIN_CELL_TEMPERATURE_C)
        && (cell->temperature_C <= EV_MAX_CELL_TEMPERATURE_C)
        && cell->sensor_valid;
}

static float ev_lut_interp(const float lut[EV_LUT_SEGMENTS + 1U], float soc, float *slope_per_soc) {
    float f_idx = ev_clamp(soc, 0.0f, 1.0f) * (float)EV_LUT_SEGMENTS;
    uint8_t idx = (uint8_t)f_idx;
    if (idx >= EV_LUT_SEGMENTS) { idx = (uint8_t)(EV_LUT_SEGMENTS - 1U); }
    float t = f_idx - (float)idx;
    float delta = lut[idx + 1U] - lut[idx];
    if (slope_per_soc != NULL) { *slope_per_soc = delta * (float)EV_LUT_SEGMENTS; }
    return lut[idx] + (t * delta);
}

void ev_full_cell_ocv_V_and_dv_dsoc(float soc, float *ocv, float *dv_dsoc) {
    float graphite_slope = 0.0f;
    float nmc_slope = 0.0f;
    float graphite_val = ev_lut_interp(LUT_GRAPHITE_OCP_101, soc, &graphite_slope);
    float nmc_val = ev_lut_interp(LUT_NMC_OCP_101, soc, &nmc_slope);

    if (ocv != NULL) { *ocv = nmc_val - graphite_val; }
    if (dv_dsoc != NULL) { *dv_dsoc = nmc_slope - graphite_slope; }
}

static float ev_predict_temperature_C(const EvCellState *cell, float charge_current_A, float dt_s) {
    float thermal_conductance_W_K = cell->heat_transfer_W_m2_K * EV_CELL_COOLING_AREA_M2;
    float heat_capacity_J_K = EV_CELL_MASS_KG * cell->heat_capacity_J_kg_K;
    float joule_heat_W = charge_current_A * charge_current_A * cell->resistance_ohm;
    float temperature_rise_steady_C = joule_heat_W / thermal_conductance_W_K;
    float exponent = -(thermal_conductance_W_K * dt_s) / heat_capacity_J_K;
    float decay = expf(ev_clamp(exponent, -20.0f, 0.0f));

    return EV_AMBIENT_TEMPERATURE_C + temperature_rise_steady_C
        + ((cell->temperature_C - EV_AMBIENT_TEMPERATURE_C - temperature_rise_steady_C) * decay);
}

bool ev_predict_cell(const EvCellState *cell, float charge_current_A, float horizon_s, EvCellPrediction *prediction) {
    const float anode_stoich_window = EV_ANODE_STOICH_AT_SOC100 - EV_ANODE_STOICH_AT_SOC0;
    const float thermal_voltage_V = EV_GAS_J_PER_MOL_K / EV_FARADAY_C_PER_MOL;

    if (prediction == NULL) {
        return false;
    }
    prediction->anode_potential_V = 0.0f;
    prediction->predicted_temperature_C = 0.0f;
    prediction->predicted_soc = 0.0f;
    prediction->safe = false;
    prediction->error_code = EV_PREDICTION_INVALID_INPUT;

    if (!ev_is_valid_cell(cell) || !ev_is_finite(charge_current_A) || !ev_is_finite(horizon_s)
        || (charge_current_A < 0.0f) || (horizon_s <= 0.0f)) {
        return false;
    }

    double effective_capacity_Ah = (double)cell->nominal_capacity_Ah * (double)cell->soh;
    prediction->predicted_soc = (float)(cell->soc + (((double)charge_current_A * (double)horizon_s) / (effective_capacity_Ah * 3600.0)));
    prediction->predicted_temperature_C = ev_predict_temperature_C(cell, charge_current_A, horizon_s);

    /* Plating risk rises as the cell gets colder: evaluate the electrochemistry at the
     * colder of the present and predicted temperatures, minus the sensor margin. */
    float plating_temperature_C = fminf(cell->temperature_C, prediction->predicted_temperature_C) - EV_TEMPERATURE_SENSOR_MARGIN_C;
    float temperature_K = ev_clamp(plating_temperature_C, EV_MIN_CELL_TEMPERATURE_C - EV_TEMPERATURE_SENSOR_MARGIN_C,
                                   EV_MAX_CELL_TEMPERATURE_C) + EV_KELVIN_OFFSET;

    /* Solid phase: steady-state polynomial (Subramanian) surface-concentration correction.
     * Charging inserts lithium into graphite, so the surface runs ahead of the bulk. */
    float solid_diffusivity_m2_s = EV_SOLID_DIFFUSIVITY_REF_M2_S * ev_arrhenius(EV_DIFFUSION_ACTIVATION_ENERGY_J_MOL, temperature_K);
    float current_density_A_m2 = charge_current_A / EV_GEOMETRIC_ELECTRODE_AREA_M2;
    float specific_surface_area_m_inv = (3.0f * cell->active_solid_fraction) / EV_PARTICLE_RADIUS_M;
    float interfacial_current_A_m2 = current_density_A_m2 / (specific_surface_area_m_inv * EV_ANODE_THICKNESS_M);
    float molar_flux_mol_m2_s = interfacial_current_A_m2 / EV_FARADAY_C_PER_MOL;   /* into the particle */
    float surface_stoich_rise = (molar_flux_mol_m2_s * EV_PARTICLE_RADIUS_M) / (5.0f * solid_diffusivity_m2_s * EV_MAX_CONCENTRATION_MOL_M3);
    float surface_soc = ev_clamp(prediction->predicted_soc + (surface_stoich_rise / anode_stoich_window), 0.0f, 1.0f);
    float surface_stoich = ev_clamp(EV_ANODE_STOICH_AT_SOC0 + (surface_soc * anode_stoich_window), 0.001f, 0.999f);

    /* Electrolyte: steady-state concentration profile across the anode with uniform
     * reaction; electrode-averaged depletion relative to the separator. */
    float porosity_factor = powf(EV_ELECTROLYTE_POROSITY, EV_BRUGGEMAN_EXPONENT);
    float effective_conductivity_S_m = EV_ELECTROLYTE_CONDUCTIVITY_S_M * porosity_factor
        * ev_arrhenius(EV_ELECTROLYTE_ACTIVATION_ENERGY_J_MOL, temperature_K);
    float effective_electrolyte_diffusivity_m2_s = EV_ELECTROLYTE_DIFFUSIVITY_REF_M2_S * porosity_factor
        * ev_arrhenius(EV_ELECTROLYTE_DIFFUSIVITY_ACTIVATION_J_MOL, temperature_K);
    float depletion_fraction = ((1.0f - EV_TRANSFERENCE_NUMBER) * current_density_A_m2 * EV_ANODE_THICKNESS_M)
        / (3.0f * EV_FARADAY_C_PER_MOL * effective_electrolyte_diffusivity_m2_s * EV_ELECTROLYTE_CONCENTRATION_MOL_M3);
    if (depletion_fraction >= 0.99f) {
        prediction->error_code = EV_PREDICTION_ELECTROLYTE_DEPLETION;
        return false;
    }
    float concentration_overpotential_V = -2.0f * thermal_voltage_V * temperature_K * (1.0f - EV_TRANSFERENCE_NUMBER)
        * log1pf(-depletion_fraction);

    /* Butler-Volmer (alpha = 0.5) inverted for the cathodic (charging) overpotential. */
    float exchange_current_A_m2 = EV_EXCHANGE_CURRENT_REF_A_M2 * sqrtf(surface_stoich * (1.0f - surface_stoich))
        * ev_arrhenius(EV_KINETIC_ACTIVATION_ENERGY_J_MOL, temperature_K);
    float kinetic_overpotential_V = 2.0f * thermal_voltage_V * temperature_K
        * asinhf(interfacial_current_A_m2 / (2.0f * exchange_current_A_m2));

    float electrolyte_ohmic_drop_V = (current_density_A_m2 * EV_ANODE_THICKNESS_M) / (2.0f * effective_conductivity_S_m);

    prediction->anode_potential_V = ev_lut_interp(LUT_GRAPHITE_OCP_101, surface_soc, NULL)
        - kinetic_overpotential_V - electrolyte_ohmic_drop_V - concentration_overpotential_V;
    prediction->safe = (prediction->predicted_soc <= EV_MAX_SAFE_SOC)
        && (prediction->predicted_temperature_C <= EV_MAX_CELL_TEMPERATURE_C)
        && (prediction->anode_potential_V >= EV_MIN_SAFE_ANODE_POTENTIAL_V);
    prediction->error_code = EV_PREDICTION_OK;
    return true;
}

static bool ev_current_is_safe(const EvCellState *cells, float current_A, float horizon_s, uint8_t *limiting_index, EvCellPrediction *limiting_prediction) {
    uint8_t index;
    EvCellPrediction worst = {FLT_MAX, -FLT_MAX, -FLT_MAX, true, EV_PREDICTION_OK};
    uint8_t worst_index = 0U;

    for (index = 0U; index < EV_NUM_CELLS; ++index) {
        EvCellPrediction prediction;
        if (!ev_predict_cell(&cells[index], current_A, horizon_s, &prediction) || !prediction.safe) {
            *limiting_index = index;
            *limiting_prediction = prediction;
            return false;
        }
        if (prediction.anode_potential_V < worst.anode_potential_V) {
            worst = prediction;
            worst_index = index;
        }
    }
    *limiting_index = worst_index;
    *limiting_prediction = worst;
    return true;
}

/*
 * Largest current in [0, requested] below the first unsafe point found by a coarse
 * upward sweep, refined by bisection. Sweeping from zero means a non-monotonic
 * constraint can only make the result more conservative, never less.
 */
static float ev_find_safe_charge_current(const EvCellState *cells, float requested_current_A, float horizon_s, uint8_t *limiting_index, EvCellPrediction *limiting_prediction) {
    float safe_A = 0.0f;
    float unsafe_A = 0.0f;
    uint32_t step;
    bool found_unsafe = false;

    if (!ev_current_is_safe(cells, 0.0f, horizon_s, limiting_index, limiting_prediction)) {
        return 0.0f;
    }
    if (requested_current_A <= 0.0f) {
        return 0.0f;
    }

    for (step = 1U; step <= EV_COARSE_SEARCH_STEPS; ++step) {
        float candidate_A = (requested_current_A * (float)step) / (float)EV_COARSE_SEARCH_STEPS;
        if (!ev_current_is_safe(cells, candidate_A, horizon_s, limiting_index, limiting_prediction)) {
            unsafe_A = candidate_A;
            found_unsafe = true;
            break;
        }
        safe_A = candidate_A;
    }

    if (found_unsafe) {
        for (step = 0U; step < EV_BISECTION_STEPS; ++step) {
            float mid_A = 0.5f * (safe_A + unsafe_A);
            if (ev_current_is_safe(cells, mid_A, horizon_s, limiting_index, limiting_prediction)) {
                safe_A = mid_A;
            } else {
                unsafe_A = mid_A;
            }
        }
        /* The binding cell is the one that fails just above the approved current;
         * report its prediction at the approved current. */
        (void)ev_current_is_safe(cells, unsafe_A, horizon_s, limiting_index, limiting_prediction);
        (void)ev_predict_cell(&cells[*limiting_index], safe_A, horizon_s, limiting_prediction);
    }
    return safe_A;
}

static void ev_ekf_and_ageing_update(EvCellState *cell, float measured_current_A, float measured_voltage_V, float measured_temperature_C, float dt_s) {
    bool voltage_plausible = ev_is_finite(measured_voltage_V)
        && (measured_voltage_V >= EV_MIN_PLAUSIBLE_VOLTAGE_V) && (measured_voltage_V <= EV_MAX_PLAUSIBLE_VOLTAGE_V);
    bool temperature_plausible = ev_is_finite(measured_temperature_C)
        && (measured_temperature_C >= EV_MIN_PLAUSIBLE_TEMPERATURE_C) && (measured_temperature_C <= EV_MAX_PLAUSIBLE_TEMPERATURE_C);

    /* Ageing: SEI growth with throughput plus calendar loss, Arrhenius in temperature. */
    float processed_Ah = (fabsf(measured_current_A) * dt_s) / 3600.0f;
    float ageing_temperature_K = ev_clamp(cell->temperature_C, EV_MIN_PLAUSIBLE_TEMPERATURE_C, EV_MAX_PLAUSIBLE_TEMPERATURE_C) + EV_KELVIN_OFFSET;
    float ageing_factor = ev_arrhenius(EV_AGING_ACTIVATION_ENERGY_J_MOL, ageing_temperature_K);
    float capacity_loss = ((EV_SEI_GROWTH_RATE_PER_AH * processed_Ah) + (EV_CALENDRIC_AGING_RATE_PER_S * dt_s)) * ageing_factor;
    float previous_soh = cell->soh;
    cell->throughput_Ah += processed_Ah;
    cell->soh = ev_clamp(cell->soh - capacity_loss, EV_MIN_SOH, 1.0f);
    cell->resistance_ohm += EV_SEI_RESISTANCE_FACTOR_OHM * (previous_soh - cell->soh);

    /* EKF time update: coulomb counting. */
    double effective_capacity_Ah = (double)cell->nominal_capacity_Ah * (double)cell->soh;
    double soc_pred = ev_clamp_d(cell->soc + (((double)measured_current_A * (double)dt_s) / (effective_capacity_Ah * 3600.0)), 0.0, 1.0);
    double variance_pred = cell->soc_variance + (EV_SOC_PROCESS_NOISE_PER_S * (double)dt_s);

    /* EKF measurement update (Joseph form), only with a plausible voltage. */
    if (voltage_plausible) {
        float ocv = 0.0f;
        float dv_dsoc = 0.0f;
        ev_full_cell_ocv_V_and_dv_dsoc((float)soc_pred, &ocv, &dv_dsoc);
        double h = (double)dv_dsoc;
        double v_pred = (double)ocv + ((double)measured_current_A * (double)cell->resistance_ohm);
        double s = (h * variance_pred * h) + EV_VOLTAGE_MEASUREMENT_NOISE_V2;
        double k = (variance_pred * h) / s;
        double one_minus_kh = 1.0 - (k * h);
        soc_pred = ev_clamp_d(soc_pred + (k * ((double)measured_voltage_V - v_pred)), 0.0, 1.0);
        variance_pred = (one_minus_kh * variance_pred * one_minus_kh) + (k * EV_VOLTAGE_MEASUREMENT_NOISE_V2 * k);
    }
    cell->soc = soc_pred;
    cell->soc_variance = variance_pred;

    /* Keep the real measured temperature (no clamping) so over/under-temperature is
     * reported as a fault; fall back to the thermal model if the sensor is implausible. */
    if (temperature_plausible) {
        cell->temperature_C = measured_temperature_C;
    } else {
        cell->temperature_C = ev_predict_temperature_C(cell, measured_current_A, dt_s);
    }
    cell->sensor_valid = voltage_plausible && temperature_plausible;
}

EvPackResult ev_step(EvCellState *cells, float measured_pack_current_A, float requested_charge_current_A,
                     const float measured_voltages_V[EV_NUM_CELLS], const float measured_temperatures_C[EV_NUM_CELLS],
                     float control_dt_s, float prediction_horizon_s) {
    EvPackResult result = {EV_DECISION_FAULT, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0U};
    EvCellPrediction limiting_prediction = {0.0f, 0.0f, 0.0f, false, EV_PREDICTION_OK};
    uint8_t limiting_index = 0U;
    uint8_t index;

    if ((cells == NULL) || (measured_voltages_V == NULL) || (measured_temperatures_C == NULL)
        || !ev_is_finite(measured_pack_current_A) || (fabsf(measured_pack_current_A) > EV_MAX_PLAUSIBLE_CURRENT_A)
        || !ev_is_finite(requested_charge_current_A)
        || !ev_is_finite(control_dt_s) || (control_dt_s <= 0.0f)
        || !ev_is_finite(prediction_horizon_s) || (prediction_horizon_s <= 0.0f)) {
        return result;
    }
    requested_charge_current_A = (requested_charge_current_A > 0.0f) ? requested_charge_current_A : 0.0f;
    result.requested_charge_current_A = requested_charge_current_A;

    for (index = 0U; index < EV_NUM_CELLS; ++index) {
        if (!ev_is_valid_config(&cells[index])) {
            result.limiting_cell_index = index;
            return result;
        }
    }

    /* Estimate first: the measurements describe the interval that has just elapsed. */
    for (index = 0U; index < EV_NUM_CELLS; ++index) {
        ev_ekf_and_ageing_update(&cells[index], measured_pack_current_A, measured_voltages_V[index], measured_temperatures_C[index], control_dt_s);
    }
    for (index = 0U; index < EV_NUM_CELLS; ++index) {
        if (!ev_is_valid_cell(&cells[index])) {
            result.limiting_cell_index = index;
            result.limiting_temperature_C = cells[index].temperature_C;
            result.limiting_soc = (float)cells[index].soc;
            return result;
        }
    }

    result.approved_charge_current_A = ev_find_safe_charge_current(cells, requested_charge_current_A, prediction_horizon_s, &limiting_index, &limiting_prediction);
    result.limiting_cell_index = limiting_index;
    result.limiting_anode_potential_V = limiting_prediction.anode_potential_V;
    result.limiting_temperature_C = limiting_prediction.predicted_temperature_C;
    result.limiting_soc = limiting_prediction.predicted_soc;
    result.decision = ((result.approved_charge_current_A + EV_CURRENT_TOLERANCE_A) < requested_charge_current_A) ? EV_DECISION_LIMIT : EV_DECISION_ALLOW;
    return result;
}
