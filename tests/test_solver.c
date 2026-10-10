/*
 * Self-test and demo for the EV MODIFY charge limiter.
 *
 *   test_solver --self-test   run all checks, exit 0 on success
 *   test_solver               print a short CSV demo run
 *
 * The plant used here is the controller's own OCV + 1-RC voltage model (optionally with
 * deliberately mismatched parameters); the validation against an independent model is
 * tools/run_pybamm_sim.py.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spme_v7_apex_solver.h"

#define DT_S 1.0f
#define HORIZON_S 30.0f

static int g_failures = 0;
static int g_checks = 0;

static void check(bool condition, const char *name) {
    ++g_checks;
    if (!condition) {
        printf("  FAIL: %s\n", name);
        ++g_failures;
    }
}

static float make_nan(void) { return strtof("nan", NULL); }
static float make_inf(void) { return strtof("inf", NULL); }

/* 1-RC parameters fitted to PyBaMM SPMe Chen2020 by tools/run_pybamm_sim.py */
static EvCellConfig lgm50_config(void) {
    EvCellConfig config;
    config.nominal_capacity_Ah = 5.0f;
    config.resistance_ohm = 0.041f;
    config.rc_resistance_ohm = 0.0135f;
    config.rc_time_constant_s = 190.0f;
    config.heat_transfer_W_m2_K = 10.0f;
    config.heat_capacity_J_kg_K = 1100.0f;
    config.coolant_temperature_C = 25.0f;
    return config;
}

typedef struct {
    EvCellConfig config[EV_NUM_CELLS];
    EvCellState state[EV_NUM_CELLS];
    double true_soc[EV_NUM_CELLS];
    double true_rc_V[EV_NUM_CELLS];
    float plant_r0_scale;            /* plant R0 = scale * configured R0 */
    float plant_r1_scale;            /* plant R1 = scale * configured R1 */
    float plant_tau_scale;           /* plant tau = scale * configured tau */
    float plant_capacity_scale;      /* plant capacity = scale * configured capacity */
    float adc_lsb_V;                 /* voltage quantisation (0 = none) */
    float voltages_V[EV_NUM_CELLS];
    float temperatures_C[EV_NUM_CELLS];
} Pack;

static void pack_init(Pack *pack, const double soc[EV_NUM_CELLS], const float temperature_C[EV_NUM_CELLS]) {
    uint8_t i;
    memset(pack, 0, sizeof(*pack));
    pack->plant_r0_scale = 1.0f;
    pack->plant_r1_scale = 1.0f;
    pack->plant_tau_scale = 1.0f;
    pack->plant_capacity_scale = 1.0f;
    for (i = 0U; i < EV_NUM_CELLS; ++i) {
        pack->config[i] = lgm50_config();
        ev_cell_init(&pack->state[i], soc[i], temperature_C[i]);
        pack->state[i].covariance[0] = 1.0e-4;   /* sigma 1 %: SOC known from a rested OCV reading */
        pack->true_soc[i] = soc[i];
        pack->temperatures_C[i] = temperature_C[i];
    }
}

static void pack_uniform(Pack *pack, double soc, float temperature_C) {
    double socs[EV_NUM_CELLS];
    float temps[EV_NUM_CELLS];
    uint8_t i;
    for (i = 0U; i < EV_NUM_CELLS; ++i) {
        socs[i] = soc;
        temps[i] = temperature_C;
    }
    pack_init(pack, socs, temps);
}

/* Resistance rise in the cold, as in the controller (Arrhenius, EV_RESISTANCE_ACTIVATION_ENERGY_J_MOL). */
static double plant_resistance_scale(float temperature_C) {
    return exp(((double)EV_RESISTANCE_ACTIVATION_ENERGY_J_MOL / (double)EV_GAS_J_PER_MOL_K)
               * ((1.0 / ((double)temperature_C + (double)EV_KELVIN_OFFSET)) - (1.0 / (double)EV_REFERENCE_TEMPERATURE_K)));
}

/* Advances the plant by dt with the current applied over the last interval. */
static void plant_step(Pack *pack, float current_A, float dt_s) {
    uint8_t i;
    for (i = 0U; i < EV_NUM_CELLS; ++i) {
        const EvCellConfig *c = &pack->config[i];
        double a = exp(-(double)dt_s / ((double)c->rc_time_constant_s * (double)pack->plant_tau_scale));
        double scale = plant_resistance_scale(pack->temperatures_C[i]);
        float ocv = 0.0f;
        pack->true_soc[i] += ((double)current_A * (double)dt_s)
            / ((double)c->nominal_capacity_Ah * (double)pack->plant_capacity_scale * pack->state[i].soh * 3600.0);
        pack->true_rc_V[i] = (a * pack->true_rc_V[i])
            + ((double)c->rc_resistance_ohm * (double)pack->plant_r1_scale * scale * (1.0 - a) * (double)current_A);
        ev_full_cell_ocv_V_and_dv_dsoc((float)pack->true_soc[i], &ocv, NULL);
        pack->voltages_V[i] = ocv + (float)pack->true_rc_V[i] + (float)((double)pack->plant_r0_scale * (double)c->resistance_ohm * scale * (double)current_A);
        if (pack->adc_lsb_V > 0.0f) {
            pack->voltages_V[i] = roundf(pack->voltages_V[i] / pack->adc_lsb_V) * pack->adc_lsb_V;
        }
    }
}

static EvPackResult pack_step(Pack *pack, float applied_A, float request_A) {
    plant_step(pack, applied_A, DT_S);
    return ev_step(pack->config, pack->state, applied_A, request_A, pack->voltages_V, pack->temperatures_C, DT_S, HORIZON_S);
}

/* ---- Tests ----------------------------------------------------------------------------- */

static void test_ocv_table(void) {
    float ocv = 0.0f;
    float slope = 0.0f;
    float soc;
    float previous = 0.0f;
    bool monotone = true;
    ev_full_cell_ocv_V_and_dv_dsoc(0.0f, &ocv, NULL);
    check(fabsf(ocv - 2.5f) < 1.0e-3f, "ocv: 2.5 V at SOC 0");
    ev_full_cell_ocv_V_and_dv_dsoc(1.0f, &ocv, NULL);
    check(fabsf(ocv - EV_MAX_CELL_VOLTAGE_V) < 1.0e-3f, "ocv: 4.2 V at SOC 1");
    for (soc = 0.0f; soc <= 1.0f; soc += 0.001f) {
        ev_full_cell_ocv_V_and_dv_dsoc(soc, &ocv, &slope);
        monotone = monotone && (ocv >= previous) && (slope > 0.0f);
        previous = ocv;
    }
    check(monotone, "ocv: strictly increasing with SOC");
    ev_full_cell_ocv_V_and_dv_dsoc(make_nan(), &ocv, &slope);
    check(isnan(ocv) && isnan(slope), "ocv: NaN SOC gives NaN");
    ev_full_cell_ocv_V_and_dv_dsoc(make_inf(), &ocv, NULL);
    check(isnan(ocv), "ocv: inf SOC gives NaN");
}

static void test_closed_loop_charge_to_full(void) {
    const double socs[EV_NUM_CELLS] = {0.20, 0.30, 0.25, 0.22};
    const float temps[EV_NUM_CELLS] = {25.0f, 15.0f, 30.0f, 22.0f};
    Pack pack;
    EvPackResult result;
    float applied_A = 0.0f;
    float max_voltage_V = 0.0f;
    bool any_fault = false;
    bool bounded = true;
    bool margin_kept = true;
    bool saw_limit = false;
    uint32_t t;
    uint8_t i;
    pack_init(&pack, socs, temps);

    for (t = 0U; t < 6U * 3600U; ++t) {
        result = pack_step(&pack, applied_A, 10.0f);   /* 2C request */
        any_fault = any_fault || (result.decision == EV_DECISION_FAULT);
        bounded = bounded && (result.approved_charge_current_A >= 0.0f) && (result.approved_charge_current_A <= 10.0f);
        if (result.decision == EV_DECISION_LIMIT) {
            saw_limit = true;
            margin_kept = margin_kept && (result.limiting_constraints != 0U) && (result.limiting_cell_index < EV_NUM_CELLS);
        }
        if (result.decision != EV_DECISION_FAULT) {
            margin_kept = margin_kept && (result.limiting_anode_potential_V >= EV_MIN_SAFE_ANODE_POTENTIAL_V - 1.0e-4f);
        }
        for (i = 0U; i < EV_NUM_CELLS; ++i) {
            if (pack.voltages_V[i] > max_voltage_V) { max_voltage_V = pack.voltages_V[i]; }
        }
        applied_A = result.approved_charge_current_A;
        if (any_fault || (applied_A < 0.05f)) { break; }
    }
    check(!any_fault, "closed loop: no fault");
    check(bounded, "closed loop: 0 <= approved <= request");
    check(saw_limit, "closed loop: 2C request is limited");
    check(margin_kept, "closed loop: anode margin kept and binding constraint reported");
    check(applied_A < 0.05f, "closed loop: charge tapers to zero");
    check(max_voltage_V <= EV_MAX_CELL_VOLTAGE_V + 0.005f, "closed loop: terminal voltage never exceeds 4.2 V");
    for (i = 0U; i < EV_NUM_CELLS; ++i) {
        check(pack.true_soc[i] <= (double)EV_MAX_SAFE_SOC + 1.0e-3, "closed loop: true SOC stays below the SOC limit");
    }
    check(pack.true_soc[1] > 0.85, "closed loop: pack charges to high SOC");
}

static void test_allow_small_request(void) {
    Pack pack;
    EvPackResult result;
    pack_uniform(&pack, 0.4, 25.0f);
    result = pack_step(&pack, 0.0f, 1.0f);
    check(result.decision == EV_DECISION_ALLOW, "allow: 0.2C at 40 % SOC, 25 degC is allowed");
    check(fabsf(result.approved_charge_current_A - 1.0f) < 1.0e-6f, "allow: full request approved");
    check(result.limiting_constraints == 0U, "allow: no binding constraint");
}

static void test_monotone_constraints(void) {
    const double socs[] = {0.05, 0.3, 0.6, 0.9};
    const float temps[] = {2.0f, 15.0f, 30.0f, 44.0f};
    bool monotone = true;
    bool colder_is_worse = true;
    size_t s;
    size_t k;
    for (s = 0U; s < sizeof(socs) / sizeof(socs[0]); ++s) {
        for (k = 0U; k < sizeof(temps) / sizeof(temps[0]); ++k) {
            EvCellConfig config = lgm50_config();
            EvCellState state;
            EvCellPrediction previous;
            float current_A;
            ev_cell_init(&state, socs[s], temps[k]);
            (void)ev_predict_cell(&config, &state, 0.0f, HORIZON_S, &previous);
            for (current_A = 0.05f; current_A <= 30.0f; current_A += 0.05f) {
                EvCellPrediction now;
                (void)ev_predict_cell(&config, &state, current_A, HORIZON_S, &now);
                monotone = monotone && (now.anode_potential_V <= previous.anode_potential_V)
                    && (now.terminal_voltage_V >= previous.terminal_voltage_V)
                    && (now.predicted_temperature_C >= previous.predicted_temperature_C)
                    && (now.predicted_soc >= previous.predicted_soc);
                previous = now;
            }
            if (k > 0U) {
                EvCellState colder;
                EvCellPrediction warm_prediction;
                EvCellPrediction cold_prediction;
                ev_cell_init(&colder, socs[s], temps[k - 1U]);
                (void)ev_predict_cell(&config, &state, 5.0f, HORIZON_S, &warm_prediction);
                (void)ev_predict_cell(&config, &colder, 5.0f, HORIZON_S, &cold_prediction);
                colder_is_worse = colder_is_worse && (cold_prediction.anode_potential_V < warm_prediction.anode_potential_V);
            }
        }
    }
    check(monotone, "physics: constraints monotone in current (search assumption)");
    check(colder_is_worse, "physics: colder cell has lower anode potential");
}

static void test_full_cell_blocks_charge(void) {
    Pack pack;
    EvPackResult result;
    pack_uniform(&pack, 0.5, 25.0f);
    pack.state[2].soc = 0.99;
    pack.true_soc[2] = 0.99;
    result = pack_step(&pack, 0.0f, 5.0f);
    check(result.decision == EV_DECISION_LIMIT, "full cell: LIMIT");
    check(result.approved_charge_current_A == 0.0f, "full cell: zero current");
    check(result.limiting_cell_index == 2U, "full cell: limiting cell identified");
    check((result.limiting_constraints & EV_CONSTRAINT_SOC) != 0U, "full cell: SOC constraint reported");
}

static void test_soc_uncertainty_derates(void) {
    EvCellConfig config = lgm50_config();
    EvCellState certain;
    EvCellState uncertain;
    EvCellPrediction p_certain;
    EvCellPrediction p_uncertain;
    ev_cell_init(&certain, 0.94, 25.0f);
    ev_cell_init(&uncertain, 0.94, 25.0f);
    certain.covariance[0] = 1.0e-6;
    uncertain.covariance[0] = 4.0e-4;   /* sigma 2 % */
    (void)ev_predict_cell(&config, &certain, 1.0f, HORIZON_S, &p_certain);
    (void)ev_predict_cell(&config, &uncertain, 1.0f, HORIZON_S, &p_uncertain);
    check((p_certain.violated & EV_CONSTRAINT_SOC) == 0U, "uncertainty: certain SOC 0.94 may charge");
    check((p_uncertain.violated & EV_CONSTRAINT_SOC) != 0U, "uncertainty: 3-sigma bound blocks SOC 0.94 +/- 2 %");
}

static void test_cold_cell_blocks_charge(void) {
    Pack pack;
    EvPackResult result;
    pack_uniform(&pack, 0.5, 25.0f);
    pack.temperatures_C[1] = -5.0f;
    result = pack_step(&pack, 0.0f, 5.0f);
    check(result.decision == EV_DECISION_LIMIT, "cold: LIMIT, not FAULT, below 0 degC");
    check(result.approved_charge_current_A == 0.0f, "cold: no charging below 0 degC");
    check((result.limiting_constraints & EV_CONSTRAINT_TEMPERATURE_LOW) != 0U, "cold: low-temperature constraint reported");
}

static void test_predicted_overtemperature_limits(void) {
    EvCellConfig config = lgm50_config();
    EvCellState state;
    EvCellPrediction prediction;
    ev_cell_init(&state, 0.3, 44.0f);
    (void)ev_predict_cell(&config, &state, 20.0f, HORIZON_S, &prediction);
    check((prediction.violated & EV_CONSTRAINT_TEMPERATURE_HIGH) != 0U, "thermal: 4C at 44 degC exceeds 45 degC");
    (void)ev_predict_cell(&config, &state, 0.5f, HORIZON_S, &prediction);
    check((prediction.violated & EV_CONSTRAINT_TEMPERATURE_HIGH) == 0U, "thermal: 0.1C at 44 degC stays below 45 degC");
}

static void test_overtemperature_fault_hysteresis_and_hold(void) {
    Pack pack;
    EvPackResult result;
    uint32_t cycle;
    pack_uniform(&pack, 0.5, 25.0f);

    pack.temperatures_C[3] = 46.0f;
    result = pack_step(&pack, 0.0f, 5.0f);
    check(result.decision == EV_DECISION_FAULT, "overtemp: FAULT");
    check(result.approved_charge_current_A == 0.0f, "overtemp: zero current");
    check(result.limiting_cell_index == 3U, "overtemp: faulted cell identified");
    check((result.fault_flags & EV_FAULT_OVER_TEMPERATURE) != 0U, "overtemp: flag set");

    pack.temperatures_C[3] = 44.0f;   /* back below 45 but inside the 3 degC hysteresis band */
    result = pack_step(&pack, 0.0f, 5.0f);
    check(result.decision == EV_DECISION_FAULT, "overtemp: hysteresis keeps the fault at 44 degC");

    pack.temperatures_C[3] = 40.0f;
    for (cycle = 1U; cycle < EV_FAULT_HOLD_CYCLES; ++cycle) {
        result = pack_step(&pack, 0.0f, 5.0f);
        check((result.decision == EV_DECISION_FAULT) && ((result.fault_flags & EV_FAULT_HOLD) != 0U), "overtemp: fault held");
    }
    result = pack_step(&pack, 0.0f, 5.0f);
    check(result.decision != EV_DECISION_FAULT, "overtemp: clears after the hold cycles");
}

static void test_invalid_inputs(void) {
    Pack pack;
    EvPackResult result;
    double soc_before;
    double variance_before;
    pack_uniform(&pack, 0.5, 25.0f);
    plant_step(&pack, 0.0f, DT_S);

    soc_before = pack.state[0].soc;
    variance_before = pack.state[0].covariance[0];
    result = ev_step(pack.config, pack.state, make_nan(), 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    check((result.decision == EV_DECISION_FAULT) && (result.fault_flags == EV_FAULT_INPUT), "invalid: NaN current faults");
    check(result.limiting_cell_index == EV_NO_CELL, "invalid: pack-level fault has no cell index");
    check(pack.state[0].soc == soc_before, "invalid: SOC mean kept when current is unknown");
    check(pack.state[0].covariance[0] > variance_before, "invalid: SOC variance inflated when current is unknown");

    result = ev_step(pack.config, pack.state, 2000.0f, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    check(result.decision == EV_DECISION_FAULT, "invalid: implausible current faults");
    result = ev_step(pack.config, pack.state, 0.0f, 5.0f, pack.voltages_V, NULL, DT_S, HORIZON_S);
    check(result.decision == EV_DECISION_FAULT, "invalid: missing temperatures fault");
    result = ev_step(NULL, pack.state, 0.0f, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    check(result.decision == EV_DECISION_FAULT, "invalid: missing config faults");

    soc_before = pack.state[0].soc;
    result = ev_step(pack.config, pack.state, 0.0f, 5.0f, pack.voltages_V, pack.temperatures_C, 4.0e6f, HORIZON_S);
    check(result.decision == EV_DECISION_FAULT, "invalid: dt above the limit faults");
    check((pack.state[0].soc == soc_before) && (pack.state[0].soh == 1.0), "invalid: glitched dt does not age or move the state");
    check(pack.state[0].covariance[0] >= 0.01 - 1.0e-12, "invalid: unknown dt widens the SOC uncertainty");
    result = ev_step(pack.config, pack.state, 0.0f, 5.0f, pack.voltages_V, pack.temperatures_C, make_nan(), HORIZON_S);
    check(result.decision == EV_DECISION_FAULT, "invalid: NaN dt faults");
    result = ev_step(pack.config, pack.state, 0.0f, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, 1.0e6f);
    check((result.decision == EV_DECISION_FAULT) && ((result.fault_flags & EV_FAULT_INPUT) != 0U), "invalid: horizon above the limit faults");
    result = ev_step(pack.config, pack.state, 0.0f, make_nan(), pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    check(result.decision == EV_DECISION_FAULT, "invalid: NaN request faults");

    pack_uniform(&pack, 0.5, 25.0f);
    plant_step(&pack, 0.0f, DT_S);
    result = ev_step(pack.config, pack.state, 0.0f, 1.0e9f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    check(result.requested_charge_current_A == EV_MAX_PLAUSIBLE_CURRENT_A, "invalid: huge request clamped");

    pack_uniform(&pack, 0.5, 25.0f);
    plant_step(&pack, 0.0f, DT_S);
    pack.voltages_V[1] = make_nan();
    result = ev_step(pack.config, pack.state, 0.0f, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    check((result.decision == EV_DECISION_FAULT) && ((result.fault_flags & EV_FAULT_VOLTAGE_SENSOR) != 0U), "invalid: NaN cell voltage faults");
    check(isfinite(pack.state[1].soc) && (pack.state[1].soc >= 0.0) && (pack.state[1].soc <= 1.0), "invalid: SOC not corrupted by NaN voltage");

    pack_uniform(&pack, 0.5, 25.0f);
    plant_step(&pack, 0.0f, DT_S);
    pack.temperatures_C[0] = make_nan();
    result = ev_step(pack.config, pack.state, 0.0f, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    check((result.fault_flags & EV_FAULT_TEMPERATURE_SENSOR) != 0U, "invalid: NaN temperature faults");
    check(isfinite(pack.state[0].temperature_C), "invalid: temperature falls back to the model");

    pack_uniform(&pack, 0.5, 25.0f);
    plant_step(&pack, 0.0f, DT_S);
    pack.config[0].nominal_capacity_Ah = make_inf();
    pack.config[1].heat_transfer_W_m2_K = 0.0f;
    result = ev_step(pack.config, pack.state, 0.0f, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    check((result.decision == EV_DECISION_FAULT) && ((result.fault_flags & EV_FAULT_CONFIG) != 0U), "invalid: inf capacity / zero heat transfer config faults");
    check(!ev_config_is_valid(&pack.config[0]) && !ev_config_is_valid(NULL), "invalid: config validator rejects inf and NULL");

    pack_uniform(&pack, 0.5, 25.0f);
    plant_step(&pack, 0.0f, DT_S);
    pack.state[2].soh = 0.3;
    pack.state[3].covariance[0] = (double)make_inf();
    result = ev_step(pack.config, pack.state, 0.0f, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    check((result.decision == EV_DECISION_FAULT) && ((result.fault_flags & EV_FAULT_STATE) != 0U), "invalid: corrupted state is reset and faulted");
    check((pack.state[2].soh >= EV_MIN_SOH) && (pack.config[2].resistance_ohm + pack.state[2].resistance_growth_ohm > 0.0f), "invalid: reset SOH in range, resistance positive");
    check(isfinite(pack.state[3].covariance[0]) && (pack.state[3].covariance[0] <= 0.25), "invalid: reset variance finite");
}

static void test_ageing(void) {
    EvCellConfig config[EV_NUM_CELLS];
    EvCellState state[EV_NUM_CELLS];
    float voltages_V[EV_NUM_CELLS];
    float temperatures_C[EV_NUM_CELLS];
    double expected_loss;
    double r0_before;
    uint32_t t;
    uint8_t i;
    float ocv = 0.0f;
    ev_full_cell_ocv_V_and_dv_dsoc(0.5f, &ocv, NULL);
    for (i = 0U; i < EV_NUM_CELLS; ++i) {
        config[i] = lgm50_config();
        ev_cell_init(&state[i], 0.5, 25.0f);
        voltages_V[i] = ocv;
        temperatures_C[i] = 25.0f;
    }
    for (t = 0U; t < 100000U; ++t) {   /* 100000 s at rest, 10 s steps = 11.6 days */
        (void)ev_step(config, state, 0.0f, 0.0f, voltages_V, temperatures_C, 10.0f, HORIZON_S);
    }
    expected_loss = EV_CALENDAR_LOSS_PER_S * 1.0e6;
    check(fabs((1.0 - state[0].soh) - expected_loss) < 0.01 * expected_loss, "ageing: calendar loss accumulates (no float swallowing)");
    r0_before = (double)state[0].resistance_growth_ohm;
    check(r0_before > 0.0, "ageing: resistance grows with capacity loss");

    for (i = 0U; i < EV_NUM_CELLS; ++i) { ev_cell_init(&state[i], 0.5, 25.0f); }
    for (t = 0U; t < 3600U; ++t) {   /* 1 h of +/-5 A square wave: 5 Ah = 0.5 full equivalent cycles */
        float current_A = ((t / 60U) % 2U == 0U) ? 5.0f : -5.0f;
        (void)ev_step(config, state, current_A, 0.0f, voltages_V, temperatures_C, 1.0f, HORIZON_S);
    }
    expected_loss = (EV_CYCLE_LOSS_PER_FEC * 0.5) + (EV_CALENDAR_LOSS_PER_S * 3600.0);
    check(fabs(state[0].throughput_Ah - 5.0) < 1.0e-6, "ageing: throughput counted in double precision");
    check(fabs((1.0 - state[0].soh) - expected_loss) < 0.01 * expected_loss, "ageing: cycle loss follows full equivalent cycles");
}

static void test_ekf_converges_with_model_mismatch(void) {
    Pack pack;
    EvPackResult result;
    float applied_A = 0.0f;
    uint32_t t;
    uint8_t i;
    bool converged = true;
    pack_uniform(&pack, 0.3, 25.0f);
    pack.plant_r0_scale = 1.2f;                      /* plant R0 20 % above the configured value */
    for (i = 0U; i < EV_NUM_CELLS; ++i) {
        pack.state[i].soc = 0.5;                     /* estimator starts 20 SOC points high */
        pack.state[i].covariance[0] = 0.04;
    }
    for (t = 0U; t < 3600U; ++t) {
        float request_A = ((t / 300U) % 2U == 0U) ? 5.0f : 0.0f;   /* 1C pulses, 5 min on / off */
        result = pack_step(&pack, applied_A, request_A);
        applied_A = result.approved_charge_current_A;
    }
    for (i = 0U; i < EV_NUM_CELLS; ++i) {
        converged = converged && (fabs(pack.state[i].soc - pack.true_soc[i]) < 0.03);
    }
    check(converged, "ekf: converges to within 3 % SOC from a 20 point error with 20 % R0 mismatch");
    check(sqrt(pack.state[0].covariance[0]) < 0.03, "ekf: SOC uncertainty shrinks");
}


static void test_terminal_voltage_bound_covers_whole_horizon(void) {
    EvCellConfig config = lgm50_config();
    EvCellState state;
    EvCellPrediction prediction;
    float ocv_now = 0.0f;
    ev_cell_init(&state, 0.869, 25.0f);
    state.covariance[0] = 1.0e-6;
    state.rc_voltage_V = 0.1f;   /* relaxing after a charge pulse: V(t) peaks at t = 0+ */
    ev_full_cell_ocv_V_and_dv_dsoc(0.869f, &ocv_now, NULL);
    (void)ev_predict_cell(&config, &state, 1.0f, 300.0f, &prediction);
    check(prediction.terminal_voltage_V >= ocv_now + 0.1f + (config.resistance_ohm * 1.0f) - 1.0e-4f,
          "voltage: bound includes the t = 0+ peak of a relaxing RC voltage");
}

static void test_over_voltage_fault(void) {
    Pack pack;
    EvPackResult result;
    uint32_t cycle;
    pack_uniform(&pack, 0.5, 25.0f);
    plant_step(&pack, 0.0f, DT_S);
    pack.voltages_V[2] = 4.30f;
    result = ev_step(pack.config, pack.state, 0.0f, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    check((result.decision == EV_DECISION_FAULT) && ((result.fault_flags & EV_FAULT_OVER_VOLTAGE) != 0U), "over-voltage: measured 4.30 V faults");
    pack.voltages_V[2] = 4.22f;   /* below the 4.25 V trip but inside the hysteresis band */
    result = ev_step(pack.config, pack.state, 0.0f, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    check((result.fault_flags & EV_FAULT_OVER_VOLTAGE) != 0U, "over-voltage: hysteresis keeps the fault at 4.22 V");
    plant_step(&pack, 0.0f, DT_S);
    for (cycle = 0U; cycle < EV_FAULT_HOLD_CYCLES; ++cycle) {
        result = ev_step(pack.config, pack.state, 0.0f, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    }
    check(result.decision != EV_DECISION_FAULT, "over-voltage: clears after the hold once the voltage is normal");
}

static void test_stuck_voltage_sensor_is_detected(void) {
    const double socs[EV_NUM_CELLS] = {0.5, 0.5, 0.75, 0.5};
    const float temps[EV_NUM_CELLS] = {25.0f, 25.0f, 25.0f, 25.0f};
    Pack pack;
    EvPackResult result;
    float applied_A = 0.0f;
    bool detected = false;
    uint32_t t;
    pack_init(&pack, socs, temps);
    for (t = 0U; t < 3U * 3600U; ++t) {
        plant_step(&pack, applied_A, DT_S);
        if (t >= 600U) {
            pack.voltages_V[2] = 3.70f;   /* sensor stuck at a plausible value */
        }
        result = ev_step(pack.config, pack.state, applied_A, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
        detected = detected || ((result.fault_flags & EV_FAULT_VOLTAGE_PLAUSIBILITY) != 0U);
        applied_A = result.approved_charge_current_A;
    }
    check(detected, "stuck sensor: voltage plausibility fault raised");
    check(pack.true_soc[2] <= (double)EV_MAX_SAFE_SOC, "stuck sensor: cell is not overcharged");
    check(result.decision == EV_DECISION_FAULT, "stuck sensor: fault stays latched while the sensor is stuck");
}

/* Round-3 finding: a sensor frozen at its own in-range reading during a charge used to drag the
 * SOC estimate down with it (estimate 0.68 while the cell reached 0.96) and was never flagged. */
static void test_frozen_in_range_sensor_is_detected(void) {
    const double socs[EV_NUM_CELLS] = {0.30, 0.30, 0.40, 0.30};   /* cell 2 is the fullest: it limits the pack */
    const float temps[EV_NUM_CELLS] = {25.0f, 25.0f, 25.0f, 25.0f};
    Pack pack;
    EvPackResult result;
    float applied_A = 0.0f;
    float frozen_V = 0.0f;
    uint32_t detected_at = 0U;
    uint32_t t;
    pack_init(&pack, socs, temps);
    for (t = 0U; t < 2U * 3600U; ++t) {
        plant_step(&pack, applied_A, DT_S);
        if (t == 300U) { frozen_V = pack.voltages_V[2]; }
        if (t >= 300U) { pack.voltages_V[2] = frozen_V; }   /* frozen at the true value */
        result = ev_step(pack.config, pack.state, applied_A, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
        if ((detected_at == 0U) && ((result.fault_flags & EV_FAULT_VOLTAGE_PLAUSIBILITY) != 0U)) { detected_at = t; }
        applied_A = result.approved_charge_current_A;
    }
    check((detected_at > 300U) && (detected_at < 300U + 600U), "frozen sensor: flagged within 10 min of a 1C charge");
    check(fabs(pack.state[2].soc - pack.true_soc[2]) < 0.01, "frozen sensor: SOC estimate not dragged by the frozen reading");
    check(result.decision == EV_DECISION_FAULT, "frozen sensor: fault latched while the reading stays frozen");
}

static void test_frozen_sensor_detected_when_charge_starts(void) {
    Pack pack;
    EvPackResult result;
    float frozen_V;
    uint32_t cycle;
    bool detected = false;
    pack_uniform(&pack, 0.5, 25.0f);
    plant_step(&pack, 0.0f, DT_S);
    for (cycle = 0U; cycle < 3U; ++cycle) {   /* at rest: the tracking window starts */
        result = ev_step(pack.config, pack.state, 0.0f, 0.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    }
    result = ev_step(pack.config, pack.state, 0.0f, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    frozen_V = pack.voltages_V[1];
    for (cycle = 0U; cycle < EV_TRACK_FAULT_CYCLES + 1U; ++cycle) {
        float applied_A = result.approved_charge_current_A;
        plant_step(&pack, applied_A, DT_S);
        pack.voltages_V[1] = frozen_V;   /* no R0 step when the current starts */
        result = ev_step(pack.config, pack.state, applied_A, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
        detected = detected || ((result.fault_flags & EV_FAULT_VOLTAGE_PLAUSIBILITY) != 0U);
    }
    check(detected, "frozen sensor: missing R0 step at charge start flagged within 4 cycles");
}

/* Round-4 finding: +/-1.5 mV of noise on a frozen reading restarted the old dead-band run every
 * cycle; the cell reached 106 % SOC without a fault. */
static void test_noisy_frozen_sensor_is_detected(void) {
    const float amplitudes_mV[] = {1.5f, 5.0f, 15.0f};
    size_t k;
    for (k = 0U; k < sizeof(amplitudes_mV) / sizeof(amplitudes_mV[0]); ++k) {
        const double socs[EV_NUM_CELLS] = {0.30, 0.30, 0.40, 0.30};
        const float temps[EV_NUM_CELLS] = {25.0f, 25.0f, 25.0f, 25.0f};
        Pack pack;
        EvPackResult result;
        float applied_A = 0.0f;
        float frozen_V = 0.0f;
        double max_true_soc = 0.0;
        bool detected = false;
        uint32_t t;
        pack_init(&pack, socs, temps);
        for (t = 0U; t < 4U * 3600U; ++t) {
            plant_step(&pack, applied_A, DT_S);
            if (t == 300U) { frozen_V = pack.voltages_V[2]; }
            if (t >= 300U) {
                pack.voltages_V[2] = frozen_V + ((((t & 1U) != 0U) ? 1.0f : -1.0f) * amplitudes_mV[k] * 1.0e-3f);
            }
            result = ev_step(pack.config, pack.state, applied_A, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
            detected = detected || ((result.fault_flags & EV_FAULT_VOLTAGE_PLAUSIBILITY) != 0U);
            if (pack.true_soc[2] > max_true_soc) { max_true_soc = pack.true_soc[2]; }
            applied_A = result.approved_charge_current_A;
        }
        check(detected, "noisy frozen sensor: plausibility fault raised");
        check(max_true_soc < 0.55, "noisy frozen sensor: charging stops soon after the freeze");
        check(fabs(pack.state[2].soc - pack.true_soc[2]) < 0.02, "noisy frozen sensor: SOC estimate not dragged");
    }
}

/* Round-4 finding: the old fault cleared at rest, charging resumed in bursts and each rest pulled
 * the SOC estimate back toward the frozen reading. */
static void test_frozen_sensor_fault_stays_latched(void) {
    Pack pack;
    EvPackResult result;
    float applied_A = 0.0f;
    float frozen_V = 0.0f;
    uint32_t toggles = 0U;
    uint32_t t;
    bool was_fault = false;
    pack_uniform(&pack, 0.80, 25.0f);
    for (t = 0U; t < 6U * 3600U; ++t) {
        plant_step(&pack, applied_A, DT_S);
        if (t == 10U) { frozen_V = pack.voltages_V[2]; }
        if (t >= 10U) { pack.voltages_V[2] = frozen_V; }   /* frozen at its rested value */
        result = ev_step(pack.config, pack.state, applied_A, (t >= 60U) ? 10.0f : 0.0f,
                         pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
        if ((result.decision == EV_DECISION_FAULT) != was_fault) { ++toggles; }
        was_fault = (result.decision == EV_DECISION_FAULT);
        applied_A = result.approved_charge_current_A;
    }
    check(was_fault && (toggles <= 1U), "frozen sensor: fault latched, no charge bursts");
    check(pack.true_soc[2] < 0.81, "frozen sensor: no blind charge after the fault");
}

/* Round-4 review: the tracking tolerance must cover realistic plant/model mismatch, including
 * an unconverged EKF, or a healthy sensor is faulted (and the EKF stops learning). */
static void test_tracking_tolerates_parameter_mismatch(void) {
    /* capacity, R1, tau scale; initial SOC error */
    const float cases[][4] = {
        {0.85f, 1.0f, 1.0f, 0.0f}, {1.1f, 1.0f, 1.0f, 0.0f}, {1.0f, 2.0f, 1.0f, 0.0f}, {1.0f, 1.0f, 0.3f, 0.0f},
        {1.0f, 1.0f, 3.0f, 0.0f}, {0.9f, 1.5f, 2.0f, 0.0f}, {1.0f, 0.5f, 0.5f, 0.0f}, {1.0f, 1.0f, 1.0f, 0.15f},
    };
    bool no_false_fault = true;
    bool estimate_ok = true;
    size_t k;
    for (k = 0U; k < sizeof(cases) / sizeof(cases[0]); ++k) {
        Pack pack;
        EvPackResult result;
        float applied_A = 0.0f;
        uint32_t t;
        uint8_t i;
        pack_uniform(&pack, 0.2, 25.0f);
        pack.plant_capacity_scale = cases[k][0];
        pack.plant_r1_scale = cases[k][1];
        pack.plant_tau_scale = cases[k][2];
        pack.adc_lsb_V = 0.001f;
        for (i = 0U; i < EV_NUM_CELLS; ++i) {
            pack.state[i].soc = 0.2 + (double)cases[k][3];
            if (cases[k][3] > 0.0f) { pack.state[i].covariance[0] = 0.04; }
        }
        for (t = 0U; t < 5U * 3600U; ++t) {
            /* CC-CV charge to the limit; with an SOC error, 1C pulses and rests first so the EKF can learn */
            float request_A = ((cases[k][3] > 0.0f) && (t < 3600U) && (((t / 600U) % 2U) == 1U)) ? 0.0f : 10.0f;
            result = pack_step(&pack, applied_A, request_A);
            no_false_fault = no_false_fault && ((result.fault_flags & EV_FAULT_VOLTAGE_PLAUSIBILITY) == 0U);
            applied_A = result.approved_charge_current_A;
        }
        estimate_ok = estimate_ok && (fabs(pack.state[0].soc - pack.true_soc[0]) < 0.05);
    }
    check(no_false_fault, "tracking: no false fault with capacity/R1/tau mismatch, 1 mV ADC or an unconverged EKF");
    check(estimate_ok, "tracking: EKF still learns SOC under those mismatches");
}

static void test_no_false_frozen_fault_with_mismatch(void) {
    const float scales[] = {0.8f, 1.2f};
    size_t k;
    for (k = 0U; k < sizeof(scales) / sizeof(scales[0]); ++k) {
        const double socs[EV_NUM_CELLS] = {0.20, 0.30, 0.25, 0.22};
        const float temps[EV_NUM_CELLS] = {25.0f, 5.0f, 40.0f, 15.0f};
        Pack pack;
        EvPackResult result;
        float applied_A = 0.0f;
        uint16_t flags = 0U;
        uint32_t t;
        pack_init(&pack, socs, temps);
        pack.plant_r0_scale = scales[k];
        for (t = 0U; t < 5U * 3600U; ++t) {
            result = pack_step(&pack, applied_A, 10.0f);
            flags |= result.fault_flags;
            applied_A = result.approved_charge_current_A;
            if (applied_A < 0.05f) { break; }
        }
        check((flags & EV_FAULT_VOLTAGE_PLAUSIBILITY) == 0U, "frozen sensor: no false fault through CC-CV taper with +/-20 % R0 error");
        check(pack.true_soc[1] > 0.85, "frozen sensor: pack still charges with R0 error");
    }
}

/* Round-3 finding: Joule-heating credit let a larger current pass the 0 degC and plating checks
 * that a smaller one failed, so the approved current fell as the request rose. */
static void test_limit_monotone_in_request(void) {
    /* soc, horizon_s, coolant degC, cell degC, heat transfer W/m2/K: cases that failed before */
    const float cases[][5] = {
        {0.05f, 10.0f, -7.5f, 11.5f, 1000.0f}, {0.05f, 30.0f, -30.0f, 2.5f, 10.0f}, {0.05f, 30.0f, -7.5f, 4.0f, 100.0f},
        {0.60f, 120.0f, -22.5f, 34.0f, 100.0f}, {0.60f, 10.0f, 0.0f, 4.0f, 1000.0f}, {0.30f, 600.0f, -7.5f, 4.0f, 10.0f},
    };
    bool monotone = true;
    bool interval = true;
    size_t k;
    for (k = 0U; k < sizeof(cases) / sizeof(cases[0]); ++k) {
        EvCellConfig config = lgm50_config();
        EvCellState state;
        float previous_A = 0.0f;
        float current_A;
        bool was_safe = true;
        config.coolant_temperature_C = cases[k][2];
        config.heat_transfer_W_m2_K = cases[k][4];
        ev_cell_init(&state, (double)cases[k][0], cases[k][3]);
        state.covariance[0] = 1.0e-4;
        for (current_A = 0.05f; current_A <= 40.0f; current_A += 0.05f) {
            EvCellPrediction p;
            float approved_A = ev_cell_max_charge_current(&config, &state, current_A, cases[k][1], NULL, NULL);
            monotone = monotone && (approved_A >= previous_A - 1.0e-4f);
            previous_A = approved_A;
            (void)ev_predict_cell(&config, &state, current_A, cases[k][1], &p);
            interval = interval && (was_safe || !p.safe);
            was_safe = p.safe;
        }
    }
    check(monotone, "search: approved current never falls when the request rises (cold coolant)");
    check(interval, "physics: once unsafe, a larger current is never safe again");
}

/* Runs ev_step with a charger that slews toward the approved current at slew_A_s (or ignores it
 * when slew_A_s < 0 and stays at stuck_A). Returns the union of fault flags. */
static uint16_t run_charger(Pack *pack, float dt_s, float slew_A_s, float stuck_A, float request_A, float stop_at_s, float end_s) {
    EvPackResult result;
    float applied_A = 0.0f;
    float command_A = 0.0f;
    float t_s;
    uint16_t flags = 0U;
    for (t_s = 0.0f; t_s < end_s; t_s += dt_s) {
        if (slew_A_s < 0.0f) {
            applied_A = (t_s >= stop_at_s) ? stuck_A : command_A;
        } else {
            float step_A = slew_A_s * dt_s;
            float target_A = command_A;
            applied_A = (applied_A > target_A) ? fmaxf(applied_A - step_A, target_A) : fminf(applied_A + step_A, target_A);
        }
        plant_step(pack, applied_A, dt_s);
        result = ev_step(pack->config, pack->state, applied_A, (t_s >= stop_at_s) ? 0.0f : request_A,
                         pack->voltages_V, pack->temperatures_C, dt_s, HORIZON_S);
        flags |= result.fault_flags;
        command_A = result.approved_charge_current_A;
    }
    return flags;
}

/* Round-3 finding: nothing checked that the measured current follows the approved command.
 * Round-4 finding: a lag counted in cycles tripped at short control periods. */
static void test_over_current_fault(void) {
    Pack pack;
    EvPackResult result;
    uint32_t cycle;
    pack_uniform(&pack, 0.5, 25.0f);
    check((run_charger(&pack, 0.1f, 5.0f, 0.0f, 3.0f, 60.0f, 70.0f) & EV_FAULT_OVER_CURRENT) == 0U,
          "over-current: 5 A/s ramp-down after a stop at dt 0.1 s is tolerated");
    pack_uniform(&pack, 0.5, 25.0f);
    check((run_charger(&pack, 0.1f, 1.6f, 0.0f, 3.0f, 60.0f, 70.0f) & EV_FAULT_OVER_CURRENT) == 0U,
          "over-current: 1.6 A/s ramp-down (1.9 s) at dt 0.1 s is tolerated");
    pack_uniform(&pack, 0.5, 25.0f);
    check((run_charger(&pack, 1.0f, 2.0f, 0.0f, 3.0f, 60.0f, 70.0f) & EV_FAULT_OVER_CURRENT) == 0U,
          "over-current: 2 A/s ramp-down at dt 1 s is tolerated");
    pack_uniform(&pack, 0.5, 25.0f);
    check((run_charger(&pack, 0.1f, -1.0f, 3.0f, 3.0f, 60.0f, 70.0f) & EV_FAULT_OVER_CURRENT) != 0U,
          "over-current: charger that keeps 3 A after a stop faults (dt 0.1 s)");
    pack_uniform(&pack, 0.5, 25.0f);
    check((run_charger(&pack, 10.0f, -1.0f, 3.0f, 3.0f, 60.0f, 100.0f) & EV_FAULT_OVER_CURRENT) != 0U,
          "over-current: charger that keeps 3 A after a stop faults (dt 10 s)");

    pack_uniform(&pack, 0.5, 25.0f);
    result = pack_step(&pack, 0.0f, 2.0f);
    check(result.approved_charge_current_A == 2.0f, "over-current: 2 A approved");
    result = pack_step(&pack, (2.0f * 1.05f) + 0.4f, 2.0f);
    check(result.decision != EV_DECISION_FAULT, "over-current: small excess tolerated");
    result = pack_step(&pack, 6.0f, 2.0f);
    check((result.decision == EV_DECISION_FAULT) && ((result.fault_flags & EV_FAULT_OVER_CURRENT) != 0U),
          "over-current: charger ignoring the limit faults");
    for (cycle = 0U; cycle < EV_FAULT_HOLD_CYCLES + 3U; ++cycle) {
        result = pack_step(&pack, 0.0f, 2.0f);
    }
    check(result.decision != EV_DECISION_FAULT, "over-current: clears once the charger stops");
}

static void test_config_fault_marks_soc_stale(void) {
    Pack pack;
    pack_uniform(&pack, 0.5, 25.0f);
    plant_step(&pack, 0.0f, DT_S);
    pack.config[1].rc_time_constant_s = 0.0f;
    (void)ev_step(pack.config, pack.state, 50.0f, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    check(pack.state[1].covariance[0] >= 0.25 - 1.0e-12, "config fault: SOC marked unknown while charge is not counted");
}

static void test_input_fault_is_held(void) {
    Pack pack;
    EvPackResult result;
    uint32_t cycle;
    bool held = true;
    pack_uniform(&pack, 0.5, 25.0f);
    plant_step(&pack, 0.0f, DT_S);
    (void)ev_step(pack.config, pack.state, make_nan(), 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    for (cycle = 1U; cycle < EV_FAULT_HOLD_CYCLES; ++cycle) {
        result = ev_step(pack.config, pack.state, 0.0f, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
        held = held && (result.decision == EV_DECISION_FAULT);
    }
    check(held, "input fault: a single bad current sample is held, not toggled");
    result = ev_step(pack.config, pack.state, 0.0f, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    check(result.decision != EV_DECISION_FAULT, "input fault: clears after the hold");
}

static void test_horizon_covers_control_period(void) {
    Pack a;
    Pack b;
    EvPackResult short_horizon;
    EvPackResult long_horizon;
    pack_uniform(&a, 0.85, 25.0f);
    pack_uniform(&b, 0.85, 25.0f);
    plant_step(&a, 0.0f, 10.0f);
    plant_step(&b, 0.0f, 10.0f);
    short_horizon = ev_step(a.config, a.state, 0.0f, 20.0f, a.voltages_V, a.temperatures_C, 10.0f, 0.001f);
    long_horizon = ev_step(b.config, b.state, 0.0f, 20.0f, b.voltages_V, b.temperatures_C, 10.0f, 10.0f);
    check(short_horizon.approved_charge_current_A == long_horizon.approved_charge_current_A,
          "horizon: never shorter than the control period");
}

static void test_cold_charge_hysteresis(void) {
    Pack pack;
    EvPackResult result;
    pack_uniform(&pack, 0.5, 25.0f);
    plant_step(&pack, 0.0f, DT_S);
    pack.temperatures_C[0] = 1.0f;   /* within the 2 degC sensor margin of 0 degC */
    result = ev_step(pack.config, pack.state, 0.0f, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    check((result.approved_charge_current_A == 0.0f) && ((result.limiting_constraints & EV_CONSTRAINT_TEMPERATURE_LOW) != 0U),
          "cold: 1 degC measured blocks charging (sensor margin)");
    pack.temperatures_C[0] = 4.0f;
    result = ev_step(pack.config, pack.state, 0.0f, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    check(result.approved_charge_current_A == 0.0f, "cold: block held at 4 degC (hysteresis)");
    pack.temperatures_C[0] = 6.0f;
    result = ev_step(pack.config, pack.state, 0.0f, 5.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    check(result.approved_charge_current_A > 0.0f, "cold: released at 6 degC");
}

static void test_end_of_life_and_growth_clamp(void) {
    Pack pack;
    EvPackResult result;
    pack_uniform(&pack, 0.5, 25.0f);
    plant_step(&pack, 0.0f, DT_S);
    pack.state[0].soh = 0.5000000001;
    result = ev_step(pack.config, pack.state, 5.0f, 5.0f, pack.voltages_V, pack.temperatures_C, 10.0f, HORIZON_S);
    check((result.fault_flags & EV_FAULT_END_OF_LIFE) != 0U, "ageing: end of life is reported");

    /* A 0.1 Ah, 1 ohm cell cycled at 1000 A ages by ~3e-3 per 10 s step: growth would pass 1 ohm. */
    pack_uniform(&pack, 0.5, 25.0f);
    plant_step(&pack, 0.0f, DT_S);
    pack.config[0].nominal_capacity_Ah = 0.1f;
    pack.config[0].resistance_ohm = 1.0f;
    pack.state[0].resistance_growth_ohm = 0.999f;
    (void)ev_step(pack.config, pack.state, 1000.0f, 0.0f, pack.voltages_V, pack.temperatures_C, 10.0f, HORIZON_S);
    check(pack.state[0].resistance_growth_ohm <= 1.0f, "ageing: resistance growth stays inside its valid range");
    result = ev_step(pack.config, pack.state, 0.0f, 0.0f, pack.voltages_V, pack.temperatures_C, 10.0f, HORIZON_S);
    check((result.fault_flags & EV_FAULT_STATE) == 0U, "ageing: no reset loop from resistance growth");
}

static void test_ageing_and_cold_reduce_the_limit(void) {
    EvCellConfig config = lgm50_config();
    EvCellState fresh;
    EvCellState faded;
    EvCellState filmed;
    EvCellState cold;
    EvCellPrediction p_fresh;
    EvCellPrediction p_faded;
    EvCellPrediction p_filmed;
    EvCellPrediction p_cold;
    uint8_t constraints = 0U;
    ev_cell_init(&fresh, 0.5, 25.0f);
    ev_cell_init(&faded, 0.5, 25.0f);
    ev_cell_init(&filmed, 0.5, 25.0f);
    ev_cell_init(&cold, 0.5, 5.0f);
    faded.soh = 0.8;
    filmed.resistance_growth_ohm = 0.018f;
    (void)ev_predict_cell(&config, &fresh, 3.0f, HORIZON_S, &p_fresh);
    (void)ev_predict_cell(&config, &faded, 3.0f, HORIZON_S, &p_faded);
    (void)ev_predict_cell(&config, &filmed, 3.0f, HORIZON_S, &p_filmed);
    check(p_faded.anode_potential_V < p_fresh.anode_potential_V - 0.002f, "ageing: active-material loss lowers the anode potential");
    check(p_filmed.anode_potential_V < p_fresh.anode_potential_V - 0.05f, "ageing: SEI film drop lowers the anode potential");
    (void)ev_predict_cell(&config, &fresh, 3.0f, HORIZON_S, &p_fresh);
    (void)ev_predict_cell(&config, &cold, 3.0f, HORIZON_S, &p_cold);
    check(p_cold.terminal_voltage_V > p_fresh.terminal_voltage_V, "thermal: resistances rise in the cold");
    check(ev_cell_max_charge_current(NULL, &fresh, 5.0f, HORIZON_S, &constraints, NULL) == 0.0f
          && (constraints == EV_CONSTRAINT_INVALID_INPUT), "api: single-cell limit rejects invalid arguments");
}

/* Round-4 physics finding: steady gradients at the new, lower current overstated the anode
 * potential right after a step-down. */
static void test_gradients_remember_previous_current(void) {
    EvCellConfig config = lgm50_config();
    EvCellState fresh;
    EvCellState after_step_down;
    EvCellPrediction p_fresh;
    EvCellPrediction p_after;
    ev_cell_init(&fresh, 0.5, 5.0f);
    ev_cell_init(&after_step_down, 0.5, 5.0f);
    after_step_down.last_current_A = 6.0f;
    (void)ev_predict_cell(&config, &fresh, 1.5f, HORIZON_S, &p_fresh);
    (void)ev_predict_cell(&config, &after_step_down, 1.5f, HORIZON_S, &p_after);
    check(p_after.anode_potential_V < p_fresh.anode_potential_V - 0.005f, "physics: gradients after a 6 A -> 1.5 A step use the 6 A history");
    after_step_down.last_current_A = -6.0f;   /* discharge history adds nothing */
    (void)ev_predict_cell(&config, &after_step_down, 1.5f, HORIZON_S, &p_after);
    check(p_after.anode_potential_V == p_fresh.anode_potential_V, "physics: discharge history does not change the charge prediction");
}

/* Round-4 physics finding: at exactly 45 degC equilibrium the old thermal formula rounded above
 * 45 degC at a small current and back to 45 degC at a larger one. */
static void test_thermal_prediction_exact_at_equilibrium(void) {
    EvCellConfig config = lgm50_config();
    EvCellState state;
    EvCellPrediction p;
    float current_A;
    float previous_C = -1.0e9f;
    bool monotone = true;
    config.coolant_temperature_C = 45.0f;
    ev_cell_init(&state, 0.3, 45.0f);
    (void)ev_predict_cell(&config, &state, 0.0f, HORIZON_S, &p);
    check(p.predicted_temperature_C == 45.0f, "thermal: no heat at equilibrium keeps exactly 45 degC");
    for (current_A = 0.0f; current_A <= 2.0f; current_A += 0.005f) {
        (void)ev_predict_cell(&config, &state, current_A, HORIZON_S, &p);
        monotone = monotone && (p.predicted_temperature_C >= previous_C);
        previous_C = p.predicted_temperature_C;
    }
    check(monotone, "thermal: prediction monotone in current at the 45 degC edge");
}

static void test_single_cell_limit_matches_pack(void) {
    Pack pack;
    EvPackResult result;
    float single_A;
    pack_uniform(&pack, 0.6, 25.0f);
    plant_step(&pack, 0.0f, DT_S);
    result = ev_step(pack.config, pack.state, 0.0f, 20.0f, pack.voltages_V, pack.temperatures_C, DT_S, HORIZON_S);
    single_A = ev_cell_max_charge_current(&pack.config[0], &pack.state[0], 20.0f, HORIZON_S, NULL, NULL);
    check(fabsf(single_A - result.approved_charge_current_A) < 1.0e-5f, "api: single-cell limit equals the pack limit for identical cells");
}

static bool self_test(void) {
    test_ocv_table();
    test_closed_loop_charge_to_full();
    test_allow_small_request();
    test_monotone_constraints();
    test_full_cell_blocks_charge();
    test_soc_uncertainty_derates();
    test_cold_cell_blocks_charge();
    test_predicted_overtemperature_limits();
    test_overtemperature_fault_hysteresis_and_hold();
    test_invalid_inputs();
    test_ageing();
    test_ekf_converges_with_model_mismatch();
    test_terminal_voltage_bound_covers_whole_horizon();
    test_over_voltage_fault();
    test_stuck_voltage_sensor_is_detected();
    test_frozen_in_range_sensor_is_detected();
    test_frozen_sensor_detected_when_charge_starts();
    test_noisy_frozen_sensor_is_detected();
    test_frozen_sensor_fault_stays_latched();
    test_no_false_frozen_fault_with_mismatch();
    test_tracking_tolerates_parameter_mismatch();
    test_limit_monotone_in_request();
    test_over_current_fault();
    test_config_fault_marks_soc_stale();
    test_input_fault_is_held();
    test_horizon_covers_control_period();
    test_cold_charge_hysteresis();
    test_end_of_life_and_growth_clamp();
    test_ageing_and_cold_reduce_the_limit();
    test_single_cell_limit_matches_pack();
    test_gradients_remember_previous_current();
    test_thermal_prediction_exact_at_equilibrium();
    return g_failures == 0;
}

static const char *decision_name(EvDecision decision) {
    switch (decision) {
        case EV_DECISION_ALLOW: return "ALLOW";
        case EV_DECISION_LIMIT: return "LIMIT";
        default: return "FAULT";
    }
}

static int demo(void) {
    const double socs[EV_NUM_CELLS] = {0.55, 0.78, 0.70, 0.65};
    const float temps[EV_NUM_CELLS] = {25.0f, 8.0f, 30.0f, 22.0f};
    Pack pack;
    EvPackResult result;
    float applied_A = 0.0f;
    uint32_t step;
    pack_init(&pack, socs, temps);
    puts("EV MODIFY - BMS charge-limit demo (4 x LG M50, 30 A request)");
    puts("step,decision,request_A,approved_A,limit_cell,constraints,anode_V,terminal_V,temp_C,soc_upper");
    for (step = 0U; step < 12U; ++step) {
        result = pack_step(&pack, applied_A, 30.0f);
        printf("%u,%s,%.3f,%.3f,%u,0x%02X,%.5f,%.4f,%.3f,%.4f\n", (unsigned)step, decision_name(result.decision),
               (double)result.requested_charge_current_A, (double)result.approved_charge_current_A,
               (unsigned)result.limiting_cell_index, (unsigned)result.limiting_constraints,
               (double)result.limiting_anode_potential_V, (double)result.limiting_terminal_voltage_V,
               (double)result.limiting_temperature_C, (double)result.limiting_soc);
        if (result.decision == EV_DECISION_FAULT) {
            return 1;
        }
        applied_A = result.approved_charge_current_A;
    }
    return 0;
}

int main(int argc, char *argv[]) {
    if ((argc == 2) && (strcmp(argv[1], "--self-test") == 0)) {
        bool passed = self_test();
        printf("SELF-TEST: %s (%d checks, %d failed)\n", passed ? "PASS" : "FAIL", g_checks, g_failures);
        return passed ? 0 : 1;
    }
    if (argc != 1) {
        fprintf(stderr, "usage: test_solver [--self-test]\n");
        return 2;
    }
    return demo();
}
