#include <math.h>
#include <stdio.h>
#include <string.h>
#include "spme_v7_apex_solver.h"

#define EV_DEMO_REQUEST_A 30.0f
#define EV_DEMO_DT_S 1.0f
#define EV_DEMO_HORIZON_S 30.0f

static const char *ev_decision_name(EvDecision decision) {
    switch (decision) {
        case EV_DECISION_ALLOW: return "ALLOW";
        case EV_DECISION_LIMIT: return "LIMIT";
        default: return "FAULT";
    }
}

static void ev_init_cells(EvCellState cells[EV_NUM_CELLS]) {
    const EvCellState initial[EV_NUM_CELLS] = {
        {0.55, 0.01, 25.0f, 1.00f, 0.0f, 5.0f, 0.020f, true, 12.5f, 950.0f, 0.75f},
        {0.78, 0.01, -10.0f, 0.93f, 0.0f, 5.0f, 0.026f, true, 12.5f, 950.0f, 0.75f},
        {0.70, 0.01, 30.0f, 0.90f, 0.0f, 5.0f, 0.028f, true, 12.5f, 950.0f, 0.75f},
        {0.65, 0.01, 22.0f, 0.97f, 0.0f, 5.0f, 0.021f, true, 12.5f, 950.0f, 0.75f}
    };
    memcpy(cells, initial, sizeof(initial));
}

/* Simple plant: true SOC is coulomb counted; terminal voltage = OCV(true SOC) + I*R. */
typedef struct {
    double soc[EV_NUM_CELLS];
    float voltages_V[EV_NUM_CELLS];
    float temperatures_C[EV_NUM_CELLS];
} EvPlant;

static void ev_plant_init(EvPlant *plant, const EvCellState cells[EV_NUM_CELLS]) {
    uint8_t i;
    for (i = 0U; i < EV_NUM_CELLS; ++i) {
        plant->soc[i] = cells[i].soc;
        plant->temperatures_C[i] = cells[i].temperature_C;
        plant->voltages_V[i] = 0.0f;
    }
}

static void ev_plant_step(EvPlant *plant, const EvCellState cells[EV_NUM_CELLS], float current_A, float dt_s) {
    uint8_t i;
    for (i = 0U; i < EV_NUM_CELLS; ++i) {
        float ocv = 0.0f;
        double capacity_Ah = (double)cells[i].nominal_capacity_Ah * (double)cells[i].soh;
        plant->soc[i] += ((double)current_A * (double)dt_s) / (capacity_Ah * 3600.0);
        ev_full_cell_ocv_V_and_dv_dsoc((float)plant->soc[i], &ocv, NULL);
        plant->voltages_V[i] = ocv + (current_A * cells[i].resistance_ohm);
    }
}

static int ev_failures = 0;

static void ev_check(bool condition, const char *name) {
    if (!condition) {
        printf("  FAIL: %s\n", name);
        ++ev_failures;
    }
}

static void test_closed_loop_charge(void) {
    EvCellState cells[EV_NUM_CELLS];
    EvPlant plant;
    EvPackResult result;
    float applied_A = 0.0f;
    uint32_t iteration;
    ev_init_cells(cells);
    ev_plant_init(&plant, cells);

    for (iteration = 0U; iteration < 120U; ++iteration) {
        ev_plant_step(&plant, cells, applied_A, EV_DEMO_DT_S);
        result = ev_step(cells, applied_A, EV_DEMO_REQUEST_A, plant.voltages_V, plant.temperatures_C, EV_DEMO_DT_S, EV_DEMO_HORIZON_S);
        ev_check(result.decision != EV_DECISION_FAULT, "closed loop: no fault");
        ev_check(result.approved_charge_current_A >= 0.0f, "closed loop: approved >= 0");
        ev_check(result.approved_charge_current_A <= EV_DEMO_REQUEST_A, "closed loop: approved <= request");
        ev_check(result.limiting_cell_index < EV_NUM_CELLS, "closed loop: limiting index in range");
        ev_check(isfinite(result.limiting_anode_potential_V), "closed loop: finite anode potential");
        if (result.decision == EV_DECISION_FAULT) {
            return;
        }
        applied_A = result.approved_charge_current_A;
    }
    /* The -10 degC cell must dominate a 6C request. */
    ev_check(result.decision == EV_DECISION_LIMIT, "closed loop: 6C request is limited");
    ev_check(result.limiting_anode_potential_V >= EV_MIN_SAFE_ANODE_POTENTIAL_V - 1.0e-4f, "closed loop: approved current keeps anode margin");
}

static void test_physics_trends(void) {
    EvCellState cells[EV_NUM_CELLS];
    EvCellPrediction low;
    EvCellPrediction high;
    EvCellPrediction cold;
    EvCellState cold_cell;
    ev_init_cells(cells);
    cold_cell = cells[0];
    cold_cell.temperature_C = -10.0f;

    ev_check(ev_predict_cell(&cells[0], 2.0f, EV_DEMO_HORIZON_S, &low), "physics: low-current prediction valid");
    ev_check(ev_predict_cell(&cells[0], 10.0f, EV_DEMO_HORIZON_S, &high), "physics: high-current prediction valid");
    ev_check(ev_predict_cell(&cold_cell, 2.0f, EV_DEMO_HORIZON_S, &cold), "physics: cold prediction valid");
    ev_check(high.anode_potential_V < low.anode_potential_V, "physics: more current lowers anode potential");
    ev_check(cold.anode_potential_V < low.anode_potential_V, "physics: colder cell lowers anode potential");
    ev_check(high.predicted_temperature_C > low.predicted_temperature_C, "physics: more current heats more");
    ev_check(high.predicted_soc > low.predicted_soc, "physics: charging raises SOC");
}

static void test_full_cell_blocks_charge(void) {
    EvCellState cells[EV_NUM_CELLS];
    EvPlant plant;
    EvPackResult result;
    ev_init_cells(cells);
    cells[2].soc = 0.99;
    ev_plant_init(&plant, cells);
    ev_plant_step(&plant, cells, 0.0f, EV_DEMO_DT_S);
    result = ev_step(cells, 0.0f, 5.0f, plant.voltages_V, plant.temperatures_C, EV_DEMO_DT_S, EV_DEMO_HORIZON_S);
    ev_check(result.decision == EV_DECISION_LIMIT, "full cell: decision LIMIT");
    ev_check(result.approved_charge_current_A == 0.0f, "full cell: zero current");
    ev_check(result.limiting_cell_index == 2U, "full cell: limiting cell identified");
}

static void test_overtemperature_faults_and_recovers(void) {
    EvCellState cells[EV_NUM_CELLS];
    EvPlant plant;
    EvPackResult result;
    ev_init_cells(cells);
    ev_plant_init(&plant, cells);
    ev_plant_step(&plant, cells, 0.0f, EV_DEMO_DT_S);

    plant.temperatures_C[3] = 50.0f;
    result = ev_step(cells, 0.0f, 5.0f, plant.voltages_V, plant.temperatures_C, EV_DEMO_DT_S, EV_DEMO_HORIZON_S);
    ev_check(result.decision == EV_DECISION_FAULT, "overtemp: FAULT");
    ev_check(result.approved_charge_current_A == 0.0f, "overtemp: zero current");
    ev_check(result.limiting_cell_index == 3U, "overtemp: faulted cell identified");

    plant.temperatures_C[3] = 30.0f;
    result = ev_step(cells, 0.0f, 5.0f, plant.voltages_V, plant.temperatures_C, EV_DEMO_DT_S, EV_DEMO_HORIZON_S);
    ev_check(result.decision != EV_DECISION_FAULT, "overtemp: recovers once temperature is back in range");
}

static void test_invalid_inputs_fault(void) {
    EvCellState cells[EV_NUM_CELLS];
    EvPlant plant;
    EvPackResult result;
    ev_init_cells(cells);
    ev_plant_init(&plant, cells);
    ev_plant_step(&plant, cells, 0.0f, EV_DEMO_DT_S);

    result = ev_step(cells, NAN, 5.0f, plant.voltages_V, plant.temperatures_C, EV_DEMO_DT_S, EV_DEMO_HORIZON_S);
    ev_check(result.decision == EV_DECISION_FAULT, "invalid: NaN current faults");
    result = ev_step(cells, 0.0f, 5.0f, plant.voltages_V, NULL, EV_DEMO_DT_S, EV_DEMO_HORIZON_S);
    ev_check(result.decision == EV_DECISION_FAULT, "invalid: missing temperatures fault");
    result = ev_step(cells, 0.0f, 5.0f, plant.voltages_V, plant.temperatures_C, NAN, EV_DEMO_HORIZON_S);
    ev_check(result.decision == EV_DECISION_FAULT, "invalid: NaN dt faults");

    plant.voltages_V[1] = NAN;
    result = ev_step(cells, 0.0f, 5.0f, plant.voltages_V, plant.temperatures_C, EV_DEMO_DT_S, EV_DEMO_HORIZON_S);
    ev_check(result.decision == EV_DECISION_FAULT, "invalid: NaN cell voltage faults");
    ev_check(isfinite(cells[1].soc) && (cells[1].soc >= 0.0) && (cells[1].soc <= 1.0), "invalid: SOC state not corrupted by NaN voltage");

    ev_init_cells(cells);
    cells[0].heat_transfer_W_m2_K = 0.0f;
    result = ev_step(cells, 0.0f, 5.0f, plant.voltages_V, plant.temperatures_C, EV_DEMO_DT_S, EV_DEMO_HORIZON_S);
    ev_check(result.decision == EV_DECISION_FAULT, "invalid: zero heat transfer config faults");
}

static bool ev_self_test(void) {
    ev_failures = 0;
    test_closed_loop_charge();
    test_physics_trends();
    test_full_cell_blocks_charge();
    test_overtemperature_faults_and_recovers();
    test_invalid_inputs_fault();
    return ev_failures == 0;
}

int main(int argc, char *argv[]) {
    EvCellState cells[EV_NUM_CELLS];
    EvPlant plant;
    EvPackResult result;
    float applied_A = 0.0f;
    uint32_t step;

    if ((argc == 2) && (strcmp(argv[1], "--self-test") == 0)) {
        if (ev_self_test()) {
            puts("SELF-TEST: PASS");
            return 0;
        }
        printf("SELF-TEST: FAIL (%d checks)\n", ev_failures);
        return 1;
    }

    ev_init_cells(cells);
    ev_plant_init(&plant, cells);
    puts("EV MODIFY - BMS Firmware Simulator");
    puts("step,decision,request_A,approved_A,limit_cell,anode_V,temp_C,soc,soh");
    for (step = 0U; step < 12U; ++step) {
        ev_plant_step(&plant, cells, applied_A, EV_DEMO_DT_S);
        result = ev_step(cells, applied_A, EV_DEMO_REQUEST_A, plant.voltages_V, plant.temperatures_C, EV_DEMO_DT_S, EV_DEMO_HORIZON_S);
        printf("%u,%s,%.3f,%.3f,%u,%.5f,%.3f,%.5f,%.5f\n", (unsigned)step,
               ev_decision_name(result.decision), (double)result.requested_charge_current_A,
               (double)result.approved_charge_current_A, (unsigned)result.limiting_cell_index,
               (double)result.limiting_anode_potential_V, (double)result.limiting_temperature_C,
               (double)result.limiting_soc, (double)cells[result.limiting_cell_index].soh);
        if (result.decision == EV_DECISION_FAULT) {
            return 1;
        }
        applied_A = result.approved_charge_current_A;
    }
    return 0;
}
