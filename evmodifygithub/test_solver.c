#include <stdio.h>
#include <string.h>
#include "spme_v7_apex_solver.h"

static const char *ev_decision_name(EvDecision decision) {
    switch (decision) {
        case EV_DECISION_ALLOW: return "ALLOW";
        case EV_DECISION_LIMIT: return "LIMIT";
        default: return "FAULT";
    }
}

static bool ev_self_test(void) {
    EvCellState cells[EV_NUM_CELLS] = {
        {0.55, 0.01, 25.0f, 1.00f, 0.0f, 5.0f, 0.020f, true, 12.5f, 950.0f, 0.75f},
        {0.78, 0.01, -10.0f, 0.93f, 0.0f, 5.0f, 0.026f, true, 12.5f, 950.0f, 0.75f},
        {0.70, 0.01, 30.0f, 0.90f, 0.0f, 5.0f, 0.028f, true, 12.5f, 950.0f, 0.75f},
        {0.65, 0.01, 22.0f, 0.97f, 0.0f, 5.0f, 0.021f, true, 12.5f, 950.0f, 0.75f}
    };
    EvPackResult result;
    const double starting_soc = cells[0].soc;
    uint32_t iteration;
    float applied_current_A = 30.0f;
    float requested_current_A = 30.0f;
    float mock_voltages_V[EV_NUM_CELLS];
    float mock_temperatures_C[EV_NUM_CELLS];

    for (iteration = 0U; iteration < 60U; ++iteration) {
        for (uint8_t i = 0U; i < EV_NUM_CELLS; i++) {
            // Fake measurement based on ohmic drop
            mock_voltages_V[i] = 3.6f + applied_current_A * cells[i].resistance_ohm; 
            mock_temperatures_C[i] = cells[i].temperature_C;
        }
        
        result = ev_step(cells, applied_current_A + (float)((int)(iteration % 3U) - 1), requested_current_A, mock_voltages_V, mock_temperatures_C, 1.0f, 30.0f);
        if (result.decision == EV_DECISION_FAULT) {
            return false;
        }
        applied_current_A = result.approved_charge_current_A;
    }
    return result.limiting_cell_index < EV_NUM_CELLS;
}

int main(int argc, char *argv[]) {
    EvCellState cells[EV_NUM_CELLS] = {
        {0.55, 0.01, 25.0f, 1.00f, 0.0f, 5.0f, 0.020f, true, 12.5f, 950.0f, 0.75f},
        {0.78, 0.01, -10.0f, 0.93f, 0.0f, 5.0f, 0.026f, true, 12.5f, 950.0f, 0.75f},
        {0.70, 0.01, 30.0f, 0.90f, 0.0f, 5.0f, 0.028f, true, 12.5f, 950.0f, 0.75f},
        {0.65, 0.01, 22.0f, 0.97f, 0.0f, 5.0f, 0.021f, true, 12.5f, 950.0f, 0.75f}
    };
    EvPackResult result;
    uint32_t step;
    float mock_voltages_V[EV_NUM_CELLS];
    float mock_temperatures_C[EV_NUM_CELLS];

    if (argc == 2 && strcmp(argv[1], "--self-test") == 0) {
        if (ev_self_test()) {
            puts("SELF-TEST: PASS");
            return 0;
        }
        puts("SELF-TEST: FAIL");
        return 1;
    }

    puts("EV MODIFY - BMS Firmware Simulator");
    puts("step,decision,request_A,approved_A,limit_cell,anode_V,temp_C,soc,soh");
    for (step = 0U; step < 12U; ++step) {
        const float noisy_measurement_A = 30.0f + (float)((int)(step % 3U) - 1);
        const float requested_current_A = 30.0f;
        for (uint8_t i = 0U; i < EV_NUM_CELLS; i++) {
            mock_voltages_V[i] = 3.6f + noisy_measurement_A * cells[i].resistance_ohm;
            mock_temperatures_C[i] = cells[i].temperature_C;
        }
        
        result = ev_step(cells, noisy_measurement_A, requested_current_A, mock_voltages_V, mock_temperatures_C, 1.0f, 30.0f);
        printf("%u,%s,%.3f,%.3f,%u,%.5f,%.3f,%.5f,%.5f\n", step,
               ev_decision_name(result.decision), result.requested_charge_current_A,
               result.approved_charge_current_A, result.limiting_cell_index,
               result.limiting_anode_potential_V, result.limiting_temperature_C,
               result.limiting_soc, cells[result.limiting_cell_index].soh);
        if (result.decision == EV_DECISION_FAULT) {
            return 1;
        }
    }
    return 0;
}
