/* Smoke test for a downstream project linking ev_modify::ev_modify. */
#include <stdio.h>

#include "ev_modify.h"

int main(void) {
    EvCellConfig config[EV_NUM_CELLS];
    EvCellState state[EV_NUM_CELLS];
    float voltages_V[EV_NUM_CELLS];
    float temperatures_C[EV_NUM_CELLS];
    EvPackResult result;
    unsigned int i;

    for (i = 0U; i < EV_NUM_CELLS; ++i) {
        config[i].nominal_capacity_Ah = 50.0f;
        config[i].resistance_ohm = 0.002f;
        config[i].rc_resistance_ohm = 0.001f;
        config[i].rc_time_constant_s = 30.0f;
        config[i].heat_transfer_W_m2_K = 10.0f;
        config[i].heat_capacity_J_kg_K = 1100.0f;
        config[i].coolant_temperature_C = 25.0f;
        if (!ev_config_is_valid(&config[i])) {
            fprintf(stderr, "configuration rejected\n");
            return 1;
        }
        ev_cell_init(&state[i], 0.5, 25.0f);
        voltages_V[i] = 3.7f;
        temperatures_C[i] = 25.0f;
    }
    result = ev_step(config, state, 0.0f, 10.0f, voltages_V, temperatures_C, 1.0f, 30.0f);
    printf("decision %u, approved %.3f A\n", (unsigned)result.decision, (double)result.approved_charge_current_A);
    return (result.decision <= EV_DECISION_FAULT) ? 0 : 1;
}
