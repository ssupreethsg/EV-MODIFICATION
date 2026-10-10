/*
 * Command-line bridge to the C model, used by the Python tools and tests.
 *
 *   ev_cli ocv
 *       stdin : one SOC in [0, 1] per line
 *       stdout: "<ocv_V> <docv_dsoc_V>" per line
 *
 *   ev_cli limit <capacity_Ah> <r0_ohm> <r1_ohm> <tau_s> <coolant_temperature_C>
 *       stdin : "<soc> <soc_sigma> <temperature_C> <rc_voltage_V> <soh> <request_A> <horizon_s>" per line
 *       stdout: "<approved_A> <anode_V> <terminal_V> <temperature_C> <constraints>" per line, from
 *               ev_cell_max_charge_current() on one cell in that state (output flushed per line)
 *
 *   ev_cli ekf <capacity_Ah> <r0_ohm> <r1_ohm> <tau_s> <initial_soc> <initial_temperature_C>
 *       stdin : "<dt_s> <current_A> <voltage_V> <temperature_C>" per line (positive = charge)
 *       stdout: "<soc> <soc_sigma> <rc_voltage_V> <decision>" per line, running ev_step() on a
 *               pack of identical cells with a 0 A request and a 30 s horizon
 *
 * Exits 2 on a usage error and 3 on malformed or out-of-range input (with a message on stderr).
 */
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spme_v7_apex_solver.h"

#define EV_CLI_MAX_FIELDS 7
#define EV_CLI_LINE_LENGTH 512

/* Parses exactly `count` finite numbers from `text`; returns 0 on success. */
static int ev_parse_numbers(const char *text, double *values, int count) {
    const char *cursor = text;
    int i;
    for (i = 0; i < count; ++i) {
        char *end = NULL;
        errno = 0;
        values[i] = strtod(cursor, &end);
        if ((end == cursor) || (errno != 0) || !isfinite(values[i])) {
            return -1;
        }
        cursor = end;
    }
    while ((*cursor == ' ') || (*cursor == '\t') || (*cursor == '\r') || (*cursor == '\n')) {
        ++cursor;
    }
    return (*cursor == '\0') ? 0 : -1;
}

static int ev_run_ocv(void) {
    char line[EV_CLI_LINE_LENGTH];
    unsigned long line_number = 0UL;
    while (fgets(line, (int)sizeof(line), stdin) != NULL) {
        double soc;
        float ocv = 0.0f;
        float slope = 0.0f;
        ++line_number;
        if ((ev_parse_numbers(line, &soc, 1) != 0) || (soc < 0.0) || (soc > 1.0)) {
            fprintf(stderr, "line %lu: expected one SOC in [0, 1]\n", line_number);
            return 3;
        }
        ev_full_cell_ocv_V_and_dv_dsoc((float)soc, &ocv, &slope);
        printf("%.6f %.6f\n", (double)ocv, (double)slope);
        fflush(stdout);
    }
    return 0;
}

/* Parses the five cell parameters shared by the limit and ekf modes. */
static int ev_parse_config(char *argv[], float coolant_temperature_C, EvCellConfig *config) {
    double parameters[4];
    unsigned int i;
    for (i = 0U; i < 4U; ++i) {
        if (ev_parse_numbers(argv[i], &parameters[i], 1) != 0) {
            fprintf(stderr, "argument %u is not a finite number\n", i + 1U);
            return 2;
        }
    }
    config->nominal_capacity_Ah = (float)parameters[0];
    config->resistance_ohm = (float)parameters[1];
    config->rc_resistance_ohm = (float)parameters[2];
    config->rc_time_constant_s = (float)parameters[3];
    config->heat_transfer_W_m2_K = 10.0f;
    config->heat_capacity_J_kg_K = 1100.0f;
    config->coolant_temperature_C = coolant_temperature_C;
    if (!ev_config_is_valid(config)) {
        fprintf(stderr, "cell parameters out of range\n");
        return 2;
    }
    return 0;
}

static int ev_run_limit(char *argv[]) {
    EvCellConfig config;
    double coolant_C;
    char line[EV_CLI_LINE_LENGTH];
    unsigned long line_number = 0UL;
    int status;

    if (ev_parse_numbers(argv[4], &coolant_C, 1) != 0) {
        fprintf(stderr, "argument 5 is not a finite number\n");
        return 2;
    }
    status = ev_parse_config(argv, (float)coolant_C, &config);
    if (status != 0) {
        return status;
    }
    while (fgets(line, (int)sizeof(line), stdin) != NULL) {
        double f[EV_CLI_MAX_FIELDS];
        EvCellState state;
        EvCellPrediction prediction;
        uint8_t constraints = 0U;
        float approved_A;
        ++line_number;
        if ((ev_parse_numbers(line, f, 7) != 0) || (f[0] < 0.0) || (f[0] > 1.0) || (f[1] < 0.0) || (f[1] > 0.5)
            || (f[4] < EV_MIN_SOH) || (f[4] > 1.0) || (f[5] < 0.0) || (f[6] <= 0.0)) {
            fprintf(stderr, "line %lu: expected '<soc> <soc_sigma> <temperature_C> <rc_voltage_V> <soh> <request_A> <horizon_s>'\n", line_number);
            return 3;
        }
        ev_cell_init(&state, f[0], (float)f[2]);
        state.covariance[0] = f[1] * f[1];
        state.rc_voltage_V = (float)f[3];
        state.soh = f[4];
        approved_A = ev_cell_max_charge_current(&config, &state, (float)f[5], (float)f[6], &constraints, &prediction);
        printf("%.6f %.6f %.6f %.4f %u\n", (double)approved_A, (double)prediction.anode_potential_V,
               (double)prediction.terminal_voltage_V, (double)prediction.predicted_temperature_C, (unsigned)constraints);
        fflush(stdout);
    }
    return 0;
}

static int ev_run_ekf(char *argv[]) {
    EvCellConfig config[EV_NUM_CELLS];
    EvCellState state[EV_NUM_CELLS];
    double initial_soc;
    double initial_temperature_C;
    char line[EV_CLI_LINE_LENGTH];
    unsigned long line_number = 0UL;
    unsigned int i;

    if ((ev_parse_numbers(argv[4], &initial_soc, 1) != 0) || (ev_parse_numbers(argv[5], &initial_temperature_C, 1) != 0)) {
        fprintf(stderr, "initial SOC and temperature must be finite numbers\n");
        return 2;
    }
    for (i = 0U; i < EV_NUM_CELLS; ++i) {
        int status = ev_parse_config(argv, (float)initial_temperature_C, &config[i]);
        if (status != 0) {
            return status;
        }
        ev_cell_init(&state[i], initial_soc, (float)initial_temperature_C);
    }

    while (fgets(line, (int)sizeof(line), stdin) != NULL) {
        double fields[EV_CLI_MAX_FIELDS];
        float voltages_V[EV_NUM_CELLS];
        float temperatures_C[EV_NUM_CELLS];
        EvPackResult result;
        ++line_number;
        if (ev_parse_numbers(line, fields, 4) != 0) {
            fprintf(stderr, "line %lu: expected '<dt_s> <current_A> <voltage_V> <temperature_C>'\n", line_number);
            return 3;
        }
        for (i = 0U; i < EV_NUM_CELLS; ++i) {
            voltages_V[i] = (float)fields[2];
            temperatures_C[i] = (float)fields[3];
        }
        result = ev_step(config, state, (float)fields[1], 0.0f, voltages_V, temperatures_C, (float)fields[0], 30.0f);
        printf("%.6f %.6f %.6f %u\n", state[0].soc, sqrt(state[0].covariance[0]), (double)state[0].rc_voltage_V,
               (unsigned)result.decision);
        fflush(stdout);
    }
    return 0;
}

int main(int argc, char *argv[]) {
    if ((argc == 2) && (strcmp(argv[1], "ocv") == 0)) {
        return ev_run_ocv();
    }
    if ((argc == 7) && (strcmp(argv[1], "limit") == 0)) {
        return ev_run_limit(&argv[2]);
    }
    if ((argc == 8) && (strcmp(argv[1], "ekf") == 0)) {
        return ev_run_ekf(&argv[2]);
    }
    fprintf(stderr, "usage: ev_cli ocv\n"
                    "       ev_cli limit <capacity_Ah> <r0_ohm> <r1_ohm> <tau_s> <coolant_temperature_C>\n"
                    "       ev_cli ekf <capacity_Ah> <r0_ohm> <r1_ohm> <tau_s> <initial_soc> <initial_temperature_C>\n");
    return 2;
}
