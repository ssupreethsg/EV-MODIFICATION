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
 *       stdout: "<soc> <soc_sigma> <rc_voltage_V> <decision> <fault_flags>" per line, running
 *               ev_step() on a pack of identical cells with a 0 A request and a 30 s horizon; the
 *               measured current is not supervised (there is no charger to follow)
 *
 *   ev_cli charge <capacity_Ah> <r0_ohm> <r1_ohm> <tau_s> <initial_soc> <initial_temperature_C>
 *       stdin : "<dt_s> <current_A> <voltage_V> <temperature_C> <request_A> <horizon_s>" per line
 *       stdout: "<approved_A> <decision> <fault_flags> <soc> <soc_sigma>" per line: the complete
 *               controller in a closed loop (ev_step() with the state carried between lines: EKF,
 *               gradient history, voltage tracking, over-current check). current_A is the current
 *               that flowed over the dt_s that just ended, voltage_V the reading at its end.
 *
 * Exits 2 on a usage error (including out-of-range arguments), 3 on malformed, out-of-range or
 * overlong (> EV_CLI_LINE_LENGTH - 2 characters) input lines or lines containing a NUL byte, and
 * 4 when writing to stdout fails (with a message on stderr). dt_s must be in
 * [EV_MIN_CONTROL_DT_S, EV_MAX_CONTROL_DT_S].
 */
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ev_modify.h"

#define EV_CLI_MAX_FIELDS 7
#define EV_CLI_LINE_LENGTH 512
/* Bounds that keep every double -> float cast defined (C11 6.3.1.5); the library judges plausibility. */
#define EV_CLI_MAX_PARAMETER 1.0e6
#define EV_CLI_MAX_VOLTAGE_V 1.0e3
#define EV_CLI_MAX_TEMPERATURE_C 1.0e3

static int ev_within(double value, double limit) {
    return (value >= -limit) && (value <= limit);
}

/* dt in [EV_MIN_CONTROL_DT_S, EV_MAX_CONTROL_DT_S], judged in float exactly as the library does. */
static int ev_dt_in_range(double dt) {
    return (dt > 0.0) && (dt <= EV_CLI_MAX_PARAMETER) && ((float)dt >= EV_MIN_CONTROL_DT_S)
           && ((float)dt <= EV_MAX_CONTROL_DT_S);
}

/* Flushes stdout; returns 4 (after a message) when any output was lost, else 0. */
static int ev_flush_status(void) {
    if ((fflush(stdout) != 0) || ferror(stdout)) {
        fprintf(stderr, "error writing to stdout\n");
        return 4;
    }
    return 0;
}

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

/* Reads one line (without the newline): 1 = line read, 0 = end of input, -1 = rejected (reported):
 * longer than the buffer or containing a NUL byte. */
static int ev_read_line(char *line, int size, unsigned long line_number) {
    int length = 0;
    int c = getc(stdin);
    if (c == EOF) {
        return 0;
    }
    while ((c != EOF) && (c != '\n')) {
        if (c == '\0') {
            fprintf(stderr, "line %lu: contains a NUL byte\n", line_number);
            return -1;
        }
        if (length >= (size - 2)) {
            fprintf(stderr, "line %lu: longer than %d characters\n", line_number, size - 2);
            return -1;
        }
        line[length++] = (char)c;
        c = getc(stdin);
    }
    line[length] = '\0';
    return 1;
}

static int ev_run_ocv(void) {
    char line[EV_CLI_LINE_LENGTH];
    unsigned long line_number = 0UL;
    int status;
    while ((status = ev_read_line(line, (int)sizeof(line), line_number + 1UL)) > 0) {
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
        if (ev_flush_status() != 0) {
            return 4;
        }
    }
    return (status < 0) ? 3 : 0;
}

/* Parses the five cell parameters shared by the limit and ekf modes. */
static int ev_parse_config(char *argv[], float coolant_temperature_C, EvCellConfig *config) {
    double parameters[4];
    unsigned int i;
    for (i = 0U; i < 4U; ++i) {
        if ((ev_parse_numbers(argv[i], &parameters[i], 1) != 0) || !ev_within(parameters[i], EV_CLI_MAX_PARAMETER)) {
            fprintf(stderr, "argument %u is not a finite number of plausible size\n", i + 1U);
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

    if ((ev_parse_numbers(argv[4], &coolant_C, 1) != 0) || !ev_within(coolant_C, EV_CLI_MAX_TEMPERATURE_C)) {
        fprintf(stderr, "argument 5 is not a finite temperature\n");
        return 2;
    }
    status = ev_parse_config(argv, (float)coolant_C, &config);
    if (status != 0) {
        return status;
    }
    while ((status = ev_read_line(line, (int)sizeof(line), line_number + 1UL)) > 0) {
        double f[EV_CLI_MAX_FIELDS];
        EvCellState state;
        EvCellPrediction prediction;
        uint8_t constraints = 0U;
        float approved_A;
        ++line_number;
        if ((ev_parse_numbers(line, f, 7) != 0) || (f[0] < 0.0) || (f[0] > 1.0) || (f[1] < 0.0) || (f[1] > 0.5)
            || (f[2] < (double)EV_MIN_PLAUSIBLE_TEMPERATURE_C) || (f[2] > (double)EV_MAX_PLAUSIBLE_TEMPERATURE_C)
            || (f[3] < -1.0) || (f[3] > 1.0) || (f[4] < EV_MIN_SOH) || (f[4] > 1.0)
            || (f[5] < 0.0) || (f[5] > (double)EV_MAX_PLAUSIBLE_CURRENT_A)
            || (f[6] <= 0.0) || (f[6] > (double)EV_MAX_PREDICTION_HORIZON_S)) {
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
        if (ev_flush_status() != 0) {
            return 4;
        }
    }
    return (status < 0) ? 3 : 0;
}

/* Replays measurements through ev_step() on a pack of identical cells. closed_loop: each line
 * also carries the request and horizon, and the charger is supervised against the approvals. */
static int ev_run_pack(char *argv[], int closed_loop) {
    EvCellConfig config[EV_NUM_CELLS];
    EvCellState state[EV_NUM_CELLS];
    double initial_soc;
    double initial_temperature_C;
    char line[EV_CLI_LINE_LENGTH];
    unsigned long line_number = 0UL;
    int field_count = closed_loop ? 6 : 4;
    unsigned int i;
    int status;

    if ((ev_parse_numbers(argv[4], &initial_soc, 1) != 0) || (ev_parse_numbers(argv[5], &initial_temperature_C, 1) != 0)) {
        fprintf(stderr, "initial SOC and temperature must be finite numbers\n");
        return 2;
    }
    if ((initial_soc < 0.0) || (initial_soc > 1.0) || (initial_temperature_C < (double)EV_MIN_PLAUSIBLE_TEMPERATURE_C)
        || (initial_temperature_C > (double)EV_MAX_PLAUSIBLE_TEMPERATURE_C)) {
        fprintf(stderr, "initial SOC must be in [0, 1] and the initial temperature plausible\n");
        return 2;
    }
    for (i = 0U; i < EV_NUM_CELLS; ++i) {
        status = ev_parse_config(argv, (float)initial_temperature_C, &config[i]);
        if (status != 0) {
            return status;
        }
        ev_cell_init(&state[i], initial_soc, (float)initial_temperature_C);
    }

    while ((status = ev_read_line(line, (int)sizeof(line), line_number + 1UL)) > 0) {
        double f[EV_CLI_MAX_FIELDS];
        float voltages_V[EV_NUM_CELLS];
        float temperatures_C[EV_NUM_CELLS];
        float request_A = 0.0f;
        float horizon_s = 30.0f;
        EvPackResult result;
        ++line_number;
        if ((ev_parse_numbers(line, f, field_count) != 0) || !ev_dt_in_range(f[0])
            || (fabs(f[1]) > (double)EV_MAX_PLAUSIBLE_CURRENT_A) || !ev_within(f[2], EV_CLI_MAX_VOLTAGE_V)
            || !ev_within(f[3], EV_CLI_MAX_TEMPERATURE_C)
            || (closed_loop && ((f[4] < 0.0) || (f[4] > (double)EV_MAX_PLAUSIBLE_CURRENT_A)
                                || (f[5] <= 0.0) || (f[5] > (double)EV_MAX_PREDICTION_HORIZON_S)))) {
            fprintf(stderr, closed_loop
                ? "line %lu: expected '<dt_s> <current_A> <voltage_V> <temperature_C> <request_A> <horizon_s>' in range\n"
                : "line %lu: expected '<dt_s> <current_A> <voltage_V> <temperature_C>' in range\n", line_number);
            return 3;
        }
        for (i = 0U; i < EV_NUM_CELLS; ++i) {
            voltages_V[i] = (float)f[2];
            temperatures_C[i] = (float)f[3];
            if (!closed_loop) {
                /* the estimator replay has no charger to supervise */
                state[i].commanded_current_A = EV_MAX_PLAUSIBLE_CURRENT_A;
                state[i].ramp_from_A = EV_MAX_PLAUSIBLE_CURRENT_A;
            }
        }
        if (closed_loop) {
            request_A = (float)f[4];
            horizon_s = (float)f[5];
        }
        result = ev_step(config, state, (float)f[1], request_A, voltages_V, temperatures_C, (float)f[0], horizon_s);
        if (closed_loop) {
            printf("%.6f %u %u %.6f %.6f\n", (double)result.approved_charge_current_A, (unsigned)result.decision,
                   (unsigned)result.fault_flags, state[0].soc, sqrt(state[0].covariance[0]));
        } else {
            printf("%.6f %.6f %.6f %u %u\n", state[0].soc, sqrt(state[0].covariance[0]), (double)state[0].rc_voltage_V,
                   (unsigned)result.decision, (unsigned)result.fault_flags);
        }
        if (ev_flush_status() != 0) {
            return 4;
        }
    }
    return (status < 0) ? 3 : 0;
}

int main(int argc, char *argv[]) {
    if ((argc == 2) && (strcmp(argv[1], "ocv") == 0)) {
        return ev_run_ocv();
    }
    if ((argc == 7) && (strcmp(argv[1], "limit") == 0)) {
        return ev_run_limit(&argv[2]);
    }
    if ((argc == 8) && (strcmp(argv[1], "ekf") == 0)) {
        return ev_run_pack(&argv[2], 0);
    }
    if ((argc == 8) && (strcmp(argv[1], "charge") == 0)) {
        return ev_run_pack(&argv[2], 1);
    }
    fprintf(stderr, "usage: ev_cli ocv\n"
                    "       ev_cli limit <capacity_Ah> <r0_ohm> <r1_ohm> <tau_s> <coolant_temperature_C>\n"
                    "       ev_cli ekf <capacity_Ah> <r0_ohm> <r1_ohm> <tau_s> <initial_soc> <initial_temperature_C>\n"
                    "       ev_cli charge <capacity_Ah> <r0_ohm> <r1_ohm> <tau_s> <initial_soc> <initial_temperature_C>\n");
    return 2;
}
