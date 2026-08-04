#ifndef SPME_V7_APEX_SOLVER_INTERNAL_H
#define SPME_V7_APEX_SOLVER_INTERNAL_H

#include "spme_v7_apex_solver.h"

#ifdef EV_TEST_MODE
#define STATIC
#else
#define STATIC static
#endif

STATIC float ev_clamp(float value, float lower, float upper);
STATIC double ev_clamp_d(double value, double lower, double upper);
STATIC bool ev_is_valid_cell(const EvCellState *cell);
STATIC float ev_graphite_ocp_V(float theta);
STATIC float ev_nmc_ocp_V(float soc);
STATIC void ev_full_cell_ocv_V_and_dv_dsoc(float soc, float *ocv, float *dv_dsoc);
STATIC float ev_diffusivity_m2_s(float temperature_K);
STATIC float ev_predict_temperature_C(const EvCellState *cell, float charge_current_A, float dt_s);
STATIC float fast_asinhf(float x);
STATIC float fast_log1pf(float x);
STATIC bool ev_predict_cell(const EvCellState *cell, float charge_current_A, float horizon_s, EvCellPrediction *prediction);
STATIC bool ev_current_is_safe(const EvCellState cells[EV_NUM_CELLS], float current_A, float horizon_s, uint8_t *limiting_index, EvCellPrediction *limiting_prediction);
STATIC float ev_find_safe_charge_current(const EvCellState cells[EV_NUM_CELLS], float requested_current_A, float horizon_s, uint8_t *limiting_index, EvCellPrediction *limiting_prediction);
STATIC void ev_ekf_and_sei_update(EvCellState *cell, float measured_current_A, float measured_voltage_V, float measured_temperature_C, float dt_s);
STATIC const char *ev_decision_name(EvDecision decision);

#endif // SPME_V7_APEX_SOLVER_INTERNAL_H
