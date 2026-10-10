/*
 * EV MODIFY - conservative charge-current limiter.
 * Copyright (c) 2026 Supreeth S G. All rights reserved. See LICENSE.
 */

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "spme_v7_apex_solver.h"

#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && (__FINITE_MATH_ONLY__ != 0)) || defined(_M_FP_FAST)
#error "Do not build with -ffast-math, -ffinite-math-only or /fp:fast: the NaN/inf checks are safety relevant."
#endif

/* The bit-pattern finiteness checks below assume IEEE-754 binary32/binary64. */
_Static_assert(sizeof(float) == 4U, "float must be IEEE-754 binary32");
_Static_assert(sizeof(double) == 8U, "double must be IEEE-754 binary64");

#define EV_COARSE_SEARCH_STEPS 20U
#define EV_BISECTION_STEPS 24U
#define EV_SOC_PROCESS_NOISE_PER_S 1.0e-7     /* SOC^2 / s */
#define EV_RC_PROCESS_NOISE_V2_PER_S 1.0e-6   /* V^2 / s */
#define EV_VOLTAGE_MEASUREMENT_NOISE_V2 4.0e-4 /* (20 mV)^2: sensor noise and OCV/hysteresis error */
#define EV_RESISTANCE_UNCERTAINTY 0.3         /* relative R0+R1 error, adds (0.3*(R0+R1)*I)^2 to R */
#define EV_MAX_SOC_VARIANCE 0.25
#define EV_UNKNOWN_DT_SOC_VARIANCE 0.01
#define EV_INITIAL_SOC_VARIANCE 0.0025
#define EV_INITIAL_RC_VARIANCE_V2 1.0e-4
#define EV_MAX_RC_VOLTAGE_V 1.0
/* Electrode-averaged depletion at which the steady parabolic profile falls to 10 % of c_e0 at
 * the current collector (c_min = c_e0 * (1 - 1.5 * depletion)). */
#define EV_ELECTROLYTE_DEPLETION_LIMIT 0.6f

/* Electrode OCPs (PyBaMM Chen2020) tabulated against cell SOC on a uniform grid using the
 * stoichiometry windows in the header. Do not edit by hand. */
/* BEGIN GENERATED OCV TABLES (tools/generate_ocv_lut.py) */
#define EV_LUT_POINTS 201U
static const float LUT_GRAPHITE_OCP[EV_LUT_POINTS] = {
    1.105430f, 0.993178f, 0.898779f, 0.819357f, 0.752490f, 0.696136f, 0.648566f, 0.608315f,
    0.574136f, 0.544959f, 0.519863f, 0.498043f, 0.478795f, 0.461492f, 0.445577f, 0.430557f,
    0.416010f, 0.401592f, 0.387060f, 0.372287f, 0.357281f, 0.342179f, 0.327229f, 0.312746f,
    0.299059f, 0.286448f, 0.275113f, 0.265154f, 0.256572f, 0.249292f, 0.243189f, 0.238108f,
    0.233890f, 0.230380f, 0.227434f, 0.224930f, 0.222760f, 0.220834f, 0.219078f, 0.217430f,
    0.215837f, 0.214257f, 0.212655f, 0.211001f, 0.209273f, 0.207449f, 0.205517f, 0.203465f,
    0.201286f, 0.198979f, 0.196543f, 0.193986f, 0.191316f, 0.188548f, 0.185698f, 0.182787f,
    0.179838f, 0.176875f, 0.173922f, 0.171004f, 0.168147f, 0.165370f, 0.162694f, 0.160135f,
    0.157706f, 0.155417f, 0.153274f, 0.151281f, 0.149437f, 0.147742f, 0.146190f, 0.144775f,
    0.143492f, 0.142333f, 0.141288f, 0.140350f, 0.139510f, 0.138759f, 0.138090f, 0.137495f,
    0.136966f, 0.136497f, 0.136082f, 0.135715f, 0.135391f, 0.135105f, 0.134852f, 0.134630f,
    0.134434f, 0.134262f, 0.134110f, 0.133976f, 0.133859f, 0.133756f, 0.133665f, 0.133585f,
    0.133515f, 0.133453f, 0.133398f, 0.133350f, 0.133307f, 0.133268f, 0.133234f, 0.133202f,
    0.133172f, 0.133144f, 0.133117f, 0.133089f, 0.133060f, 0.133028f, 0.132992f, 0.132949f,
    0.132897f, 0.132833f, 0.132753f, 0.132652f, 0.132522f, 0.132355f, 0.132141f, 0.131866f,
    0.131512f, 0.131060f, 0.130484f, 0.129757f, 0.128846f, 0.127719f, 0.126344f, 0.124697f,
    0.122767f, 0.120560f, 0.118109f, 0.115473f, 0.112733f, 0.109986f, 0.107329f, 0.104846f,
    0.102600f, 0.100627f, 0.098938f, 0.097525f, 0.096363f, 0.095423f, 0.094672f, 0.094078f,
    0.093611f, 0.093247f, 0.092964f, 0.092745f, 0.092576f, 0.092447f, 0.092347f, 0.092270f,
    0.092211f, 0.092166f, 0.092132f, 0.092106f, 0.092085f, 0.092070f, 0.092058f, 0.092049f,
    0.092042f, 0.092037f, 0.092033f, 0.092030f, 0.092028f, 0.092026f, 0.092024f, 0.092023f,
    0.092023f, 0.092022f, 0.092022f, 0.092021f, 0.092021f, 0.092021f, 0.092021f, 0.092020f,
    0.092020f, 0.092020f, 0.092020f, 0.092020f, 0.092020f, 0.092020f, 0.092020f, 0.092020f,
    0.092020f, 0.092020f, 0.092020f, 0.092020f, 0.092020f, 0.092020f, 0.092020f, 0.092020f,
    0.092020f, 0.092020f, 0.092020f, 0.092020f, 0.092020f, 0.092020f, 0.092020f, 0.092020f,
    0.092020f
};
static const float LUT_NMC_OCP[EV_LUT_POINTS] = {
    3.605436f, 3.607823f, 3.610210f, 3.612597f, 3.614985f, 3.617372f, 3.619759f, 3.622147f,
    3.624534f, 3.626922f, 3.629309f, 3.631697f, 3.634084f, 3.636472f, 3.638860f, 3.641248f,
    3.643636f, 3.646024f, 3.648412f, 3.650800f, 3.653188f, 3.655577f, 3.657965f, 3.660354f,
    3.662743f, 3.665133f, 3.667522f, 3.669912f, 3.672302f, 3.674692f, 3.677083f, 3.679475f,
    3.681866f, 3.684259f, 3.686651f, 3.689045f, 3.691439f, 3.693834f, 3.696230f, 3.698627f,
    3.701025f, 3.703425f, 3.705826f, 3.708228f, 3.710632f, 3.713038f, 3.715446f, 3.717857f,
    3.720270f, 3.722687f, 3.725106f, 3.727530f, 3.729957f, 3.732389f, 3.734826f, 3.737269f,
    3.739718f, 3.742174f, 3.744638f, 3.747110f, 3.749592f, 3.752085f, 3.754590f, 3.757108f,
    3.759640f, 3.762188f, 3.764754f, 3.767340f, 3.769948f, 3.772579f, 3.775237f, 3.777924f,
    3.780644f, 3.783398f, 3.786191f, 3.789026f, 3.791906f, 3.794837f, 3.797822f, 3.800866f,
    3.803973f, 3.807147f, 3.810395f, 3.813720f, 3.817127f, 3.820620f, 3.824204f, 3.827882f,
    3.831657f, 3.835532f, 3.839507f, 3.843583f, 3.847759f, 3.852033f, 3.856401f, 3.860857f,
    3.865395f, 3.870007f, 3.874683f, 3.879411f, 3.884180f, 3.888978f, 3.893791f, 3.898605f,
    3.903408f, 3.908186f, 3.912928f, 3.917624f, 3.922262f, 3.926836f, 3.931339f, 3.935766f,
    3.940113f, 3.944379f, 3.948563f, 3.952667f, 3.956693f, 3.960643f, 3.964522f, 3.968336f,
    3.972089f, 3.975787f, 3.979437f, 3.983045f, 3.986618f, 3.990163f, 3.993686f, 3.997194f,
    4.000694f, 4.004192f, 4.007695f, 4.011209f, 4.014740f, 4.018294f, 4.021877f, 4.025494f,
    4.029151f, 4.032852f, 4.036603f, 4.040407f, 4.044270f, 4.048194f, 4.052184f, 4.056242f,
    4.060370f, 4.064570f, 4.068843f, 4.073189f, 4.077607f, 4.082096f, 4.086651f, 4.091270f,
    4.095945f, 4.100670f, 4.105435f, 4.110230f, 4.115041f, 4.119854f, 4.124651f, 4.129414f,
    4.134123f, 4.138754f, 4.143284f, 4.147688f, 4.151940f, 4.156016f, 4.159889f, 4.163536f,
    4.166937f, 4.170073f, 4.172933f, 4.175510f, 4.177801f, 4.179817f, 4.181571f, 4.183090f,
    4.184408f, 4.185569f, 4.186626f, 4.187639f, 4.188676f, 4.189809f, 4.191112f, 4.192659f,
    4.194523f, 4.196773f, 4.199468f, 4.202661f, 4.206392f, 4.210691f, 4.215573f, 4.221043f,
    4.227091f, 4.233697f, 4.240830f, 4.248449f, 4.256507f, 4.264950f, 4.273723f, 4.282767f,
    4.292021f
};
/* END GENERATED OCV TABLES */

/* ---- Numeric helpers --------------------------------------------------------------------- */

/* Bit-pattern test: independent of compiler floating-point optimisation settings. */
static bool ev_is_finite(float value) {
    uint32_t bits;
    memcpy(&bits, &value, sizeof(bits));
    return (bits & 0x7F800000UL) != 0x7F800000UL;
}

static bool ev_is_finite_d(double value) {
    uint64_t bits;
    memcpy(&bits, &value, sizeof(bits));
    return (bits & 0x7FF0000000000000ULL) != 0x7FF0000000000000ULL;
}

static bool ev_in_range(float value, float lower, float upper) {
    return ev_is_finite(value) && (value >= lower) && (value <= upper);
}

static bool ev_in_range_d(double value, double lower, double upper) {
    return ev_is_finite_d(value) && (value >= lower) && (value <= upper);
}

/* Callers guarantee a finite value. */
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

/* exp(-Ea/R * (1/T - 1/Tref)) with a bounded exponent */
static float ev_arrhenius(float activation_energy_J_mol, float temperature_K) {
    float exponent = (-activation_energy_J_mol / EV_GAS_J_PER_MOL_K) * ((1.0f / temperature_K) - (1.0f / EV_REFERENCE_TEMPERATURE_K));
    return expf(ev_clamp(exponent, -20.0f, 20.0f));
}

/* soc must be finite */
static float ev_lut_interp(const float lut[EV_LUT_POINTS], float soc, float *slope_per_soc) {
    const float segments = (float)(EV_LUT_POINTS - 1U);
    float f_idx = ev_clamp(soc, 0.0f, 1.0f) * segments;
    uint32_t idx = (uint32_t)f_idx;
    if (idx >= (EV_LUT_POINTS - 1U)) { idx = EV_LUT_POINTS - 2U; }
    float t = f_idx - (float)idx;
    float delta = lut[idx + 1U] - lut[idx];
    if (slope_per_soc != NULL) { *slope_per_soc = delta * segments; }
    return lut[idx] + (t * delta);
}

void ev_full_cell_ocv_V_and_dv_dsoc(float soc, float *ocv, float *dv_dsoc) {
    float graphite_slope = 0.0f;
    float nmc_slope = 0.0f;
    float graphite_val;
    float nmc_val;

    if (!ev_is_finite(soc)) {
        float not_a_number = soc - soc;   /* NaN for both NaN and +/-inf input */
        if (ocv != NULL) { *ocv = not_a_number; }
        if (dv_dsoc != NULL) { *dv_dsoc = not_a_number; }
        return;
    }
    graphite_val = ev_lut_interp(LUT_GRAPHITE_OCP, soc, &graphite_slope);
    nmc_val = ev_lut_interp(LUT_NMC_OCP, soc, &nmc_slope);
    if (ocv != NULL) { *ocv = nmc_val - graphite_val; }
    if (dv_dsoc != NULL) { *dv_dsoc = nmc_slope - graphite_slope; }
}

/* ---- Validation -------------------------------------------------------------------------- */

bool ev_config_is_valid(const EvCellConfig *config) {
    return (config != NULL)
        && ev_in_range(config->nominal_capacity_Ah, 0.1f, 1000.0f)
        && ev_in_range(config->resistance_ohm, 0.0f, 1.0f)
        && ev_in_range(config->rc_resistance_ohm, 0.0f, 1.0f)
        && ev_in_range(config->rc_time_constant_s, 1.0f, 1.0e4f)
        && ev_in_range(config->heat_transfer_W_m2_K, 1.0f, 1000.0f)
        && ev_in_range(config->heat_capacity_J_kg_K, 100.0f, 5000.0f)
        && ev_in_range(config->coolant_temperature_C, EV_MIN_PLAUSIBLE_TEMPERATURE_C, EV_MAX_PLAUSIBLE_TEMPERATURE_C);
}

static bool ev_state_is_valid(const EvCellState *state) {
    return (state != NULL)
        && ev_in_range_d(state->soc, 0.0, 1.0)
        && ev_in_range_d(state->covariance[0], 0.0, EV_MAX_SOC_VARIANCE)
        && ev_is_finite_d(state->covariance[1])
        && ev_in_range_d(state->covariance[2], 0.0, 1.0)
        && ev_in_range_d(state->soh, EV_MIN_SOH, 1.0)
        && ev_in_range_d(state->throughput_Ah, 0.0, DBL_MAX)
        && ev_in_range(state->rc_voltage_V, (float)-EV_MAX_RC_VOLTAGE_V, (float)EV_MAX_RC_VOLTAGE_V)
        && ev_in_range(state->temperature_C, EV_MIN_PLAUSIBLE_TEMPERATURE_C, EV_MAX_PLAUSIBLE_TEMPERATURE_C)
        && ev_in_range(state->resistance_growth_ohm, 0.0f, 1.0f)
        && ev_in_range(state->last_voltage_V, 0.0f, EV_MAX_PLAUSIBLE_VOLTAGE_V)
        && ev_in_range(state->last_current_A, -EV_MAX_PLAUSIBLE_CURRENT_A, EV_MAX_PLAUSIBLE_CURRENT_A)
        && ev_in_range(state->track_start_V, 0.0f, EV_MAX_PLAUSIBLE_VOLTAGE_V)
        && ev_is_finite(state->track_start_model_V)
        && ev_in_range(state->track_start_ocv_V, 0.0f, EV_MAX_PLAUSIBLE_VOLTAGE_V)
        && ev_in_range_d(state->track_start_soc, 0.0, 1.0)
        && ev_in_range(state->track_start_current_A, -EV_MAX_PLAUSIBLE_CURRENT_A, EV_MAX_PLAUSIBLE_CURRENT_A)
        && ev_in_range(state->track_max_current_step_A, 0.0f, 2.0f * EV_MAX_PLAUSIBLE_CURRENT_A)
        && ev_in_range(state->track_rc_voltage_V, (float)-EV_MAX_RC_VOLTAGE_V, (float)EV_MAX_RC_VOLTAGE_V)
        && ev_in_range_d(state->track_soc, 0.0, 1.0)
        && ev_in_range(state->commanded_current_A, 0.0f, EV_MAX_PLAUSIBLE_CURRENT_A)
        && ev_in_range(state->ramp_from_A, 0.0f, EV_MAX_PLAUSIBLE_CURRENT_A)
        && ev_in_range(state->ramp_age_s, 0.0f, FLT_MAX)
        && ev_in_range(state->over_current_time_s, 0.0f, FLT_MAX);
}

void ev_cell_init(EvCellState *state, double initial_soc, float initial_temperature_C) {
    if (state == NULL) {
        return;
    }
    memset(state, 0, sizeof(*state));
    state->soc = ev_is_finite_d(initial_soc) ? ev_clamp_d(initial_soc, 0.0, 1.0) : 0.5;
    state->covariance[0] = ev_is_finite_d(initial_soc) ? EV_INITIAL_SOC_VARIANCE : EV_MAX_SOC_VARIANCE;
    state->covariance[1] = 0.0;
    state->covariance[2] = EV_INITIAL_RC_VARIANCE_V2;
    state->soh = 1.0;
    state->throughput_Ah = 0.0;
    state->rc_voltage_V = 0.0f;
    state->temperature_C = ev_is_finite(initial_temperature_C) ? initial_temperature_C : EV_MAX_PLAUSIBLE_TEMPERATURE_C;
    state->resistance_growth_ohm = 0.0f;
    state->last_voltage_V = 0.0f;
    state->last_current_A = 0.0f;
    state->innovation_fail_cycles = 0U;
    state->track_soc = state->soc;
    state->fault_hold_cycles = 0U;
    state->fault_flags = 0U;
    state->cold_charge_block = false;
}

/* Recovers a corrupted estimator: SOC unknown (sigma 0.5), polarisation zero. */
static void ev_reset_estimator(EvCellState *state) {
    if (!ev_in_range_d(state->soc, 0.0, 1.0)) { state->soc = 0.5; }
    state->covariance[0] = EV_MAX_SOC_VARIANCE;
    state->covariance[1] = 0.0;
    state->covariance[2] = EV_INITIAL_RC_VARIANCE_V2;
    if (!ev_in_range_d(state->soh, EV_MIN_SOH, 1.0)) { state->soh = EV_MIN_SOH; }
    if (!ev_in_range_d(state->throughput_Ah, 0.0, DBL_MAX)) { state->throughput_Ah = 0.0; }
    state->rc_voltage_V = 0.0f;
    if (!ev_in_range(state->temperature_C, EV_MIN_PLAUSIBLE_TEMPERATURE_C, EV_MAX_PLAUSIBLE_TEMPERATURE_C)) {
        state->temperature_C = EV_MAX_PLAUSIBLE_TEMPERATURE_C;   /* forces a temperature fault */
    }
    if (!ev_in_range(state->resistance_growth_ohm, 0.0f, 1.0f)) { state->resistance_growth_ohm = 1.0f; }
    state->last_voltage_V = 0.0f;
    state->last_current_A = 0.0f;
    state->track_start_V = 0.0f;   /* restart the voltage-tracking window */
    state->track_start_model_V = 0.0f;
    state->track_start_ocv_V = 0.0f;
    state->track_rc_voltage_V = 0.0f;
    state->track_soc = state->soc;
    state->track_start_soc = state->soc;
    state->track_start_current_A = 0.0f;
    state->track_max_current_step_A = 0.0f;
    state->track_fail_cycles = 0U;
    state->track_suspect = false;
    /* Unknown command history: allow nothing beyond the tolerance from now on. */
    state->commanded_current_A = 0.0f;
    state->ramp_from_A = 0.0f;
    state->ramp_age_s = 0.0f;
    if (!ev_in_range(state->over_current_time_s, 0.0f, FLT_MAX)) { state->over_current_time_s = 0.0f; }
}

/* ---- Models ------------------------------------------------------------------------------ */

/* Resistances rise as the cell cools: R(T) = R(25 degC) * exp(Ea/R * (1/T - 1/Tref)). */
static float ev_resistance_scale(float temperature_C) {
    float temperature_K = ev_clamp(temperature_C, EV_MIN_PLAUSIBLE_TEMPERATURE_C, EV_MAX_PLAUSIBLE_TEMPERATURE_C) + EV_KELVIN_OFFSET;
    return 1.0f / ev_arrhenius(EV_RESISTANCE_ACTIVATION_ENERGY_J_MOL, temperature_K);
}

static float ev_effective_r0(const EvCellConfig *config, const EvCellState *state, float temperature_C) {
    return (config->resistance_ohm + state->resistance_growth_ohm) * ev_resistance_scale(temperature_C);
}

static float ev_effective_r1(const EvCellConfig *config, float temperature_C) {
    return config->rc_resistance_ohm * ev_resistance_scale(temperature_C);
}

/* Mean of ln(c/c_e0) over the anode for the steady parabolic profile
 * c/c_e0 = 1 - a*(1 - u^2), a = 1.5 * depletion, u in [0, 1] (closed form). */
static float ev_mean_log_concentration(float depletion) {
    double a = 1.5 * (double)depletion;
    double b;
    if (a < 1.0e-4) {
        return -depletion;   /* first-order expansion */
    }
    b = sqrt((1.0 - a) / a);
    return (float)(log(a) + log(1.0 + (b * b)) - 2.0 + (2.0 * b * atan(1.0 / b)));
}

/* Lumped thermal model, exact solution for constant heat input over dt_s. */
static float ev_predict_temperature_C(const EvCellConfig *config, float start_temperature_C, float heat_W, float dt_s) {
    float thermal_conductance_W_K = config->heat_transfer_W_m2_K * EV_CELL_COOLING_AREA_M2;
    float heat_capacity_J_K = EV_CELL_MASS_KG * config->heat_capacity_J_kg_K;
    float steady_C = config->coolant_temperature_C + (heat_W / thermal_conductance_W_K);
    /* T0 + (T_inf - T0) * (1 - exp(-t/tau)): exactly T0 at zero heat and equilibrium, and monotone in heat_W */
    float approach = -expm1f(ev_clamp(-(thermal_conductance_W_K * dt_s) / heat_capacity_J_K, -50.0f, 0.0f));
    return start_temperature_C + ((steady_C - start_temperature_C) * approach);
}

bool ev_predict_cell(const EvCellConfig *config, const EvCellState *state, float charge_current_A,
                     float horizon_s, EvCellPrediction *prediction) {
    const float anode_stoich_window = EV_ANODE_STOICH_AT_SOC100 - EV_ANODE_STOICH_AT_SOC0;
    const float r_over_f = EV_GAS_J_PER_MOL_K / EV_FARADAY_C_PER_MOL;

    if (prediction == NULL) {
        return false;
    }
    prediction->anode_potential_V = 0.0f;
    prediction->terminal_voltage_V = 0.0f;
    prediction->predicted_temperature_C = 0.0f;
    prediction->predicted_soc = 0.0f;
    prediction->violated = EV_CONSTRAINT_INVALID_INPUT;
    prediction->safe = false;

    if (!ev_config_is_valid(config) || !ev_state_is_valid(state)
        || !ev_in_range(charge_current_A, 0.0f, EV_MAX_PLAUSIBLE_CURRENT_A)
        || !ev_in_range(horizon_s, FLT_MIN, EV_MAX_PREDICTION_HORIZON_S)) {
        return false;
    }

    /* Thermal prediction first. The over-temperature check uses the prediction with Joule heat.
     * Resistances, electrochemistry and the 0 degC check use the colder of the present
     * temperature and the prediction WITHOUT Joule heat: taking credit for self-heating would
     * let a larger current pass a check that a smaller one fails, so the safe set would not be
     * an interval [0, I_max] and the approved current would not rise monotonically with the
     * request. Without it every constraint below is monotone in the current. */
    float r0_now_ohm = ev_effective_r0(config, state, state->temperature_C);
    float r1_now_ohm = ev_effective_r1(config, state->temperature_C);
    prediction->predicted_temperature_C = ev_predict_temperature_C(config, state->temperature_C,
        charge_current_A * charge_current_A * (r0_now_ohm + r1_now_ohm), horizon_s);
    float coldest_C = fminf(state->temperature_C, ev_predict_temperature_C(config, state->temperature_C, 0.0f, horizon_s));
    float r0_ohm = ev_effective_r0(config, state, coldest_C);
    float r1_ohm = ev_effective_r1(config, coldest_C);

    /* SOC at the end of the horizon, as a 3-sigma upper bound. */
    double capacity_Ah = (double)config->nominal_capacity_Ah * state->soh;
    double soc_upper = state->soc + (EV_SOC_SIGMA_MULTIPLIER * sqrt(state->covariance[0]))
        + (((double)charge_current_A * (double)horizon_s) / (capacity_Ah * 3600.0));
    prediction->predicted_soc = (float)soc_upper;

    /* Terminal voltage bound over the whole horizon: OCV rises monotonically with SOC and the RC
     * voltage moves exponentially between its end points, so OCV(end) + max(v_rc(0), v_rc(end))
     * + R0*I bounds V(t). It is also never below the latest measurement plus the R0 step. */
    float rc_decay = expf(-horizon_s / config->rc_time_constant_s);
    float rc_end_V = (rc_decay * state->rc_voltage_V) + (r1_ohm * (1.0f - rc_decay) * charge_current_A);
    float ocv_end_V = 0.0f;
    ev_full_cell_ocv_V_and_dv_dsoc((float)ev_clamp_d(soc_upper, 0.0, 1.0), &ocv_end_V, NULL);
    float model_bound_V = ocv_end_V + fmaxf(state->rc_voltage_V, rc_end_V) + (r0_ohm * charge_current_A);
    float measured_bound_V = (state->last_voltage_V > 0.0f)
        ? (state->last_voltage_V + (r0_ohm * (charge_current_A - state->last_current_A))) : 0.0f;
    prediction->terminal_voltage_V = fmaxf(model_bound_V, measured_bound_V);

    float temperature_K = ev_clamp(coldest_C - EV_TEMPERATURE_SENSOR_MARGIN_C,
                                   EV_MIN_PLAUSIBLE_TEMPERATURE_C, EV_MAX_PLAUSIBLE_TEMPERATURE_C) + EV_KELVIN_OFFSET;

    /* Solid phase: steady-state polynomial (Subramanian) surface correction. Charging inserts
     * lithium, so the surface stoichiometry runs ahead of the bulk. Capacity fade is treated as
     * loss of active material, which raises the interfacial current density by 1/SOH.
     * The concentration gradients (solid surface excess and electrolyte depletion) relax with time
     * constants of tens of seconds to minutes, so right after a current step-down they still
     * belong to the previous current: they are evaluated at max(I, previous charge current). The
     * kinetic and ohmic terms follow the current instantly and use I. */
    float area_m2 = EV_GEOMETRIC_ELECTRODE_AREA_M2 * (config->nominal_capacity_Ah / EV_REFERENCE_CAPACITY_AH);
    float current_density_A_m2 = charge_current_A / area_m2;
    float gradient_current_density_A_m2 = fmaxf(charge_current_A, state->last_current_A) / area_m2;
    float specific_area_m_inv = (3.0f * EV_ANODE_ACTIVE_FRACTION) / EV_PARTICLE_RADIUS_M;
    float interfacial_current_A_m2 = current_density_A_m2 / (specific_area_m_inv * EV_ANODE_THICKNESS_M * (float)state->soh);
    float gradient_interfacial_A_m2 = gradient_current_density_A_m2 / (specific_area_m_inv * EV_ANODE_THICKNESS_M * (float)state->soh);
    float solid_diffusivity_m2_s = EV_SOLID_DIFFUSIVITY_REF_M2_S * ev_arrhenius(EV_DIFFUSION_ACTIVATION_ENERGY_J_MOL, temperature_K);
    float surface_stoich_rise = (gradient_interfacial_A_m2 * EV_PARTICLE_RADIUS_M)
        / (5.0f * EV_FARADAY_C_PER_MOL * solid_diffusivity_m2_s * EV_MAX_CONCENTRATION_MOL_M3);
    float surface_soc = ev_clamp((float)ev_clamp_d(soc_upper, 0.0, 1.0) + (surface_stoich_rise / anode_stoich_window), 0.0f, 1.0f);
    float surface_stoich = ev_clamp(EV_ANODE_STOICH_AT_SOC0 + (surface_soc * anode_stoich_window), 0.001f, 0.999f);

    /* Electrolyte: steady parabolic profile across the anode for a uniform reaction. The diffusion
     * time constant 4*eps*L^2/(pi^2*De_eff) is about 33 s, so steady state is reached during any
     * sustained charge and is the appropriate (conservative) bound for the horizon. */
    float porosity_factor = powf(EV_ELECTROLYTE_POROSITY, EV_BRUGGEMAN_EXPONENT);
    float conductivity_S_m = EV_ELECTROLYTE_CONDUCTIVITY_S_M * porosity_factor
        * ev_arrhenius(EV_ELECTROLYTE_ACTIVATION_ENERGY_J_MOL, temperature_K);
    float electrolyte_diffusivity_m2_s = EV_ELECTROLYTE_DIFFUSIVITY_M2_S * porosity_factor
        * ev_arrhenius(EV_ELECTROLYTE_DIFFUSIVITY_ACTIVATION_J_MOL, temperature_K);
    float depletion = ((1.0f - EV_TRANSFERENCE_NUMBER) * gradient_current_density_A_m2 * EV_ANODE_THICKNESS_M)
        / (3.0f * EV_FARADAY_C_PER_MOL * electrolyte_diffusivity_m2_s * EV_ELECTROLYTE_CONCENTRATION_MOL_M3);
    uint8_t violated = 0U;
    if (depletion >= EV_ELECTROLYTE_DEPLETION_LIMIT) {
        violated |= EV_CONSTRAINT_ELECTROLYTE_DEPLETION;
        depletion = EV_ELECTROLYTE_DEPLETION_LIMIT;
    }
    float concentration_overpotential_V = -2.0f * r_over_f * temperature_K * EV_DIFFUSION_POTENTIAL_FACTOR
        * ev_mean_log_concentration(depletion);

    /* Butler-Volmer (alpha = 0.5) inverted for the cathodic overpotential, i0 at the depleted c_e. */
    float exchange_current_A_m2 = EV_EXCHANGE_CURRENT_REF_A_M2
        * sqrtf(surface_stoich * (1.0f - surface_stoich) * (1.0f - depletion))
        * ev_arrhenius(EV_KINETIC_ACTIVATION_ENERGY_J_MOL, temperature_K);
    float kinetic_overpotential_V = 2.0f * r_over_f * temperature_K * asinhf(interfacial_current_A_m2 / (2.0f * exchange_current_A_m2));

    /* Electrolyte ohmic drop for a uniform reaction: phi_e rises by i*L/(2*kappa) from the
     * current collector to the separator, and the separator side (where plating starts) lies
     * i*L/(3*kappa) above the electrode average. The full i*L/(2*kappa) is used as a margin. */
    float electrolyte_ohmic_drop_V = (current_density_A_m2 * EV_ANODE_THICKNESS_M) / (2.0f * conductivity_S_m);

    /* Ageing resistance growth is assumed to be anode SEI film (the conservative attribution). */
    float film_drop_V = state->resistance_growth_ohm * ev_resistance_scale(coldest_C) * charge_current_A;

    prediction->anode_potential_V = ev_lut_interp(LUT_GRAPHITE_OCP, surface_soc, NULL)
        - kinetic_overpotential_V - electrolyte_ohmic_drop_V - concentration_overpotential_V - film_drop_V;

    if (!ev_is_finite(prediction->anode_potential_V) || !ev_is_finite(prediction->terminal_voltage_V)
        || !ev_is_finite(prediction->predicted_temperature_C) || !ev_is_finite_d(soc_upper)) {
        violated |= EV_CONSTRAINT_INVALID_INPUT;   /* never treat a non-finite prediction as safe */
    }
    if (soc_upper > (double)EV_MAX_SAFE_SOC) { violated |= EV_CONSTRAINT_SOC; }
    if (prediction->predicted_temperature_C > EV_MAX_CELL_TEMPERATURE_C) { violated |= EV_CONSTRAINT_TEMPERATURE_HIGH; }
    if ((charge_current_A > 0.0f)
        && (state->cold_charge_block || ((coldest_C - EV_TEMPERATURE_SENSOR_MARGIN_C) < EV_MIN_CHARGE_TEMPERATURE_C))) {
        violated |= EV_CONSTRAINT_TEMPERATURE_LOW;
    }
    if (prediction->anode_potential_V < EV_MIN_SAFE_ANODE_POTENTIAL_V) { violated |= EV_CONSTRAINT_ANODE_POTENTIAL; }
    if (prediction->terminal_voltage_V > EV_MAX_CELL_VOLTAGE_V) { violated |= EV_CONSTRAINT_CELL_VOLTAGE; }
    prediction->violated = violated;
    prediction->safe = (violated == 0U);
    return true;
}

/* ---- Current search ---------------------------------------------------------------------- */

static bool ev_current_is_safe(const EvCellConfig *config, const EvCellState *state, uint8_t cell_count, float current_A,
                               float horizon_s, uint8_t *limiting_index, EvCellPrediction *limiting_prediction) {
    uint8_t index;
    uint8_t worst_index = 0U;
    EvCellPrediction worst;
    memset(&worst, 0, sizeof(worst));
    worst.anode_potential_V = FLT_MAX;

    for (index = 0U; index < cell_count; ++index) {
        EvCellPrediction prediction;
        if (!ev_predict_cell(&config[index], &state[index], current_A, horizon_s, &prediction) || !prediction.safe) {
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
 * Largest current in [0, requested] below the first unsafe point of a 20-step upward sweep,
 * refined by 24 bisection steps. The constraints are monotone in current for this model
 * (checked by the self-test); the sweep from zero additionally stops at any non-monotone
 * unsafe band wider than 5 % of the request. Result resolution: requested / (20 * 2^24).
 */
static float ev_find_safe_charge_current(const EvCellConfig *config, const EvCellState *state, uint8_t cell_count,
                                         float requested_A, float horizon_s, uint8_t *limiting_index,
                                         uint8_t *limiting_constraints, EvCellPrediction *limiting_prediction) {
    float safe_A = 0.0f;
    float unsafe_A = 0.0f;
    uint32_t step;

    *limiting_constraints = 0U;
    if (!ev_current_is_safe(config, state, cell_count, 0.0f, horizon_s, limiting_index, limiting_prediction)) {
        *limiting_constraints = limiting_prediction->violated;
        return 0.0f;
    }
    for (step = 1U; step <= EV_COARSE_SEARCH_STEPS; ++step) {
        float candidate_A = requested_A * ((float)step / (float)EV_COARSE_SEARCH_STEPS);
        if (!ev_current_is_safe(config, state, cell_count, candidate_A, horizon_s, limiting_index, limiting_prediction)) {
            unsafe_A = candidate_A;
            break;
        }
        safe_A = candidate_A;
    }
    if (step > EV_COARSE_SEARCH_STEPS) {
        return requested_A;   /* limiting_* already describe the worst cell at the request */
    }

    for (step = 0U; step < EV_BISECTION_STEPS; ++step) {
        float mid_A = 0.5f * (safe_A + unsafe_A);
        if (ev_current_is_safe(config, state, cell_count, mid_A, horizon_s, limiting_index, limiting_prediction)) {
            safe_A = mid_A;
        } else {
            unsafe_A = mid_A;
        }
    }
    /* Binding cell and constraints just above the approved current; its prediction at the approved current. */
    (void)ev_current_is_safe(config, state, cell_count, unsafe_A, horizon_s, limiting_index, limiting_prediction);
    *limiting_constraints = limiting_prediction->violated;
    (void)ev_predict_cell(&config[*limiting_index], &state[*limiting_index], safe_A, horizon_s, limiting_prediction);
    return safe_A;
}

float ev_cell_max_charge_current(const EvCellConfig *config, const EvCellState *state, float requested_A,
                                 float horizon_s, uint8_t *binding_constraints, EvCellPrediction *at_limit) {
    EvCellPrediction prediction;
    uint8_t index = 0U;
    uint8_t constraints = EV_CONSTRAINT_INVALID_INPUT;
    float approved_A = 0.0f;

    memset(&prediction, 0, sizeof(prediction));
    prediction.violated = EV_CONSTRAINT_INVALID_INPUT;
    if (ev_config_is_valid(config) && ev_state_is_valid(state)
        && ev_in_range(requested_A, 0.0f, EV_MAX_PLAUSIBLE_CURRENT_A)
        && ev_in_range(horizon_s, FLT_MIN, EV_MAX_PREDICTION_HORIZON_S)) {
        approved_A = ev_find_safe_charge_current(config, state, 1U, requested_A, horizon_s, &index, &constraints, &prediction);
    }
    if (binding_constraints != NULL) { *binding_constraints = constraints; }
    if (at_limit != NULL) { *at_limit = prediction; }
    return approved_A;
}

/* ---- Estimation -------------------------------------------------------------------------- */

static void ev_update_ageing(const EvCellConfig *config, EvCellState *state, float current_A, float dt_s) {
    double processed_Ah = ((double)fabsf(current_A) * (double)dt_s) / 3600.0;
    double full_equivalent_cycles = processed_Ah / (2.0 * (double)config->nominal_capacity_Ah);
    float temperature_K = ev_clamp(state->temperature_C, EV_MIN_PLAUSIBLE_TEMPERATURE_C, EV_MAX_PLAUSIBLE_TEMPERATURE_C) + EV_KELVIN_OFFSET;
    double loss = ((EV_CALENDAR_LOSS_PER_S * (double)dt_s) + (EV_CYCLE_LOSS_PER_FEC * full_equivalent_cycles))
        * (double)ev_arrhenius(EV_AGING_ACTIVATION_ENERGY_J_MOL, temperature_K);
    double new_soh = ev_clamp_d(state->soh - loss, EV_MIN_SOH, state->soh);   /* never increases */
    float growth_ohm = state->resistance_growth_ohm
        + (float)(EV_RESISTANCE_GROWTH_PER_LOSS * (double)config->resistance_ohm * (state->soh - new_soh));

    state->throughput_Ah += processed_Ah;
    state->resistance_growth_ohm = ev_clamp(growth_ohm, 0.0f, 1.0f);
    state->soh = new_soh;
}

/* Largest change of OCV(to) - OCV(from) when both SOCs are shifted by up to +/- shift. */
static float ev_ocv_change_spread_V(double from_soc, double to_soc, double shift) {
    float base_from = 0.0f;
    float base_to = 0.0f;
    float spread = 0.0f;
    int side;
    ev_full_cell_ocv_V_and_dv_dsoc((float)from_soc, &base_from, NULL);
    ev_full_cell_ocv_V_and_dv_dsoc((float)to_soc, &base_to, NULL);
    for (side = -1; side <= 1; side += 2) {
        float from = 0.0f;
        float to = 0.0f;
        ev_full_cell_ocv_V_and_dv_dsoc((float)ev_clamp_d(from_soc + ((double)side * shift), 0.0, 1.0), &from, NULL);
        ev_full_cell_ocv_V_and_dv_dsoc((float)ev_clamp_d(to_soc + ((double)side * shift), 0.0, 1.0), &to, NULL);
        spread = fmaxf(spread, fabsf((to - from) - (base_to - base_from)));
    }
    return spread;
}

/* Voltage-tracking check, run before the EKF update. Over a window the open-loop model is
 * propagated from the EKF state at the window start (coulomb counting, 1-RC, R0), and the change
 * of the reading is compared with the change of the model. What parameter errors can explain is
 * subtracted: 30 % of R0 times the present current step, 100 % of R1 times the largest step in
 * the window, 20 % of the open-loop OCV change (capacity and OCV-slope error) and the change in
 * that OCV difference if the SOC were 3 sigma higher or lower (an unconverged EKF must still be
 * able to learn from the reading). A frozen reading, with or without noise, does not follow the R0 step when the current
 * changes or the OCV rise during a charge, so the mismatch grows until it faults.
 *  - mismatch > EV_TRACK_SUSPECT_V: the EKF stops using the reading for the rest of the window,
 *    otherwise a frozen in-range value drags the SOC estimate along and the cell is overcharged;
 *  - mismatch > EV_TRACK_FAULT_V for EV_TRACK_FAULT_CYCLES cycles: plausibility fault, latched
 *    for the rest of the window (stopping the charge removes the R0 step a frozen reading never
 *    showed, so the mismatch alone would clear and charging would resume in bursts);
 *  - the reading has moved by EV_TRACK_VALIDATE_V more than the tolerance, in agreement with the
 *    model (for example the R0 step when the charger stops): the sensor is shown to work, the fault
 *    and suspicion are cleared and a new window starts after the EKF update;
 *  - at rest (|I| < C/50) with no latched fault, the window restarts once two consecutive
 *    readings agree, so slow relaxation cannot build up a mismatch. At rest the model predicts
 *    no change, so nothing is lost; a sensor frozen at rest is caught when the current changes.
 * Returns true when the EKF must not use the reading. */
static bool ev_voltage_tracking_suspect(const EvCellConfig *config, EvCellState *state, float current_A,
                                        float measured_voltage_V, float dt_s) {
    double a;
    float ocv = 0.0f;
    float r0_ohm = ev_effective_r0(config, state, state->temperature_C);
    float r1_ohm = ev_effective_r1(config, state->temperature_C);
    float expected_change_V;
    float measured_change_V;
    float tolerance_V;
    float mismatch_V;

    if (state->track_start_V <= 0.0f) {
        return false;   /* no window yet: one starts after the EKF update */
    }
    a = exp(-(double)dt_s / (double)config->rc_time_constant_s);
    state->track_soc = ev_clamp_d(state->track_soc
        + (((double)current_A * (double)dt_s) / ((double)config->nominal_capacity_Ah * state->soh * 3600.0)), 0.0, 1.0);
    state->track_rc_voltage_V = (float)ev_clamp_d((a * (double)state->track_rc_voltage_V) + ((double)r1_ohm * (1.0 - a) * (double)current_A),
                                                  -EV_MAX_RC_VOLTAGE_V, EV_MAX_RC_VOLTAGE_V);
    state->track_max_current_step_A = fmaxf(state->track_max_current_step_A, fabsf(current_A - state->track_start_current_A));
    ev_full_cell_ocv_V_and_dv_dsoc((float)state->track_soc, &ocv, NULL);
    expected_change_V = (ocv + state->track_rc_voltage_V + (r0_ohm * current_A)) - state->track_start_model_V;
    measured_change_V = measured_voltage_V - state->track_start_V;
    tolerance_V = ((float)EV_RESISTANCE_UNCERTAINTY * r0_ohm * fabsf(current_A - state->track_start_current_A))
        + (EV_TRACK_R1_UNCERTAINTY * r1_ohm * state->track_max_current_step_A)
        + (EV_TRACK_CAPACITY_UNCERTAINTY * fabsf(ocv - state->track_start_ocv_V))
        + ev_ocv_change_spread_V(state->track_start_soc, state->track_soc, EV_SOC_SIGMA_MULTIPLIER * sqrt(state->covariance[0]));
    mismatch_V = fabsf(measured_change_V - expected_change_V) - tolerance_V;

    if (mismatch_V > EV_TRACK_FAULT_V) {
        if (state->track_fail_cycles < UINT16_MAX) { state->track_fail_cycles++; }
    } else if (state->track_fail_cycles < EV_TRACK_FAULT_CYCLES) {
        state->track_fail_cycles = 0U;   /* not latched yet */
    }
    if (mismatch_V > EV_TRACK_SUSPECT_V) {
        state->track_suspect = true;
    } else if ((fabsf(measured_change_V) >= (EV_TRACK_VALIDATE_V + tolerance_V))
               || ((fabsf(current_A) < (EV_TRACK_REST_C_RATE * config->nominal_capacity_Ah))
                   && (state->track_fail_cycles < EV_TRACK_FAULT_CYCLES)
                   && (fabsf(measured_voltage_V - state->last_voltage_V) <= EV_TRACK_SPIKE_V))) {
        state->track_start_V = 0.0f;   /* validated or at rest: restart from the updated EKF state */
        state->track_fail_cycles = 0U;
        state->track_suspect = false;
        return false;
    }
    return state->track_suspect;
}

/* Starts a voltage-tracking window at this reading, from the updated EKF state. */
static void ev_start_tracking_window(const EvCellConfig *config, EvCellState *state, float current_A, float measured_voltage_V) {
    float ocv = 0.0f;
    ev_full_cell_ocv_V_and_dv_dsoc((float)state->soc, &ocv, NULL);
    state->track_start_V = measured_voltage_V;
    state->track_start_ocv_V = ocv;
    state->track_start_soc = state->soc;
    state->track_soc = state->soc;
    state->track_rc_voltage_V = state->rc_voltage_V;
    state->track_start_model_V = ocv + state->rc_voltage_V + (ev_effective_r0(config, state, state->temperature_C) * current_A);
    state->track_start_current_A = current_A;
    state->track_max_current_step_A = 0.0f;
    state->track_fail_cycles = 0U;
    state->track_suspect = false;
}

/* Allowed measured current for the interval that just ended: the approved command plus what is
 * left of a linear ramp-down from a higher earlier command, averaged over the interval. */
static float ev_allowed_current_A(const EvCellState *state, float dt_s) {
    const float ramp_s = EV_OVER_CURRENT_RAMP_S;
    float drop_A = state->ramp_from_A - state->commanded_current_A;
    float start_s = fminf(state->ramp_age_s, ramp_s);
    float end_s = fminf(state->ramp_age_s + dt_s, ramp_s);
    /* integral of max(0, 1 - t / ramp) over [age, age + dt] */
    float remaining_s = (end_s - ((end_s * end_s) / (2.0f * ramp_s))) - (start_s - ((start_s * start_s) / (2.0f * ramp_s)));
    float allowed_A = state->commanded_current_A + ((drop_A > 0.0f) ? (drop_A * (remaining_s / dt_s)) : 0.0f);
    return (allowed_A * (1.0f + EV_OVER_CURRENT_TOLERANCE_REL)) + EV_OVER_CURRENT_TOLERANCE_A;
}

/* Records the new command; a lower one starts a ramp from the level the charger may still be at. */
static void ev_record_command(EvCellState *state, float approved_A) {
    float progress = fminf(state->ramp_age_s / EV_OVER_CURRENT_RAMP_S, 1.0f);
    float level_A = state->commanded_current_A + (fmaxf(state->ramp_from_A - state->commanded_current_A, 0.0f) * (1.0f - progress));
    if (approved_A < level_A) {
        state->ramp_from_A = level_A;
        state->ramp_age_s = 0.0f;
    } else {
        state->ramp_from_A = approved_A;
        state->ramp_age_s = 0.0f;
    }
    state->commanded_current_A = approved_A;
}

/* Two-state EKF, x = [soc, v_rc], model V = OCV(soc) + v_rc + R0(T)*I, Joseph-form update.
 * Innovations beyond EV_INNOVATION_GATE_SIGMA are rejected (time update only) and counted. */
static void ev_update_estimator(const EvCellConfig *config, EvCellState *state, float current_A,
                                float measured_voltage_V, bool voltage_valid, float dt_s) {
    double capacity_Ah = (double)config->nominal_capacity_Ah * state->soh;
    double a = exp(-(double)dt_s / (double)config->rc_time_constant_s);
    double soc = state->soc + (((double)current_A * (double)dt_s) / (capacity_Ah * 3600.0));
    double v_rc = (a * (double)state->rc_voltage_V) + ((double)ev_effective_r1(config, state->temperature_C) * (1.0 - a) * (double)current_A);
    double p00 = state->covariance[0] + (EV_SOC_PROCESS_NOISE_PER_S * (double)dt_s);
    double p01 = a * state->covariance[1];
    double p11 = (a * a * state->covariance[2]) + (EV_RC_PROCESS_NOISE_V2_PER_S * (double)dt_s);

    if (voltage_valid) {
        float ocv = 0.0f;
        float slope = 0.0f;
        ev_full_cell_ocv_V_and_dv_dsoc((float)ev_clamp_d(soc, 0.0, 1.0), &ocv, &slope);
        double h0 = (double)slope;
        double predicted_V = (double)ocv + v_rc + ((double)ev_effective_r0(config, state, state->temperature_C) * (double)current_A);
        /* Resistance errors bias the voltage in proportion to the current, so the measurement
         * noise grows with it: the filter learns SOC mainly from low-current samples. */
        double model_error_V = EV_RESISTANCE_UNCERTAINTY * (double)(ev_effective_r0(config, state, state->temperature_C)
            + ev_effective_r1(config, state->temperature_C)) * (double)current_A;
        double r = EV_VOLTAGE_MEASUREMENT_NOISE_V2 + (model_error_V * model_error_V);
        double s = (h0 * h0 * p00) + (2.0 * h0 * p01) + p11 + r;
        double innovation = (double)measured_voltage_V - predicted_V;

        if ((innovation * innovation) > (EV_INNOVATION_GATE_SIGMA * EV_INNOVATION_GATE_SIGMA * s)) {
            if (state->innovation_fail_cycles < UINT16_MAX) { state->innovation_fail_cycles++; }
        } else {
            double k0 = ((p00 * h0) + p01) / s;
            double k1 = ((p01 * h0) + p11) / s;
            /* A = I - K*H */
            double a00 = 1.0 - (k0 * h0);
            double a01 = -k0;
            double a10 = -k1 * h0;
            double a11 = 1.0 - k1;
            /* B = A * P */
            double b00 = (a00 * p00) + (a01 * p01);
            double b01 = (a00 * p01) + (a01 * p11);
            double b10 = (a10 * p00) + (a11 * p01);
            double b11 = (a10 * p01) + (a11 * p11);

            state->innovation_fail_cycles = 0U;
            soc += k0 * innovation;
            v_rc += k1 * innovation;
            p00 = (b00 * a00) + (b01 * a01) + (k0 * r * k0);
            p01 = (b00 * a10) + (b01 * a11) + (k0 * r * k1);
            p11 = (b10 * a10) + (b11 * a11) + (k1 * r * k1);
        }
    }

    state->soc = ev_clamp_d(soc, 0.0, 1.0);
    state->rc_voltage_V = (float)ev_clamp_d(v_rc, -EV_MAX_RC_VOLTAGE_V, EV_MAX_RC_VOLTAGE_V);
    state->covariance[0] = ev_clamp_d(p00, 0.0, EV_MAX_SOC_VARIANCE);
    state->covariance[1] = p01;
    state->covariance[2] = ev_clamp_d(p11, 0.0, 1.0);
}

/* Updates one cell from this cycle's measurements and returns its EV_FAULT_* flags. */
static uint16_t ev_update_cell(const EvCellConfig *config, EvCellState *state, float current_A,
                               float voltage_V, float temperature_C, float dt_s) {
    uint16_t flags = 0U;
    bool holding = (state->fault_hold_cycles > 0U);
    bool voltage_valid = ev_in_range(voltage_V, EV_MIN_PLAUSIBLE_VOLTAGE_V, EV_MAX_PLAUSIBLE_VOLTAGE_V);
    bool temperature_valid = ev_in_range(temperature_C, EV_MIN_PLAUSIBLE_TEMPERATURE_C, EV_MAX_PLAUSIBLE_TEMPERATURE_C);
    float hysteresis_C = holding ? EV_TEMPERATURE_HYSTERESIS_C : 0.0f;
    float over_voltage_V = EV_MAX_CELL_VOLTAGE_V + EV_OVER_VOLTAGE_MARGIN_V - (holding ? EV_VOLTAGE_HYSTERESIS_V : 0.0f);
    float over_current_limit_A;

    if (!ev_config_is_valid(config)) {
        flags |= EV_FAULT_CONFIG;
        if (ev_state_is_valid(state)) {
            state->covariance[0] = EV_MAX_SOC_VARIANCE;   /* charge is not being counted: SOC is stale */
        }
        state->track_start_V = 0.0f;
    } else {
        bool suspect = false;
        if (!ev_state_is_valid(state)) {
            ev_reset_estimator(state);
            flags |= EV_FAULT_STATE;
        }
        ev_update_ageing(config, state, current_A, dt_s);
        if (voltage_valid) {
            suspect = ev_voltage_tracking_suspect(config, state, current_A, voltage_V, dt_s);
        }
        ev_update_estimator(config, state, current_A, voltage_V, voltage_valid && !suspect, dt_s);
        if (!voltage_valid) {
            state->track_start_V = 0.0f;
            state->track_fail_cycles = 0U;
            state->track_suspect = false;
        } else if ((state->track_start_V <= 0.0f) && (state->last_voltage_V > 0.0f)
                   && (fabsf(voltage_V - state->last_voltage_V) <= EV_TRACK_SPIKE_V)) {
            /* a window starts only on two consistent readings, never on an isolated spike */
            ev_start_tracking_window(config, state, current_A, voltage_V);
        }
        if (state->track_fail_cycles >= EV_TRACK_FAULT_CYCLES) { flags |= EV_FAULT_VOLTAGE_PLAUSIBILITY; }
        if (temperature_valid) {
            state->temperature_C = temperature_C;
        } else {
            state->temperature_C = ev_predict_temperature_C(config, state->temperature_C, current_A * current_A
                * (ev_effective_r0(config, state, state->temperature_C) + ev_effective_r1(config, state->temperature_C)), dt_s);
        }
        if (state->innovation_fail_cycles >= EV_INNOVATION_FAULT_CYCLES) { flags |= EV_FAULT_VOLTAGE_PLAUSIBILITY; }
        if (state->soh <= EV_MIN_SOH) { flags |= EV_FAULT_END_OF_LIFE; }
    }
    state->last_voltage_V = voltage_valid ? voltage_V : 0.0f;
    state->last_current_A = current_A;

    /* The charger must follow the previous approval, allowing for its ramp-down time. */
    over_current_limit_A = ev_allowed_current_A(state, dt_s);
    state->ramp_age_s = fminf(state->ramp_age_s + dt_s, EV_OVER_CURRENT_RAMP_S);
    if (current_A > over_current_limit_A) {
        state->over_current_time_s = fminf(state->over_current_time_s + dt_s, FLT_MAX);
    } else {
        state->over_current_time_s = 0.0f;
    }
    if (state->over_current_time_s >= EV_OVER_CURRENT_CONFIRM_S) { flags |= EV_FAULT_OVER_CURRENT; }

    if ((state->temperature_C - EV_TEMPERATURE_SENSOR_MARGIN_C) < EV_MIN_CHARGE_TEMPERATURE_C) {
        state->cold_charge_block = true;
    } else if ((state->temperature_C - EV_TEMPERATURE_SENSOR_MARGIN_C) >= (EV_MIN_CHARGE_TEMPERATURE_C + EV_TEMPERATURE_HYSTERESIS_C)) {
        state->cold_charge_block = false;
    }

    if (!voltage_valid) { flags |= EV_FAULT_VOLTAGE_SENSOR; }
    if (voltage_valid && (voltage_V > over_voltage_V)) { flags |= EV_FAULT_OVER_VOLTAGE; }
    if (!temperature_valid) { flags |= EV_FAULT_TEMPERATURE_SENSOR; }
    if (state->temperature_C > (EV_MAX_CELL_TEMPERATURE_C - hysteresis_C)) { flags |= EV_FAULT_OVER_TEMPERATURE; }
    if (state->temperature_C < (EV_MIN_CELL_TEMPERATURE_C + hysteresis_C)) { flags |= EV_FAULT_UNDER_TEMPERATURE; }

    if (flags != 0U) {
        state->fault_hold_cycles = (uint16_t)EV_FAULT_HOLD_CYCLES;
    } else if (state->fault_hold_cycles > 0U) {
        state->fault_hold_cycles--;
        if (state->fault_hold_cycles > 0U) { flags |= EV_FAULT_HOLD; }
    }
    state->fault_flags = flags;
    return flags;
}

/* Current or dt could not be trusted: keep the SOC mean, widen its uncertainty. */
static void ev_inflate_uncertainty(const EvCellConfig *config, EvCellState *state, bool dt_valid, float dt_s) {
    if (!ev_config_is_valid(config) || !ev_state_is_valid(state)) {
        return;
    }
    if (dt_valid) {
        double max_soc_change = ((double)EV_MAX_PLAUSIBLE_CURRENT_A * (double)dt_s)
            / ((double)config->nominal_capacity_Ah * state->soh * 3600.0);
        state->covariance[0] = ev_clamp_d(state->covariance[0] + (max_soc_change * max_soc_change), 0.0, EV_MAX_SOC_VARIANCE);
    } else if (state->covariance[0] < EV_UNKNOWN_DT_SOC_VARIANCE) {
        state->covariance[0] = EV_UNKNOWN_DT_SOC_VARIANCE;
    }
}

/* A pack-level input fault is held like a cell fault so that a flapping signal cannot toggle
 * the charger every cycle. */
static void ev_hold_input_fault(EvCellState state[EV_NUM_CELLS]) {
    uint8_t index;
    for (index = 0U; index < EV_NUM_CELLS; ++index) {
        state[index].fault_hold_cycles = (uint16_t)EV_FAULT_HOLD_CYCLES;
        state[index].fault_flags |= EV_FAULT_INPUT;
    }
}

/* ---- Control cycle ----------------------------------------------------------------------- */

static EvPackResult ev_step_unrecorded(const EvCellConfig config[EV_NUM_CELLS], EvCellState state[EV_NUM_CELLS],
                                       float measured_pack_current_A, float requested_charge_current_A,
                                       const float measured_voltages_V[EV_NUM_CELLS], const float measured_temperatures_C[EV_NUM_CELLS],
                                       float control_dt_s, float prediction_horizon_s) {
    EvPackResult result;
    EvCellPrediction limiting_prediction;
    uint8_t limiting_index = 0U;
    uint8_t limiting_constraints = 0U;
    uint8_t index;
    bool dt_valid = ev_in_range(control_dt_s, FLT_MIN, EV_MAX_CONTROL_DT_S);
    bool current_valid = ev_in_range(measured_pack_current_A, -EV_MAX_PLAUSIBLE_CURRENT_A, EV_MAX_PLAUSIBLE_CURRENT_A);

    memset(&result, 0, sizeof(result));
    result.decision = EV_DECISION_FAULT;
    result.limiting_cell_index = EV_NO_CELL;
    result.fault_flags = EV_FAULT_INPUT;

    if ((config == NULL) || (state == NULL) || (measured_voltages_V == NULL) || (measured_temperatures_C == NULL)) {
        return result;
    }
    if (!dt_valid || !current_valid) {
        for (index = 0U; index < EV_NUM_CELLS; ++index) {
            ev_inflate_uncertainty(&config[index], &state[index], dt_valid, control_dt_s);
        }
        ev_hold_input_fault(state);
        return result;
    }

    result.fault_flags = 0U;
    for (index = 0U; index < EV_NUM_CELLS; ++index) {
        uint16_t flags = ev_update_cell(&config[index], &state[index], measured_pack_current_A,
                                        measured_voltages_V[index], measured_temperatures_C[index], control_dt_s);
        if ((flags != 0U) && (result.limiting_cell_index == EV_NO_CELL)) {
            result.limiting_cell_index = index;
            result.limiting_temperature_C = state[index].temperature_C;
            result.limiting_soc = (float)state[index].soc;
        }
        result.fault_flags |= flags;
    }

    if (!ev_is_finite(requested_charge_current_A) || !ev_in_range(prediction_horizon_s, FLT_MIN, EV_MAX_PREDICTION_HORIZON_S)) {
        result.fault_flags |= EV_FAULT_INPUT;
        ev_hold_input_fault(state);
        return result;
    }
    result.requested_charge_current_A = ev_clamp(requested_charge_current_A, 0.0f, EV_MAX_PLAUSIBLE_CURRENT_A);
    if (result.fault_flags != 0U) {
        return result;
    }

    /* The approved current flows for at least one control period. */
    result.approved_charge_current_A = ev_find_safe_charge_current(config, state, (uint8_t)EV_NUM_CELLS,
        result.requested_charge_current_A, fmaxf(prediction_horizon_s, control_dt_s),
        &limiting_index, &limiting_constraints, &limiting_prediction);
    result.limiting_cell_index = limiting_index;
    result.limiting_anode_potential_V = limiting_prediction.anode_potential_V;
    result.limiting_terminal_voltage_V = limiting_prediction.terminal_voltage_V;
    result.limiting_temperature_C = limiting_prediction.predicted_temperature_C;
    result.limiting_soc = limiting_prediction.predicted_soc;
    if ((result.approved_charge_current_A + EV_CURRENT_TOLERANCE_A) < result.requested_charge_current_A) {
        result.decision = EV_DECISION_LIMIT;
        result.limiting_constraints = limiting_constraints;
    } else {
        result.decision = EV_DECISION_ALLOW;
    }
    return result;
}

EvPackResult ev_step(const EvCellConfig config[EV_NUM_CELLS], EvCellState state[EV_NUM_CELLS],
                     float measured_pack_current_A, float requested_charge_current_A,
                     const float measured_voltages_V[EV_NUM_CELLS], const float measured_temperatures_C[EV_NUM_CELLS],
                     float control_dt_s, float prediction_horizon_s) {
    EvPackResult result = ev_step_unrecorded(config, state, measured_pack_current_A, requested_charge_current_A,
                                             measured_voltages_V, measured_temperatures_C, control_dt_s, prediction_horizon_s);
    uint8_t index;
    if (state != NULL) {
        /* Command for the next interval, checked against the next measured current. */
        for (index = 0U; index < EV_NUM_CELLS; ++index) {
            ev_record_command(&state[index], result.approved_charge_current_A);
        }
    }
    return result;
}
