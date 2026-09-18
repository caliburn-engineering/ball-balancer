// src/app_state.h
#pragma once

#include "linear_system.h"
#include "auto_balance.h"
#include "analysis/controller_builders.h"
#include "analysis/frequency_response.h"
#include "analysis/loop_diagnostics.h"
#include "analysis/lqr.h"
#include "analysis/model_library.h"
#include "analysis/pole_zero.h"
#include "analysis/system_connect.h"
#include "analysis/system_properties.h"
#include "analysis/time_response.h"
#include "imgui.h"

#include <algorithm>
#include <cstring>
#include <string>

namespace caliburn {

constexpr int NUM_SYSTEMS = 4;

inline const char* system_names[NUM_SYSTEMS] = {
    "Plant", "Controller", "Open-Loop", "Closed-Loop"};

inline ImVec4 system_colors[NUM_SYSTEMS] = {
    {0.22f, 0.74f, 0.97f, 1.0f},  // Plant — blue  #38bdf8
    {0.66f, 0.33f, 0.97f, 1.0f},  // Controller — purple  #a855f7
    {0.20f, 0.83f, 0.60f, 1.0f},  // Open-Loop — green  #34d399
    {0.98f, 0.45f, 0.09f, 1.0f},  // Closed-Loop — orange  #f97316
};

inline ImU32 system_colors_u32[NUM_SYSTEMS] = {
    IM_COL32(56, 189, 248, 255),
    IM_COL32(168, 85, 247, 255),
    IM_COL32(52, 211, 153, 255),
    IM_COL32(249, 115, 22, 255),
};

// LQR is appended, never inserted.  The combo-box label array in the model
// panel is index-coupled to this enum, and #16 proposed inserting the new
// member in the middle — which would silently shift every label after it.
// Appending shifts nothing; the combo now sizes itself from the array.
enum class ControllerType { None, PID, LeadLag, StateSpace, GainMatrix, LQR };
enum class InputType { Step, Impulse, Ramp };

// Flat 4-way; two entries are conditionally valid.  PlantLocus is the old
// UnityFB, unchanged and honestly relabelled — it sweeps a proportional gain on
// the bare plant channel and does not put the controller in the path.
// LoopLocus sweeps kappa on the whole compensator.  See issue #10.
enum class PZMode { PoleZero, PlantLocus, LoopLocus, StateFB };

struct AppState {
    // --- Model ---
    int preset_index = 0;
    LinearSystem plant;
    ControllerType ctrl_type = ControllerType::None;
    LinearSystem ctrl_ss;
    Eigen::MatrixXd ctrl_K;

    // --- LQR design weights ---
    // Diagonals only.  A full Q and R would be 49 + 9 editable cells on the
    // cascade plant, and the off-diagonal terms have no physical reading for
    // a user picking a tuning: the diagonal says "how much do I care about
    // this state" and "how expensive is this leg", which is the whole
    // vocabulary the design surface needs.  Resized by the recompute guard
    // whenever the plant's dimensions change.
    Eigen::VectorXd lqr_q;
    Eigen::VectorXd lqr_r;
    LqrResult lqr_result;

    // --- Derived systems ---
    LinearSystem systems[NUM_SYSTEMS];
    bool system_valid[NUM_SYSTEMS] = {true, false, false, false};

    // --- Controller: loop list ---
    // Maintained for the current plant regardless of ctrl_type, so
    // None -> PID -> None -> PID round-trips without losing tuning and picking
    // PID needs no seeding path of its own.  See issue #2.
    std::vector<Loop> loops;
    int cached_p = -1;  // plant dimensions the loop list was seeded for
    int cached_m = -1;
    int selected_loop = -1;  // last cell clicked in the pairing grid; -1 = none

    // Per-output scale: the maximum allowed deviation for that output.  Channel
    // Share is not scale-invariant and S&P's scaling is engineering intent, not
    // a model property, so it comes from the user.  Reseeded to 1.0 by the same
    // (p, m) guard as `loops`.  See issue #8.
    std::vector<float> output_scales;

    // --- Channel selection (rule C, issue #9) ---
    int output_i = 0;  // [0, p) — plant output, controller input, loop output
    int input_j = 0;   // [0, m) — plant input, controller output
    int ref_r = 0;     // [0, p) — reference channel of the open/closed loop
    bool show_all_channels = false;

    // Set by the recompute path: false when the selected channel has no loop on
    // it.  Dead channels are suppressed with a stated reason rather than drawn
    // with a fabricated flat 0 phase.
    bool channel_live[NUM_SYSTEMS] = {true, true, true, true};
    std::string channel_reason[NUM_SYSTEMS];

    // --- Analysis results ---
    FrequencyResponse bode[NUM_SYSTEMS];
    PoleZeroResult pole_zero[NUM_SYSTEMS];
    TimeResponse time_resp[NUM_SYSTEMS];
    PropertyResult controllability{};
    PropertyResult observability{};
    PropertyResult cl_controllability{};
    PropertyResult cl_observability{};
    std::vector<RootLocusPoint> root_locus;

    // Cached, not recomputed per frame: the structurally-dead-channel test
    // sweeps the whole Bode grid once per loop.
    LoopDiagnostics diagnostics;

    // --- Physical parameters (mutable copy of current preset's params) ---
    std::vector<PhysicalParam> current_params;

    // All-channel results (populated when show_all_channels is true)
    std::vector<FrequencyResponse> bode_grid[NUM_SYSTEMS];

    // --- Recompute flag ---
    bool needs_recompute = true;

    // --- Panel visibility ---
    bool show_pole_zero = true;
    bool show_bode = true;
    bool show_nyquist = true;
    bool show_time_response = true;

    // --- Trace visibility ---
    bool trace_visible[NUM_SYSTEMS] = {true, false, false, false};

    // --- Bode settings ---
    float freq_min_hz = 0.01f;
    float freq_max_hz = 100.0f;
    int num_freq_points = 500;

    // --- Pairing diagnostic (RGA / Channel Share) ---
    // One omega, plus a collapsed RGA-number-vs-omega sweep across the existing
    // Bode grid.  S&P rule 2 gets no separate steady-state frequency: it is a
    // DIC theorem requiring a stable plant, and neither Ball-Balancer (poles at
    // the origin) nor Inverted Pendulum qualifies.  See issues #4, #6.
    float diag_omega_hz = 1.0f;
    bool diag_sweep = false;

    // --- Time response settings ---
    InputType input_type = InputType::Step;
    float amplitude = 1.0f;
    float slope = 1.0f;
    float duration = 10.0f;
    float dt_sim = 0.01f;

    // --- Pole-zero / root locus settings ---
    PZMode pz_mode = PZMode::PoleZero;
    float rl_k_min = 0.0f;
    float rl_k_max = 100.0f;
    int rl_num_points = 200;
    float rl_current_k = 1.0f;
    float rl_alpha_min = 0.0f;
    float rl_alpha_max = 2.0f;
    float rl_current_alpha = 1.0f;

    // Loop Locus.  A third explicit triple beside rl_k_* and rl_alpha_*: the
    // three scales are genuinely different (kappa 0.01-100 log, K 0-100 linear,
    // alpha 0-2 linear), so a shared triple would re-default on nearly every
    // mode switch.  kappa_min > 0 is a correctness constraint: kappa = 0 would
    // collapse a loop's states and make the closed-loop dimension ragged
    // mid-sweep, which matchPoles reads out of bounds.
    float rl_kappa_min = 0.01f;
    float rl_kappa_max = 100.0f;
    float rl_kappa_current = 1.0f;
    double rl_loop_gain_margin = -1.0;  // kappa*, or < 0 for no crossing

    // --- Transfer function parameterization (1st/2nd order SISO) ---
    bool plant_is_1st_order = false;
    bool plant_is_2nd_order = false;
    float tf_K = 1.0f;
    float tf_tau = 1.0f;
    float tf_wn = 1.0f;
    float tf_zeta = 0.7f;

    // --- Text fields for matrix entry ---
    char A_text[1024] = "";
    char B_text[1024] = "";
    char C_text[1024] = "";
    char D_text[1024] = "";
    char Ac_text[1024] = "";
    char Bc_text[1024] = "";
    char Cc_text[1024] = "";
    char Dc_text[1024] = "";
    char K_text[1024] = "";

    // --- Cursors (for synchronized vertical lines + value readout) ---
    double time_cursor_x = -1.0;
    double bode_cursor_x = -1.0;
    double time_x_min = 0.0;
    double time_x_max = 10.0;
};

// Seed the identity pairing y[k] -> u[k] for k < min(p, m), and reset the
// per-output scales.  Applied unconditionally on a dimension change.
//
// On Ball-Balancer this seeds a structurally DEAD pairing, deliberately: B(2,1)
// and B(3,0) mean input 0 drives output 1's chain, so G is exactly
// anti-diagonal and identity pairing leaves both loops open at every frequency.
// That is the pairing diagnostic's demonstration case — RGA of an anti-diagonal
// G is [[0,1],[1,0]], the loudest possible signal — and smart seeding would
// delete the lesson.  See issue #2.
inline void seedIdentityLoops(AppState& state) {
    state.loops.clear();
    const int p = state.plant.outputs();
    const int m = state.plant.inputs();
    for (int k = 0; k < std::min(p, m); ++k) {
        state.loops.push_back(Loop{k, k, {}, {}});
    }
    state.output_scales.assign(static_cast<std::size_t>(p), 1.0f);
    state.selected_loop = state.loops.empty() ? -1 : 0;
}

// Helper: write matrix to text buffer for ImGui display
inline void matrixToTextBuf(const Eigen::MatrixXd& mat, char* buf,
                            size_t buf_size) {
    std::string s = matrixToString(mat);
    std::strncpy(buf, s.c_str(), buf_size - 1);
    buf[buf_size - 1] = '\0';
}

}  // namespace caliburn
