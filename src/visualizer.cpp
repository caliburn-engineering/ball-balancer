// src/visualizer.cpp
#include "app_state.h"
#include "analysis/loop_locus.h"

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#include <emscripten/html5.h>  // emscripten_get_device_pixel_ratio — see uiScale()
#else
#include <glad/glad.h>
#endif
#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_internal.h"  // DockBuilder — the first-run layout
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "implot.h"

#include <algorithm>
#include <cstdio>
#include <string>

#include "panels/model_panel.h"
#include "panels/properties_panel.h"
#include "panels/pole_zero_panel.h"
#include "panels/bode_panel.h"
#include "panels/nyquist_panel.h"
#include "panels/time_response_panel.h"

#include "plate_view.h"
#include "sim_step.h"

// The design surface's current answer, in the form the plate needs it.
//
// Assembled here rather than inside PlateView because only this side knows
// which preset is loaded, which controller type is selected and whether the
// Riccati solve landed.  The plate adds the two checks it alone can make — the
// gain's shape and whether the plant described here is the one it simulates.
//
// The servo lag and the home leg angle are handed over whether or not a gain
// is on offer: they are properties of the plate, and the model panel's sliders
// are the one place they live.
static void handDesignToPlate(caliburn::AppState& state,
                              const std::vector<caliburn::ModelEntry>& presets,
                              caliburn::PlateView& plate) {
    const bool is_cascade =
        state.preset_index >= 0 &&
        state.preset_index < static_cast<int>(presets.size()) &&
        caliburn::isCascadeModel(presets[state.preset_index]);

    // Assembled by `cascadeDesign`, which every test harness also calls: the
    // plant the gain was designed against and the plant the plate checks
    // against cannot be assembled two different ways, and neither can the
    // plant a test claims to be measuring.
    caliburn::AutoBalanceDesign d;
    if (is_cascade) d = caliburn::cascadeDesign(state.current_params);

    std::string reason;
    bool offered = false;
    if (!is_cascade) {
        reason = "plant is not the Ball-Balancer Cascade";
    } else if (state.ctrl_type != caliburn::ControllerType::LQR) {
        reason = "select LQR as the controller type";
    } else if (!state.lqr_result.success) {
        reason = "the LQR solve failed";
    } else {
        d.K = state.lqr_result.K;
        offered = true;
    }
    plate.setDesign(d, offered, reason);
}

struct FrameContext {
    GLFWwindow* window;
    caliburn::AppState* state;
    std::vector<caliburn::ModelEntry>* presets;
    caliburn::PlateView* plate;
};

// How many framebuffer pixels one logical UI pixel should occupy.
//
// ImGui lays out in framebuffer pixels, and the web shell sizes the canvas
// backing store to innerWidth * devicePixelRatio.  On a 390 px phone at dpr 3
// that is a 1170 px framebuffer displayed across 390 px of glass, so a 16 px
// glyph arrives on screen about five pixels tall — working, animating, and
// completely unreadable (#25).  Scaling the font atlas and every style metric
// by the same ratio puts a 16 px glyph back at 16 px *as seen*, and because the
// dock layout is proportional rather than absolute the panels keep their places.
//
// Desktop returns 1.0 deliberately.  The desktop window's framebuffer and its
// logical size already agree on the machines this is developed on, and #25 asks
// for desktop rendering to be untouched at 1440x900.
static float uiScale() {
#ifdef __EMSCRIPTEN__
    const double dpr = emscripten_get_device_pixel_ratio();
    // A dpr below 1 would shrink text that is already the right size, and an
    // absurd one would blow the font atlas up for no gain.  Neither is a real
    // device; both are cheap to refuse.
    if (!(dpr > 0.0)) return 1.0f;
    return static_cast<float>(std::clamp(dpr, 1.0, 4.0));
#else
    return 1.0f;
#endif
}

// One dockspace for both panel sets.  Built once, and only when ImGui has no
// layout of its own — a saved imgui.ini always wins, so a user's arrangement
// survives every restart.
static void buildDefaultLayout(ImGuiID dock_id, const ImVec2& size) {
    ImGui::DockBuilderRemoveNode(dock_id);
    ImGui::DockBuilderAddNode(dock_id, ImGuiDockNodeFlags_DockSpace |
                                       ImGuiDockNodeFlags_PassthruCentralNode);
    ImGui::DockBuilderSetNodeSize(dock_id, size);

    // Ratios read off a layout arrived at by using the thing, rather than
    // guessed: 0.22 / 0.46 / 0.34 were already what a session settled on, and
    // the two that were not are corrected here.
    //
    // The central node is deliberately left empty: it is where the 3D plate
    // shows through the passthru dockspace.
    ImGuiID centre = dock_id;
    ImGuiID left   = ImGui::DockBuilderSplitNode(centre, ImGuiDir_Left,  0.22f, nullptr, &centre);
    ImGuiID right  = ImGui::DockBuilderSplitNode(centre, ImGuiDir_Right, 0.46f, nullptr, &centre);
    ImGuiID bottom = ImGui::DockBuilderSplitNode(centre, ImGuiDir_Down,  0.34f, nullptr, &centre);
    // System Properties is a readout, not a workspace: it needs the height of
    // its own text and no more.  At 0.34 it took a third of the left column
    // away from the panel that actually gets used.
    ImGuiID left_bottom =
        ImGui::DockBuilderSplitNode(left, ImGuiDir_Down, 0.125f, nullptr, &left);
    // Plate Plots and Model Comparison side by side rather than tabbed.  They
    // are meant to be read against each other — the whole point of the
    // comparison panel is the divergence between the linear prediction and the
    // nonlinear ball — and a tab makes that a memory test.
    ImGuiID bottom_right =
        ImGui::DockBuilderSplitNode(bottom, ImGuiDir_Right, 0.5f, nullptr, &bottom);

    ImGui::DockBuilderDockWindow("Model Configuration", left);
    ImGui::DockBuilderDockWindow("Plate Control", left);
    ImGui::DockBuilderDockWindow("System Properties", left_bottom);

    ImGui::DockBuilderDockWindow("Bode Plot", right);
    ImGui::DockBuilderDockWindow("Nyquist Plot", right);
    ImGui::DockBuilderDockWindow("Pole-Zero / Root Locus", right);
    ImGui::DockBuilderDockWindow("Time Response", right);

    ImGui::DockBuilderDockWindow("Plate Plots", bottom);
    ImGui::DockBuilderDockWindow("Model Comparison", bottom_right);

    ImGui::DockBuilderFinish(dock_id);
}

static void render_frame(void* arg) {
    auto* ctx = static_cast<FrameContext*>(arg);
    auto& state = *ctx->state;
    auto& presets = *ctx->presets;
    GLFWwindow* window = ctx->window;

    glfwPollEvents();

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    // After NewFrame, which is what computes io.WantCaptureMouse: the plate
    // arbitrates its orbit drag on that flag, and stepping earlier read the
    // previous frame's answer.
    //
    // Fixed step, as the plate view has always used: vsync on the desktop and
    // requestAnimationFrame on the web both hold the frame rate at 60 Hz, and a
    // fixed step keeps the ball's trajectory reproducible.
    //
    // The design is handed over before the step, so the plate drives on the
    // gain the previous frame's recompute produced.  One frame of staleness at
    // 60 Hz against a servo lag of 0.05 s; the alternative is reordering the
    // whole frame around a dependency the loop does not have.
    handDesignToPlate(state, presets, *ctx->plate);
    ctx->plate->step(window, 1.0f / 60.0f);

    // --- Full-viewport dockspace ---
    ImGuiViewport* vp = ImGui::GetMainViewport();
    const ImGuiID dock_id = ImGui::GetID("MainDockSpace");
    static bool layout_checked = false;
    static int focus_defaults_frames_left = 0;
    static int place_toggles_frames_left = 0;
    if (!layout_checked) {
        layout_checked = true;
        if (ImGui::DockBuilderGetNode(dock_id) == nullptr) {
            buildDefaultLayout(dock_id, vp->WorkSize);
            // DockBuilder decides which tab of a fresh node is on top, and it
            // does not pick this one.  Windows have to exist before they can be
            // focused, so the pick is made from the frames that follow.
            focus_defaults_frames_left = 2;
            // The central node has no geometry until DockSpace has laid the
            // fresh layout out, so the toggle bar is placed over the first few
            // frames rather than pinned to a position that is still (0, 0).
            place_toggles_frames_left = 4;
        }
    }

    ImGui::SetNextWindowPos(vp->WorkPos);
    ImGui::SetNextWindowSize(vp->WorkSize);
    ImGui::SetNextWindowViewport(vp->ID);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    ImGui::Begin("DockSpace", nullptr,
        ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoTitleBar |
        ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize |
        ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoBringToFrontOnFocus |
        ImGuiWindowFlags_NoNavFocus | ImGuiWindowFlags_NoBackground);
    ImGui::PopStyleVar(3);
    ImGui::DockSpace(dock_id, ImVec2(0, 0),
                     ImGuiDockNodeFlags_PassthruCentralNode);
    ImGui::End();

    // --- Recompute analysis if needed ---
    if (state.needs_recompute) {
        state.needs_recompute = false;

        // Plant dimensions changed: reseed the loop list wholesale.
        //
        // Reseed, not clamp and not partial-drop.  Clamping is the dangerous
        // option — a loop clamped from y1 to y0 controls a different physical
        // quantity while still carrying gains tuned for the old one, with
        // nothing on screen saying so.  Gains are lost, but a Kp tuned for a
        // quarter-car is meaningless on a pendulum.
        //
        // Keyed on dimensions, not plant equality, so the physical-parameter
        // sliders never disturb tuning mid-sweep.  Located here rather than in
        // the preset handler because FOUR sites mutate state.plant — the preset
        // combo, Apply Plant, the physical-param sliders and applyTFToSS — and
        // only a central guard covers all of them.  See issue #2.
        if (state.plant.outputs() != state.cached_p ||
            state.plant.inputs() != state.cached_m) {
            state.cached_p = state.plant.outputs();
            state.cached_m = state.plant.inputs();
            caliburn::seedIdentityLoops(state);
            state.output_i = 0;
            state.input_j = 0;
            state.ref_r = 0;
        }

        // LQR weights follow the STATE count, which the (p, m) guard above
        // does not track — a plant can keep its shape and change order.
        if (state.lqr_q.size() != state.plant.states())
            state.lqr_q = caliburn::defaultLqrStateWeights(state.plant.states());
        if (state.lqr_r.size() != state.plant.inputs())
            state.lqr_r = caliburn::defaultLqrInputWeights(state.plant.inputs());

        state.systems[0] = state.plant;
        state.system_valid[0] = true;

        const bool loop_backed =
            state.ctrl_type == caliburn::ControllerType::PID ||
            state.ctrl_type == caliburn::ControllerType::LeadLag;
        const caliburn::LoopKind loop_kind =
            (state.ctrl_type == caliburn::ControllerType::LeadLag)
                ? caliburn::LoopKind::LeadLag
                : caliburn::LoopKind::PID;

        if (loop_backed && !state.loops.empty()) {
            // The controller is m-out/p-in by construction, so seriesConnect's
            // precondition holds and the open loop is p x p — feedbackConnect
            // needs no shape test.  systems[1] is written directly rather than
            // through state.ctrl_ss, so switching to PID does not clobber a
            // hand-entered State-Space controller.
            state.systems[1] = caliburn::buildLoopController(
                state.loops, loop_kind, state.plant);
            state.system_valid[1] = true;
            state.systems[2] =
                caliburn::seriesConnect(state.systems[1], state.plant);
            state.system_valid[2] = true;
            state.systems[3] = caliburn::feedbackConnect(state.systems[2]);
            state.system_valid[3] = true;
        } else if (state.ctrl_type == caliburn::ControllerType::StateSpace &&
                   state.ctrl_ss.D.size() > 0) {
            // Was `ctrl_ss.A.size() > 0`, which conflated "is a controller
            // configured" with "is it safe to analyze".  The n == 0 guards make
            // the second question moot, and a 0-state controller is a
            // legitimate controller — D is the matrix that is always present
            // once Apply Controller has run.  See issue #5.
            state.systems[1] = state.ctrl_ss;
            state.system_valid[1] = true;
            if (state.ctrl_ss.outputs() == state.plant.inputs()) {
                state.systems[2] = caliburn::seriesConnect(state.ctrl_ss, state.plant);
                state.system_valid[2] = true;
                if (state.systems[2].outputs() == state.systems[2].inputs()) {
                    state.systems[3] = caliburn::feedbackConnect(state.systems[2]);
                    state.system_valid[3] = true;
                } else {
                    state.system_valid[3] = false;
                }
            } else {
                state.system_valid[2] = false;
                state.system_valid[3] = false;
            }
        } else if (state.ctrl_type == caliburn::ControllerType::LQR) {
            // Solved every recompute, so the gain tracks the weight sliders
            // and the physical-parameter sliders without a Design button.
            // The plant is re-linearised upstream of here, so a change to the
            // leg geometry moves K as well as the poles.
            state.lqr_result = caliburn::computeLQR(
                state.plant,
                state.lqr_q.asDiagonal().toDenseMatrix(),
                state.lqr_r.asDiagonal().toDenseMatrix());
            state.system_valid[1] = false;
            state.system_valid[2] = false;
            if (state.lqr_result.success) {
                state.ctrl_K = state.lqr_result.K;
                state.systems[3] =
                    caliburn::stateFeedbackClose(state.plant, state.ctrl_K);
                state.system_valid[3] = true;
            } else {
                state.system_valid[3] = false;
            }
        } else if (state.ctrl_type == caliburn::ControllerType::GainMatrix &&
                   state.ctrl_K.size() > 0) {
            state.system_valid[1] = false;
            state.system_valid[2] = false;
            state.systems[3] = caliburn::stateFeedbackClose(state.plant, state.ctrl_K);
            state.system_valid[3] = true;
        } else {
            // Includes a loop-backed type with an empty loop list, which is
            // equivalent to ControllerType::None (issue #2).  An all-zero m x p
            // controller connects cleanly but gives an open loop of -inf dB at
            // every frequency, which is not worth handing to Bode and Nyquist.
            state.system_valid[1] = false;
            state.system_valid[2] = false;
            state.system_valid[3] = false;
        }

        // Rule C.  output_i in [0,p), input_j in [0,m), ref_r in [0,p).
        //   plant           p x m -> (i, j)
        //   controller      m x p -> (j, i)   transposed
        //   open/closed     p x p -> (i, r)
        // The output index is genuinely shared: plant output i, controller
        // input i and loop-system output i are the same physical signal — so
        // this is three indices, not six, and every system is in range by
        // construction.  See issue #9.
        auto channelFor = [&](int s, int& out, int& in) {
            switch (s) {
                case 1:  out = state.input_j;  in = state.output_i; break;
                case 2:
                case 3:  out = state.output_i; in = state.ref_r;    break;
                default: out = state.output_i; in = state.input_j;  break;
            }
        };

        // Liveness is read off the loop list, never by inspecting matrices.
        auto loopOn = [&](int out, int in) {
            for (const auto& l : state.loops)
                if (l.out == out && l.in == in) return true;
            return false;
        };
        auto refLive = [&](int r) {
            for (const auto& l : state.loops)
                if (l.out == r) return true;
            return false;
        };

        for (int s = 0; s < caliburn::NUM_SYSTEMS; ++s) {
            state.channel_live[s] = true;
            state.channel_reason[s].clear();
        }
        if (loop_backed && !state.loops.empty()) {
            char reason[96];
            if (!loopOn(state.output_i, state.input_j)) {
                state.channel_live[1] = false;
                std::snprintf(reason, sizeof(reason),
                              "no loop pairs y%d -> u%d",
                              state.output_i, state.input_j);
                state.channel_reason[1] = reason;
            }
            if (!refLive(state.ref_r)) {
                state.channel_live[2] = false;
                state.channel_live[3] = false;
                std::snprintf(reason, sizeof(reason),
                              "reference r%d drives no loop", state.ref_r);
                state.channel_reason[2] = reason;
                state.channel_reason[3] = reason;
            }
        }

        state.diagnostics = caliburn::computeLoopDiagnostics(
            state.plant, state.loops, state.output_scales,
            state.diag_omega_hz,
            state.freq_min_hz, state.freq_max_hz, state.num_freq_points);

        for (int s = 0; s < caliburn::NUM_SYSTEMS; ++s) {
            // s == 2 (open-loop combined) is computed for feedbackConnect only.
            if (!state.system_valid[s] || s == 2 || !state.channel_live[s]) {
                state.bode[s] = {};
                state.pole_zero[s] = {};
                continue;
            }
            int out = 0, in = 0;
            channelFor(s, out, in);
            state.bode[s] = caliburn::computeBode(
                state.systems[s], out, in,
                state.freq_min_hz, state.freq_max_hz, state.num_freq_points);
            // The controller's pole plot stays unfiltered: poles are eig(A) for
            // the whole system regardless of channel, only zeros are
            // channel-specific.  That is the correct convention and what the
            // plant already does.
            state.pole_zero[s] = caliburn::computePoleZero(
                state.systems[s], out, in);
        }

        if (state.show_all_channels) {
            const int p = state.plant.outputs();
            const int m = state.plant.inputs();
            // Two grids, by role.  Rule C applied to the grid, so the grid
            // needs no rule of its own and each table is indexed by the shape
            // it actually holds — the plant-dimensioned read in the panel
            // cannot recur.
            //   open grid  p x m : plant at (i,j), controller at (j,i)
            //   loop grid  p x p : closed loop at (i,r)
            for (int s = 0; s < caliburn::NUM_SYSTEMS; ++s) {
                if (!state.system_valid[s] || s == 2) {
                    state.bode_grid[s].clear();
                    continue;
                }
                const int rows = p;
                const int cols = (s == 3) ? p : m;
                state.bode_grid[s].assign(
                    static_cast<std::size_t>(rows * cols), {});
                for (int i = 0; i < rows; ++i) {
                    for (int j = 0; j < cols; ++j) {
                        // Cell (i,j) of the open grid means plant channel (i,j)
                        // and controller channel (j,i) — the same physical
                        // path, read from each end.
                        int out = i, in = j;
                        if (s == 1) { out = j; in = i; }
                        state.bode_grid[s][i * cols + j] = caliburn::computeBode(
                            state.systems[s], out, in,
                            state.freq_min_hz, state.freq_max_hz,
                            state.num_freq_points);
                    }
                }
            }
        }

        for (int s : {0, 3}) {
            if (!state.system_valid[s]) continue;
            // Plant is p x m, closed loop is p x p — the same shape mismatch,
            // in the time domain.  No parallel time-domain index pair.
            const int u_idx = (s == 0) ? state.input_j : state.ref_r;
            switch (state.input_type) {
                case caliburn::InputType::Step:
                    state.time_resp[s] = caliburn::computeStepResponse(
                        state.systems[s], u_idx,
                        state.amplitude, state.duration, state.dt_sim);
                    break;
                case caliburn::InputType::Impulse:
                    state.time_resp[s] = caliburn::computeImpulseResponse(
                        state.systems[s], u_idx,
                        state.amplitude, state.duration, state.dt_sim);
                    break;
                case caliburn::InputType::Ramp:
                    state.time_resp[s] = caliburn::computeRampResponse(
                        state.systems[s], u_idx,
                        state.slope, state.duration, state.dt_sim);
                    break;
            }
        }

        state.controllability = caliburn::checkControllability(state.plant);
        state.observability = caliburn::checkObservability(state.plant);
        if (state.system_valid[3]) {
            state.cl_controllability =
                caliburn::checkControllability(state.systems[3]);
            state.cl_observability =
                caliburn::checkObservability(state.systems[3]);
        }

        // Clear on EVERY non-computing path.  Before this, when pz_mode was
        // StateFB and ctrl_type was not GainMatrix, neither branch ran and
        // root_locus was never cleared — so the panel kept drawing the PREVIOUS
        // mode's locus labelled as State FB, with its marker read off
        // rl_current_alpha.  Clearing here kills that by construction rather
        // than by a guard.  See issue #10.
        state.root_locus.clear();
        state.rl_loop_gain_margin = -1.0;
        switch (state.pz_mode) {
            case caliburn::PZMode::PlantLocus:
                state.root_locus = caliburn::computeRootLocus(
                    state.plant, state.output_i, state.input_j,
                    state.rl_k_min, state.rl_k_max, state.rl_num_points);
                break;
            case caliburn::PZMode::LoopLocus:
                if (loop_backed && !state.loops.empty()) {
                    state.root_locus = caliburn::computeLoopLocus(
                        state.plant, state.loops, loop_kind,
                        state.output_i, state.input_j,
                        state.rl_kappa_min, state.rl_kappa_max,
                        state.rl_num_points);
                    state.rl_loop_gain_margin =
                        caliburn::loopGainMargin(state.root_locus);
                }
                break;
            case caliburn::PZMode::StateFB:
                if ((state.ctrl_type == caliburn::ControllerType::GainMatrix ||
                     state.ctrl_type == caliburn::ControllerType::LQR) &&
                    state.ctrl_K.size() > 0) {
                    state.root_locus = caliburn::computeStateFeedbackLocus(
                        state.plant, state.ctrl_K,
                        state.rl_alpha_min, state.rl_alpha_max,
                        state.rl_num_points);
                }
                break;
            case caliburn::PZMode::PoleZero:
                break;
        }
    }

    // --- Panel toggle bar ---
    // Floating, and parked in the corner of the central node: anywhere else and
    // its default position lands on top of the model panel.
    if (place_toggles_frames_left > 0) {
        const ImGuiDockNode* centre = ImGui::DockBuilderGetCentralNode(dock_id);
        if (centre && centre->Size.x > 0.0f) {
            --place_toggles_frames_left;
            ImGui::SetNextWindowPos(
                ImVec2(centre->Pos.x + 12.0f, centre->Pos.y + 12.0f),
                ImGuiCond_Always);
        }
    }
    ImGui::Begin("##toggles", nullptr,
                 ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoScrollbar |
                 ImGuiWindowFlags_AlwaysAutoResize);
    auto toggleBtn = [](const char* label, bool& flag) {
        bool active = flag;
        if (active) {
            ImGui::PushStyleColor(ImGuiCol_Button,
                                  ImVec4(0.2f, 0.5f, 0.8f, 1.0f));
        }
        if (ImGui::Button(label)) flag = !flag;
        if (active) ImGui::PopStyleColor();
        ImGui::SameLine();
    };
    toggleBtn("Pole-Zero", state.show_pole_zero);
    toggleBtn("Bode", state.show_bode);
    toggleBtn("Nyquist", state.show_nyquist);
    toggleBtn("Time Resp", state.show_time_response);
    ImGui::End();

    // --- Draw panels ---
    caliburn::drawModelPanel(state, presets);
    caliburn::drawPropertiesPanel(state);
    caliburn::drawPoleZeroPanel(state);
    caliburn::drawBodePanel(state);
    caliburn::drawNyquistPanel(state);
    caliburn::drawTimeResponsePanel(state);
    ctx->plate->drawPanels();

    if (focus_defaults_frames_left > 0) {
        --focus_defaults_frames_left;
        ImGui::SetWindowFocus("Plate Plots");
    }

    // --- Render ---
    ImGui::Render();
    int fb_w, fb_h;
    glfwGetFramebufferSize(window, &fb_w, &fb_h);
    glViewport(0, 0, fb_w, fb_h);
    glClearColor(0.08f, 0.08f, 0.10f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    // The 3D scene is drawn into the central dock node rather than the whole
    // framebuffer, so the plate stays centred in the space the panels leave
    // rather than hiding behind them.
    if (const ImGuiDockNode* centre = ImGui::DockBuilderGetCentralNode(dock_id)) {
        const float sx = fb_w / std::max(vp->Size.x, 1.0f);
        const float sy = fb_h / std::max(vp->Size.y, 1.0f);
        const int cx = static_cast<int>((centre->Pos.x - vp->Pos.x) * sx);
        const int cy = static_cast<int>(
            (vp->Pos.y + vp->Size.y - (centre->Pos.y + centre->Size.y)) * sy);
        const int cw = static_cast<int>(centre->Size.x * sx);
        const int ch = static_cast<int>(centre->Size.y * sy);
        if (cw > 0 && ch > 0) {
            glViewport(cx, cy, cw, ch);
            glClear(GL_DEPTH_BUFFER_BIT);
            glEnable(GL_DEPTH_TEST);
            ctx->plate->drawScene(static_cast<float>(cw) / static_cast<float>(ch));
            glDisable(GL_DEPTH_TEST);
            glViewport(0, 0, fb_w, fb_h);
        }
    }

    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    glfwSwapBuffers(window);
}

int main() {
    // --- Init GLFW ---
    glfwSetErrorCallback([](int err, const char* desc) {
        std::fprintf(stderr, "GLFW error %d: %s\n", err, desc);
    });
    if (!glfwInit()) return 1;

#ifdef __EMSCRIPTEN__
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_ES_API);
#else
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);
#endif

    GLFWwindow* window = glfwCreateWindow(
        1600, 1000, "Ball-Balancer \xe2\x80\x94 Caliburn",
        nullptr, nullptr);
    if (!window) { glfwTerminate(); return 1; }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

#ifndef __EMSCRIPTEN__
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::fprintf(stderr, "Failed to initialize GLAD\n");
        return 1;
    }
#endif

    // --- Init ImGui + ImPlot ---
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();

    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
    io.IniFilename = "imgui.ini";

    // --- Load Unicode font ---
    ImFontConfig font_cfg;
    font_cfg.OversampleH = 2;
    font_cfg.OversampleV = 2;

    // Only ranges the bundled NotoSans subset actually covers.  A range the
    // font lacks costs nothing but buys nothing either: the glyph still draws
    // as a hollow substitution box.  Measured against this exact TTF:
    //
    //   General Punctuation  111/112   em dash, bullet          <- was missing
    //   Phonetic Extensions  128/128   the subscript letters    <- was missing
    //   Arrows                 0/112   nothing at all           <- was listed
    //   Mathematical Ops       1/256   U+2212 MINUS, and no more
    //
    // So the arrow range is dropped, and anything from Mathematical Operators
    // — infinity, angle, approximately-equal — has to be written in ASCII no
    // matter what this list says.
    static const ImWchar glyph_ranges[] = {
        0x0020, 0x00FF,  // Basic Latin + Latin Supplement
        0x0370, 0x03FF,  // Greek and Coptic
        0x1D00, 0x1D7F,  // Phonetic Extensions — subscript letters
        0x2000, 0x206F,  // General Punctuation — em dash, bullet
        0x2070, 0x209F,  // Superscripts and Subscripts
        0x2212, 0x2212,  // MINUS SIGN, the one Mathematical Operator present
        0,
    };

    // Baked once at startup rather than tracked across resizes: rebuilding the
    // font atlas costs a texture upload, and a phone does not change its pixel
    // ratio mid-session.  A desktop browser dragged between monitors of
    // different density is the one case this misses, and it costs a reload.
    const float ui_scale = uiScale();

    const char* font_path = "vendor/fonts/NotoSans-Regular.ttf";
    if (FILE* f = std::fopen(font_path, "rb")) {
        std::fclose(f);
        io.Fonts->AddFontFromFileTTF(font_path, 16.0f * ui_scale, &font_cfg, glyph_ranges);
    } else {
        std::fprintf(stderr, "Warning: %s not found, using default font\n", font_path);
        io.Fonts->AddFontDefault();
        // The built-in font is a fixed-size bitmap, so it cannot be baked
        // larger.  Stretching it is ugly but legible, which beats correct and
        // unreadable on the one path where the bundled TTF is missing.
        io.FontGlobalScale = ui_scale;
    }

    ImGui::StyleColorsDark();
    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding = 6.0f;
    style.FrameRounding = 4.0f;
    style.GrabRounding = 4.0f;
    // After the rounding values, so they scale with everything else: padding,
    // scrollbar and grab sizes, and the dock separator you have to hit with a
    // finger rather than a mouse.
    style.ScaleAllSizes(ui_scale);

    // ImPlot keeps its own metrics, and ScaleAllSizes does not reach them. Left
    // unscaled, axis labels sit hard against plot edges and the tick text runs
    // into the frame once the font grows.
    if (ui_scale != 1.0f) {
        ImPlotStyle& plot = ImPlot::GetStyle();
        auto scale2 = [ui_scale](ImVec2 v) { return ImVec2(v.x * ui_scale, v.y * ui_scale); };

        // Padding: without this the tick labels grow into the plot frame.
        plot.PlotPadding        = scale2(plot.PlotPadding);
        plot.LabelPadding       = scale2(plot.LabelPadding);
        plot.LegendPadding      = scale2(plot.LegendPadding);
        plot.LegendInnerPadding = scale2(plot.LegendInnerPadding);
        plot.LegendSpacing      = scale2(plot.LegendSpacing);
        plot.MousePosPadding    = scale2(plot.MousePosPadding);
        plot.AnnotationPadding  = scale2(plot.AnnotationPadding);

        // Ticks and rules are specified in pixels, so at dpr 3 an unscaled
        // 1 px grid line is a third of a pixel on the glass — present in the
        // buffer, invisible on the screen.
        plot.MajorTickLen  = scale2(plot.MajorTickLen);
        plot.MinorTickLen  = scale2(plot.MinorTickLen);
        plot.MajorTickSize = scale2(plot.MajorTickSize);
        plot.MinorTickSize = scale2(plot.MinorTickSize);
        plot.MajorGridSize = scale2(plot.MajorGridSize);
        plot.MinorGridSize = scale2(plot.MinorGridSize);
        plot.PlotBorderSize *= ui_scale;

        // Minimum and default plot sizes are a floor in pixels. Scaling them
        // keeps the floor at a constant *logical* size, so a plot cannot
        // collapse to something smaller than it would be on a desktop.
        plot.PlotMinSize     = scale2(plot.PlotMinSize);
        plot.PlotDefaultSize = scale2(plot.PlotDefaultSize);
    }

    // Constructed before the ImGui backend so that its scroll callback is the
    // *previous* one, which the backend chains to instead of replacing.
    //
    // static, like everything else the frame callback reaches through: see the
    // main loop below.
    static caliburn::PlateView plate;
    plate.attach(window);

    ImGui_ImplGlfw_InitForOpenGL(window, true);
#ifdef __EMSCRIPTEN__
    ImGui_ImplOpenGL3_Init("#version 300 es");
#else
    ImGui_ImplOpenGL3_Init("#version 330");
#endif

    // --- Init app state ---
    static caliburn::AppState state;
    static auto presets = caliburn::getBuiltinModels();
    state.preset_index = caliburn::defaultModelIndex(presets);
    state.plant = presets[state.preset_index].system;
    state.current_params = presets[state.preset_index].params;
    caliburn::matrixToTextBuf(state.plant.A, state.A_text, sizeof(state.A_text));
    caliburn::matrixToTextBuf(state.plant.B, state.B_text, sizeof(state.B_text));
    caliburn::matrixToTextBuf(state.plant.C, state.C_text, sizeof(state.C_text));
    caliburn::matrixToTextBuf(state.plant.D, state.D_text, sizeof(state.D_text));
    caliburn::extractTFFromSS(state);

    // The application opens with a controller already designed, because the
    // demo opens with one already running: a visitor gets about ten seconds
    // and will not go looking for a combo box to make the page do something.
    // The default preset is the cascade, so LQR here means the plate is handed
    // a usable gain on the second frame and the opening closes the loop —
    // see PlateView's `auto_engaged_` and issue #17.
    //
    // The closed-loop trace comes on with it.  Under LQR the plant is the only
    // trace enabled by default, so the pole-zero map would open showing the
    // OPEN-loop poles of a page whose whole claim is that the loop is closed.
    // Written as the model panel's own rule rather than as a bare `true`: the
    // combo applies exactly this line when a visitor changes the type, and
    // selecting LQR here has to mean what selecting it there means.
    state.ctrl_type = caliburn::ControllerType::LQR;
    state.trace_visible[3] = state.ctrl_type != caliburn::ControllerType::None;

#ifndef __EMSCRIPTEN__
    // Neither enum exists in ES 3.0; WebGL2 multisamples the default
    // framebuffer on its own and has no line-smoothing knob at all.
    glEnable(GL_MULTISAMPLE);
    glEnable(GL_LINE_SMOOTH);
#endif
    plate.initGL();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // --- Main loop ---
    //
    // static, not automatic.  emscripten_set_main_loop_arg with
    // simulate_infinite_loop escapes main by throwing, and every frame after
    // that reads the app through this pointer — so ctx and everything it points
    // at have to outlive main's frame.  Whether the runtime leaves the shadow
    // stack pointer where main left it is an implementation detail of the
    // unwind, and one static keyword is cheaper than depending on the answer.
    static FrameContext ctx{window, &state, &presets, &plate};
#ifdef __EMSCRIPTEN__
    emscripten_set_main_loop_arg(render_frame, &ctx, 0, true);
#else
    while (!glfwWindowShouldClose(window)) {
        render_frame(&ctx);
    }
#endif

    // --- Cleanup ---
    // Unreachable on the web rather than excluded from it: the main loop is
    // installed with simulate_infinite_loop, so main never returns there and
    // the browser tears the context down for us.
    //
    // Before the context goes away: ~LineRenderer calls glDeleteBuffers, and
    // running that after glfwTerminate is undefined.
    plate.shutdownGL();
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
