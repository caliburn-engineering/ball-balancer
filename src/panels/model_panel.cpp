// src/panels/model_panel.cpp
#include "model_panel.h"
#include "../analysis/loop_diagnostics.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

namespace caliburn {

void extractTFFromSS(AppState& state) {
    const auto& plant = state.plant;
    state.plant_is_1st_order = false;
    state.plant_is_2nd_order = false;

    if (plant.states() == 1 && plant.inputs() == 1 && plant.outputs() == 1) {
        double a = -plant.A(0, 0);  // a = 1/tau
        if (a > 1e-10) {
            state.plant_is_1st_order = true;
            state.tf_tau = static_cast<float>(1.0 / a);
            state.tf_K = static_cast<float>(plant.C(0, 0) * plant.B(0, 0) / a);
        }
    } else if (plant.states() == 2 && plant.inputs() == 1 && plant.outputs() == 1) {
        double neg_wn2 = plant.A(1, 0);  // = -wn^2
        if (neg_wn2 < -1e-10) {
            double wn = std::sqrt(-neg_wn2);
            double zeta = -plant.A(1, 1) / (2.0 * wn);
            double K = plant.C(0, 0) / (wn * wn);
            if (wn > 0 && zeta >= 0) {
                state.plant_is_2nd_order = true;
                state.tf_wn = static_cast<float>(wn);
                state.tf_zeta = static_cast<float>(zeta);
                state.tf_K = static_cast<float>(K);
            }
        }
    }
}

void applyTFToSS(AppState& state) {
    if (state.plant_is_1st_order) {
        double tau = state.tf_tau;
        double K = state.tf_K;
        state.plant.A = (Eigen::MatrixXd(1, 1) << -1.0 / tau).finished();
        state.plant.B = (Eigen::MatrixXd(1, 1) << 1.0 / tau).finished();
        state.plant.C = (Eigen::MatrixXd(1, 1) << K).finished();
        state.plant.D = Eigen::MatrixXd::Zero(1, 1);
    } else if (state.plant_is_2nd_order) {
        double wn = state.tf_wn;
        double zeta = state.tf_zeta;
        double K = state.tf_K;
        state.plant.A = (Eigen::MatrixXd(2, 2) << 0, 1, -wn * wn, -2.0 * zeta * wn).finished();
        state.plant.B = (Eigen::MatrixXd(2, 1) << 0, 1).finished();
        state.plant.C = (Eigen::MatrixXd(1, 2) << K * wn * wn, 0).finished();
        state.plant.D = Eigen::MatrixXd::Zero(1, 1);
    }
    matrixToTextBuf(state.plant.A, state.A_text, sizeof(state.A_text));
    matrixToTextBuf(state.plant.B, state.B_text, sizeof(state.B_text));
    matrixToTextBuf(state.plant.C, state.C_text, sizeof(state.C_text));
    matrixToTextBuf(state.plant.D, state.D_text, sizeof(state.D_text));
}

namespace {

bool drawMatrixSliders(const char* name, Eigen::MatrixXd& mat,
                       char* text_buf, size_t buf_size) {
    bool changed = false;
    int rows = static_cast<int>(mat.rows());
    int cols = static_cast<int>(mat.cols());

    char table_id[64];
    std::snprintf(table_id, sizeof(table_id), "##sliders_%s", name);

    if (ImGui::BeginTable(table_id, cols)) {
        for (int i = 0; i < rows; ++i) {
            ImGui::TableNextRow();
            for (int j = 0; j < cols; ++j) {
                ImGui::TableNextColumn();
                char slider_id[64];
                std::snprintf(slider_id, sizeof(slider_id), "##%s_%d_%d", name, i, j);
                float val = static_cast<float>(mat(i, j));
                float absv = std::abs(val);
                float range = std::max(10.0f, absv * 3.0f);
                ImGui::SetNextItemWidth(-FLT_MIN);
                if (ImGui::SliderFloat(slider_id, &val, -range, range, "%.3f")) {
                    mat(i, j) = static_cast<double>(val);
                    matrixToTextBuf(mat, text_buf, buf_size);
                    changed = true;
                }
            }
        }
        ImGui::EndTable();
    }
    return changed;
}


// Ki and Kd must reach exactly 0 and still resolve small values, and the spec
// asks for "linear 0-1, then log 1-hi".  ImGuiSliderFlags_Logarithmic DOES work
// with v_min = 0 — verified in the prototype's bench, no assert, no crash,
// ImGui's zero epsilon handles it — but it is rejected on emphasis: a value of
// 0.5 sits at ~62% of travel, compressing the 1-hi decades into the last third,
// the inverse of the split the spec asks for.  This puts 0.5 at exactly 25%.
bool gainSliderPiecewise(const char* id, float* v, float hi) {
    const float lin_hi = 1.0f;
    float t = (*v <= lin_hi)
                  ? 0.5f * (*v / lin_hi)
                  : 0.5f + 0.5f * std::log(*v / lin_hi) / std::log(hi / lin_hi);
    char fmt[32];
    std::snprintf(fmt, sizeof(fmt), "%.4g", *v);  // display the true value
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (ImGui::SliderFloat(id, &t, 0.0f, 1.0f, fmt)) {
        *v = (t <= 0.5f) ? 2.0f * t * lin_hi
                         : lin_hi * std::pow(hi / lin_hi, (t - 0.5f) * 2.0f);
        return true;
    }
    return false;
}

// Severity is graded from the RGA number, never from |lambda_ii|, and applies
// to RGA cells only — Channel Share is never coloured.
ImVec4 severityFill(double rga_number) {
    if (rga_number < 0.5) return ImVec4(0.12f, 0.32f, 0.16f, 1.0f);
    if (rga_number < 2.0) return ImVec4(0.36f, 0.30f, 0.10f, 1.0f);
    return ImVec4(0.40f, 0.14f, 0.14f, 1.0f);
}

}  // anonymous namespace

// The controller section for a loop-backed controller: the p x m matrix IS the
// editor.  Click a cell to pair or unpair; the same widget that grades the
// pairing is the one that creates it.
//
// 354 px — the docked panel's width in the checked-in imgui.ini — is the
// binding constraint, and it is what eliminated the ledger layout the spec
// sketched.  There are deliberately NO add / remove / reorder controls: cells
// are the affordance and grid order is loop order.  One consequence is that an
// exact duplicate pairing is unreachable here, so issue #2's duplicate caution
// needs no UI surface — the builder still sums duplicates, which is correct for
// a loop list built programmatically.  All fan-out shapes stay expressible:
// many-to-one is several cells in one row, one-to-many several in one column.
static void drawPairingGrid(AppState& state) {
    const int gp = state.plant.outputs();
    const int gm = state.plant.inputs();
    const auto& diag = state.diagnostics;

    ImGui::TextDisabled("left-click to pair / select \xe2\x80\xa2 right-click to unpair");

    // Columns: per-output scale, row label, then one cell per plant input.
    // The scale column supersedes a collapsed "scales: default (1, 1)" status
    // line — the values are on screen, so nothing has to state whether they are
    // still at default, and the row header is where per-output metadata would
    // naturally go if that fog ever clears.
    if (ImGui::BeginTable("##pairgrid", gm + 2,
                          ImGuiTableFlags_Borders |
                          ImGuiTableFlags_SizingFixedFit)) {
        ImGui::TableSetupColumn("scale", ImGuiTableColumnFlags_WidthFixed, 56);
        ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed, 24);
        for (int j = 0; j < gm; ++j) {
            char h[8];
            std::snprintf(h, sizeof(h), "u%d", j);
            ImGui::TableSetupColumn(h, ImGuiTableColumnFlags_WidthFixed, 74);
        }
        ImGui::TableHeadersRow();

        for (int i = 0; i < gp; ++i) {
            ImGui::PushID(i);
            ImGui::TableNextRow();

            ImGui::TableNextColumn();
            if (i < static_cast<int>(state.output_scales.size())) {
                ImGui::SetNextItemWidth(-FLT_MIN);
                if (ImGui::DragFloat("##sc", &state.output_scales[i], 0.01f,
                                     0.001f, 1000.0f, "%.3g"))
                    state.needs_recompute = true;
            }

            ImGui::TableNextColumn();
            ImGui::Text("y%d", i);

            for (int j = 0; j < gm; ++j) {
                ImGui::TableNextColumn();
                ImGui::PushID(j);

                int found = -1;
                for (std::size_t k = 0; k < state.loops.size(); ++k)
                    if (state.loops[k].out == i && state.loops[k].in == j) {
                        found = static_cast<int>(k);
                        break;
                    }
                const bool paired = found >= 0;

                // Two lines: the diagnostic over |g_ij| in dB, dimmed.  The dB
                // line is REQUIRED under Channel Share, which cannot tell
                // "drives both strongly" from "drives neither"; under the RGA
                // it is free and still worth reading.
                char top[24] = "-";
                if (diag.has_rga) {
                    const auto ro = std::find(diag.sub_out.begin(),
                                              diag.sub_out.end(), i);
                    const auto co = std::find(diag.sub_in.begin(),
                                              diag.sub_in.end(), j);
                    if (ro != diag.sub_out.end() && co != diag.sub_in.end())
                        std::snprintf(top, sizeof(top), "%.2f",
                            diag.lambda(ro - diag.sub_out.begin(),
                                        co - diag.sub_in.begin()));
                } else if (diag.has_share && j == diag.share_in &&
                           i < static_cast<int>(diag.share.size())) {
                    std::snprintf(top, sizeof(top), "%.0f%%",
                                  100.0 * diag.share[i]);
                }

                char cell[64];
                if (i < static_cast<int>(diag.mag_db.size()) &&
                    diag.has_share && j == diag.share_in) {
                    std::snprintf(cell, sizeof(cell), "%s\n%.1f dB##c",
                                  top, diag.mag_db[i]);
                } else {
                    std::snprintf(cell, sizeof(cell), "%s\n##c", top);
                }

                // Paired-ness is the BORDER; severity is the FILL.  One
                // encoding cannot carry both.
                const bool is_selected = (found >= 0 && found == state.selected_loop);
                const bool colour_cell = paired && diag.has_rga;
                ImGui::PushStyleColor(ImGuiCol_Button,
                    colour_cell ? severityFill(diag.rga_number)
                                : ImVec4(0.16f, 0.16f, 0.18f, 1.0f));
                ImGui::PushStyleColor(ImGuiCol_Border,
                    is_selected ? ImVec4(1.0f, 1.0f, 1.0f, 1.0f)
                                : ImVec4(0.35f, 0.72f, 1.0f, 1.0f));
                ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize,
                                    paired ? 2.0f : 0.0f);

                // LEFT click is always non-destructive: it pairs an unpaired
                // cell (and selects it), or selects an already-paired one so
                // its gains appear below.  RIGHT click unpairs.
                //
                // Issue #6 specified a pure toggle — click to pair, click again
                // to unpair — which leaves an existing loop's gains unreachable
                // without destroying the loop first.  Overloading the same
                // button with "unpair if this one is already selected" fixes
                // that but makes whether a click destroys depend on hidden
                // state; it was mis-predicted repeatedly in testing.  Giving
                // destruction its own gesture keeps "click a cell to pair"
                // true and makes no left-click ever lose a tuning.
                if (ImGui::Button(cell, ImVec2(-FLT_MIN, 34))) {
                    if (!paired) {
                        state.loops.push_back(Loop{i, j, {}, {}});
                        state.selected_loop =
                            static_cast<int>(state.loops.size()) - 1;
                        state.needs_recompute = true;
                    } else {
                        state.selected_loop = found;
                    }
                }
                if (paired && ImGui::IsItemHovered() &&
                    ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
                    state.loops.erase(state.loops.begin() + found);
                    state.selected_loop = -1;
                    state.needs_recompute = true;
                }

                ImGui::PopStyleVar();
                ImGui::PopStyleColor(2);

                // Dead-channel marker, drawn INTO the cell's top-right rather
                // than with SameLine: the button is full column width, so a
                // SameLine marker lands past the column edge and is clipped by
                // the next column — visible only in the last column.
                if (paired &&
                    found < static_cast<int>(diag.loop_dead.size()) &&
                    diag.loop_dead[found]) {
                    const ImVec2 mn = ImGui::GetItemRectMin();
                    const ImVec2 mx = ImGui::GetItemRectMax();
                    ImGui::GetWindowDrawList()->AddText(
                        ImVec2(mx.x - 11.0f, mn.y + 1.0f),
                        IM_COL32(255, 90, 90, 255), "!");
                    if (ImGui::IsItemHovered())
                        ImGui::SetTooltip(
                            "u%d does not affect y%d at any frequency", j, i);
                }
                ImGui::PopID();
            }
            ImGui::PopID();
        }
        ImGui::EndTable();
    }

    // One omega, plus a collapsed sweep.  No separate steady-state frequency:
    // S&P rule 2 is a DIC theorem requiring a stable plant, and neither
    // Ball-Balancer nor Inverted Pendulum qualifies, so a permanently
    // inapplicable readout would earn its space on no preset.
    ImGui::TextUnformatted("\xcf\x89 [Hz]"); ImGui::SameLine();
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (ImGui::SliderFloat("##diagw", &state.diag_omega_hz, 0.01f, 100.0f,
                           "%.3g", ImGuiSliderFlags_Logarithmic))
        state.needs_recompute = true;

    // Only offered where an RGA exists — a flat zero strip under Channel Share
    // would imply a number that is not defined.
    if (state.diagnostics.has_rga) {
        ImGui::Checkbox("RGA number vs \xcf\x89", &state.diag_sweep);
    }
    if (state.diag_sweep && state.diagnostics.has_rga) {
        float vals[40];
        float vmin = FLT_MAX, vmax = -FLT_MAX;
        for (int k = 0; k < 40; ++k) {
            const double f = state.freq_min_hz *
                std::pow(state.freq_max_hz / state.freq_min_hz, k / 39.0);
            const double v = rgaNumberAt(state.plant, state.loops, f);
            vals[k] = std::isnan(v) ? 0.0f : static_cast<float>(v);
            vmin = std::min(vmin, vals[k]);
            vmax = std::max(vmax, vals[k]);
        }
        // Explicit headroom.  Auto-scaling puts a CONSTANT sweep flush against
        // the frame's top edge, where it is invisible — and a constant RGA
        // number is exactly the Ball-Balancer demo case (4.00 at every omega,
        // by scaling invariance).
        char overlay[64];
        std::snprintf(overlay, sizeof(overlay), "%.2f - %.2f over %.3g-%.3g Hz",
                      vmin, vmax, state.freq_min_hz, state.freq_max_hz);
        ImGui::PlotLines("##sweep", vals, 40, 0, overlay,
                         0.0f, vmax * 1.25f + 0.1f, ImVec2(-FLT_MIN, 46));
    }

    // The subtitle carries the definition AND the explicit negative.  The
    // failure mode guarded against is a user reading a number in the pairing
    // area and assuming it grades their pairing.
    ImGui::TextWrapped("%s", state.diagnostics.headline.c_str());
    if (state.diagnostics.has_share) {
        // Wrapped, not TextDisabled: at the panel's real 354 px this sentence
        // truncates mid-word, and it is the load-bearing negative — the whole
        // point is that a user must not read it as a pairing grade.
        ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyle().Colors[ImGuiCol_TextDisabled]);
        ImGui::TextWrapped(
            "Channel Share is not the RGA and does not grade the pairing.");
        ImGui::PopStyleColor();
    }

    // Gain block for the selected loop, full width, below the grid.
    if (state.selected_loop >= 0 &&
        state.selected_loop < static_cast<int>(state.loops.size())) {
        Loop& l = state.loops[state.selected_loop];
        char hdr[64];
        std::snprintf(hdr, sizeof(hdr), "Gains  y%d <- u%d",
                      l.out, l.in);
        ImGui::SeparatorText(hdr);

        bool g_changed = false;
        if (state.ctrl_type == ControllerType::PID) {
            ImGui::TextUnformatted("Kp"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            g_changed |= ImGui::SliderFloat("##Kp", &l.pid.Kp, 0.01f, 100.0f,
                                            "%.4g", ImGuiSliderFlags_Logarithmic);
            ImGui::TextUnformatted("Ki"); ImGui::SameLine();
            g_changed |= gainSliderPiecewise("##Ki", &l.pid.Ki, 50.0f);
            ImGui::TextUnformatted("Kd"); ImGui::SameLine();
            g_changed |= gainSliderPiecewise("##Kd", &l.pid.Kd, 20.0f);
            ImGui::TextUnformatted("\xcf\x84" "f"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            g_changed |= ImGui::SliderFloat("##tauf", &l.pid.tau_f,
                                            0.001f, 1.0f, "%.4g",
                                            ImGuiSliderFlags_Logarithmic);
        } else {
            // Kind is one top-level combo for the whole list; Lead/Lag MODE is
            // per loop, inside the gain block.
            const char* modes[] = {"Lead", "Lag", "Lead-Lag"};
            int mi = static_cast<int>(l.leadlag.mode);
            if (ImGui::Combo("Mode", &mi, modes, 3)) {
                l.leadlag.mode = static_cast<CompensatorMode>(mi);
                g_changed = true;
            }
            ImGui::TextUnformatted("Kc"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            g_changed |= ImGui::SliderFloat("##Kc", &l.leadlag.Kc, 0.1f, 100.0f,
                                            "%.4g", ImGuiSliderFlags_Logarithmic);
            if (l.leadlag.mode != CompensatorMode::Lag) {
                ImGui::TextUnformatted("\xce\xb1 lead"); ImGui::SameLine();
                ImGui::SetNextItemWidth(-FLT_MIN);
                g_changed |= ImGui::SliderFloat("##al", &l.leadlag.alpha_lead,
                                                0.01f, 0.99f, "%.4g",
                                                ImGuiSliderFlags_Logarithmic);
                ImGui::TextUnformatted("T lead"); ImGui::SameLine();
                ImGui::SetNextItemWidth(-FLT_MIN);
                g_changed |= ImGui::SliderFloat("##Tl", &l.leadlag.T_lead,
                                                0.001f, 100.0f, "%.4g",
                                                ImGuiSliderFlags_Logarithmic);
            }
            if (l.leadlag.mode != CompensatorMode::Lead) {
                ImGui::TextUnformatted("\xce\xb1 lag"); ImGui::SameLine();
                ImGui::SetNextItemWidth(-FLT_MIN);
                g_changed |= ImGui::SliderFloat("##ag", &l.leadlag.alpha_lag,
                                                1.01f, 100.0f, "%.4g",
                                                ImGuiSliderFlags_Logarithmic);
                ImGui::TextUnformatted("T lag"); ImGui::SameLine();
                ImGui::SetNextItemWidth(-FLT_MIN);
                g_changed |= ImGui::SliderFloat("##Tg", &l.leadlag.T_lag,
                                                0.001f, 100.0f, "%.4g",
                                                ImGuiSliderFlags_Logarithmic);
            }
        }
        if (g_changed) state.needs_recompute = true;
    } else {
        ImGui::TextDisabled("select a paired cell to edit its gains");
    }
}

// The LQR design surface: state weights, input weights, the resulting gain,
// and where it put the closed-loop poles.  No Design button — the gain is
// re-solved on every recompute, so dragging a weight moves K and the poles
// together, which is the point of designing this way rather than by hand.
//
// LQR needs the full state vector.  The simulator has it; a physical plate
// would need an observer.  That is a documented limitation, not a defect.
//
// `is_cascade` is passed in rather than inferred, because a plant is
// identified by its NAME and this function cannot see one.  `7 states, 3
// inputs` is a shape, and shapes coincide — see `isCascadeModel`.  The state
// labels below still key on shape, which is a weaker claim they were always
// making; the presets are a claim about a particular ball and need the
// stronger one.
static void drawLqrDesigner(AppState& state, bool is_cascade) {
    const int n = state.plant.states();
    const int m = state.plant.inputs();
    if (state.lqr_q.size() != n || state.lqr_r.size() != m) {
        ImGui::TextDisabled("Sizing weights to the plant...");
        state.needs_recompute = true;
        return;
    }

    // State names, where the plant is one we know the physical reading of.
    // Anonymous x0..xn-1 otherwise: inventing labels for an arbitrary plant
    // would be worse than admitting we do not know them.
    static const char* kCascadeStates[7] = {
        "\xce\xb4\xce\xb1\xe2\x82\x81 leg 1", "\xce\xb4\xce\xb1\xe2\x82\x82 leg 2",
        "\xce\xb4\xce\xb1\xe2\x82\x83 leg 3", "x  ball", "y  ball",
        "x' ball", "y' ball"};
    const bool named = (n == 7 && m == 3);

    // --- Preset tunings (#19) ---
    // Three named tunings, first, because they are the fastest route from "I
    // have never seen a weighting matrix" to "controller design changes what
    // the ball does".  Flipping between Detuned and Aggressive teaches more in
    // five seconds than the sliders below teach in five minutes; the sliders
    // are what somebody does next, not what they do first.
    //
    // Offered only for the cascade, by name.  A preset is a claim about how a
    // particular ball behaves, and there is no such claim to make about some
    // other plant that happens to have seven states and three inputs — which
    // is all the `named` flag below ever established.
    if (is_cascade) {
        ImGui::TextDisabled("Tunings - one plant, three designs");
        const auto& presets = lqrPresets();
        const int active = activePreset(state.lqr_q, state.lqr_r);

        for (int i = 0; i < static_cast<int>(presets.size()); ++i) {
            if (i > 0) ImGui::SameLine();
            const bool on = (i == active);
            if (on) ImGui::PushStyleColor(ImGuiCol_Button,
                                          ImGui::GetStyleColorVec4(ImGuiCol_ButtonActive));
            if (ImGui::Button(presets[i].name)) {
                state.lqr_q = presetStateWeights(presets[i], n);
                state.lqr_r = presetInputWeights(presets[i], m);
                state.needs_recompute = true;
            }
            if (on) ImGui::PopStyleColor();
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Q on ball position %.0f, on ball velocity "
                                  "%.0f, R %.0f\n\n%s",
                                  presets[i].q_position, presets[i].q_velocity,
                                  presets[i].r, presets[i].blurb);
        }

        // Fixed height, wrapped: the blurbs are different lengths and this
        // panel is narrow, so an unreserved line would move every slider below
        // it by a row at the exact moment somebody was reaching for one.  Two
        // lines is what the longest blurb takes at the docked width; the
        // tooltip above carries the whole text whatever the panel is doing.
        ImGui::BeginChild("lqr_preset_blurb", ImVec2(0.0f,
                          ImGui::GetTextLineHeightWithSpacing() * 2.0f), false,
                          ImGuiWindowFlags_NoScrollbar |
                          ImGuiWindowFlags_NoScrollWithMouse);
        // No preset matches once a slider has been dragged, which is the
        // honest reading: the weights below are still whatever the user made
        // them, and the buttons are an offer rather than a mode.
        if (active >= 0) ImGui::TextWrapped("%s", presets[active].blurb);
        else ImGui::TextDisabled("custom weights");
        ImGui::EndChild();
        ImGui::Spacing();
    }

    ImGui::TextDisabled("Q - how much each state deviation costs");
    bool changed = false;
    char label[64];
    for (int i = 0; i < n; ++i) {
        double v = state.lqr_q(i);
        float f = static_cast<float>(v);
        if (named) std::snprintf(label, sizeof(label), "%s##q%d", kCascadeStates[i], i);
        else       std::snprintf(label, sizeof(label), "x%d##q%d", i, i);
        if (ImGui::SliderFloat(label, &f, 0.01f, 1000.0f, "%.2f",
                               ImGuiSliderFlags_Logarithmic)) {
            state.lqr_q(i) = f;
            changed = true;
        }
    }

    // What the last two sliders ALSO do, said where they are rather than left
    // for a visitor to deduce.  The trajectory feedforward has no gain of its
    // own — the path's velocity enters as a reference velocity and is
    // multiplied by K's velocity columns, which LQR designs from exactly these
    // two weights.  So the strength of the feedforward is set here, from
    // controls that would otherwise say nothing about it, and the alternative
    // — a feedforward slider — would be a second opinion about a number K
    // already owns.  Asserted in `test_auto_balance`; see #24.
    if (is_cascade) {
        ImGui::PushStyleColor(ImGuiCol_Text,
                              ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
        ImGui::TextWrapped("x' and y' also set how hard the loop chases a "
                           "MOVING setpoint: a path's own velocity enters the "
                           "reference and is multiplied by these columns of K. "
                           "The trajectory feedforward has no separate gain.");
        ImGui::PopStyleColor();
    }

    ImGui::Spacing();
    ImGui::TextDisabled("R - how expensive each input is");
    for (int j = 0; j < m; ++j) {
        float f = static_cast<float>(state.lqr_r(j));
        if (named) std::snprintf(label, sizeof(label), "leg %d cmd##r%d", j + 1, j);
        else       std::snprintf(label, sizeof(label), "u%d##r%d", j, j);
        if (ImGui::SliderFloat(label, &f, 0.001f, 1000.0f, "%.3f",
                               ImGuiSliderFlags_Logarithmic)) {
            state.lqr_r(j) = f;
            changed = true;
        }
    }
    if (changed) state.needs_recompute = true;

    // Only where there are no presets to reset TO.  On the cascade, Nominal
    // above IS this button — a second control doing the same thing would only
    // raise the question of how they differ.
    if (!is_cascade) {
        ImGui::Spacing();
        if (ImGui::Button("Reset weights")) {
            // Back to the tuning the application opens on, not to unit
            // weights: on this plant those park the ball in the
            // rolling-friction dead band, 37 mm off centre.  See
            // defaultLqrStateWeights.
            state.lqr_q = defaultLqrStateWeights(n);
            state.lqr_r = defaultLqrInputWeights(m);
            state.needs_recompute = true;
        }
    }

    // --- Result ---
    ImGui::SeparatorText("Gain");
    if (!state.lqr_result.success) {
        ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.35f, 1.0f), "LQR failed");
        ImGui::TextWrapped("%s", state.lqr_result.error.c_str());
        return;
    }

    ImGui::Text("K (%d x %d), u = -Kx",
                (int)state.lqr_result.K.rows(), (int)state.lqr_result.K.cols());
    if (ImGui::BeginTable("lqr_K", (int)state.lqr_result.K.cols() + 1,
                          ImGuiTableFlags_Borders | ImGuiTableFlags_SizingFixedFit)) {
        ImGui::TableNextColumn();
        ImGui::TextDisabled("row");
        for (int j = 0; j < state.lqr_result.K.cols(); ++j) {
            ImGui::TableNextColumn();
            ImGui::TextDisabled("x%d", j);
        }
        for (int i = 0; i < state.lqr_result.K.rows(); ++i) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextDisabled("u%d", i);
            for (int j = 0; j < state.lqr_result.K.cols(); ++j) {
                ImGui::TableNextColumn();
                ImGui::Text("%+.3f", state.lqr_result.K(i, j));
            }
        }
        ImGui::EndTable();
    }

    ImGui::SeparatorText("Closed-loop poles");
    double slowest = 0.0;
    bool stable = true;
    for (const auto& pole : state.lqr_result.closed_loop_poles) {
        if (pole.real() >= 0.0) stable = false;
        if (slowest == 0.0 || pole.real() > slowest) slowest = pole.real();
    }
    for (const auto& pole : state.lqr_result.closed_loop_poles) {
        if (std::abs(pole.imag()) < 1e-9)
            ImGui::Text("  %+.3f", pole.real());
        else
            ImGui::Text("  %+.3f %s %.3fj", pole.real(),
                        pole.imag() < 0 ? "-" : "+", std::abs(pole.imag()));
    }
    if (stable) {
        // The slowest pole sets how long the ball takes to settle, which is
        // the number the 3D view makes visible.
        ImGui::TextColored(ImVec4(0.3f, 1.0f, 0.3f, 1.0f),
                           "Stable - slowest pole %.3f (settles ~%.1f s)",
                           slowest, 4.0 / std::abs(slowest));
    } else {
        ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.35f, 1.0f),
                           "Closed loop is NOT stable");
    }
}

void drawModelPanel(AppState& state, const std::vector<ModelEntry>& presets) {
    ImGui::Begin("Model Configuration");

    // --- Plant section ---
    ImGui::SeparatorText("Plant G(s)");

    // Preset dropdown
    auto getter = [](void* data, int idx) -> const char* {
        auto* v = static_cast<const std::vector<ModelEntry>*>(data);
        return (*v)[idx].name.c_str();
    };
    if (ImGui::Combo("Preset", &state.preset_index, getter,
                     const_cast<void*>(
                         static_cast<const void*>(&presets)),
                     static_cast<int>(presets.size()))) {
        const auto& model = presets[state.preset_index];
        state.plant = model.system;
        state.current_params = model.params;
        matrixToTextBuf(state.plant.A, state.A_text, sizeof(state.A_text));
        matrixToTextBuf(state.plant.B, state.B_text, sizeof(state.B_text));
        matrixToTextBuf(state.plant.C, state.C_text, sizeof(state.C_text));
        matrixToTextBuf(state.plant.D, state.D_text, sizeof(state.D_text));
        state.output_i = 0;
        state.input_j = 0;
        extractTFFromSS(state);
        state.needs_recompute = true;
    }

    if (!presets.empty()) {
        ImGui::TextDisabled("%s", presets[state.preset_index].description.c_str());
    }

    // --- Physical parameters ---
    if (!state.current_params.empty()) {
        const auto& preset = presets[state.preset_index];
        if (ImGui::CollapsingHeader("Physical Parameters",
                                     ImGuiTreeNodeFlags_DefaultOpen)) {
            bool param_changed = false;
            for (size_t i = 0; i < state.current_params.size(); ++i) {
                auto& pp = state.current_params[i];
                char label[64];
                std::snprintf(label, sizeof(label), "%s [%s]",
                              pp.name.c_str(), pp.unit.c_str());
                ImGui::TextUnformatted(label);
                ImGui::SameLine();
                char slider_id[32];
                std::snprintf(slider_id, sizeof(slider_id), "##pp_%d",
                              static_cast<int>(i));
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGuiSliderFlags flags = pp.logarithmic
                    ? ImGuiSliderFlags_Logarithmic : 0;
                param_changed |= ImGui::SliderFloat(
                    slider_id, &pp.value, pp.min_val, pp.max_val,
                    "%.4g", flags);
            }
            if (param_changed && preset.builder) {
                state.plant = preset.builder(state.current_params);
                matrixToTextBuf(state.plant.A, state.A_text, sizeof(state.A_text));
                matrixToTextBuf(state.plant.B, state.B_text, sizeof(state.B_text));
                matrixToTextBuf(state.plant.C, state.C_text, sizeof(state.C_text));
                matrixToTextBuf(state.plant.D, state.D_text, sizeof(state.D_text));
                extractTFFromSS(state);
                state.needs_recompute = true;
            }
        }
    }

    // --- Derivation steps ---
    if (!presets.empty() && !presets[state.preset_index].derivation.empty()) {
        if (ImGui::CollapsingHeader("Derivation")) {
            for (const auto& step : presets[state.preset_index].derivation) {
                if (ImGui::TreeNode(step.title.c_str())) {
                    ImGui::TextWrapped("%s", step.content.c_str());
                    ImGui::TreePop();
                }
            }
        }
    }

    // Matrix text fields
    ImGui::InputText("A", state.A_text, sizeof(state.A_text));
    ImGui::InputText("B", state.B_text, sizeof(state.B_text));
    ImGui::InputText("C", state.C_text, sizeof(state.C_text));
    ImGui::InputText("D", state.D_text, sizeof(state.D_text));

    if (ImGui::Button("Apply Plant")) {
        auto A = parseMatrix(state.A_text);
        auto B = parseMatrix(state.B_text);
        auto C = parseMatrix(state.C_text);
        auto D = parseMatrix(state.D_text);
        if (A && B && C && D) {
            state.plant = {*A, *B, *C, *D};
            state.output_i = 0;
            state.input_j = 0;
            extractTFFromSS(state);
            state.needs_recompute = true;
        }
    }

    // System info
    ImGui::Text("n=%d  m=%d  p=%d",
                state.plant.states(), state.plant.inputs(),
                state.plant.outputs());
    ImGui::SameLine();

    // Stability indicator helper lambda
    auto showStability = [](const PoleZeroResult& pz) {
        bool has_positive = false;
        bool has_zero = false;
        for (const auto& p : pz.poles) {
            if (p.real() > 1e-10) has_positive = true;
            else if (std::abs(p.real()) <= 1e-10) has_zero = true;
        }
        if (has_positive) {
            ImGui::TextColored(ImVec4(1, 0.3f, 0.3f, 1), "UNSTABLE");
        } else if (has_zero) {
            ImGui::TextColored(ImVec4(1, 0.8f, 0.2f, 1), "MARGINAL");
        } else {
            ImGui::TextColored(ImVec4(0.3f, 1, 0.3f, 1), "STABLE");
        }
    };

    if (state.system_valid[0]) {
        showStability(state.pole_zero[0]);
    }
    if (state.system_valid[3]) {
        ImGui::SameLine();
        ImGui::Text("CL:");
        ImGui::SameLine();
        showStability(state.pole_zero[3]);
    }

    // --- Transfer function display (1st/2nd order SISO only) ---
    if (state.plant_is_1st_order) {
        ImGui::SeparatorText("Transfer Function");
        ImGui::TextDisabled("G(s) = K / (\xcf\x84s + 1)");
        bool tf_changed = false;
        ImGui::TextUnformatted("K"); ImGui::SameLine();
        ImGui::SetNextItemWidth(-FLT_MIN);
        tf_changed |= ImGui::SliderFloat("##tf_K", &state.tf_K, -10.0f, 10.0f);
        ImGui::TextUnformatted("\xcf\x84"); ImGui::SameLine();
        ImGui::SetNextItemWidth(-FLT_MIN);
        tf_changed |= ImGui::SliderFloat("##tf_tau", &state.tf_tau, 0.01f, 100.0f,
                                          "%.3f", ImGuiSliderFlags_Logarithmic);
        if (tf_changed) {
            applyTFToSS(state);
            extractTFFromSS(state);
            state.needs_recompute = true;
        }
    } else if (state.plant_is_2nd_order) {
        ImGui::SeparatorText("Transfer Function");
        ImGui::TextDisabled("G(s) = K\xc2\xb7\xcf\x89n\xc2\xb2 / (s\xc2\xb2 + 2\xce\xb6\xcf\x89n\xc2\xb7s + \xcf\x89n\xc2\xb2)");
        bool tf_changed = false;
        ImGui::TextUnformatted("K"); ImGui::SameLine();
        ImGui::SetNextItemWidth(-FLT_MIN);
        tf_changed |= ImGui::SliderFloat("##tf_K", &state.tf_K, -10.0f, 10.0f);
        ImGui::TextUnformatted("\xcf\x89n"); ImGui::SameLine();
        ImGui::SetNextItemWidth(-FLT_MIN);
        tf_changed |= ImGui::SliderFloat("##tf_wn", &state.tf_wn, 0.01f, 100.0f,
                                          "%.3f", ImGuiSliderFlags_Logarithmic);
        ImGui::TextUnformatted("\xce\xb6"); ImGui::SameLine();
        ImGui::SetNextItemWidth(-FLT_MIN);
        tf_changed |= ImGui::SliderFloat("##tf_zeta", &state.tf_zeta, 0.0f, 5.0f);
        if (tf_changed) {
            applyTFToSS(state);
            extractTFFromSS(state);
            state.needs_recompute = true;
        }
    }

    // --- State-space matrix sliders ---
    ImGui::SeparatorText("State-Space Matrices");
    bool slider_changed = false;
    ImGui::TextUnformatted("A:");
    slider_changed |= drawMatrixSliders("A", state.plant.A, state.A_text, sizeof(state.A_text));
    ImGui::TextUnformatted("B:");
    slider_changed |= drawMatrixSliders("B", state.plant.B, state.B_text, sizeof(state.B_text));
    ImGui::TextUnformatted("C:");
    slider_changed |= drawMatrixSliders("C", state.plant.C, state.C_text, sizeof(state.C_text));
    ImGui::TextUnformatted("D:");
    slider_changed |= drawMatrixSliders("D", state.plant.D, state.D_text, sizeof(state.D_text));
    if (slider_changed) {
        extractTFFromSS(state);
        state.needs_recompute = true;
    }

    // --- Controller section ---
    ImGui::SeparatorText("Controller");
    // Index-coupled to ControllerType.  IM_ARRAYSIZE rather than a literal:
    // the previous hardcoded 5 would have hidden any member added here.
    const char* ctrl_types[] = {"None", "PID", "Lead/Lag",
                                "State-Space", "Gain Matrix K", "LQR"};
    static_assert(IM_ARRAYSIZE(ctrl_types) ==
                      static_cast<int>(ControllerType::LQR) + 1,
                  "ctrl_types is index-coupled to ControllerType");
    int ctrl_idx = static_cast<int>(state.ctrl_type);
    if (ImGui::Combo("Type", &ctrl_idx, ctrl_types, IM_ARRAYSIZE(ctrl_types))) {
        state.ctrl_type = static_cast<ControllerType>(ctrl_idx);
        state.needs_recompute = true;
        const bool lb = state.ctrl_type == ControllerType::PID ||
                        state.ctrl_type == ControllerType::LeadLag;
        state.trace_visible[3] = state.ctrl_type != ControllerType::None;
        state.trace_visible[1] =
            lb || state.ctrl_type == ControllerType::StateSpace;
    }

    if (state.ctrl_type == ControllerType::PID ||
        state.ctrl_type == ControllerType::LeadLag) {
        drawPairingGrid(state);
    } else if (state.ctrl_type == ControllerType::StateSpace) {
        ImGui::InputText("Ac", state.Ac_text, sizeof(state.Ac_text));
        ImGui::InputText("Bc", state.Bc_text, sizeof(state.Bc_text));
        ImGui::InputText("Cc", state.Cc_text, sizeof(state.Cc_text));
        ImGui::InputText("Dc", state.Dc_text, sizeof(state.Dc_text));
        if (ImGui::Button("Apply Controller")) {
            auto Ac = parseMatrix(state.Ac_text);
            auto Bc = parseMatrix(state.Bc_text);
            auto Cc = parseMatrix(state.Cc_text);
            auto Dc = parseMatrix(state.Dc_text);
            if (Ac && Bc && Cc && Dc &&
                Ac->rows() == Ac->cols() &&
                Bc->rows() == Ac->rows() &&
                Cc->cols() == Ac->cols() &&
                Dc->rows() == Cc->rows() &&
                Dc->cols() == Bc->cols()) {
                state.ctrl_ss = {*Ac, *Bc, *Cc, *Dc};
                state.needs_recompute = true;
            }
        }
        if (state.ctrl_ss.A.size() > 0) {
            bool ctrl_slider_changed = false;
            ImGui::TextUnformatted("Ac:");
            ctrl_slider_changed |= drawMatrixSliders("Ac", state.ctrl_ss.A, state.Ac_text, sizeof(state.Ac_text));
            ImGui::TextUnformatted("Bc:");
            ctrl_slider_changed |= drawMatrixSliders("Bc", state.ctrl_ss.B, state.Bc_text, sizeof(state.Bc_text));
            ImGui::TextUnformatted("Cc:");
            ctrl_slider_changed |= drawMatrixSliders("Cc", state.ctrl_ss.C, state.Cc_text, sizeof(state.Cc_text));
            ImGui::TextUnformatted("Dc:");
            ctrl_slider_changed |= drawMatrixSliders("Dc", state.ctrl_ss.D, state.Dc_text, sizeof(state.Dc_text));
            if (ctrl_slider_changed) state.needs_recompute = true;
        }
    } else if (state.ctrl_type == ControllerType::GainMatrix) {
        ImGui::TextDisabled("K must be %d x %d (inputs x states)",
                            state.plant.inputs(), state.plant.states());
        ImGui::InputText("K", state.K_text, sizeof(state.K_text));
        if (ImGui::Button("Apply K")) {
            auto K = parseMatrix(state.K_text);
            if (K && K->rows() == state.plant.inputs() &&
                K->cols() == state.plant.states()) {
                state.ctrl_K = *K;
                state.needs_recompute = true;
            } else if (K) {
                ImGui::OpenPopup("K_dim_err");
            }
        }
        if (ImGui::BeginPopup("K_dim_err")) {
            ImGui::TextColored(ImVec4(1, 0.3f, 0.3f, 1),
                "K must be %d x %d, got %d x %d",
                state.plant.inputs(), state.plant.states(),
                (int)state.ctrl_K.rows(), (int)state.ctrl_K.cols());
            ImGui::EndPopup();
        }
        if (state.ctrl_K.size() > 0) {
            ImGui::TextUnformatted("K:");
            if (drawMatrixSliders("K", state.ctrl_K, state.K_text, sizeof(state.K_text)))
                state.needs_recompute = true;
        }
    } else if (state.ctrl_type == ControllerType::LQR) {
        drawLqrDesigner(state, state.preset_index >= 0 &&
                               state.preset_index < (int)presets.size() &&
                               isCascadeModel(presets[state.preset_index]));
    }

    // --- Channel selector ---
    ImGui::SeparatorText("Channel");
    int p = state.plant.outputs();
    int m = state.plant.inputs();
    if (ImGui::SliderInt("Output i", &state.output_i, 0, std::max(p - 1, 0))) {
        state.needs_recompute = true;
    }
    if (ImGui::SliderInt("Input j", &state.input_j, 0, std::max(m - 1, 0))) {
        state.needs_recompute = true;
    }
    if (ImGui::SliderInt("Reference r", &state.ref_r, 0, std::max(p - 1, 0))) {
        state.needs_recompute = true;
    }
    if (p <= 3 && m <= 3) {
        if (ImGui::Checkbox("Show All Channels", &state.show_all_channels)) {
            state.needs_recompute = true;
        }
    }

    // --- Trace toggles (skip system 2: open-loop combined) ---
    ImGui::SeparatorText("Traces");
    for (int s = 0; s < NUM_SYSTEMS; ++s) {
        if (s == 2) continue;
        if (!state.system_valid[s]) continue;
        ImGui::PushStyleColor(ImGuiCol_CheckMark, system_colors[s]);
        ImGui::Checkbox(system_names[s], &state.trace_visible[s]);
        ImGui::PopStyleColor();
        // Find next visible system for SameLine
        bool has_next = false;
        for (int ns = s + 1; ns < NUM_SYSTEMS; ++ns) {
            if (ns == 2) continue;
            if (state.system_valid[ns]) { has_next = true; break; }
        }
        if (has_next) ImGui::SameLine();
    }

    // --- Reset ---
    ImGui::Spacing();
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.7f, 0.2f, 0.2f, 1.0f));
    if (ImGui::Button("Reset All", ImVec2(-1, 0))) {
        state = AppState{};
        state.preset_index = defaultModelIndex(presets);
        state.plant = presets[state.preset_index].system;
        state.current_params = presets[state.preset_index].params;
        matrixToTextBuf(state.plant.A, state.A_text, sizeof(state.A_text));
        matrixToTextBuf(state.plant.B, state.B_text, sizeof(state.B_text));
        matrixToTextBuf(state.plant.C, state.C_text, sizeof(state.C_text));
        matrixToTextBuf(state.plant.D, state.D_text, sizeof(state.D_text));
        extractTFFromSS(state);
        state.needs_recompute = true;
    }
    ImGui::PopStyleColor();

    ImGui::End();
}

}  // namespace caliburn
