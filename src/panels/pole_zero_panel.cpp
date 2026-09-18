// src/panels/pole_zero_panel.cpp
#include "pole_zero_panel.h"
#include "panel_utils.h"
#include "implot.h"
#include <cmath>

namespace caliburn {

void drawPoleZeroPanel(AppState& state) {
    if (!state.show_pole_zero) return;
    ImGui::Begin("Pole-Zero / Root Locus", &state.show_pole_zero);
    drawHelpMarker(
        "Pole-Zero / Root Locus Plot\n\n"
        "Poles (x) are roots of det(sI - A) = 0.\n"
        "Zeros (o) are transmission zeros.\n\n"
        "Stability: All poles must be in the left half-plane (Re < 0).\n"
        "Green shading = stable region.\n\n"
        "Pole location:\n"
        "  Real part: decay/growth rate\n"
        "  Imag part: oscillation frequency\n"
        "  Distance from origin: speed of response\n"
        "  Angle from neg. real axis: damping ratio\n\n"
        "Pole-Zero: poles and zeros of each visible system.\n"
        "Plant Locus: proportional gain K on the bare plant channel (i,j).\n"
        "  The controller is NOT in the path - this is the pre-tuning tool.\n"
        "Loop Locus: kappa scaling the whole compensator on channel (i,j),\n"
        "  a true locus of 1 + kappa*L(s) = 0. kappa = 1 is the live design.\n"
        "State FB: alpha scaling the gain matrix K.\n\n"
        "No mode draws zeros.\n"
        "Unstable when a locus crosses into the right half-plane.");
    ImGui::SameLine();

    // Mode selector
    int mode = static_cast<int>(state.pz_mode);
    // All four stay selectable.  Disabling them plus an auto-fallback would
    // silently yank the user out of their chosen mode on a controller-type
    // switch; hiding them would change the control row's width as the
    // controller is edited, with nothing on screen saying the mode exists.
    ImGui::RadioButton("Pole-Zero", &mode, 0); ImGui::SameLine();
    ImGui::RadioButton("Plant Locus", &mode, 1); ImGui::SameLine();
    ImGui::RadioButton("Loop Locus", &mode, 2); ImGui::SameLine();
    ImGui::RadioButton("State FB", &mode, 3);
    if (static_cast<PZMode>(mode) != state.pz_mode) {
        state.pz_mode = static_cast<PZMode>(mode);
        state.needs_recompute = true;
    }

    // Root locus controls
    if (state.pz_mode == PZMode::PlantLocus) {
        if (ImGui::BeginTable("##rl_ctrls", 3)) {
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("K"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderFloat("##K", &state.rl_current_k,
                               state.rl_k_min, state.rl_k_max, "%.2f",
                               ImGuiSliderFlags_Logarithmic);
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("K min"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::SliderFloat("##Kmin", &state.rl_k_min, 0.0f, 1.0f))
                state.needs_recompute = true;
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("K max"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::SliderFloat("##Kmax", &state.rl_k_max, 1.0f, 1000.0f))
                state.needs_recompute = true;
            ImGui::EndTable();
        }
    } else if (state.pz_mode == PZMode::LoopLocus) {
        if (ImGui::BeginTable("##ll_ctrls", 3)) {
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("\xce\xba"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderFloat("##kappa", &state.rl_kappa_current,
                               state.rl_kappa_min, state.rl_kappa_max, "%.3g",
                               ImGuiSliderFlags_Logarithmic);
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("\xce\xba min"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            // Clamped away from zero: kappa = 0 would collapse the loop's
            // states and make the closed-loop dimension ragged mid-sweep.
            if (ImGui::SliderFloat("##kmin", &state.rl_kappa_min, 0.001f, 1.0f,
                                   "%.3g", ImGuiSliderFlags_Logarithmic))
                state.needs_recompute = true;
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("\xce\xba max"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::SliderFloat("##kmax", &state.rl_kappa_max, 1.0f, 1000.0f,
                                   "%.3g", ImGuiSliderFlags_Logarithmic))
                state.needs_recompute = true;
            ImGui::EndTable();
        }
        // The selected loop's gain margin, with every other loop held at its
        // tuning — read off the sweep already computed, so it costs nothing.
        if (state.rl_loop_gain_margin > 0.0) {
            const double ratio = state.rl_loop_gain_margin /
                                 std::max(state.rl_kappa_current, 1e-6f);
            ImGui::TextColored(ImVec4(0.3f, 1, 0.3f, 1),
                "\xce\xba* = %.3g   margin %.2f\xc3\x97 / %.1f dB",
                state.rl_loop_gain_margin, ratio, 20.0 * std::log10(ratio));
        } else if (!state.root_locus.empty()) {
            ImGui::TextDisabled("no crossing in [%.3g, %.3g]",
                                state.rl_kappa_min, state.rl_kappa_max);
        }
    } else if (state.pz_mode == PZMode::StateFB) {
        if (ImGui::BeginTable("##sfb_ctrls", 2)) {
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("\xce\xb1"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderFloat("##alpha", &state.rl_current_alpha,
                               state.rl_alpha_min, state.rl_alpha_max);
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("\xce\xb1 max"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::SliderFloat("##alphamax", &state.rl_alpha_max, 0.1f, 10.0f))
                state.needs_recompute = true;
            ImGui::EndTable();
        }
    }

    // A locus mode with an unmet precondition draws the s-plane furniture and
    // states why, rather than the previous mode's data under the wrong label.
    if (state.pz_mode != PZMode::PoleZero && state.root_locus.empty()) {
        const char* why =
            (state.pz_mode == PZMode::LoopLocus)
                ? "Loop Locus needs a PID or Lead/Lag controller with at least "
                  "one loop on the selected channel"
          : (state.pz_mode == PZMode::StateFB)
                ? "State FB needs a Gain Matrix K controller"
                : "no locus data";
        ImGui::TextDisabled("%s", why);
    }

    // Drawn BEFORE the plots — they consume GetContentRegionAvail(), and
    // anything after them is clipped out of the panel entirely.  Which is also
    // why this is one fixed line rather than none, one or two: a line that came
    // and went here would resize the plot underneath it.
    drawSuppressedTraces(state.channel_reason, system_names, NUM_SYSTEMS, 2);

    // Compute axis limits from visible data (mode-dependent)
    double pz_xmin = 0, pz_xmax = 0, pz_ymin = 0, pz_ymax = 0;
    bool pz_has_data = false;
    auto pz_extend = [&](double re, double im) {
        if (!pz_has_data) { pz_xmin = pz_xmax = re; pz_ymin = pz_ymax = im; pz_has_data = true; }
        else { pz_xmin = std::min(pz_xmin, re); pz_xmax = std::max(pz_xmax, re);
               pz_ymin = std::min(pz_ymin, im); pz_ymax = std::max(pz_ymax, im); }
    };
    if (state.pz_mode == PZMode::PoleZero) {
        for (int s = 0; s < NUM_SYSTEMS; ++s) {
            if (!state.trace_visible[s] || !state.system_valid[s]) continue;
            for (const auto& p : state.pole_zero[s].poles) pz_extend(p.real(), p.imag());
            for (const auto& z : state.pole_zero[s].zeros) pz_extend(z.real(), z.imag());
        }
    } else {
        for (const auto& step : state.root_locus)
            for (const auto& p : step.poles) pz_extend(p.real(), p.imag());
    }
    if (!pz_has_data) { pz_xmin = -2; pz_xmax = 2; pz_ymin = -2; pz_ymax = 2; }
    {
        double mx = std::max((pz_xmax - pz_xmin) * 0.15, 0.5);
        double my = std::max((pz_ymax - pz_ymin) * 0.15, 0.5);
        pz_xmin -= mx; pz_xmax += mx; pz_ymin -= my; pz_ymax += my;
        pz_xmin = std::min(pz_xmin, -1.0);
        pz_xmax = std::max(pz_xmax, 1.0);
    }

    // Plot
    ImVec2 avail = ImGui::GetContentRegionAvail();
    if (ImPlot::BeginPlot("##pz_plot", avail, ImPlotFlags_Equal)) {
        ImPlot::SetupAxes("\xcf\x83 (Real)", "j\xcf\x89 (Imag)",
                          ImPlotAxisFlags_NoGridLines,
                          ImPlotAxisFlags_NoGridLines);
        ImPlot::SetupAxisLimits(ImAxis_X1, pz_xmin, pz_xmax, ImPlotCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, pz_ymin, pz_ymax, ImPlotCond_Always);

        // Stable region shading (Re < 0)
        {
            auto limits = ImPlot::GetPlotLimits();
            double shade_x[4] = {limits.X.Min, 0, 0, limits.X.Min};
            double shade_y[4] = {limits.Y.Min, limits.Y.Min, limits.Y.Max, limits.Y.Max};
            ImPlot::PlotShaded("##stable", shade_x, shade_y, 4, 0.0,
                               ImPlotSpec(ImPlotProp_FillColor,
                                          ImVec4(0.20f, 0.78f, 0.20f, 0.08f)));
        }

        // Stability boundary (imaginary axis)
        double zero = 0.0;
        ImPlot::PlotInfLines("##jw_axis", &zero, 1,
                             ImPlotSpec(ImPlotProp_LineColor,
                                        ImVec4(1, 1, 1, 0.16f)));

        // Real axis
        ImPlot::PlotInfLines("##sigma_axis", &zero, 1,
                             ImPlotSpec(ImPlotProp_LineColor,
                                        ImVec4(1, 1, 1, 0.16f),
                                        ImPlotProp_Flags,
                                        (int)ImPlotInfLinesFlags_Horizontal));

        // S-plane grid (wn circles + zeta lines)
        drawSPlaneGrid();

        if (state.pz_mode == PZMode::PoleZero) {
            for (int s = 0; s < NUM_SYSTEMS; ++s) {
                if (!state.trace_visible[s] || !state.system_valid[s]) continue;
                const auto& pz = state.pole_zero[s];

                std::vector<double> pre, pim;
                for (const auto& p : pz.poles) {
                    pre.push_back(p.real());
                    pim.push_back(p.imag());
                }
                if (!pre.empty()) {
                    char label_p[32];
                    std::snprintf(label_p, sizeof(label_p), "P(%s)", system_names[s]);
                    ImPlot::PlotScatter(label_p, pre.data(), pim.data(),
                                       static_cast<int>(pre.size()),
                                       ImPlotSpec(ImPlotProp_Marker, ImPlotMarker_Cross,
                                                  ImPlotProp_MarkerSize, 8.0f,
                                                  ImPlotProp_MarkerLineColor, system_colors[s],
                                                  ImPlotProp_LineWeight, 2.0f));
                }

                std::vector<double> zre, zim;
                for (const auto& z : pz.zeros) {
                    zre.push_back(z.real());
                    zim.push_back(z.imag());
                }
                if (!zre.empty()) {
                    char label_z[32];
                    std::snprintf(label_z, sizeof(label_z), "Z(%s)", system_names[s]);
                    ImPlot::PlotScatter(label_z, zre.data(), zim.data(),
                                       static_cast<int>(zre.size()),
                                       ImPlotSpec(ImPlotProp_Marker, ImPlotMarker_Circle,
                                                  ImPlotProp_MarkerSize, 8.0f,
                                                  ImPlotProp_MarkerFillColor, ImVec4(0, 0, 0, 0),
                                                  ImPlotProp_MarkerLineColor, system_colors[s],
                                                  ImPlotProp_LineWeight, 2.0f));
                }
            }
        } else {
            if (!state.root_locus.empty()) {
                int n_poles = static_cast<int>(state.root_locus[0].poles.size());
                int n_steps = static_cast<int>(state.root_locus.size());

                for (int p = 0; p < n_poles; ++p) {
                    std::vector<double> re(n_steps), im(n_steps);
                    for (int k = 0; k < n_steps; ++k) {
                        re[k] = state.root_locus[k].poles[p].real();
                        im[k] = state.root_locus[k].poles[p].imag();
                    }
                    char lbl[16];
                    std::snprintf(lbl, sizeof(lbl), "##rl%d", p);
                    ImPlot::PlotLine(lbl, re.data(), im.data(), n_steps,
                                     ImPlotSpec(ImPlotProp_LineColor, system_colors[0],
                                                ImPlotProp_LineWeight, 1.5f));
                }

                double current_gain = state.rl_current_k;
                if (state.pz_mode == PZMode::LoopLocus)
                    current_gain = state.rl_kappa_current;
                else if (state.pz_mode == PZMode::StateFB)
                    current_gain = state.rl_current_alpha;
                int best = 0;
                double best_d = 1e30;
                for (int k = 0; k < n_steps; ++k) {
                    double d = std::abs(state.root_locus[k].gain - current_gain);
                    if (d < best_d) { best_d = d; best = k; }
                }
                std::vector<double> cre, cim;
                for (const auto& pole : state.root_locus[best].poles) {
                    cre.push_back(pole.real());
                    cim.push_back(pole.imag());
                }
                ImPlot::PlotScatter("##current", cre.data(), cim.data(),
                                   static_cast<int>(cre.size()),
                                   ImPlotSpec(ImPlotProp_Marker, ImPlotMarker_Diamond,
                                              ImPlotProp_MarkerSize, 10.0f,
                                              ImPlotProp_MarkerLineColor, system_colors[3],
                                              ImPlotProp_LineWeight, 2.5f));
            }

            // One rule for both locus modes: the reference is where the sweep
            // actually starts.  For Plant Locus with k_min = 0 that IS the
            // plant's poles, so nothing changes there.  For Loop Locus the
            // plant's poles are simply the wrong reference — the branches start
            // near the poles of plant-plus-other-loops-plus-this-compensator.
            // It also deletes a plant-shaped pole_zero[0] read, the kind of
            // assumption rule C exists to remove.
            if (!state.root_locus.empty()) {
                std::vector<double> pre, pim;
                for (const auto& pole : state.root_locus.front().poles) {
                    pre.push_back(pole.real());
                    pim.push_back(pole.imag());
                }
                const char* sym = (state.pz_mode == PZMode::LoopLocus)
                                      ? "\xce\xba"
                                  : (state.pz_mode == PZMode::StateFB)
                                      ? "\xce\xb1"
                                      : "K";
                char lbl[48];
                std::snprintf(lbl, sizeof(lbl), "Start (%s = %.3g)", sym,
                              state.root_locus.front().gain);
                ImPlot::PlotScatter(lbl, pre.data(), pim.data(),
                                   static_cast<int>(pre.size()),
                                   ImPlotSpec(ImPlotProp_Marker, ImPlotMarker_Cross,
                                              ImPlotProp_MarkerSize, 6.0f,
                                              ImPlotProp_MarkerLineColor,
                                              ImVec4(0.5f, 0.5f, 0.5f, 0.7f),
                                              ImPlotProp_LineWeight, 1.5f));
            }
        }

        ImPlot::EndPlot();
    }

    ImGui::End();
}

}  // namespace caliburn
