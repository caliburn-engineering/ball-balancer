// src/panels/bode_panel.cpp
#include "bode_panel.h"
#include "panel_utils.h"
#include "implot.h"
#include <cmath>
#include <cstdio>

namespace caliburn {

namespace {

// Interpolate Bode data at a given frequency (log-space)
void interpBode(const FrequencyResponse& bode, double freq_hz,
                double& mag_db, double& phase_deg) {
    mag_db = 0; phase_deg = 0;
    if (bode.points.empty()) return;
    int n = static_cast<int>(bode.points.size());
    if (freq_hz <= bode.points[0].freq_hz) {
        mag_db = bode.points[0].magnitude_db;
        phase_deg = bode.points[0].phase_deg;
        return;
    }
    if (freq_hz >= bode.points[n - 1].freq_hz) {
        mag_db = bode.points[n - 1].magnitude_db;
        phase_deg = bode.points[n - 1].phase_deg;
        return;
    }
    for (int k = 0; k < n - 1; ++k) {
        if (freq_hz >= bode.points[k].freq_hz && freq_hz <= bode.points[k + 1].freq_hz) {
            double t = (freq_hz - bode.points[k].freq_hz) /
                       (bode.points[k + 1].freq_hz - bode.points[k].freq_hz);
            mag_db = bode.points[k].magnitude_db * (1 - t) + bode.points[k + 1].magnitude_db * t;
            phase_deg = bode.points[k].phase_deg * (1 - t) + bode.points[k + 1].phase_deg * t;
            return;
        }
    }
}

void plotBodeSingle(const char* suffix, AppState& state, float plot_h, int margin_sys) {

    // Precompute y-limits for magnitude and phase
    double mag_lo = 1e30, mag_hi = -1e30, ph_lo = 1e30, ph_hi = -1e30;
    bool has_bode = false;
    for (int s = 0; s < NUM_SYSTEMS; ++s) {
        if (!state.trace_visible[s] || !state.system_valid[s]) continue;
        for (const auto& pt : state.bode[s].points) {
            mag_lo = std::min(mag_lo, pt.magnitude_db);
            mag_hi = std::max(mag_hi, pt.magnitude_db);
            ph_lo = std::min(ph_lo, pt.phase_deg);
            ph_hi = std::max(ph_hi, pt.phase_deg);
            has_bode = true;
        }
    }
    if (has_bode) {
        double mm = std::max((mag_hi - mag_lo) * 0.1, 5.0);
        mag_lo -= mm; mag_hi += mm;
        double pm = std::max((ph_hi - ph_lo) * 0.1, 10.0);
        ph_lo -= pm; ph_hi += pm;
    } else {
        mag_lo = -60; mag_hi = 10; ph_lo = -200; ph_hi = 10;
    }

    // Build plot titles with margin info
    char mag_id[128];
    char phase_id[128];
    if (margin_sys >= 0) {
        const auto& ms = state.bode[margin_sys];
        if (ms.phase_crossover_hz > 0)
            std::snprintf(mag_id, sizeof(mag_id), "Magnitude [dB]  |  GM = %.1f dB @ %.2f Hz%s",
                          ms.gain_margin_db, ms.phase_crossover_hz, suffix);
        else
            std::snprintf(mag_id, sizeof(mag_id), "Magnitude [dB]  |  GM = inf%s", suffix);

        if (ms.gain_crossover_hz > 0)
            std::snprintf(phase_id, sizeof(phase_id), "Phase [deg]  |  PM = %.1f\xc2\xb0 @ %.2f Hz%s",
                          ms.phase_margin_deg, ms.gain_crossover_hz, suffix);
        else
            std::snprintf(phase_id, sizeof(phase_id), "Phase [deg]  |  PM = inf%s", suffix);
    } else {
        std::snprintf(mag_id, sizeof(mag_id), "Magnitude [dB]%s", suffix);
        std::snprintf(phase_id, sizeof(phase_id), "Phase [deg]%s", suffix);
    }

    // --- Magnitude ---
    if (ImPlot::BeginPlot(mag_id, ImVec2(-1, plot_h))) {
        ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
        ImPlot::SetupAxes("Frequency [Hz]", "|G| [dB]");
        ImPlot::SetupAxisLimits(ImAxis_X1, state.freq_min_hz, state.freq_max_hz,
                                ImPlotCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, mag_lo, mag_hi, ImPlotCond_Always);

        for (int s = 0; s < NUM_SYSTEMS; ++s) {
            if (!state.trace_visible[s] || !state.system_valid[s]) continue;
            const auto& bode = state.bode[s];
            if (bode.points.empty()) continue;
            int n = static_cast<int>(bode.points.size());
            std::vector<double> f(n), m(n);
            for (int k = 0; k < n; ++k) {
                f[k] = bode.points[k].freq_hz;
                m[k] = bode.points[k].magnitude_db;
            }
            ImPlot::PlotLine(system_names[s], f.data(), m.data(), n,
                             ImPlotSpec(ImPlotProp_LineColor, system_colors[s],
                                        ImPlotProp_LineWeight, 2.0f));
        }

        double zero_db = 0.0;
        ImPlot::PlotInfLines("##0dB", &zero_db, 1,
                             ImPlotSpec(ImPlotProp_LineColor,
                                        ImVec4(1, 1, 1, 0.16f),
                                        ImPlotProp_Flags,
                                        (int)ImPlotInfLinesFlags_Horizontal));

        // GM vertical line (red) on magnitude plot
        if (margin_sys >= 0 && state.bode[margin_sys].phase_crossover_hz > 0) {
            double x = state.bode[margin_sys].phase_crossover_hz;
            double gm = state.bode[margin_sys].gain_margin_db;
            ImPlot::PlotInfLines("##gm_line", &x, 1,
                                 ImPlotSpec(ImPlotProp_LineColor,
                                            ImVec4(1, 0.39f, 0.39f, 0.47f),
                                            ImPlotProp_LineWeight, 1.0f));
            ImPlot::Annotation(x, 0.0, ImVec4(1, 0.4f, 0.4f, 1),
                               ImVec2(5, -10), true, "GM=%.1f dB", gm);
        }

        // PM vertical line (green) on magnitude plot too
        if (margin_sys >= 0 && state.bode[margin_sys].gain_crossover_hz > 0) {
            double x = state.bode[margin_sys].gain_crossover_hz;
            ImPlot::PlotInfLines("##pm_line_mag", &x, 1,
                                 ImPlotSpec(ImPlotProp_LineColor,
                                            ImVec4(0.39f, 0.78f, 0.39f, 0.47f),
                                            ImPlotProp_LineWeight, 1.0f));
        }

        // Cursor
        if (ImPlot::IsPlotHovered()) {
            double mx = ImPlot::GetPlotMousePos().x;
            ImPlot::PlotInfLines("##cursor", &mx, 1,
                                 ImPlotSpec(ImPlotProp_LineColor,
                                            ImVec4(1, 1, 1, 0.31f)));
            state.bode_cursor_x = mx;
        }
        if (state.bode_cursor_x > 0) {
            double cx = state.bode_cursor_x;
            if (!ImPlot::IsPlotHovered()) {
                ImPlot::PlotInfLines("##cursor", &cx, 1,
                                     ImPlotSpec(ImPlotProp_LineColor,
                                                ImVec4(1, 1, 1, 0.31f)));
            }
        }

        ImPlot::EndPlot();
    }

    // --- Phase ---
    if (ImPlot::BeginPlot(phase_id, ImVec2(-1, plot_h))) {
        ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
        ImPlot::SetupAxes("Frequency [Hz]", "arg G [\xc2\xb0]");
        ImPlot::SetupAxisLimits(ImAxis_X1, state.freq_min_hz, state.freq_max_hz,
                                ImPlotCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, ph_lo, ph_hi, ImPlotCond_Always);

        for (int s = 0; s < NUM_SYSTEMS; ++s) {
            if (!state.trace_visible[s] || !state.system_valid[s]) continue;
            const auto& bode = state.bode[s];
            if (bode.points.empty()) continue;
            int n = static_cast<int>(bode.points.size());
            std::vector<double> f(n), p(n);
            for (int k = 0; k < n; ++k) {
                f[k] = bode.points[k].freq_hz;
                p[k] = bode.points[k].phase_deg;
            }
            ImPlot::PlotLine(system_names[s], f.data(), p.data(), n,
                             ImPlotSpec(ImPlotProp_LineColor, system_colors[s],
                                        ImPlotProp_LineWeight, 2.0f));
        }

        double neg180 = -180.0;
        ImPlot::PlotInfLines("##-180", &neg180, 1,
                             ImPlotSpec(ImPlotProp_LineColor,
                                        ImVec4(1, 1, 1, 0.16f),
                                        ImPlotProp_Flags,
                                        (int)ImPlotInfLinesFlags_Horizontal));

        // PM vertical line (green) on phase plot
        if (margin_sys >= 0 && state.bode[margin_sys].gain_crossover_hz > 0) {
            double x = state.bode[margin_sys].gain_crossover_hz;
            double pm_val = state.bode[margin_sys].phase_margin_deg;
            ImPlot::PlotInfLines("##pm_line", &x, 1,
                                 ImPlotSpec(ImPlotProp_LineColor,
                                            ImVec4(0.39f, 0.78f, 0.39f, 0.47f),
                                            ImPlotProp_LineWeight, 1.0f));
            ImPlot::Annotation(x, -180.0, ImVec4(0.4f, 0.8f, 0.4f, 1),
                               ImVec2(5, -10), true, "PM=%.1f\xc2\xb0", pm_val);
        }

        // GM vertical line (red) on phase plot too
        if (margin_sys >= 0 && state.bode[margin_sys].phase_crossover_hz > 0) {
            double x = state.bode[margin_sys].phase_crossover_hz;
            ImPlot::PlotInfLines("##gm_line_ph", &x, 1,
                                 ImPlotSpec(ImPlotProp_LineColor,
                                            ImVec4(1, 0.39f, 0.39f, 0.47f),
                                            ImPlotProp_LineWeight, 1.0f));
        }

        // Cursor
        if (ImPlot::IsPlotHovered()) {
            double mx = ImPlot::GetPlotMousePos().x;
            ImPlot::PlotInfLines("##cursor", &mx, 1,
                                 ImPlotSpec(ImPlotProp_LineColor,
                                            ImVec4(1, 1, 1, 0.31f)));
            state.bode_cursor_x = mx;
        }
        if (state.bode_cursor_x > 0) {
            double cx = state.bode_cursor_x;
            if (!ImPlot::IsPlotHovered()) {
                ImPlot::PlotInfLines("##cursor", &cx, 1,
                                     ImPlotSpec(ImPlotProp_LineColor,
                                                ImVec4(1, 1, 1, 0.31f)));
            }
        }

        ImPlot::EndPlot();
    }

    // Cursor readout tooltip
    if (state.bode_cursor_x > 0) {
        ImGui::BeginTooltip();
        ImGui::Text("f = %.3f Hz", state.bode_cursor_x);
        for (int s = 0; s < NUM_SYSTEMS; ++s) {
            if (!state.trace_visible[s] || !state.system_valid[s]) continue;
            if (state.bode[s].points.empty()) continue;
            double mag, ph;
            interpBode(state.bode[s], state.bode_cursor_x, mag, ph);
            ImGui::TextColored(system_colors[s], "%s: %.1f dB, %.1f\xc2\xb0",
                               system_names[s], mag, ph);
        }
        ImGui::EndTooltip();
        state.bode_cursor_x = -1.0;
    }
}
}  // anonymous namespace

void drawBodePanel(AppState& state) {
    if (!state.show_bode) return;
    ImGui::Begin("Bode Plot", &state.show_bode);
    drawHelpMarker(
        "Bode Plot\n\n"
        "Magnitude (top) and phase (bottom) vs frequency.\n\n"
        "Gain Margin (GM): How much gain can increase before\n"
        "instability. Measured at phase crossover (-180 deg).\n"
        "GM > 6 dB is typical for robust design.\n\n"
        "Phase Margin (PM): How much phase lag can increase\n"
        "before instability. Measured at gain crossover (0 dB).\n"
        "PM > 30-60 deg is typical for good design.\n\n"
        "Bandwidth: Freq where gain drops to -3 dB.\n"
        "Higher bandwidth = faster response.\n\n"
        "Roll-off: -20 dB/dec per pole, +20 dB/dec per zero.\n\n"
        "0 dB line: unity gain reference.\n"
        "-180 deg line: critical phase for stability.");
    ImGui::SameLine();

    bool freq_changed = false;
    if (ImGui::BeginTable("##bode_ctrls", 2)) {
        ImGui::TableNextColumn();
        ImGui::TextUnformatted("f min [Hz]"); ImGui::SameLine();
        ImGui::SetNextItemWidth(-FLT_MIN);
        freq_changed |= ImGui::SliderFloat(
            "##fmin", &state.freq_min_hz, 0.001f, 1.0f, "%.3f",
            ImGuiSliderFlags_Logarithmic);
        ImGui::TableNextColumn();
        ImGui::TextUnformatted("f max [Hz]"); ImGui::SameLine();
        ImGui::SetNextItemWidth(-FLT_MIN);
        freq_changed |= ImGui::SliderFloat(
            "##fmax", &state.freq_max_hz, 1.0f, 10000.0f, "%.0f",
            ImGuiSliderFlags_Logarithmic);
        ImGui::EndTable();
    }
    if (freq_changed) state.needs_recompute = true;

    // Drawn BEFORE the plots — they consume GetContentRegionAvail(), and
    // anything after them is clipped out of the panel entirely.  Which is also
    // why this is one fixed line rather than none, one or two: a line that came
    // and went here would resize the plot underneath it.
    drawSuppressedTraces(state.channel_reason, system_names, NUM_SYSTEMS, 2);


    if (state.show_all_channels) {
        const int p = state.plant.outputs();
        const int m = state.plant.inputs();

        // Two grids, by role.  The old single grid sized its table from the
        // PLANT's (p, m) and indexed idx = i*m + j with the plant's m, while
        // the writer filled each bode_grid[s] from that system's own shape —
        // so for any non-square plant the panel read cells belonging to other
        // channels, silently, behind an idx < size() guard.  Each table is now
        // sized by the shape it actually holds.  See issue #9.
        auto drawGrid = [&](const char* id, const char* title, int rows,
                            int cols, const int* systems, int n_systems) {
            bool any = false;
            for (int k = 0; k < n_systems; ++k) {
                const int s = systems[k];
                if (state.trace_visible[s] && state.system_valid[s] &&
                    !state.bode_grid[s].empty()) { any = true; break; }
            }
            if (!any) return;
            ImGui::SeparatorText(title);
            if (!ImGui::BeginTable(id, cols, ImGuiTableFlags_Borders)) return;
            for (int i = 0; i < rows; ++i) {
                ImGui::TableNextRow();
                for (int j = 0; j < cols; ++j) {
                    ImGui::TableNextColumn();
                    ImGui::Text("(%d,%d)", i, j);
                    for (int k = 0; k < n_systems; ++k) {
                        const int s = systems[k];
                        if (!state.trace_visible[s] || !state.system_valid[s])
                            continue;
                        const int idx = i * cols + j;
                        if (idx >= (int)state.bode_grid[s].size()) continue;
                        const auto& bode = state.bode_grid[s][idx];
                        if (bode.points.empty()) continue;
                        const int nn = static_cast<int>(bode.points.size());
                        std::vector<double> ff(nn), mm(nn);
                        for (int q = 0; q < nn; ++q) {
                            ff[q] = bode.points[q].freq_hz;
                            mm[q] = bode.points[q].magnitude_db;
                        }
                        char pid[64];
                        std::snprintf(pid, sizeof(pid), "##bode_%s_%d_%d_%d",
                                      id, i, j, s);
                        if (ImPlot::BeginPlot(pid, ImVec2(-1, 120))) {
                            ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
                            ImPlot::PlotLine(system_names[s], ff.data(),
                                             mm.data(), nn,
                                             ImPlotSpec(ImPlotProp_LineColor,
                                                        system_colors[s]));
                            ImPlot::EndPlot();
                        }
                    }
                }
            }
            ImGui::EndTable();
        };

        // The controller's cells are mostly blank by construction — one nonzero
        // per loop — so the open grid becomes a picture of which channels are
        // actually closed.  systems[2] stays excluded from Bode, as today.
        const int open_systems[2] = {0, 1};
        const int loop_systems[1] = {3};
        drawGrid("open", "Open channels  (y_i <- u_j)", p, m,
                 open_systems, 2);
        drawGrid("loop", "Loop channels  (y_i <- r_j)", p, p,
                 loop_systems, 1);
    } else {
        // Determine which system to show margins for
        int margin_sys = -1;
        if (state.system_valid[3] && state.trace_visible[3]) {
            margin_sys = 3;
        } else if (state.system_valid[0] && state.trace_visible[0]) {
            margin_sys = 0;
        }

        float avail_h = ImGui::GetContentRegionAvail().y;
        float plot_h = avail_h / 2.0f;
        plotBodeSingle("##single", state, plot_h, margin_sys);
    }

    ImGui::End();
}

}  // namespace caliburn
