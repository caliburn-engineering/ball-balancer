#include "plot_panel.h"
#include <cstdio>

static constexpr float MIN_PLOT_HEIGHT = 120.0f;
static constexpr float Y_MARGIN = 0.05f;  // 5% margin on Y-limits
static const ImVec4 CURSOR_COLOR  = ImVec4(1.0f, 1.0f, 1.0f, 0.5f);
static const ImVec4 MARKER_COLOR  = ImVec4(1.0f, 0.8f, 0.2f, 0.7f);

void draw_time_series_panel(
    const char* window_title,
    PlotState& state,
    std::vector<PlotConfig>& plots)
{
    ImGui::Begin(window_title);

    // --- Controls row ---
    ImGui::SliderFloat("Time window [s]", &state.time_window, 2.0f, 60.0f, "%.0f");

    if (ImGui::Button(state.paused ? "Resume" : "Pause")) {
        state.paused = !state.paused;
        if (state.paused) {
            // Initialize zoom range to current visible window
            state.pause_t_max = state.latest_time();
            state.pause_t_min = state.pause_t_max - state.time_window;
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Clear Data")) {
        state.clear();
        for (auto& pc : plots) {
            for (auto* s : pc.series) {
                s->data.clear();
            }
        }
    }
    ImGui::SameLine();
    if (!state.markers.empty()) {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "Clear Markers (%d)", (int)state.markers.size());
        if (ImGui::Button(buf)) {
            state.markers.clear();
        }
    }

    // --- Time-zoom sliders (when paused) ---
    float t_min, t_max;
    if (state.paused && state.count > 0) {
        float earliest = state.earliest_time();
        float latest = state.latest_time();
        ImGui::SliderFloat("Zoom start", &state.pause_t_min, earliest, latest, "%.1f s");
        ImGui::SliderFloat("Zoom end",   &state.pause_t_max, earliest, latest, "%.1f s");
        if (state.pause_t_min > state.pause_t_max)
            state.pause_t_min = state.pause_t_max - 0.1f;
        t_min = state.pause_t_min;
        t_max = state.pause_t_max;
    } else {
        t_max = state.latest_time();
        t_min = t_max - state.time_window;
    }

    // --- Plot visibility toggles ---
    for (int pi = 0; pi < (int)plots.size(); ++pi) {
        if (pi > 0) ImGui::SameLine();
        // The toggle carries the plot's title, and so does its ImPlot::BeginPlot
        // below — same label, same window, same ID.  ImGui flags that as a
        // conflict and the two items fight over hover and activation state.
        ImGui::PushID(pi);
        bool active = plots[pi].visible;
        if (!active) {
            ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.2f, 0.2f, 0.2f, 1.0f));
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.5f, 0.5f, 0.5f, 1.0f));
        }
        if (ImGui::SmallButton(plots[pi].title.c_str())) {
            plots[pi].visible = !plots[pi].visible;
        }
        if (!active) {
            ImGui::PopStyleColor(2);
        }
        ImGui::PopID();
    }
    ImGui::Separator();

    // --- Compute plot heights (visible plots only) ---
    int n_visible = 0;
    for (auto& pc : plots) if (pc.visible) n_visible++;
    float avail = ImGui::GetContentRegionAvail().y;
    float plot_h = std::max(MIN_PLOT_HEIGHT, (avail - 10.0f) / std::max(n_visible, 1));

    // If plots won't fit at minimum height, enable scrolling
    if (n_visible * MIN_PLOT_HEIGHT > avail) {
        plot_h = MIN_PLOT_HEIGHT;
        ImGui::BeginChild("PlotScroll", ImVec2(0, 0), false, ImGuiWindowFlags_None);
    }

    // --- Cursor tracking ---
    // Reset cursor; it'll be set by whichever plot is hovered
    state.cursor_active = false;

    ImPlotSpec scroll_spec;
    scroll_spec.Offset = state.offset;

    int n_plots = (int)plots.size();
    for (int pi = 0; pi < n_plots; ++pi) {
        auto& pc = plots[pi];
        if (!pc.visible) continue;

        // Compute Y-limits from visible data (only visible series)
        float y_lo = std::numeric_limits<float>::max();
        float y_hi = std::numeric_limits<float>::lowest();
        for (auto* s : pc.series) {
            // Check if this series is hidden in ImPlot
            // We compute range for all and let ImPlot's legend toggle handle visual hiding
            // But we also track per-series to do manual Y-fit
            auto [lo, hi] = series_range_in_window(*s, state, t_min, t_max);
            y_lo = std::min(y_lo, lo);
            y_hi = std::max(y_hi, hi);
        }
        if (y_lo >= y_hi) { y_lo = -1; y_hi = 1; }
        float margin = (y_hi - y_lo) * Y_MARGIN;
        if (margin < 0.01f) margin = 0.5f;

        if (ImPlot::BeginPlot(pc.title.c_str(), ImVec2(-1, plot_h))) {
            // No axis titles.  Every plot in this panel is a time series and
            // every one of them said "t [s]" underneath, which is a label the
            // reader has already worked out from the first plot — and it costs
            // a row of height on each of eight stacked plots.  The tick
            // numbers stay; it is the title that is redundant.
            ImPlot::SetupAxes(nullptr, nullptr);
            ImPlot::SetupAxisLimits(ImAxis_X1, t_min, t_max, ImPlotCond_Always);
            ImPlot::SetupAxisLimits(ImAxis_Y1, y_lo - margin, y_hi + margin, ImPlotCond_Always);
            ImPlot::SetupLegend(ImPlotLocation_East, ImPlotLegendFlags_Outside);

            // --- Plot data ---
            for (auto* s : pc.series) {
                if (state.count > 0 && (int)s->data.size() >= state.count) {
                    ImPlot::PlotLine(s->label.c_str(),
                                     state.time.data(), s->data.data(),
                                     state.count, scroll_spec);
                }
            }

            // --- Synchronized cursor ---
            if (ImPlot::IsPlotHovered()) {
                ImPlotPoint mouse = ImPlot::GetPlotMousePos();
                state.cursor_active = true;
                state.cursor_time = mouse.x;

                // Show tooltip with values
                ImGui::BeginTooltip();
                ImGui::Text("t = %.3f s", mouse.x);
                for (auto* s : pc.series) {
                    float val = interpolate_at_time(*s, state, mouse.x);
                    ImGui::Text("%s: %.4f", s->label.c_str(), val);
                }
                ImGui::EndTooltip();

                // Click to place marker
                if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                    PlotMarker m;
                    m.time = mouse.x;
                    m.source_plot = pi;
                    for (auto* s : pc.series) {
                        m.values.push_back(interpolate_at_time(*s, state, mouse.x));
                        m.labels.push_back(s->label);
                    }
                    state.markers.push_back(m);
                }
            }

            // Draw cursor line on this plot (if any plot is hovered)
            if (state.cursor_active) {
                double cx = state.cursor_time;
                ImPlot::PlotInfLines("##cursor", &cx, 1,
                    ImPlotSpec(ImPlotProp_LineColor, CURSOR_COLOR,
                               ImPlotProp_LineWeight, 1.0f,
                               ImPlotProp_Flags, (int)ImPlotItemFlags_NoLegend));
            }

            // Draw persistent markers
            for (const auto& mk : state.markers) {
                double mx = mk.time;
                ImPlot::PlotInfLines("##marker", &mx, 1,
                    ImPlotSpec(ImPlotProp_LineColor, MARKER_COLOR,
                               ImPlotProp_LineWeight, 1.5f,
                               ImPlotProp_Flags, (int)ImPlotItemFlags_NoLegend));

                // Show pinned tooltip on the source plot
                if (mk.source_plot == pi) {
                    // Find Y position for annotation (use first series value)
                    float y_pos = mk.values.empty() ? 0.0f : mk.values[0];
                    char annot_buf[128] = "";
                    int pos = 0;
                    for (size_t si = 0; si < mk.values.size() && pos < 120; ++si) {
                        pos += std::snprintf(annot_buf + pos, sizeof(annot_buf) - pos,
                            "%s%.4f", si > 0 ? "\n" : "", mk.values[si]);
                    }
                    ImPlot::Annotation(mx, y_pos, MARKER_COLOR,
                                       ImVec2(5, -5), true, "%s", annot_buf);
                }
            }

            ImPlot::EndPlot();
        }
    }

    // End scroll region if active
    if (n_visible * MIN_PLOT_HEIGHT > avail) {
        ImGui::EndChild();
    }

    ImGui::End();
}
