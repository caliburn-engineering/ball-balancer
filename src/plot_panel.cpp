#include "plot_panel.h"

#include <cstdio>
#include <string>
#include <utility>

static constexpr float MIN_PLOT_HEIGHT = 120.0f;
static constexpr float Y_MARGIN = 0.05f;  // 5% margin on Y-limits
static const ImVec4 CURSOR_COLOR  = ImVec4(1.0f, 1.0f, 1.0f, 0.5f);
static const ImVec4 MARKER_COLOR  = ImVec4(1.0f, 0.8f, 0.2f, 0.7f);
// How near the cursor has to be, in PIXELS, for a marker's note to appear in
// the hover.
//
// Pixels rather than seconds or a fraction of the time window, because "near
// the marker" is a judgement the reader makes with their eyes and their mouse.
// Measured in the browser: this panel docks to about 270 px of plot width, and
// a 10 s window there is 23 px per second — so a tolerance of 2% of the window
// was +/-4 px, which is narrower than the marker's own line is easy to aim at,
// and it got narrower still every time somebody zoomed out.
static constexpr float MARKER_HOVER_PIXELS = 10.0f;

// The most markers kept, oldest dropped first.
//
// Markers used to be one-per-click, so the list was bounded by how long
// somebody was willing to keep clicking.  A lost ball places one too, and with
// Auto-reset on nobody is clicking: a plate left tilted loses the ball, has it
// put back, and loses it again, for as long as the tab is open.  A cap is what
// stops an unattended demo growing a list forever, and dropping the oldest is
// right for both producers — the interesting marker is the recent one.
static constexpr size_t MAX_MARKERS = 64;

void place_marker(PlotState& state, const std::vector<PlotConfig>& plots,
                  int pi, double t, std::string note)
{
    if (pi < 0 || pi >= (int)plots.size()) return;
    PlotMarker m;
    m.time = t;
    m.source_plot = pi;
    m.note = std::move(note);
    for (auto* s : plots[pi].series) {
        m.values.push_back(interpolate_at_time(*s, state, t));
        m.labels.push_back(s->label);
    }
    state.markers.push_back(m);
    if (state.markers.size() > MAX_MARKERS)
        state.markers.erase(state.markers.begin());
}

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
                // And what a marker near the cursor has to say, on any plot
                // rather than only on the one it was placed from — its line is
                // drawn on all of them, so a reader can meet it on any of them.
                //
                // The annotation on the plot carries the first line only,
                // because it is painted over the data whether anyone wants it
                // or not.  The rest is here, where it costs nothing until
                // somebody asks.  This is the tooltip #33 puts the raw flags in.
                const float mouse_px = ImGui::GetMousePos().x;
                for (const auto& mk : state.markers) {
                    if (mk.note.empty()) continue;
                    const float marker_px =
                        ImPlot::PlotToPixels(mk.time, 0.0).x;
                    if (std::abs(marker_px - mouse_px) > MARKER_HOVER_PIXELS)
                        continue;
                    ImGui::Separator();
                    ImGui::PushTextWrapPos(ImGui::GetFontSize() * 28.0f);
                    ImGui::TextUnformatted(mk.note.c_str());
                    ImGui::PopTextWrapPos();
                }
                ImGui::EndTooltip();

                // Click to place marker
                if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                    place_marker(state, plots, pi, mouse.x);
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
                    // A std::string rather than a 128-byte buffer: a note can
                    // be a sentence, and the old one silently stopped at 120
                    // characters — which for an explanation is the same as
                    // being wrong, since the half that gets cut is the half
                    // that says what it means.
                    // A marker with a note shows the note, and a marker
                    // without one shows the numbers.  Not both: these plots
                    // dock to about 90 px each, which is four lines of text
                    // including the axis, so an annotation carrying a two-line
                    // note AND two values overflows the plot rect and ImPlot
                    // clamps it — cutting off the top line, which is the line
                    // that says what happened.  Measured in the browser at
                    // 1440x900, where the panel is at its roomiest.
                    //
                    // The note wins because it is the part that cannot be got
                    // any other way: the values at that instant are still on
                    // the hover, and they are also just where the curve is.
                    // Up to the blank line, per `PlotMarker::note`.
                    std::string annot;
                    if (!mk.note.empty()) {
                        annot = mk.note.substr(0, mk.note.find("\n\n"));
                    } else {
                        char num[32];
                        for (size_t si = 0; si < mk.values.size(); ++si) {
                            std::snprintf(num, sizeof(num), "%.4f", mk.values[si]);
                            if (!annot.empty()) annot += '\n';
                            annot += num;
                        }
                    }
                    ImPlot::Annotation(mx, y_pos, MARKER_COLOR,
                                       ImVec2(5, -5), true, "%s", annot.c_str());
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
