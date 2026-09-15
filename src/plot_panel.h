#pragma once

#include "implot.h"
#include "imgui.h"
#include <vector>
#include <string>
#include <array>
#include <cmath>
#include <algorithm>
#include <limits>

/// A single data series in a scrolling circular buffer
struct TimeSeries {
    std::string label;       // Max 2 chars for legend display
    std::vector<float> data;
    int max_size;

    explicit TimeSeries(const char* label, int max_size = 4000)
        : label(label), max_size(max_size) {
        data.reserve(max_size);
    }
};

/// A persistent marker: placed by clicking on a plot, or by the simulation
/// when something happened that is worth being able to look back at.
struct PlotMarker {
    double time;
    int source_plot;                // Index of the plot the marker was placed on
    std::vector<float> values;      // Values at marker time for the source plot's series
    std::vector<std::string> labels; // Series labels

    /// Why this marker is here, when something other than a click put it there.
    ///
    /// Everything up to the first BLANK LINE is painted on the plot as an
    /// annotation, whether the reader wants it or not, so it has to be short —
    /// what happened and the evidence, a couple of labels' worth.  What follows
    /// the blank line is hover-only, where there is room for a sentence.
    ///
    /// **The split is where it is because a phone has no hover.**  Putting the
    /// evidence behind one would mean the full picture survives for a reader
    /// with a mouse and for nobody else, which is not what "survives for anyone
    /// who looks" means.  So the annotation carries what a reader must be able
    /// to check, and the hover carries what merely helps.
    ///
    /// Empty for a marker the visitor placed themselves, which needs no
    /// explaining to the person who placed it.
    std::string note;
};

/// Shared time buffer + cursor state across all time-series plots
struct PlotState {
    // Time buffer (shared across all series)
    std::vector<float> time;
    int max_size = 4000;
    int offset = 0;
    int count = 0;

    // Cursor synchronization
    bool cursor_active = false;
    double cursor_time = 0.0;

    // Persistent markers
    std::vector<PlotMarker> markers;

    // Simulation pause
    bool paused = false;
    float pause_t_min = 0.0f;   // Time-zoom range (when paused)
    float pause_t_max = 0.0f;

    // Settings
    float time_window = 10.0f;  // Visible window [s]

    PlotState(int max_size = 4000) : max_size(max_size) {
        time.reserve(max_size);
    }

    void push_time(float t) {
        if (paused) return;
        if (count < max_size) {
            time.push_back(t);
            count++;
        } else {
            time[offset] = t;
            offset = (offset + 1) % max_size;
        }
    }

    float latest_time() const {
        if (count == 0) return 0;
        return time[(offset + count - 1) % max_size];
    }

    float earliest_time() const {
        if (count == 0) return 0;
        return time[offset % max_size];
    }

    void clear() {
        time.clear();
        offset = 0;
        count = 0;
    }
};

/// Push a value to a series (must be called AFTER PlotState::push_time for the same frame)
inline void push_series(TimeSeries& s, float value, const PlotState& state) {
    if (state.paused) return;
    // After push_time: count was incremented (or offset advanced).
    // The newest slot is at index (offset + count - 1) % max_size.
    int idx = (state.offset + state.count - 1) % state.max_size;
    if ((int)s.data.size() < state.max_size) {
        s.data.push_back(value);
    } else {
        s.data[idx] = value;
    }
}

/// A single plot containing multiple time series
struct PlotConfig {
    std::string title;
    std::vector<TimeSeries*> series;  // Non-owning pointers
    int plot_index;                    // For marker tracking
    bool visible = true;               // Toggle via UI button
};

/// Compute min/max of a series within a time range, respecting circular buffer
inline std::pair<float, float> series_range_in_window(
    const TimeSeries& s, const PlotState& state, float t_min, float t_max) {
    float lo = std::numeric_limits<float>::max();
    float hi = std::numeric_limits<float>::lowest();
    for (int i = 0; i < state.count; ++i) {
        int idx = (state.offset + i) % state.max_size;
        float t = state.time[idx];
        if (t >= t_min && t <= t_max) {
            float v = s.data[idx];
            if (std::isfinite(v)) {
                lo = std::min(lo, v);
                hi = std::max(hi, v);
            }
        }
    }
    if (lo > hi) { lo = -1; hi = 1; }  // No data in range
    return {lo, hi};
}

/// Interpolate a series value at a given time
inline float interpolate_at_time(const TimeSeries& s, const PlotState& state, double t) {
    // Find the two samples bracketing t
    int best = -1;
    float best_dt = std::numeric_limits<float>::max();
    for (int i = 0; i < state.count; ++i) {
        int idx = (state.offset + i) % state.max_size;
        float dt = std::abs(state.time[idx] - (float)t);
        if (dt < best_dt) {
            best_dt = dt;
            best = idx;
        }
    }
    if (best < 0) return 0;
    return s.data[best];
}

/// Put a marker at time `t` on plot `pi`, capturing that plot's series values.
///
/// The one place a marker is constructed.  The click path and the simulation
/// both go through it, because a marker built at a second site is a marker that
/// can be built with its labels and its values out of step — and a marker is
/// read long after the frame that made it, so nothing would ever say so.
void place_marker(PlotState& state, const std::vector<PlotConfig>& plots,
                  int pi, double t, std::string note = {});

/// Draw the plots panel with all time-series invariants
/// Returns true if any plot was clicked (marker placed)
void draw_time_series_panel(
    const char* window_title,
    PlotState& state,
    std::vector<PlotConfig>& plots);
