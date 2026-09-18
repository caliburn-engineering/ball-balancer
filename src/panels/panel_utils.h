// src/panels/panel_utils.h
#pragma once
#include "imgui.h"
#include "implot.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

namespace caliburn {

/// The colour a warning badge is drawn in.  One definition, because a warning
/// that is a different yellow in two panels reads as two different kinds of
/// warning.
inline const ImVec4 kBadgeWarn(1.0f, 0.8f, 0.2f, 1.0f);
inline const ImVec4 kBadgeBad(1.0f, 0.35f, 0.35f, 1.0f);
inline const ImVec4 kBadgeGood(0.3f, 1.0f, 0.3f, 1.0f);

/// A condition-dependent message, appended to the line just drawn rather than
/// given a line of its own.
///
/// **The rule: text that comes and goes must never change the layout's
/// height.**  A warning on its own line pushes every control below it down
/// while it is showing and lets them snap back when it stops — and these
/// conditions flicker at frame rate, so a slider a visitor is reaching for
/// moves out from under the cursor.  Two warnings that can appear
/// independently make it worse: the panel below them has three different
/// resting positions.
///
/// So a badge attaches to a line that is always there.  It costs no vertical
/// space at all, it cannot reflow anything, and the line it rides on is
/// usually the reading the warning qualifies — `error: 178.2 mm  clipped` says
/// in one line that the error is large AND why the loop is not closing it.
///
/// Keep `label` to a word or two: it shares a line, and a badge that wraps has
/// reintroduced exactly the problem it exists to solve.  The sentence explaining
/// what the condition means goes in `explain`, which is shown on hover.
///
/// The other half of the rule, which no helper can enforce: where a message is
/// one of several alternatives, draw ALL branches, so the line exists whichever
/// is true.  A trailing `else` with a disabled "nothing to report" is not
/// clutter, it is what keeps the panel still.  See `drawBallControls`.
inline void statusBadge(const char* label, const ImVec4& colour,
                        const char* explain) {
    ImGui::SameLine();
    ImGui::TextColored(colour, "%s", label);
    if (explain && ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 30.0f);
        ImGui::TextUnformatted(explain);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

/// "Which traces are missing from this plot, and on purpose."  One line,
/// always, whatever is suppressed.
///
/// This one matters more than an ordinary badge line, because of where it
/// sits.  The reasons have to be drawn BEFORE the plot — a plot consumes
/// `GetContentRegionAvail()`, so anything after it is clipped out of the panel
/// entirely — which means a line that comes and goes here does not merely
/// shift the controls below it.  It resizes the plot: the axes rescale and the
/// whole curve moves, for a reason the reader cannot see.
///
/// A suppressed trace must still be distinguishable from a broken one, so the
/// trace's name is shown in the warning colour rather than hidden.  Only the
/// sentence saying WHY moves to the hover, where there is room for it to be a
/// sentence rather than a truncated fragment.
///
/// `skip` is the index of a system that never reports a reason; pass -1 for
/// none.
inline void drawSuppressedTraces(const std::string* reasons,
                                 const char* const* names,
                                 int n, int skip) {
    ImGui::TextDisabled("suppressed:");
    bool any = false;
    for (int s = 0; s < n; ++s) {
        if (s == skip || reasons[s].empty()) continue;
        statusBadge(names[s], ImVec4(1.0f, 0.75f, 0.35f, 1.0f),
                    reasons[s].c_str());
        any = true;
    }
    if (!any) {
        ImGui::SameLine();
        ImGui::TextDisabled("none");
    }
}

inline void drawHelpMarker(const char* desc) {
    ImGui::SameLine();
    ImGui::TextDisabled("(?)");
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort)) {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 30.0f);
        ImGui::TextUnformatted(desc);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

// Compute a nice grid step for a given range, targeting ~target intervals
inline double gridNiceStep(double range, int target = 5) {
    if (range <= 0) return 1.0;
    double rough = range / target;
    double e = std::floor(std::log10(rough));
    double b = std::pow(10.0, e);
    double f = rough / b;
    if (f < 1.5) return b;
    if (f < 3.5) return 2.0 * b;
    if (f < 7.5) return 5.0 * b;
    return 10.0 * b;
}

// S-plane grid: wn circles + zeta radial lines (for pole-zero / root locus)
// Call after SetupAxes/SetupAxisLimits, before data plotting
inline void drawSPlaneGrid() {
    auto lim = ImPlot::GetPlotLimits();
    double max_r = std::max({std::abs(lim.X.Min), std::abs(lim.X.Max),
                              std::abs(lim.Y.Min), std::abs(lim.Y.Max)});

    ImVec4 gc(1, 1, 1, 0.07f);
    ImU32 lc = IM_COL32(255, 255, 255, 76);
    ImDrawList* dl = ImPlot::GetPlotDrawList();

    // wn concentric circles
    double step = gridNiceStep(max_r);
    constexpr int N = 128;
    int cid = 0;
    for (double wn = step; wn < max_r + step * 0.01 && cid < 10; wn += step) {
        double cx[N + 1], cy[N + 1];
        for (int i = 0; i <= N; ++i) {
            double th = 2.0 * M_PI * i / N;
            cx[i] = wn * std::cos(th);
            cy[i] = wn * std::sin(th);
        }
        char l[32];
        std::snprintf(l, sizeof(l), "##sgwn%d", cid++);
        ImPlot::PlotLine(l, cx, cy, N + 1,
                         ImPlotSpec(ImPlotProp_LineColor, gc));
        // wn label on negative real axis
        char txt[16];
        std::snprintf(txt, sizeof(txt), "%.4g", wn);
        ImVec2 p = ImPlot::PlotToPixels(ImPlotPoint(-wn, 0));
        dl->AddText(ImVec2(p.x + 4, p.y - 14), lc, txt);
    }

    // zeta radial lines (every 10 deg from negative real axis, LHP only)
    for (int deg = 10; deg < 90; deg += 10) {
        double th = deg * M_PI / 180.0;
        double zeta = std::cos(th);
        double dx = -std::cos(th);
        double dy = std::sin(th);

        // Upper line (Q2)
        double xu[2] = {0, dx * max_r};
        double yu[2] = {0, dy * max_r};
        char l[32];
        std::snprintf(l, sizeof(l), "##sgzu%d", deg);
        ImPlot::PlotLine(l, xu, yu, 2,
                         ImPlotSpec(ImPlotProp_LineColor, gc));

        // Lower line (Q3, mirror)
        double yl[2] = {0, -dy * max_r};
        std::snprintf(l, sizeof(l), "##sgzl%d", deg);
        ImPlot::PlotLine(l, xu, yl, 2,
                         ImPlotSpec(ImPlotProp_LineColor, gc));

        // zeta label near top edge
        double t_y = (dy > 0.01) ? std::abs(lim.Y.Max / dy) : 1e30;
        double t_x = (dx < -0.01) ? std::abs(lim.X.Min / dx) : 1e30;
        double t = std::min(t_y, t_x) * 0.92;
        double lx = dx * t, ly = dy * t;
        if (lx >= lim.X.Min && lx <= lim.X.Max &&
            ly >= lim.Y.Min && ly <= lim.Y.Max) {
            char txt[8];
            std::snprintf(txt, sizeof(txt), "%.2f", zeta);
            ImVec2 p = ImPlot::PlotToPixels(ImPlotPoint(lx, ly));
            dl->AddText(ImVec2(p.x + 3, p.y - 12), lc, txt);
        }
    }
}

// Polar grid: concentric circles + radial lines centered at origin
// For Nyquist and other complex-plane plots
inline void drawPolarGrid() {
    auto lim = ImPlot::GetPlotLimits();
    double max_r = std::max({std::abs(lim.X.Min), std::abs(lim.X.Max),
                              std::abs(lim.Y.Min), std::abs(lim.Y.Max)});

    ImVec4 gc(1, 1, 1, 0.07f);

    // Concentric circles
    double step = gridNiceStep(max_r);
    constexpr int N = 128;
    int cid = 0;
    for (double r = step; r < max_r + step * 0.01 && cid < 10; r += step) {
        double cx[N + 1], cy[N + 1];
        for (int i = 0; i <= N; ++i) {
            double th = 2.0 * M_PI * i / N;
            cx[i] = r * std::cos(th);
            cy[i] = r * std::sin(th);
        }
        char l[32];
        std::snprintf(l, sizeof(l), "##pgc%d", cid++);
        ImPlot::PlotLine(l, cx, cy, N + 1,
                         ImPlotSpec(ImPlotProp_LineColor, gc));
    }

    // Radial lines every 30 deg (full lines through origin)
    int lid = 0;
    for (int deg = 0; deg < 180; deg += 30) {
        double th = deg * M_PI / 180.0;
        double dx = std::cos(th) * max_r;
        double dy = std::sin(th) * max_r;
        double x[2] = {-dx, dx};
        double y[2] = {-dy, dy};
        char l[32];
        std::snprintf(l, sizeof(l), "##pgr%d", lid++);
        ImPlot::PlotLine(l, x, y, 2,
                         ImPlotSpec(ImPlotProp_LineColor, gc));
    }
}

}  // namespace caliburn
