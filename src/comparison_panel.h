#pragma once

#include "comparison_sim.h"
#include "ball_plant_linear.h"
#include "plot_panel.h"

namespace caliburn {

class ComparisonPanel {
public:
    ComparisonPanel();

    /// Draw the comparison panel as an ImGui window.
    /// Call once per frame from the main loop.
    void draw();

private:
    // Input configuration
    int input_type_ = 0;  // 0=Step, 1=Ramp
    float alpha_amplitude_ = 0.05f;
    float beta_amplitude_ = 0.0f;
    float duration_ = 5.0f;

    // Results
    bool has_result_ = false;
    ComparisonResult result_;
    ValidationResult validation_;

    // Plotting
    PlotState plot_state_;

    // Position traces: xN, yN, xL, yL
    TimeSeries ts_xN_{"xN"}, ts_yN_{"yN"}, ts_xL_{"xL"}, ts_yL_{"yL"};
    // Velocity traces: vxN, vyN, vxL, vyL
    TimeSeries ts_vxN_{"vN"}, ts_vyN_{"wN"}, ts_vxL_{"vL"}, ts_vyL_{"wL"};
    // Input traces
    TimeSeries ts_alpha_{"al"}, ts_beta_{"be"};
    // Error traces
    TimeSeries ts_ex_{"ex"}, ts_ey_{"ey"};

    std::vector<PlotConfig> plots_;

    // Divergence threshold
    float diverge_threshold_ = 0.01f;

    void runComparison();
    void populatePlots();
    void setupPlots();
    void drawValidationTable();
};

}  // namespace caliburn
