#include "comparison_panel.h"

#include "imgui.h"
#include "implot.h"

#include <cmath>
#include <cstdio>

namespace caliburn {

ComparisonPanel::ComparisonPanel() {
    setupPlots();
}

void ComparisonPanel::setupPlots() {
    plots_.clear();
    plots_.push_back({"Position [m]", {&ts_xN_, &ts_yN_, &ts_xL_, &ts_yL_}, 0, true});
    plots_.push_back({"Velocity [m/s]", {&ts_vxN_, &ts_vyN_, &ts_vxL_, &ts_vyL_}, 1, true});
    plots_.push_back({"Input [rad]", {&ts_alpha_, &ts_beta_}, 2, true});
    plots_.push_back({"Error [m]", {&ts_ex_, &ts_ey_}, 3, false});  // collapsed by default
}

void ComparisonPanel::runComparison() {
    auto lin = ballBalancerLinearModel();
    auto f_nl = ballBalancerNonlinearFn();

    // Build input signal
    Eigen::VectorXd amp(2);
    amp << static_cast<double>(alpha_amplitude_),
           static_cast<double>(beta_amplitude_);

    InputSignal input;
    if (input_type_ == 0) {
        input = StepInput{amp, 0.5};  // step at t=0.5s
    } else {
        input = RampInput{amp, 0.5};
    }

    double dt = 0.001;
    result_ = compareModels(f_nl, lin, Eigen::VectorXd::Zero(4),
                            input, dt, static_cast<double>(duration_));
    has_result_ = true;

    // Run validation
    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(4);
    Eigen::VectorXd u0 = Eigen::VectorXd::Zero(2);
    validation_ = validate(lin, f_nl, x0, u0);

    populatePlots();
}

void ComparisonPanel::populatePlots() {
    // Clear existing data
    plot_state_.clear();
    auto clear = [](TimeSeries& s) { s.data.clear(); };
    clear(ts_xN_); clear(ts_yN_); clear(ts_xL_); clear(ts_yL_);
    clear(ts_vxN_); clear(ts_vyN_); clear(ts_vxL_); clear(ts_vyL_);
    clear(ts_alpha_); clear(ts_beta_);
    clear(ts_ex_); clear(ts_ey_);

    // Downsample for plotting: take every Nth sample to keep buffer manageable
    int total = static_cast<int>(result_.time.size());
    int stride = std::max(1, total / 3600);

    // Temporarily unpause to allow push
    plot_state_.paused = false;

    for (int i = 0; i < total; i += stride) {
        float t = static_cast<float>(result_.time[i]);
        plot_state_.push_time(t);

        const auto& xnl = result_.nonlinear_states[i];
        const auto& xln = result_.linear_states[i];
        const auto& u = result_.inputs[i];

        push_series(ts_xN_, static_cast<float>(xnl(0)), plot_state_);
        push_series(ts_yN_, static_cast<float>(xnl(1)), plot_state_);
        push_series(ts_xL_, static_cast<float>(xln(0)), plot_state_);
        push_series(ts_yL_, static_cast<float>(xln(1)), plot_state_);

        push_series(ts_vxN_, static_cast<float>(xnl(2)), plot_state_);
        push_series(ts_vyN_, static_cast<float>(xnl(3)), plot_state_);
        push_series(ts_vxL_, static_cast<float>(xln(2)), plot_state_);
        push_series(ts_vyL_, static_cast<float>(xln(3)), plot_state_);

        push_series(ts_alpha_, static_cast<float>(u(0)), plot_state_);
        push_series(ts_beta_, static_cast<float>(u(1)), plot_state_);

        push_series(ts_ex_, static_cast<float>(std::abs(xnl(0) - xln(0))), plot_state_);
        push_series(ts_ey_, static_cast<float>(std::abs(xnl(1) - xln(1))), plot_state_);
    }

    // Pause after populating — this is static comparison data
    plot_state_.paused = true;
    plot_state_.pause_t_min = static_cast<float>(result_.time.front());
    plot_state_.pause_t_max = static_cast<float>(result_.time.back());
    plot_state_.time_window = plot_state_.pause_t_max - plot_state_.pause_t_min;
}

void ComparisonPanel::drawValidationTable() {
    if (!has_result_) return;

    ImGui::Separator();
    ImGui::Text("Linearization Validation");

    ImVec4 color = validation_.pass
        ? ImVec4(0.2f, 0.9f, 0.2f, 1.0f)
        : ImVec4(0.9f, 0.2f, 0.2f, 1.0f);
    ImGui::TextColored(color, validation_.pass ? "PASS" : "FAIL");
    ImGui::SameLine();
    ImGui::Text("max|A err|=%.2e  max|B err|=%.2e",
                validation_.max_A_error, validation_.max_B_error);

    if (ImGui::TreeNode("A error (element-wise)")) {
        for (int r = 0; r < validation_.A_error.rows(); ++r) {
            for (int c = 0; c < validation_.A_error.cols(); ++c) {
                if (c > 0) ImGui::SameLine();
                ImGui::Text("%.2e", validation_.A_error(r, c));
            }
        }
        ImGui::TreePop();
    }

    if (ImGui::TreeNode("B error (element-wise)")) {
        for (int r = 0; r < validation_.B_error.rows(); ++r) {
            for (int c = 0; c < validation_.B_error.cols(); ++c) {
                if (c > 0) ImGui::SameLine();
                ImGui::Text("%.2e", validation_.B_error(r, c));
            }
        }
        ImGui::TreePop();
    }
}

void ComparisonPanel::draw() {
    if (!ImGui::Begin("Model Comparison")) {
        ImGui::End();
        return;
    }

    // --- Controls ---
    ImGui::Text("Input Signal");
    ImGui::RadioButton("Step", &input_type_, 0);
    ImGui::SameLine();
    ImGui::RadioButton("Ramp", &input_type_, 1);

    ImGui::SliderFloat(u8"\u03b1 amplitude [rad]", &alpha_amplitude_, 0.0f, 1.0f, "%.3f");
    ImGui::SliderFloat(u8"\u03b2 amplitude [rad]", &beta_amplitude_, 0.0f, 1.0f, "%.3f");
    ImGui::SliderFloat("Duration [s]", &duration_, 1.0f, 30.0f, "%.1f");
    ImGui::SliderFloat("Diverge threshold", &diverge_threshold_, 0.001f, 0.1f, "%.3f");

    if (ImGui::Button("Run Comparison")) {
        runComparison();
    }
    ImGui::SameLine();
    if (ImGui::Button("Reset")) {
        has_result_ = false;
        plot_state_.clear();
        setupPlots();
    }

    // --- Validation results ---
    drawValidationTable();

    // --- Plots ---
    if (has_result_) {
        ImGui::Separator();
        draw_time_series_panel("##comparison_plots", plot_state_, plots_);
    }

    ImGui::End();
}

}  // namespace caliburn
