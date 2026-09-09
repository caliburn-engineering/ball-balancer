// src/plate_view.cpp
#include "plate_view.h"

// No GL header here: this file makes no GL calls of its own, and renderer.h
// (via plate_view.h) is the one place that decides between GLAD and Emscripten's
// ES 3.0 headers.
#include <GLFW/glfw3.h>
#include "imgui.h"
#include "panels/panel_utils.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace caliburn {
namespace {

constexpr double kDeg = M_PI / 180.0;
constexpr int kBufSize = 4000;

constexpr double kBallRadius = kPlateBall.radius;  // [m], used when drawing

namespace col {
    constexpr std::array<float,4> ground    = {0.5f, 0.5f, 0.5f, 1.0f};
    constexpr std::array<float,4> ground_f  = {0.3f, 0.3f, 0.3f, 0.12f};
    constexpr std::array<float,4> table_c   = {0.3f, 0.7f, 1.0f, 1.0f};
    constexpr std::array<float,4> table_f   = {0.2f, 0.5f, 0.8f, 0.12f};
    constexpr std::array<float,4> leg_L1    = {1.0f, 0.8f, 0.2f, 1.0f};
    constexpr std::array<float,4> leg_L2    = {0.2f, 1.0f, 0.4f, 1.0f};
    constexpr std::array<float,4> joint_sph = {1.0f, 0.3f, 0.3f, 1.0f};
    constexpr std::array<float,4> joint_rev = {0.3f, 1.0f, 0.3f, 1.0f};
    constexpr std::array<float,4> grid      = {0.25f, 0.25f, 0.25f, 0.5f};
    constexpr std::array<float,4> ball      = {0.98f, 0.45f, 0.09f, 1.0f};
    constexpr std::array<float,4> ball_f    = {0.98f, 0.45f, 0.09f, 0.45f};
    constexpr std::array<float,4> ball_off  = {0.9f, 0.2f, 0.2f, 1.0f};
    // Airborne: a colour of its own, because line weight carries no
    // information under WebGL2 (see CONTEXT.md, "Plate view") and this is a
    // state the visitor is meant to notice.
    constexpr std::array<float,4> ball_air  = {1.0f, 0.85f, 0.25f, 1.0f};
    constexpr std::array<float,4> path      = {0.40f, 0.85f, 0.55f, 0.7f};
    constexpr std::array<float,4> setpoint  = {0.40f, 0.95f, 0.60f, 1.0f};
}

TableParams defaultTableParams() {
    TableParams p;
    p.R_ground  = 0.300;
    p.R_table   = 0.300;
    p.L1        = 0.150;
    p.L2        = 0.150;
    p.alpha_min = 10.0 * kDeg;
    p.alpha_max = 80.0 * kDeg;
    return p;
}

OrbitCamera defaultCamera() {
    OrbitCamera cam;
    cam.azimuth = 45.0f;
    cam.elevation = 30.0f;
    cam.distance = 1.0f;
    cam.target = Eigen::Vector3f(0, 0, 0.12f);
    return cam;
}

// The plate's own gravity.  Fixed here and compared against the model's `g`
// rather than followed: the dynamics object is built once, and a design solved
// on the moon must be refused, not quietly run on Earth.
constexpr double kGravity = 9.81;

// Where this plate's legs sit at rest: what a reset returns to, and what the
// Home button commands.  Not yet the only 45 in this file — the animation
// swings about it and the command array is initialised to it, both of which
// are literals still.
constexpr double kHomeLegDeg = 45.0;
constexpr double kHomeLegRad = kHomeLegDeg * kDeg;

// The scroll wheel is the one input this class cannot read by polling: GLFW
// only delivers it as an event.  One instance drives the app, so a file-scope
// pointer is enough and keeps the callback free of captures.
PlateView* g_plate_view = nullptr;

void scrollCallback(GLFWwindow*, double, double yoffset) {
    if (!g_plate_view || ImGui::GetIO().WantCaptureMouse) return;
    OrbitCamera& cam = g_plate_view->camera();
    cam.distance *= (1.0f - 0.1f * static_cast<float>(yoffset));
    cam.distance = std::clamp(cam.distance, 0.2f, 3.0f);
}

}  // namespace

PlateView::PlateView()
    // plate_ is declared first, so its params are live here: the ball's plate
      // and the drawn plate are one object, not two that have to be kept in
      // step.  It is also exactly what every test harness runs against.
    : plate_(defaultTableParams(), kGravity),
      plot_state_(kBufSize),
      s_phi_("\xcf\x86", kBufSize),                 // phi
      s_theta_("\xce\xb8", kBufSize),               // theta
      s_a0_("\xce\xb1\xe2\x82\x80", kBufSize),      // alpha_0
      s_a1_("\xce\xb1\xe2\x82\x81", kBufSize),      // alpha_1
      s_a2_("\xce\xb1\xe2\x82\x82", kBufSize),      // alpha_2
      s_cond_("\xce\xba", kBufSize),                // kappa
      s_zc_("z", kBufSize),
      s_bx_("x", kBufSize),
      s_by_("y", kBufSize),
      s_bz_("z", kBufSize),
      s_ex_("x", kBufSize),
      s_ey_("y", kBufSize) {
    camera_ = defaultCamera();

    // The legs at home, the plate assembled there, and BOTH frames of plate
    // motion filled in — so the first contact test differences two real
    // instants rather than one real one and a default-constructed zero, which
    // would read as the plate having just been dropped.  `simStart` is where
    // that rule lives now, so a harness cannot start from a plate the
    // application never starts from.
    //
    // The opening frame has the ball already on a circle and already moving
    // along it, so the first frame that draws is a frame of the demo working
    // rather than a frame of it starting up.  Both errors are zero here, which
    // is what keeps the legs still — see `attract_mode.h`.
    path_ = openingPath();
    path_radius_mm_ = static_cast<float>(path_.radius_m * 1000.0);
    path_period_s_ = static_cast<float>(path_.period_s);
    sim_ = simStart(plate_, kHomeLegRad, attractStart(path_));

    plots_ = {
        {"Ball Position [mm]", {&s_bx_, &s_by_}, 0},
        // Height above the surface, on its own axis: it is zero almost all the
        // time and millimetres when it is not, so sharing a plot with the
        // horizontal position would flatten it to a line on the axis.
        {"Ball Height [mm]", {&s_bz_}, 5},
        // Tracking error, per axis and signed — see `s_ex_`.  Two series
        // rather than their distance, because the distance is flat whenever
        // the error is merely rotating, which under a path is most of the
        // time.  Flat under a fixed setpoint, and the whole story under a path.
        {"Tracking Error [mm]", {&s_ex_, &s_ey_}, 6},
        {"Table Angles [deg]", {&s_phi_, &s_theta_}, 1},
        {"Servo Angles [deg]", {&s_a0_, &s_a1_, &s_a2_}, 2},
        {"Condition Number",   {&s_cond_}, 3},
        {"Table Height [mm]",  {&s_zc_}, 4},
    };
}

void PlateView::attach(GLFWwindow* window) {
    g_plate_view = this;
    glfwSetScrollCallback(window, scrollCallback);
}

PlateView::~PlateView() {
    // The scroll callback outlives this object otherwise, and GLFW would call
    // it with a dangling pointer.
    if (g_plate_view == this) g_plate_view = nullptr;
}

void PlateView::initGL() {
    renderer_ = std::make_unique<LineRenderer>();
}

void PlateView::shutdownGL() {
    renderer_.reset();
}

void PlateView::commandAllServos(float degrees) {
    for (int i = 0; i < 3; ++i) alpha_cmd_deg_[i] = degrees;
}

float PlateView::legDeg(int i) const {
    return static_cast<float>(sim_.alpha_rad[i] / kDeg);
}

bool PlateView::designUsable() const {
    return design_offered_ && gainFitsCascade(design_) &&
           samePlant(design_.mechanism, design_.gravity, kinematics().params(), kGravity);
}

bool PlateView::loopDriving() const {
    // The ball is part of the precondition, not a detail: with the simulation
    // off there is nothing to balance, and the loop would hold the plate at
    // whatever tilt the frozen ball state asks for, forever.
    return balance_engaged_ && ball_enabled_ && designUsable();
}

void PlateView::setDesign(const AutoBalanceDesign& d, bool offered,
                          const std::string& reason) {
    design_ = d;
    design_offered_ = offered;
    design_reason_ = reason;

    if (offered && !gainFitsCascade(d)) {
        design_reason_ = "the gain is not 3 x 7 - this is not the cascade plant";
    } else if (offered && !designUsable()) {
        // The physical sliders move the plant the gain is designed against;
        // the simulated plate keeps the geometry and gravity it was built
        // with.  Refusing is the honest answer — engaging would drive one
        // plate with a gain solved for another, and nothing on screen would
        // say so.
        design_reason_ = "plant geometry or gravity differs from the plate";
    }

    // Losing the design mid-run drops the loop rather than freezing the last
    // command: a stale gain is not a controller.  The one place this happens.
    if (balance_engaged_ && !designUsable()) balance_engaged_ = false;

    // And the one place it is engaged without being asked.  The demo cannot
    // open with `balance_engaged_` simply set true: on the first frame the LQR
    // solve has not run yet, so the drop above would clear the flag and
    // nothing would ever set it again.  Engaging on the first usable design
    // instead means the demo starts balancing the moment it CAN, which is a
    // frame later and is what "already stabilising at load" amounts to.
    //
    // Once, and then never again.  Without the latch this would re-engage a
    // loop the visitor had deliberately dropped, on the very next frame, which
    // is the demo arguing with the person using it.
    if (!auto_engaged_ && !balance_engaged_ && ball_enabled_ && designUsable()) {
        balance_engaged_ = true;
        auto_engaged_ = true;
    }
}

void PlateView::resetBall() {
    sim_.ball = BallState{};
    ball_on_plate_ = true;
}

void PlateView::resetAll() {
    // The whole simulation back to its start state, through the one function
    // that knows what a start state is — legs at home, the plate assembled
    // there, and both frames of plate motion seeded.  Reset used to leave the
    // previous frame's plate motion behind, so the frame after a reset
    // differenced a moving plate against a still one.
    sim_ = simStart(plate_, kHomeLegRad);
    commandAllServos(static_cast<float>(kHomeLegDeg));
    animate_ = false;
    balance_engaged_ = false;
    balance_saturated_ = false;
    balance_clipped_ = false;
    report_ = SimReport{};
    sp_x_mm_ = 0.0f;
    sp_y_mm_ = 0.0f;
    anim_time_ = 0.0f;
    camera_ = defaultCamera();
    plot_state_.clear();
    plot_state_.paused = false;
    for (auto& pc : plots_)
        for (auto* s : pc.series)
            s->data.clear();
    plot_state_.markers.clear();
    sim_time_ = 0.0f;
    ball_on_plate_ = true;
}

void PlateView::step(GLFWwindow* window, float dt) {
    ImGuiIO& io = ImGui::GetIO();

    // --- Mouse orbit (the central dock node is empty, so a drag there is ours) ---
    if (!io.WantCaptureMouse) {
        double mx, my;
        glfwGetCursorPos(window, &mx, &my);
        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
            if (dragging_) {
                camera_.azimuth -= static_cast<float>(mx - last_mx_) * 0.3f;
                camera_.elevation += static_cast<float>(my - last_my_) * 0.3f;
                camera_.elevation = std::clamp(camera_.elevation, 2.0f, 89.0f);
            }
            dragging_ = true;
        } else {
            dragging_ = false;
        }
        last_mx_ = mx;
        last_my_ = my;
    }

    if (plot_state_.paused) return;

    // --- One step, the same one every harness runs ---
    //
    // Everything causal happens inside `stepSim`: the setpoint moves, the loop
    // reads the legs where they ARE and the ball where it IS, the servos
    // follow, the mechanism is solved, the plate's motion is assembled and the
    // ball decides whether it is still on it.  This function's remaining job
    // is to say what the plate is being ASKED for and to draw what came back.
    //
    // It used to carry its own copy of that order, and so did four test
    // harnesses.  See `sim_step.h` and #30.
    SimInput in;
    in.dt = dt;
    in.design = design_;
    in.closed_loop = loopDriving();
    in.ball_enabled = ball_enabled_ && ball_on_plate_;

    // --- Where the leg commands come from, while the loop is not driving ---
    //
    // Three writers, one array, in strict precedence: the closed loop, then
    // the animation, then whatever the sliders last left there.  The first is
    // `in.closed_loop`; the other two write the open-loop command below.
    if (!in.closed_loop && animate_) {
        anim_time_ += dt * anim_speed_;
        const float amp = anim_amplitude_;
        const float w = 2.0f * static_cast<float>(M_PI) * 0.5f * anim_time_;
        alpha_cmd_deg_[0] = 45.0f + amp * std::sin(w);
        alpha_cmd_deg_[1] = 45.0f + amp * std::sin(w + 2.0f * static_cast<float>(M_PI) / 3.0f);
        alpha_cmd_deg_[2] = 45.0f + amp * std::sin(w + 4.0f * static_cast<float>(M_PI) / 3.0f);
    }
    for (int i = 0; i < 3; ++i)
        in.open_loop_cmd_rad[i] = alpha_cmd_deg_[i] * kDeg;

    // --- The setpoint, and who owns it ---
    //
    // A path writes `sp_x_mm_` / `sp_y_mm_` rather than going round them, so
    // the sliders keep working as the readout of where the ball is being sent.
    // A held setpoint is the sliders' own value handed back down.  Which of
    // the two it is, and the read-then-advance rule that goes with a path, are
    // both inside the step.
    if (path_.shape != PathShape::Fixed) {
        path_.radius_m = path_radius_mm_ * 1e-3;
        path_.period_s = path_period_s_;
    }
    in.path = path_;
    in.held_setpoint = Eigen::Vector2d(sp_x_mm_ * 1e-3, sp_y_mm_ * 1e-3);

    report_ = stepSim(plate_, in, sim_);
    sim_time_ += dt;

    // --- What came back ---
    //
    // Both badges and the slider read-out belong to the loop, so they are set
    // together.  `stepSim` already reports neither complaint while the loop is
    // not driving — nobody asked the gain for anything — so this does not have
    // to clear them, only decline to echo a command the sliders themselves
    // wrote.  Echoing that back would hand them their own value through two
    // float conversions, every frame, forever.
    if (in.closed_loop) {
        balance_saturated_ = report_.saturated;
        balance_clipped_ = report_.clipped;
        for (int i = 0; i < 3; ++i)
            alpha_cmd_deg_[i] = static_cast<float>(report_.cmd_rad[i] / kDeg);
    }

    // And the setpoint, where a path owns it.  A held setpoint is the
    // sliders' and is left alone for the same reason.
    if (path_.shape != PathShape::Fixed) {
        sp_x_mm_ = static_cast<float>(report_.setpoint(0) * 1000.0);
        sp_y_mm_ = static_cast<float>(report_.setpoint(1) * 1000.0);
    }

    condition_num_ = kinematics().condition_number(sim_.alpha_rad, sim_.pose);
    manipulability_ = kinematics().manipulability(sim_.alpha_rad, sim_.pose);

    // Only while the ball is being simulated, as before: with the simulation
    // off the last reading is what the panel goes on showing, and the contact
    // line is one of the readings.
    if (in.ball_enabled) {
        if (report_.airborne) airborne_flash_s_ = 1.5f;
        else airborne_flash_s_ = std::max(0.0f, airborne_flash_s_ - dt);
    }

    if (report_.left_plate) {
        ball_on_plate_ = false;
        if (ball_auto_reset_) resetBall();
    }

    // --- Plots ---
    // `report_.ball_plate` rather than a second `plateFrame` call: this is the
    // ball the step just finished with, and asking again is how a plot comes to
    // be drawn against a different frame's plate motion from the one the
    // contact test used.  The draw pass below does call it, and has to — it
    // runs whether or not a step happened this frame.
    plot_state_.push_time(sim_time_);
    const Eigen::Matrix<double, 6, 1>& bp = report_.ball_plate;
    push_series(s_bx_,    static_cast<float>(bp(0) * 1000), plot_state_);
    push_series(s_by_,    static_cast<float>(bp(1) * 1000), plot_state_);
    // Height above the surface, not above the table centre: zero means resting
    // on it, which is what a reader of this plot wants the line to mean.
    push_series(s_bz_,    static_cast<float>((bp(2) - kBallRadius) * 1000), plot_state_);
    // Signed, and in the plate frame the setpoint is already expressed in, so
    // a positive x error means the ball is further along +x than it was told.
    push_series(s_ex_,    static_cast<float>(bp(0) * 1000) - sp_x_mm_, plot_state_);
    push_series(s_ey_,    static_cast<float>(bp(1) * 1000) - sp_y_mm_, plot_state_);
    push_series(s_phi_,   static_cast<float>(sim_.pose.phi / kDeg), plot_state_);
    push_series(s_theta_, static_cast<float>(sim_.pose.theta / kDeg), plot_state_);
    push_series(s_a0_,    legDeg(0), plot_state_);
    push_series(s_a1_,    legDeg(1), plot_state_);
    push_series(s_a2_,    legDeg(2), plot_state_);
    push_series(s_cond_,  static_cast<float>(condition_num_), plot_state_);
    push_series(s_zc_,    static_cast<float>(sim_.pose.z_c * 1000), plot_state_);
}

void PlateView::drawPanels() {
    drawControls();
    draw_time_series_panel("Plate Plots", plot_state_, plots_);
    comparison_panel_.draw();
}

void PlateView::drawControls() {
    ImGui::Begin("Plate Control");

    // --- Ball ---
    ImGui::SeparatorText("Ball");
    ImGui::Checkbox("Simulate ball", &ball_enabled_);
    ImGui::SameLine();
    ImGui::Checkbox("Auto-reset", &ball_auto_reset_);

    // Three lines, always, whatever the ball is doing.  A ball that has gone
    // off the edge used to collapse this whole block to a single line, which
    // moved every control below it up by two rows at the exact moment the
    // visitor was reaching for Reset Ball.  When the sim has stopped, the
    // readings are simply the last ones taken — which is the useful thing to
    // show anyway, since they say where it went.
    const Eigen::Matrix<double, 6, 1> bp =
        plateFrame(sim_.ball, sim_.motion, kBallRadius);
    ImGui::Text("Position: %+7.1f, %+7.1f mm", bp(0) * 1000, bp(1) * 1000);
    ImGui::Text("Velocity: %+7.1f, %+7.1f mm/s", bp(3) * 1000, bp(4) * 1000);

    // Exactly one contact line, every frame.  The `else` is not clutter: it is
    // what stops the panel below moving when the ball settles.  The flash
    // exists because a separation lasts a handful of frames and would
    // otherwise be a line of text nobody is quick enough to read.
    if (!ball_on_plate_) {
        ImGui::TextColored(ImVec4(1.0f, 0.3f, 0.3f, 1.0f),
                           "the ball has left the plate");
    } else if (sim_.ball.airborne) {
        ImGui::TextColored(ImVec4(1.0f, 0.55f, 0.15f, 1.0f),
                           "AIRBORNE  %+.1f mm, %+.0f mm/s",
                           (bp(2) - kBallRadius) * 1000, bp(5) * 1000);
    } else if (airborne_flash_s_ > 0.0f) {
        ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.2f, 1.0f),
                           "the ball left the plate");
    } else {
        ImGui::TextDisabled("in contact");
    }

    if (ImGui::Button("Reset Ball")) resetBall();
    ImGui::SameLine();
    if (ImGui::Button("Nudge +x") && !sim_.ball.airborne) sim_.ball.rolling(2) += ball_nudge_;
    ImGui::SameLine();
    if (ImGui::Button("Nudge +y") && !sim_.ball.airborne) sim_.ball.rolling(3) += ball_nudge_;
    ImGui::SameLine();
    // The disturbance the loop is watched against: put the ball a fifth of the
    // plate out and let go.  A nudge tests recovery from a kick; this tests
    // recovery from a position, which is what the design surface is aimed at.
    if (ImGui::Button("Displace")) {
        sim_.ball = BallState{};
        sim_.ball.rolling << 0.06, -0.04, 0.0, 0.0;
        ball_on_plate_ = true;
    }
    // The top of this slider is the hardest shove the interface can offer, and
    // the preset tunings are tested against exactly that number.  See
    // `kMaxNudgeSpeed`.
    ImGui::SliderFloat("Nudge [m/s]", &ball_nudge_, 0.02f,
                       static_cast<float>(kMaxNudgeSpeed), "%.2f");

    drawBalanceControls();

    // --- Servo Angles ---
    ImGui::SeparatorText("Servo Angles");

    const bool driven = loopDriving();
    // Both branches, so the sliders below do not move by a row when the loop
    // is engaged.  Saying who is driving is worth a line; saying it only half
    // the time is worth a line that jumps.
    if (driven) {
        // Not hidden — the sliders become the readout of what u = -Kx is
        // asking for, which is the most direct view of the loop working.
        ImGui::TextDisabled("commanded by the LQR gain");
    } else {
        ImGui::TextDisabled("drag to command the legs");
    }
    ImGui::BeginDisabled(driven);

    ImGui::Checkbox("Link all servos", &link_servos_);

    const float a_min = static_cast<float>(kinematics().params().alpha_min / kDeg);
    const float a_max = static_cast<float>(kinematics().params().alpha_max / kDeg);

    if (link_servos_) {
        if (ImGui::SliderFloat("All##servo", &alpha_cmd_deg_[0], a_min, a_max, "%.1f deg")) {
            alpha_cmd_deg_[1] = alpha_cmd_deg_[0];
            alpha_cmd_deg_[2] = alpha_cmd_deg_[0];
        }
    } else {
        ImGui::SliderFloat("\xce\xb1\xe2\x82\x80 [deg]", &alpha_cmd_deg_[0], a_min, a_max, "%.1f");
        ImGui::SliderFloat("\xce\xb1\xe2\x82\x81 [deg]", &alpha_cmd_deg_[1], a_min, a_max, "%.1f");
        ImGui::SliderFloat("\xce\xb1\xe2\x82\x82 [deg]", &alpha_cmd_deg_[2], a_min, a_max, "%.1f");
    }

    if (ImGui::Button("Home (45)")) commandAllServos(45.0f);
    ImGui::SameLine();
    if (ImGui::Button("Low (20)")) commandAllServos(20.0f);
    ImGui::SameLine();
    if (ImGui::Button("High (70)")) commandAllServos(70.0f);
    ImGui::EndDisabled();

    // The sliders are the COMMAND; the legs lag behind it.  Without this line
    // the two readouts disagree on screen with nothing saying why.
    ImGui::Text("legs at %.1f, %.1f, %.1f deg (lag \xcf\x84 = %.3f s)",
                legDeg(0), legDeg(1), legDeg(2), design_.servo_tau);

    ImGui::Spacing();
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.7f, 0.2f, 0.2f, 1.0f));
    if (ImGui::Button("Reset All", ImVec2(-1, 0))) resetAll();
    ImGui::PopStyleColor();

    // --- Animation ---
    ImGui::SeparatorText("Animation");
    ImGui::BeginDisabled(driven);
    ImGui::Checkbox("Animate", &animate_);
    if (animate_) {
        ImGui::SliderFloat("Speed", &anim_speed_, 0.1f, 5.0f, "%.1f");
        ImGui::SliderFloat("Amplitude [deg]", &anim_amplitude_, 1.0f, 20.0f, "%.1f");
    }
    ImGui::EndDisabled();

    // --- Pause ---
    ImGui::SeparatorText("Simulation");
    if (ImGui::Button(plot_state_.paused ? "Resume" : "Pause", ImVec2(-1, 0))) {
        plot_state_.paused = !plot_state_.paused;
        if (plot_state_.paused) {
            plot_state_.pause_t_max = plot_state_.latest_time();
            plot_state_.pause_t_min = plot_state_.pause_t_max - plot_state_.time_window;
        }
    }

    // --- Table Pose ---
    ImGui::SeparatorText("Table Pose (FK)");

    auto ok_col = [](bool ok) {
        return ok ? ImVec4(0.3f, 1.0f, 0.3f, 1.0f) : ImVec4(1.0f, 0.3f, 0.3f, 1.0f);
    };
    ImGui::TextColored(ok_col(report_.fk.converged),
                       "FK: %s (%d iter)", report_.fk.converged ? "OK" : "FAIL",
                       report_.fk.iterations);
    ImGui::Text("Roll:   %+.2f deg", sim_.pose.phi / kDeg);
    ImGui::Text("Pitch:  %+.2f deg", sim_.pose.theta / kDeg);
    ImGui::Text("Heave:  %.1f mm",   sim_.pose.z_c * 1000);

    // --- Jacobian ---
    ImGui::SeparatorText("Jacobian Analysis");
    auto cond_col = [](double c) -> ImVec4 {
        if (c < 10) return {0.3f, 1.0f, 0.3f, 1.0f};
        if (c < 20) return {1.0f, 0.8f, 0.2f, 1.0f};
        return {1.0f, 0.3f, 0.3f, 1.0f};
    };
    ImGui::TextColored(cond_col(condition_num_),
                       "Condition #:   %.2f", condition_num_);
    ImGui::Text("Manipulability: %.4f", manipulability_);

    const float cond_frac = std::min(static_cast<float>(condition_num_ / 50.0), 1.0f);
    ImGui::ProgressBar(cond_frac, ImVec2(-1, 0),
                       condition_num_ < 10 ? "Good" :
                       condition_num_ < 20 ? "Degraded" : "Poor");

    // --- Display ---
    ImGui::SeparatorText("Display");
    ImGui::Checkbox("Grid", &show_grid_);
    ImGui::SameLine();
    ImGui::Checkbox("Axes", &show_axes_);
    ImGui::SameLine();
    ImGui::Checkbox("Joints", &show_joints_);

    // --- Camera ---
    ImGui::SeparatorText("Camera");
    ImGui::SliderFloat("Azimuth",   &camera_.azimuth,   -180.0f, 180.0f, "%.0f deg");
    ImGui::SliderFloat("Elevation", &camera_.elevation,  5.0f, 89.0f, "%.0f deg");
    ImGui::SliderFloat("Distance",  &camera_.distance,   0.3f, 2.0f, "%.2f m");
    if (ImGui::Button("Reset Camera")) camera_ = defaultCamera();

    // --- Velocity Jacobian (collapsed) ---
    if (ImGui::CollapsingHeader("Velocity Jacobian")) {
        const Eigen::Matrix3d Jv = kinematics().velocity_jacobian(legsRad(), sim_.pose);
        ImGui::Text("        servo0    servo1    servo2");
        ImGui::Text("roll   %+.4f   %+.4f   %+.4f", Jv(0,0), Jv(0,1), Jv(0,2));
        ImGui::Text("pitch  %+.4f   %+.4f   %+.4f", Jv(1,0), Jv(1,1), Jv(1,2));
        ImGui::Text("heave  %+.4f   %+.4f   %+.4f", Jv(2,0), Jv(2,1), Jv(2,2));
    }

    ImGui::End();
}

// The one control the whole ticket is about: hand the plate over to the gain
// the model panel just designed.
//
// The three servo sliders and u = -Kx want to write the same array, so they
// cannot both be live.  The loop wins while it is engaged and the sliders go
// read-only rather than disappearing — a disabled slider that keeps moving is
// the clearest possible statement of what the controller is doing, and there
// is nothing left for a manual command to mean.  The setpoint takes over the
// steering job, and it is a BALL position rather than a plate tilt, because
// that is the quantity the design is regulating.
void PlateView::drawBalanceControls() {
    ImGui::SeparatorText("Balance Loop");

    // Read-only here.  `setDesign` is the one place the loop is dropped, so a
    // draw pass cannot disagree with it about whether the loop is running.
    if (!designUsable()) {
        ImGui::BeginDisabled(true);
        bool off = false;
        ImGui::Checkbox("Engage", &off);
        ImGui::EndDisabled();
        ImGui::TextWrapped("unavailable: %s", design_reason_.c_str());
        return;
    }

    ImGui::BeginDisabled(!ball_enabled_);
    ImGui::Checkbox("Engage", &balance_engaged_);
    ImGui::EndDisabled();
    if (!ball_enabled_) {
        ImGui::TextWrapped("unavailable: there is no ball to balance");
        return;
    }
    if (!balance_engaged_) {
        ImGui::TextDisabled("gain ready - %d x %d, legs from %.1f deg",
                            (int)design_.K.rows(), (int)design_.K.cols(),
                            design_.home_leg_rad / kDeg);
        return;
    }

    // A ball at rest anywhere on a flat plate is an equilibrium, so a HELD
    // setpoint needs no feedforward: it enters as a reference state that is
    // reachable with the plate flat.  A moving one does need one, and gets it
    // — see `BallReference`.  Bounded well inside the plate, since the edge is
    // where the ball leaves.
    const float lim = static_cast<float>((kinematics().params().R_table - 0.05) * 1000.0);

    // --- Trajectory ---
    // Index-coupled to PathShape, in the order it declares them.
    const char* shapes[] = {"Hold a point", "Circle", "Square", "Triangle"};
    static_assert(IM_ARRAYSIZE(shapes) == static_cast<int>(PathShape::Triangle) + 1,
                  "shapes is index-coupled to PathShape");
    int shape_idx = static_cast<int>(path_.shape);
    if (ImGui::Combo("Trajectory", &shape_idx, shapes, IM_ARRAYSIZE(shapes))) {
        path_.shape = static_cast<PathShape>(shape_idx);
        // The new shape picks up nearest to where the setpoint already is,
        // rather than at whatever the carried phase means on it.  Phase is not
        // comparable across shapes — the circle's zero is at +x, a polygon's
        // first corner is at the top — so carrying it moved the target 170 mm
        // to the far side of the path and the loop hauled the ball across
        // after it, into the workspace clip.  See `phaseNearest`.
        //
        // From the radius the slider holds, not the one `step` last copied:
        // the same one-frame skew the lap floor had, and the polygon corners
        // this reads are a function of it.
        path_.radius_m = path_radius_mm_ * 1e-3;
        sim_.path_phase = phaseNearest(
            path_, Eigen::Vector2d(sp_x_mm_ * 1e-3, sp_y_mm_ * 1e-3));
    }

    const bool on_a_path = path_.shape != PathShape::Fixed;
    if (on_a_path) {
        // Size and speed both, because they trade against each other and the
        // interesting settings are at both ends: large and slow makes the
        // shape unmistakable, small and fast makes the tracking lag and the
        // rounded corners unmistakable instead.
        ImGui::SliderFloat("size [mm]", &path_radius_mm_, 20.0f,
                           static_cast<float>(kMaxPathRadius * 1000.0), "%.0f");
        // The lap slider's floor moves with the size, because what loses the
        // ball is the setpoint's SPEED and a fixed floor would either forbid
        // fast small paths that are safe or allow fast large ones that are not.
        //
        // Taken from the size just dragged, not from the one `step` last
        // copied: the two are the same only until someone moves the slider,
        // and the frame where they differ is exactly the frame where a
        // just-enlarged path would keep the smaller path's floor and run at a
        // speed this bound exists to forbid.
        path_.radius_m = path_radius_mm_ * 1e-3;
        const float lap_min = static_cast<float>(minPeriod(path_));
        ImGui::SliderFloat("lap [s]", &path_period_s_, lap_min, 30.0f, "%.1f");
        path_period_s_ = static_cast<float>(clampPeriod(path_, path_period_s_));
        ImGui::TextDisabled("setpoint speed %.0f mm/s (max %.0f)",
                            pathLength(path_) / std::max(0.1, path_.period_s) * 1000.0,
                            kMaxSetpointSpeed * 1000.0);
    }

    // The setpoint sliders stay visible and go read-only under a path, for the
    // same reason the servo sliders do under the loop: a disabled control that
    // keeps moving is the clearest statement of what is driving it.
    ImGui::BeginDisabled(on_a_path);
    ImGui::SliderFloat("set x [mm]", &sp_x_mm_, -lim, lim, "%.0f");
    ImGui::SliderFloat("set y [mm]", &sp_y_mm_, -lim, lim, "%.0f");
    if (ImGui::Button("Centre setpoint")) { sp_x_mm_ = 0.0f; sp_y_mm_ = 0.0f; }
    ImGui::EndDisabled();

    const Eigen::Matrix<double, 6, 1> bp =
        plateFrame(sim_.ball, sim_.motion, kBallRadius);
    const double err = std::hypot(bp(0) - sp_x_mm_ * 1e-3,
                                  bp(1) - sp_y_mm_ * 1e-3);
    // Both warnings ride on the error line rather than taking lines of their
    // own.  They are independent and they flicker at frame rate, so on their
    // own lines the whole panel below has three resting positions and a slider
    // the visitor is reaching for moves out from under the cursor.  See
    // `statusBadge`.
    //
    // The error reading is also the right line to hang them on: it is the
    // number they explain.  "error: 178.2 mm  clipped" says both that the loop
    // is a long way off and why it is not closing the gap.
    ImGui::Text("error: %6.1f mm", err * 1000.0);

    if (balance_clipped_) {
        // A different thing from saturation: the legs are inside their travel
        // and the controller is still not getting what it asked for, because
        // the pose it wants is not one this mechanism has.
        statusBadge("clipped", kBadgeWarn,
                    "The command was pulled back toward the home pose until "
                    "the mechanism could actually assemble it.\n\n"
                    "Distinct from saturation: every leg is inside its travel "
                    "limits, and the pose they were asked for still does not "
                    "exist. The servo limits are a box; the workspace is not, "
                    "and barely half the box has an assembly at all.\n\n"
                    "See CONTEXT.md, \"Workspace vs. servo box\".");
    }

    if (balance_saturated_) {
        // Not a failure — the legs really do stop at 10 and 80 degrees.
        statusBadge("saturated", kBadgeWarn,
                    "At least one leg command hit a travel limit and was "
                    "clamped there.\n\n"
                    "Not a failure: the servos really do stop at 10 and 80 "
                    "degrees. But a tuning that lives against the stops is no "
                    "longer the tuning that was designed, and the closed-loop "
                    "poles on screen stop describing it.");
    }
}

void PlateView::drawScene(float aspect) {
    if (!renderer_) return;

    const Eigen::Matrix4f proj = perspective(45.0f, aspect, 0.01f, 10.0f);
    const Eigen::Matrix4f view = camera_.view_matrix();
    renderer_->set_vp(proj * view);
    renderer_->begin();
    drawMechanism();
    renderer_->flush();
}

void PlateView::drawMechanism() {
    LineRenderer& lr = *renderer_;
    const TableParams& p = kinematics().params();

    // Ground grid
    if (show_grid_) {
        const double s = 0.5;
        for (double x = -s; x <= s + 1e-9; x += 0.05) {
            lr.line({x, -s, 0}, {x, s, 0}, col::grid, 1.0f);
            lr.line({-s, x, 0}, {s, x, 0}, col::grid, 1.0f);
        }
    }

    if (show_axes_) lr.axes({0, 0, 0}, 0.08);

    // Ground: filled disc + edge
    lr.disc({0, 0, 0}, p.R_ground, {0, 0, 1}, col::ground_f, 64);
    lr.circle({0, 0, 0}, p.R_ground, {0, 0, 1}, col::ground, 64, 2.0f);

    // Table: filled disc + edge
    const Eigen::Matrix3d R = kinematics().table_rotation(sim_.pose.phi, sim_.pose.theta);
    const Eigen::Vector3d tc(0, 0, sim_.pose.z_c);
    const Eigen::Vector3d tn = R * Eigen::Vector3d(0, 0, 1);

    lr.disc(tc, p.R_table, tn, col::table_f, 64);
    lr.circle(tc, p.R_table, tn, col::table_c, 64, 2.5f);

    // Table cross-hairs
    const Eigen::Vector3d tx = R * Eigen::Vector3d(p.R_table * 0.8, 0, 0);
    const Eigen::Vector3d ty = R * Eigen::Vector3d(0, p.R_table * 0.8, 0);
    lr.line(tc - tx, tc + tx, {0.3f, 0.6f, 0.9f, 0.5f}, 1.0f);
    lr.line(tc - ty, tc + ty, {0.3f, 0.6f, 0.9f, 0.5f}, 1.0f);

    // Legs
    for (int i = 0; i < 3; ++i) {
        const double alpha_i = sim_.alpha_rad[i];
        const Eigen::Vector3d G = kinematics().ground_point(i);
        const Eigen::Vector3d K = kinematics().knee_position(i, alpha_i);
        const Eigen::Vector3d P = kinematics().table_point_world(i, sim_.pose);

        lr.line(G, K, col::leg_L1, 3.0f);
        lr.line(K, P, col::leg_L2, 3.0f);

        if (show_joints_) {
            lr.point(G, col::joint_sph, 0.008f, 4.0f);
            lr.point(K, col::joint_rev, 0.008f, 4.0f);
            lr.point(P, col::joint_sph, 0.008f, 4.0f);
        }
        if (show_axes_) {
            const Eigen::Vector3d t = kinematics().tangent_dir(i);
            lr.line(G - 0.02 * t, G + 0.02 * t, {0.7f, 0.7f, 0.0f, 0.6f}, 1.5f);
        }
    }

    if (show_joints_) lr.point(tc, col::table_c, 0.006f, 3.0f);

    // The trajectory, and the point on it the ball is chasing.  Drawn on the
    // plate's surface in the plate's own frame, so it tilts with the plate —
    // it is a target expressed in plate coordinates, and drawing it flat in
    // the world would put it somewhere the ball is not being sent.
    if (path_.shape != PathShape::Fixed) {
        Eigen::Matrix<double, 2, Eigen::Dynamic> outline;
        pathOutline(path_, 64, &outline);
        for (int i = 0; i + 1 < outline.cols(); ++i) {
            const Eigen::Vector3d a =
                tc + R * Eigen::Vector3d(outline(0, i), outline(1, i), 0.001);
            const Eigen::Vector3d b =
                tc + R * Eigen::Vector3d(outline(0, i + 1), outline(1, i + 1), 0.001);
            lr.line(a, b, col::path, 1.5f);
        }
    }
    {
        const Eigen::Vector3d sp =
            tc + R * Eigen::Vector3d(sp_x_mm_ * 1e-3, sp_y_mm_ * 1e-3, 0.001);
        lr.circle(sp, 0.012, tn, col::setpoint, 20, 2.0f);
        lr.point(sp, col::setpoint, 0.006f, 3.0f);
    }

    // Ball.  Its state is in the plate's own frame, so the same rotation that
    // places the leg attachment points places the ball — lifted off the surface
    // by its radius along the plate normal.
    if (ball_enabled_) {
        // The ball's own z, not a constant radius: since #23 it can leave, and
        // a hop that the physics performs but the renderer flattens would be
        // the same failure as the folded plate in #22 — a model doing
        // something the picture denies.
        const Eigen::Matrix<double, 6, 1> bpv =
            plateFrame(sim_.ball, sim_.motion, kBallRadius);
        const Eigen::Vector3d bc =
            tc + R * Eigen::Vector3d(bpv(0), bpv(1), bpv(2));
        const auto& edge = !ball_on_plate_ ? col::ball_off
                         : sim_.ball.airborne  ? col::ball_air
                                           : col::ball;

        // A wire sphere: a filled disc face-on to the plate, plus three great
        // circles.  The line renderer draws no triangulated solids, and three
        // circles read as a sphere where one circle reads as a hole.
        lr.disc(bc, kBallRadius, tn, col::ball_f, 24);
        lr.circle(bc, kBallRadius, tn, edge, 24, 2.0f);
        lr.circle(bc, kBallRadius, R * Eigen::Vector3d(1, 0, 0), edge, 24, 1.5f);
        lr.circle(bc, kBallRadius, R * Eigen::Vector3d(0, 1, 0), edge, 24, 1.5f);

        // A dropline to the plate centre-plane, so the ball's offset from the
        // middle of the plate is readable from any camera angle.
        const Eigen::Vector3d foot =
            tc + R * Eigen::Vector3d(bpv(0), bpv(1), 0.0);
        lr.line(tc, foot, {0.98f, 0.45f, 0.09f, 0.35f}, 1.0f);

        // While it is off the surface, the gap itself: a line from the ball
        // down to where it would be resting, and a ring on the plate under it.
        // Without them a hop of a few millimetres is invisible at this camera
        // distance, and the whole point is that it is happening.
        if (sim_.ball.airborne) {
            const Eigen::Vector3d rest =
                tc + R * Eigen::Vector3d(bpv(0), bpv(1), kBallRadius);
            lr.line(rest, bc, col::ball_air, 1.5f);
            lr.circle(foot, kBallRadius, tn, col::ball_air, 20, 1.5f);
        }
    }
}

}  // namespace caliburn
