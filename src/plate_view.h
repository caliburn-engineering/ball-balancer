// src/plate_view.h
#pragma once

#include "attract_mode.h"
#include "ball_contact.h"
#include "auto_balance.h"
#include "ball_sim.h"
#include "comparison_panel.h"
#include "plot_panel.h"
#include "renderer.h"
#include "setpoint_path.h"
#include "sim_step.h"
#include "table_kinematics.h"

#include <Eigen/Core>
#include <array>
#include <memory>
#include <string>
#include <vector>

struct GLFWwindow;

namespace caliburn {

/// The ball-balancer half of the merged application: the 3-RRS table, the ball
/// rolling on it, the 3D scene, and the panels that drive and plot them.
///
/// It owns no GLFW window, no ImGui context and no main loop — the analyzer's
/// entry point owns all three, and this is what the merge amounted to.  Frame
/// order is fixed by the caller: attach() once before the ImGui backend is
/// initialised, initGL() once after the GL loader is up, then step() ->
/// drawPanels() -> drawScene() every frame.
class PlateView {
public:
    PlateView();
    ~PlateView();

    // Non-copyable, non-movable: the plot configs hold pointers to the series
    // members, so relocating the object would dangle them.
    PlateView(const PlateView&) = delete;
    PlateView& operator=(const PlateView&) = delete;

    /// Install the scroll handler.  Must run BEFORE ImGui_ImplGlfw_InitForOpenGL,
    /// so that ImGui's backend chains to it rather than replacing it — the
    /// other order silently kills scroll-wheel zoom in every ImPlot panel.
    void attach(GLFWwindow* window);

    /// Build GL resources.  Requires a current context with GLAD loaded.
    void initGL();

    /// Release GL resources.  Must run while the context is still current —
    /// leaving it to the destructor deletes buffers after glfwTerminate().
    void shutdownGL();

    /// Advance kinematics, ball and plot buffers by one frame.
    /// Call between glfwPollEvents() and ImGui::NewFrame().
    void step(GLFWwindow* window, float dt);

    /// Draw the "Plate Control", "Plate Plots" and "Model Comparison" panels.
    void drawPanels();

    /// Draw the 3D scene into the currently-set GL viewport.
    void drawScene(float aspect);

    /// Hand over the design surface's current answer, once per frame.
    ///
    /// `offered` is the caller's claim — that LQR is the selected controller,
    /// that the solve succeeded, and that the plant it was solved against is
    /// the cascade.  Only the model panel can know any of that.  The plate
    /// adds the two checks it alone can make, the gain's shape and whether the
    /// mechanism in `d` is the one being simulated, and states its own reason
    /// when either fails.
    ///
    /// `d.servo_tau` is honoured whether or not the design is offered: the
    /// legs have first-order lag in manual driving too, and the model panel's
    /// tau slider is the one place that number lives.
    ///
    /// This is also where the opening closes the loop for the first time —
    /// see `auto_engaged_`.  Engaging and dropping are the same decision read
    /// in two directions, and they belong at the same seam.
    void setDesign(const AutoBalanceDesign& d, bool offered,
                   const std::string& reason);

    /// Camera orbiting, for the scroll callback.
    OrbitCamera& camera() { return camera_; }

private:
    void drawControls();
    void drawBalanceControls();
    void drawMechanism();
    void resetAll();
    void resetBall();

    /// Command all three legs to one angle.  A command, not a teleport: the
    /// Home/Low/High buttons ask, and the legs arrive one lag later like every
    /// other command.  The reset path does not go through here — `simStart`
    /// puts the legs where they belong, since a reset is not driving anything
    /// and has no lag to respect.
    void commandAllServos(float degrees);

    /// Where the legs actually are, in radians.  Now simply the simulation's
    /// own leg state — the panels and the 3D view read it here rather than
    /// converting a display float back to an angle, which is what they used to
    /// do and which quietly quantised the plate the application simulated
    /// against the one every harness did.
    const std::array<double, 3>& legsRad() const { return sim_.alpha_rad; }

    /// One leg, in degrees, for the readouts and the plots.
    float legDeg(int i) const;

    /// The mechanism the panels and the 3D view solve against — the plate's,
    /// not a second one built to draw with.
    const TableKinematics& kinematics() const { return plate_.kinematics(); }

    /// True when the design on hand can actually drive this plate.
    bool designUsable() const;

    /// True while the balance loop, and not the sliders, owns the leg command.
    /// The single authority: `step` acts on it and `drawControls` greys the
    /// manual controls on it, so the two cannot disagree about who is driving.
    bool loopDriving() const;

    /// The plate, the ball on it and the gravity they share, in the one object
    /// `stepSim` takes.  Both halves used to be built here and again in every
    /// harness; a test that measured a different ball from the one that ships
    /// would have looked exactly like a test that measured this one.
    SimPlate plate_;

    /// Where the legs are, where the plate is, where the ball is, and how far
    /// round the lap the setpoint has got.  Everything a frame hands to the
    /// next, and the only thing `stepSim` writes.
    SimState sim_;

    std::unique_ptr<LineRenderer> renderer_;

    // --- Servos ---
    // Two arrays, since the loop was closed: `alpha_cmd_deg_` is what the
    // sliders, the animation or the controller ASK for, `alpha_deg_` is where
    // the legs actually are.  Before the split the slider *was* the leg angle,
    // and u = -Kx around that is an algebraic loop on the leg states — the
    // servo lag the plant model claims had to become real.
    //
    // While the loop is closed the controller writes the command array, so the
    // (disabled) sliders read out what it is doing.
    //
    // Only the COMMAND is a member now.  Where the legs are is `sim_.alpha_rad`,
    // because that is a fact about the simulation rather than about this panel
    // — and holding it as a display float was quietly rounding the plate's own
    // state to seven digits every frame.
    float alpha_cmd_deg_[3] = {45.0f, 45.0f, 45.0f};
    bool link_servos_ = false;

    // --- Balance loop ---
    AutoBalanceDesign design_{};
    bool design_offered_ = false;
    std::string design_reason_ = "select LQR as the controller type";
    bool balance_engaged_ = false;
    bool balance_saturated_ = false;
    bool balance_clipped_ = false;   ///< command was not holdable, pulled back
    float sp_x_mm_ = 0.0f;  // ball setpoint, plate frame
    float sp_y_mm_ = 0.0f;

    // --- Trajectory tracking (#24) ---
    // A moving setpoint.  It writes `sp_x_mm_` / `sp_y_mm_` rather than going
    // round them, so the control law is untouched and the sliders keep working
    // as the readout of where the ball is being sent — the same arrangement
    // the servo sliders have under the balance loop.
    SetpointPath path_{};
    float path_radius_mm_ = 120.0f;
    float path_period_s_ = 10.0f;

    // --- Camera ---
    OrbitCamera camera_;
    bool dragging_ = false;
    double last_mx_ = 0.0, last_my_ = 0.0;

    // --- Display ---
    bool show_axes_ = true;
    bool show_grid_ = true;
    bool show_joints_ = true;

    // --- Animation ---
    bool animate_ = false;
    float anim_time_ = 0.0f;
    float anim_speed_ = 1.0f;
    float anim_amplitude_ = 5.0f;

    // --- Computed each frame ---
    // What the last step reported, kept only because the panels draw it.  The
    // pose itself is not here: it is `sim_.pose`, which is both the plate's
    // assembly and the seed the next solve starts from — one field, because
    // they were always the same field.
    SimReport report_{};
    double condition_num_ = 0.0;
    double manipulability_ = 0.0;
    float sim_time_ = 0.0f;

    // --- The opening ---
    // The demo running itself until somebody turns up: the loop engaged as
    // soon as a gain exists, and the ball already tracing a circle.  Without
    // it the page opens on a balanced ball sitting still, which is
    // indistinguishable from a broken build.  See issue #17 and
    // `attract_mode.h`.
    //
    // There is no schedule and no running flag any more.  The opening is a
    // STATE, not a performance: it is set once in the constructor and then the
    // visitor owns it — change the trajectory, drag the setpoint, engage or
    // drop the loop, and nothing here will argue.  The old attract mode had to
    // watch for a visitor arriving so it could stand its kicks down; a circle
    // has nothing to stand down.
    //
    // One latch survives, below, and only because the loop cannot be engaged
    // on frame zero.
    //
    // The other half of the opening state is not here and cannot be: selecting
    // LQR is a fact about the model panel's controller type, which this class
    // deliberately cannot see — `handDesignToPlate` exists for exactly that
    // reason.  It is set beside the rest of the app's opening state, in
    // `visualizer.cpp`'s main().

    /// Whether the loop has yet been engaged for the visitor, once, without
    /// being asked.  A one-shot latch rather than a mode: after it fires the
    /// checkbox is the visitor's, and a loop they drop stays dropped.
    bool auto_engaged_ = false;

    // --- Ball ---
    // Six states now, in two phases: the plate can lose contact and the ball
    // can fly.  Rarely, on the shipped tuning: measured against the step this
    // class drives, it never separates at the disturbance the Nudge buttons
    // offer, and above that in narrow slivers of direction by fractions of a
    // millimetre.  See issue #23, #30 and `ball_contact.h`.
    //
    // The ball itself and both frames of plate motion live in `sim_`; what is
    // left here is what the PANEL knows about it rather than what the physics
    // does.
    float airborne_flash_s_ = 0.0f;  ///< keeps a brief hop legible in the panel
    bool ball_enabled_ = true;
    bool ball_on_plate_ = true;
    bool ball_auto_reset_ = true;
    float ball_nudge_ = 0.15f;  // [m/s]

    // --- Plots ---
    PlotState plot_state_;
    TimeSeries s_phi_, s_theta_;
    TimeSeries s_a0_, s_a1_, s_a2_;
    TimeSeries s_cond_, s_zc_;
    TimeSeries s_bx_, s_by_, s_bz_;
    /// Tracking error, as its two SIGNED components rather than one distance.
    ///
    /// `hypot` of the pair was what this plotted, and it rectifies: it cannot
    /// go negative and it throws away which way the ball is off.  On a circle
    /// the loop trails the setpoint by a roughly fixed angle, so the distance
    /// is nearly constant while each component swings the full amplitude a
    /// quarter-lap out of phase with the other — measured, |e| held between
    /// 5.2 and 5.7 mm while ex and ey each swept +/-5 mm.  The plot drew that
    /// as a flat line, which reads as a steady offset and is the wrong
    /// conclusion: the error is a vector of near-constant length going round
    /// once per lap.
    ///
    /// The magnitude is not lost — it is the `error:` readout beside the
    /// setpoint sliders.  Plot shows the structure, number shows the size.
    TimeSeries s_ex_, s_ey_;
    std::vector<PlotConfig> plots_;

    ComparisonPanel comparison_panel_;
};

}  // namespace caliburn
