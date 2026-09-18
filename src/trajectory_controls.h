// src/trajectory_controls.h
#pragma once

#include "setpoint_path.h"
#include "sim_step.h"

#include <Eigen/Core>

namespace caliburn {

/// The trajectory half of the plate panel: where the ball is being sent, and
/// the four controls that change it.
///
/// **It exists because the panel had an invariant and kept it by convention.**
/// The path and the setpoint are one piece of state with one promise over it —
/// *the setpoint never jumps* — and every rule that keeps that promise is a
/// rule about what happens BETWEEN two of its fields: the lap floor applied
/// against the radius just dragged, the phase re-seeded from where the setpoint
/// already is, the setpoint sliders being a readout under a path and the input
/// under a held point.  Those rules were enforced at four call sites in
/// `plate_view.cpp`, each of which had to remember to refresh the path's radius
/// first, and two of them had already been got wrong once each — the second
/// time in front of a user (see #24 round three).
///
/// So the fields are held together and the rules are the methods that move
/// them.  What is left in the panel is ImGui plus forwarding, and what a test
/// drives is the object the panel drives, through the same entry points.
///
/// See [#28](https://github.com/caliburn-engineering/caliburn/issues/28) for the
/// seam, and `plate_view.h` for what deliberately stays untested on the far
/// side of it.
///
/// **The two slider floats are gone.**  `path_.radius_m` and `path_.period_s`
/// are the one representation of the size and the lap; the panel converts to a
/// `float` at the widget and hands the answer straight back.  That is what
/// makes the stale-radius skew structural rather than remembered: the panel
/// used to hold a float the visitor dragged and a double the last frame copied,
/// and every rule that read the path had to refresh the second from the first
/// or answer about the wrong path.  Two fields that must agree cannot disagree
/// if there is one field.
///
/// **What is NOT here is the phase.**  Since [#30](https://github.com/caliburn-engineering/caliburn/issues/30)
/// it lives in `SimState::path_phase` and `stepSim` advances it through
/// `stepPath`, which is where a fact one frame hands to the next belongs.  So
/// this type does not drive the setpoint — the step does — and `setShape` takes
/// the phase by reference: re-seeding it is the combo's job, owning it is the
/// step's.  That is also why it is not called `TrajectoryDriver`.
class TrajectoryControls {
public:
    /// The path the panel opens on.
    ///
    /// `accel_max` is the plate's and arrives with `opening` — pass
    /// `SimPlate::feasible(...)`, so that the questions this object answers for
    /// itself (the lap floor, the outline the panel draws, the phase a shape
    /// change re-seeds from) are asked about the same filleted path the step
    /// will run.  A default-constructed one is the sharp-cornered path, which
    /// is a reference no plate can follow; it is the geometry tests' case, not
    /// the application's.
    explicit TrajectoryControls(const SetpointPath& opening = SetpointPath{})
        : path_(opening) {}

    // ---------------------------------------------------------------------
    // The four controls
    // ---------------------------------------------------------------------

    /// The trajectory combo.  The new shape picks up nearest to where the
    /// setpoint already is, and `phase` is re-seeded to say so.
    ///
    /// Phase is not comparable across shapes — the circle's phase zero is at +x
    /// and a polygon's first corner is at the top, so equal phase is a quarter
    /// of a lap apart.  Carrying it threw the target **169.7 mm** to the far
    /// side of a 120 mm path and the loop hauled the ball across after it, into
    /// the workspace clip.  It shipped, and it was found by driving the browser.
    /// See `phaseNearest` and #24.
    ///
    /// The lap floor is applied first and the nearest phase found afterwards,
    /// in that order: a different shape is a different perimeter and so a
    /// different floor, and on a polygon the nearest point depends on the
    /// fillet, which depends on the lap.  Re-seeding against the old lap would
    /// be finding the nearest point on a path about to change.
    void setShape(PathShape s, double& phase);

    /// The size slider, in millimetres.
    ///
    /// Raising the size raises the lap floor with it — a bigger path at the
    /// same lap is a faster setpoint — so the lap comes up to meet it in the
    /// same call.  Against the NEW radius, which is the whole of why the two
    /// lines are here rather than at the widget: the frame where the dragged
    /// radius and the copied one differ is exactly the frame a just-enlarged
    /// path would keep the smaller path's floor and run at a speed the bound
    /// exists to forbid.
    void setSizeMm(double mm);

    /// The lap slider, in seconds, held off the floor by `clampPeriod`.
    void setLapS(double seconds);

    /// The setpoint sliders, in millimetres, under a held point.  Under a path
    /// they are disabled and this is not called — the path owns the setpoint,
    /// and `fromReport` is where it gets it.
    void setHeldSetpointMm(double x_mm, double y_mm);

    /// The "Centre setpoint" button.  Also what a reset leaves behind.
    void centreSetpoint() { setpoint_m_.setZero(); }

    // ---------------------------------------------------------------------
    // The frame
    // ---------------------------------------------------------------------

    /// What this frame asks of the setpoint: the path to run, and the point to
    /// hold when there is no path.
    ///
    /// Both together because the step chooses between them, and a caller that
    /// set one and forgot the other would be holding last frame's target.
    void toInput(SimInput& in) const;

    /// What came back.  Under a path the setpoint IS the path's, and the
    /// sliders become the readout of where the ball is being sent — the same
    /// arrangement the servo sliders have under the balance loop.  A held
    /// setpoint is the visitor's and is left alone.
    void fromReport(const SimReport& r);

    // ---------------------------------------------------------------------
    // Readouts
    // ---------------------------------------------------------------------

    /// The path itself, for the free functions that ask it questions — the lap
    /// floor, its length, its fillet, the outline the panel draws.  Those are
    /// pure functions of a path and stay where they are; what lives in this
    /// class is the state they are asked about and the rules that change it.
    const SetpointPath& path() const { return path_; }

    /// True while a shape, rather than the visitor, owns the setpoint.
    bool onAPath() const { return path_.shape != PathShape::Fixed; }

    /// The size and the lap, in the units their sliders are labelled in.
    ///
    /// `double`, like everything else here: ImGui wants a `float` and the panel
    /// casts to one at the widget.  A control's units are the panel's, but a
    /// control's precision is not — the two floats this class replaced are how
    /// the panel came to hold a radius the path disagreed with.
    double sizeMm() const { return path_.radius_m * 1000.0; }
    double lapS() const { return path_.period_s; }

    /// Where the ball is being sent, in metres, in the plate's own frame.
    const Eigen::Vector2d& setpoint() const { return setpoint_m_; }

    double setpointXmm() const { return setpoint_m_(0) * 1000.0; }
    double setpointYmm() const { return setpoint_m_(1) * 1000.0; }

private:
    SetpointPath path_{};
    Eigen::Vector2d setpoint_m_{Eigen::Vector2d::Zero()};
};

}  // namespace caliburn
