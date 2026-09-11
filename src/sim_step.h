// src/sim_step.h
#pragma once

#include "analysis/model_library.h"
#include "auto_balance.h"
#include "ball_contact.h"
#include "ball_sim.h"
#include "setpoint_path.h"
#include "table_kinematics.h"

#include <Eigen/Core>
#include <array>
#include <cmath>

namespace caliburn {

/// One frame of the simulation: setpoint -> loop -> servos -> kinematics ->
/// plate motion -> ball.  The application drives it and so does every test
/// harness, which is the whole of the point.
///
/// **This used to be five copies and it cost real bugs.**  `plate_view`,
/// `test_attract_mode` (three times over), `test_auto_balance` and
/// `test_trajectory` each carried the same eight calls in the same order, and
/// each copy was free to differ:
///
///   - #22's divergence was a guard the harness had and the application did
///     not, which is why an attract kick thrown away mid-flight could hide in
///     a green suite.
///   - Making the plate's accelerations analytic (#23) meant editing the same
///     `plateMotion` call in five places.  Every copy got the edit; three
///     suites went red; nobody could say which of the two normal-force
///     estimators the numbers in CONTEXT.md had been measured against, because
///     no harness pinned a baseline.
///   - And the copies had **already** drifted apart on a quantity nobody was
///     looking at.  The application computed the servo rate from the legs
///     AFTER the frame's step, the harnesses from where they were BEFORE it.
///     Those differ by `exp(-dt/tau)`, which at 60 Hz against a 0.05 s lag is
///     a factor of **1.40** — so every harness had been simulating a plate
///     whose rates, and therefore whose separation behaviour, were forty per
///     cent larger than the shipped application's.  The rate belongs at the
///     instant the pose does, which is after the step; that is what the
///     application did and what this does.
///
/// `test_ball_contact` is deliberately not one of the callers, though #30 lists
/// it as a sixth copy.  It drives `plateMotion` and `stepBallContact` against
/// plates it builds by hand — a plate held still at a tilt, a plate dropped at
/// 3g, a plate whose rates nobody believes — and there is no setpoint, no gain
/// and no servo in any of them.  That is the layer this one is built ON, and
/// building the pieces up from cases where the answer is known by hand is the
/// opposite of carrying a copy of the loop.
///
/// A harness that carries its own copy of the causal order is not evidence
/// about the application, it is evidence about itself.
///
/// See [#30](https://github.com/caliburn-engineering/caliburn/issues/30).

/// The plate a step runs against: the mechanism, the ball rolling on it, and
/// the gravity they both feel.  Built once; a step never changes it.
///
/// The rolling dynamics are built HERE rather than by each caller.  Three
/// call sites each wrote `RollingBallDynamics(kPlateBall, {R_table, g})` out
/// in full, which is three chances for a test to measure a different ball
/// from the one that ships — the same failure mode `cascade_fixture` exists to
/// prevent for the gain.
class SimPlate {
public:
    SimPlate(const TableParams& table, double gravity, double home_leg_rad);

    const TableKinematics& kinematics() const { return tk_; }
    const TableParams& params() const { return tk_.params(); }
    const RollingBallDynamics& rolling() const { return rolling_; }
    double gravity() const { return gravity_; }

    /// The most acceleration this plate can give the ball, in m/s^2.
    ///
    /// Measured once, in the constructor, because it costs a 36-direction tilt
    /// sweep — 2.3 ms on the shipped geometry — and because it cannot change
    /// without the plate changing.  A `SimPlate` IS its geometry.
    double maxBallAccel() const { return max_ball_accel_; }

    /// `p` with its corner fillet sized for THIS plate.
    ///
    /// The one expression of the rule, because it has three callers and a rule
    /// stated three times is a rule that can be forgotten a fourth: the step
    /// below, the application, and any harness that reads the path for itself.
    /// How tightly the setpoint may turn is a fact about the mechanism rather
    /// than about what the visitor asked for, so it is not the caller's to
    /// supply.  See `SetpointPath::accel_max` and #31.
    SetpointPath feasible(SetpointPath p) const {
        p.accel_max = max_ball_accel_;
        return p;
    }

    /// The simulated ball's radius.  One ball — see `kPlateBall`.
    static constexpr double ballRadius() { return kPlateBall.radius; }

private:
    TableKinematics tk_;
    RollingBallDynamics rolling_;
    double gravity_;
    double max_ball_accel_;
};

/// The most acceleration a plate can give a ball, in m/s^2, and so the
/// tightest turn a reference may ask the ball to make.
///
/// `(5/7) * g * sin(theta_max)`: a ball rolling on a plane tilted by
/// `theta_max`, with the 5/7 a solid sphere's rotational inertia takes out of
/// it.  A reference asking for more than this is asking for something no tilt
/// of this plate can produce — infeasible by construction, which is the same
/// principle `kMaxSetpointSpeed` and `kMinLapSeconds` already express and the
/// reason the corner fillet is sized by it (#31).
///
/// **`theta_max` is derived, not chosen.**  It is the largest tilt the plate can
/// actually hold whose velocity Jacobian stays under `kRatesUntrustworthyAbove`
/// in every direction, swept — the threshold this repository has already argued
/// for in `ball_contact.h` and shows in the plate panel as its "Poor" line,
/// rather than a second opinion about when the mechanism is in trouble.  So
/// `a_max` moves with the geometry: 1.52 m/s^2 at 120 mm legs, 1.89 at the
/// shipped 150, 2.25 at 180, 2.47 at 210.
///
/// **On the shipped plate it is the travel that binds, not the conditioning**,
/// which is not what #31 expected — see `max_conditioned_tilt`, where the
/// measurement and the expectation it overturns are both recorded.  The
/// condition gate is still what keeps this honest on a plate with longer legs,
/// where reaching further over means reaching into rates nobody should believe.
///
/// Measured on the shipped plate: `theta_max` = 15.63 degrees, `a_max` = 1.888
/// m/s^2, which fillets the 180 mm square's corners with a 33 mm radius at its
/// fastest offered lap and a 0.6 mm one at its slowest.
double maxBallAccel(const TableKinematics& tk, double gravity,
                    double home_leg_rad);

/// The plate a cascade parameter list describes.
///
/// One transcription, for the application's harnesses and for anything else
/// that wants to run the shipped plant.  Each of the five loop copies used to
/// spell out `cascadeMechanism` -> `TableKinematics` -> `RollingBallDynamics`
/// itself, which is five chances to measure a plate or a ball that is not the
/// one that ships.
SimPlate cascadePlate(const std::vector<PhysicalParam>& params);

/// Everything about the operating point a cascade parameter list fixes —
/// every field of `AutoBalanceDesign` except the gain, which is the designer's
/// answer rather than the plant's.
///
/// `mechanism` and `gravity` are filled in too, so a design built here passes
/// `samePlant` against the plate `cascadePlate` builds from the same list.
/// The application then adds `K` when it has one; a harness adds the K it is
/// measuring.
AutoBalanceDesign cascadeDesign(const std::vector<PhysicalParam>& params);

/// Everything one frame hands to the next.
///
/// The pose is both the plate's assembly and the seed the next solve starts
/// from, because they are the same thing: `solve_pose` has a second root with
/// the table folded flat on the base (#22), and marching the seed forward is
/// what keeps the solver on the branch the plate is actually on.
struct SimState {
    /// Where the legs ARE, in radians — not where they were told to go.
    std::array<double, 3> alpha_rad = {M_PI / 4.0, M_PI / 4.0, M_PI / 4.0};

    TablePose pose{};
    PlateMotion motion{};
    PlateMotion motion_prev{};
    BallState ball{};

    /// How far round the lap the setpoint is, in [0, 1), ACCUMULATED — never
    /// derived from elapsed time.  See `advancePhase` and #24.
    double path_phase = 0.0;
};

/// A plate standing at home with its ball at `ball0`, ready for the first step.
///
/// Both plate-motion frames are filled in, so that the first contact test
/// differences two real instants rather than one real one and a
/// default-constructed zero — which reads as the plate having just been
/// dropped.
SimState simStart(const SimPlate& plate,
                  double home_leg_rad,
                  const Eigen::Vector4d& ball0 = Eigen::Vector4d::Zero());

/// What the plate is being asked for this frame.
struct SimInput {
    double dt = 1.0 / 60.0;

    /// The gain, the operating point it was designed about, and the servo
    /// travel.  `servo_tau` and the travel limits are honoured whether or not
    /// the loop is closed: the legs lag under the sliders too.
    AutoBalanceDesign design{};

    /// True while `design.K` owns the leg command.  False leaves
    /// `open_loop_cmd_rad` in charge — the manual sliders, or the animation.
    bool closed_loop = true;
    std::array<double, 3> open_loop_cmd_rad = {M_PI / 4.0, M_PI / 4.0, M_PI / 4.0};

    /// The path the setpoint runs round.  `PathShape::Fixed` is the held
    /// setpoint the loop has always had, and holds `held_setpoint` — the phase
    /// does not advance under one, so switching to a path resumes where it
    /// left off rather than wherever the clock had got to.
    ///
    /// `path.accel_max` is the exception: the step fills it in from the plate
    /// with `SimPlate::feasible` and ignores whatever is here, because the
    /// corner fillet is what the PLANT can do rather than what the visitor
    /// asked for.  See `SetpointPath::accel_max`.
    SetpointPath path{};
    Eigen::Vector2d held_setpoint{Eigen::Vector2d::Zero()};

    /// Whether the ball is simulated at all.  The plate still moves when it is
    /// not: the pose and the plate motion are what the 3D view draws.
    bool ball_enabled = true;
};

/// What the step did, and what it saw on the way.
///
/// Everything here was computed inside the step anyway.  Handing it back is
/// what stops a caller recomputing it, and recomputing is not merely wasteful:
/// `plateFrame` called again afterwards can be called against the wrong
/// frame's motion, and `normalAccel` called again afterwards is a second
/// opinion about which ball state and which guard it was evaluated under.
struct SimReport {
    /// Where the setpoint was this frame, and how fast it was travelling.
    /// The loop was given these; a caller measuring tracking error must use
    /// them rather than re-reading the path, which would straddle the phase
    /// advance and answer with the wrong edge on a polygon.
    Eigen::Vector2d setpoint{Eigen::Vector2d::Zero()};
    Eigen::Vector2d setpoint_velocity{Eigen::Vector2d::Zero()};

    /// What the legs were asked for, whoever asked.
    std::array<double, 3> cmd_rad = {M_PI / 4.0, M_PI / 4.0, M_PI / 4.0};

    /// The loop's own two complaints, both false while it is not driving.
    bool saturated = false;   ///< a leg command hit a travel limit
    bool clipped = false;     ///< ...or asked for a pose with no assembly

    /// The pose solve, in full: the application prints the iteration count.
    FKResult fk{};

    /// The ball's centre and velocity in the plate frame AFTER the step,
    /// `[x, y, z, vx, vy, vz]`.
    Eigen::Matrix<double, 6, 1> ball_plate{Eigen::Matrix<double, 6, 1>::Zero()};

    bool airborne = false;

    /// The `N/m` this frame's contact was actually decided by, in m/s^2 —
    /// reported by `stepBallContact` rather than reconstructed here.  Positive
    /// held contact, negative ended it, zero means the ball was already flying.
    /// See `stepBallContact` for what it is under an untrusted plate.
    ///
    /// This is the margin the plate is holding the ball by, and a separation
    /// COUNT cannot stand in for it: "never lets go" and "never lets go, with a
    /// sixth of a g to spare" are different claims about the same demo.
    double normal_accel = 0.0;

    /// The ball's contact patch has left the disc.  Reported, never acted on:
    /// the application resets or freezes according to a checkbox, and a
    /// harness wants to record the frame it happened on.  Always false while
    /// `ball_enabled` is false — a ball nobody is simulating cannot fall off.
    bool left_plate = false;
};

/// Advance the plate and the ball by one frame.
///
/// The order is the whole content of this function, and every line of it is
/// answering a bug:
///
///   1. **The setpoint moves first**, so the loop sees this frame's target
///      rather than last frame's.  `stepPath` owns the read-then-advance rule
///      (#24), so it is not re-derived here either.
///   2. **The loop reads the legs where they ARE and the ball where it IS.**
///      While the ball is airborne it is shown where it will LAND instead,
///      with a zero reference velocity — the plate cannot touch a flying ball,
///      so regulating on where it is steers a quantity nothing can move.  An
///      airborne ball is the state every one of these copies had to remember
///      to ask about, and #22 is what forgetting costs: the attract kick was
///      thrown away mid-flight, and it hid because the harness asked and the
///      application did not.
///   3. **Then the servos move**, clipped at the workspace boundary and not
///      merely at each leg's travel: the servo limits are a box and the
///      workspace is not — and clipped short of its singular fringe too, since
///      a second assembly waits on the far side of one (#29).
///   4. **Then the pose**, seeded from the last one and adopted only on
///      success — a leg triple with no assembly leaves the plate where it was,
///      which is what a mechanism does when it binds.
///   5. **Then the plate's motion**, from the servo rate at the instant the
///      pose is for, and with the lag's own analytic second derivative rather
///      than a difference of a stepping command (#23).
///   6. **Then the ball**, which is the only thing here that can decide to
///      leave.
SimReport stepSim(const SimPlate& plate, const SimInput& in, SimState& s);

}  // namespace caliburn
