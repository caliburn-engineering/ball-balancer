// src/sim_step.cpp
#include "sim_step.h"

namespace caliburn {

SimPlate::SimPlate(const TableParams& table, double gravity)
    : tk_(table),
      rolling_(kPlateBall, PlateParams{table.R_table, gravity}),
      gravity_(gravity) {}

SimPlate cascadePlate(const std::vector<PhysicalParam>& params) {
    return SimPlate(cascadeMechanism(params), cascadeGravity(params));
}

AutoBalanceDesign cascadeDesign(const std::vector<PhysicalParam>& params) {
    AutoBalanceDesign d;
    d.home_leg_rad = cascadeHomeLegAngle(params);
    d.servo_tau = cascadeServoTau(params);
    d.mechanism = cascadeMechanism(params);
    d.gravity = cascadeGravity(params);
    d.alpha_min_rad = d.mechanism.alpha_min;
    d.alpha_max_rad = d.mechanism.alpha_max;
    return d;
}

SimState simStart(const SimPlate& plate,
                  double home_leg_rad,
                  const Eigen::Vector4d& ball0) {
    SimState s;
    s.alpha_rad = {home_leg_rad, home_leg_rad, home_leg_rad};
    s.pose = plate.kinematics().home_pose(home_leg_rad);
    s.motion = plateMotion(plate.kinematics(), s.pose, s.alpha_rad,
                           {0.0, 0.0, 0.0});
    s.motion_prev = s.motion;
    s.ball = BallState{};
    s.ball.rolling = ball0;
    s.path_phase = 0.0;
    return s;
}

SimReport stepSim(const SimPlate& plate, const SimInput& in, SimState& s) {
    const TableKinematics& tk = plate.kinematics();
    const double r_ball = SimPlate::ballRadius();
    const double g = plate.gravity();

    SimReport out;

    // --- 1. The setpoint ---
    //
    // Driven before the loop reads it, so the gain sees this frame's target.
    // `stepPath` owns the read-then-advance rule, so the position and the
    // velocity handed over are the same instant — on a polygon a phase
    // advanced in between would straddle a corner and answer with the wrong
    // edge.
    //
    // A held setpoint does not go through `stepPath` at all: its point is the
    // one the visitor dragged, not the origin `pathPoint` returns for `Fixed`,
    // and its phase must not move — a lap the setpoint is not running is a lap
    // that should still be where it was left when a shape is chosen.
    PathStep path_step;
    if (in.path.shape == PathShape::Fixed) {
        path_step.point = in.held_setpoint;
        path_step.next_phase = s.path_phase;
    } else {
        path_step = stepPath(in.path, s.path_phase, in.dt);
    }
    s.path_phase = path_step.next_phase;
    out.setpoint = path_step.point;
    out.setpoint_velocity = path_step.velocity;

    // --- 2. The leg command ---
    //
    // Three writers, one array, in strict precedence: the closed loop, then
    // whatever the caller last left in `open_loop_cmd_rad` — the animation, or
    // the sliders.  The loop reads the legs where they ARE and the ball where
    // it IS, and the servos move afterwards.
    std::array<double, 3> cmd_rad = in.open_loop_cmd_rad;
    if (in.closed_loop) {
        // What the gain is shown.  While the ball is on the plate this is just
        // the ball.  While it is in the air the plate cannot touch it, so
        // regulating where it IS steers on a quantity nothing can move — the
        // loop is given where it will LAND instead, and spends the flight
        // getting underneath it.
        //
        // The reference velocity is zeroed with it, for the same reason: the
        // ball's present velocity is what carries it to that landing point,
        // and feeding both would ask the plate to cancel a motion it has
        // already accounted for.
        const Eigen::Matrix<double, 6, 1> bp =
            plateFrame(s.ball, s.motion, r_ball);
        Eigen::Vector4d seen(bp(0), bp(1), bp(3), bp(4));

        BallReference ref;
        ref.position = path_step.point;
        if (s.ball.airborne) {
            const Eigen::Vector2d land = predictedLanding(bp, r_ball, g);
            seen << land(0), land(1), 0.0, 0.0;
        } else {
            ref.velocity = path_step.velocity;
        }

        const LegCommand c = legCommand(tk, in.design, s.alpha_rad, seen, ref);
        cmd_rad = c.alpha_rad;
        out.saturated = c.saturated;
        out.clipped = c.clipped_to_workspace;
    }
    out.cmd_rad = cmd_rad;

    // --- 3. The servos ---
    //
    // Stopped at the workspace boundary, not merely at each leg's travel:
    // `stepServos` lags each leg independently, so its result sits on the
    // straight line between two triples that can both be assemblable while
    // points along it are not.  See #22.
    s.alpha_rad = stepServosOnPlate(tk, s.alpha_rad, cmd_rad,
                                    in.design.servo_tau, in.dt);

    // --- 4. The pose ---
    //
    // `solve_pose`, not `forward_kinematics`: the constraint equations have a
    // second root with the table folded flat on the base, and a convergence
    // test alone passes it — with a residual of 1e-11, in one iteration.  Once
    // the pose fell onto that root the next frame was seeded from it and the
    // plate stayed flat for the rest of the session (#22).
    //
    // Adopted only on success.  A leg triple with no assembly leaves the plate
    // where it was, which is what a mechanism does when it is driven into a
    // bind: it stops.
    out.fk = tk.solve_pose(s.alpha_rad, s.pose);
    if (out.fk.converged) s.pose = out.fk.pose;

    // --- 5. The plate's motion ---
    //
    // The legs' rate is the servo lag's own derivative rather than a
    // difference: `alpha_dot = (cmd - alpha) / tau` is exactly what the model
    // says they are doing, and it costs nothing to ask it.
    //
    // Taken at the legs' position AFTER the step, because that is the instant
    // `pose` is for.  The five copies of this loop disagreed about that, and
    // the difference is not small: `(cmd - alpha)` shrinks by `exp(-dt/tau)`
    // across the frame, a factor of 1.40 at 60 Hz against a 0.05 s lag, and it
    // multiplies straight through to `omega`, to `omega_dot`, and so to
    // whether the ball leaves the plate at all.
    std::array<double, 3> alpha_dot{};
    if (in.design.servo_tau > 0.0) {
        for (int i = 0; i < 3; ++i)
            alpha_dot[i] = (cmd_rad[i] - s.alpha_rad[i]) / in.design.servo_tau;
    }
    s.motion_prev = s.motion;
    s.motion = plateMotion(tk, s.pose, s.alpha_rad, alpha_dot,
                           servoAccel(alpha_dot, in.design.servo_tau),
                           &s.motion_prev, in.dt);

    // --- 6. The ball ---
    if (in.ball_enabled) {
        s.ball = stepBallContact(plate.rolling(), s.ball, s.motion,
                                 s.motion_prev, s.pose, r_ball, g, in.dt,
                                 &out.normal_accel);
    }
    out.airborne = s.ball.airborne;
    out.ball_plate = plateFrame(s.ball, s.motion, r_ball);

    // Reported, never acted on.  A ball nobody is simulating cannot fall off,
    // so the question is not even asked of one.
    if (in.ball_enabled) {
        out.left_plate = !ballOnPlate(
            Eigen::Vector4d(out.ball_plate(0), out.ball_plate(1),
                            out.ball_plate(3), out.ball_plate(4)),
            plate.params().R_table, r_ball);
    }

    return out;
}

}  // namespace caliburn
