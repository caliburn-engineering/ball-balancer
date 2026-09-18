// src/hop_drive.cpp
#include "hop_drive.h"

#include <algorithm>

namespace caliburn {

HopSizing hopSizing(double gravity, double servo_tau, double room_m) {
    HopSizing s;
    s.rise_rate = kHopRiseMargin * gravity * servo_tau;
    s.stroke_s = kHopStrokeLags * servo_tau;
    // The strokes RAMP, so what the charge spends is the area under a triangle
    // and not under a rectangle.  Multiplying it out here rather than at each
    // reader is what keeps the two from disagreeing about the factor of two.
    s.drop_m = 0.5 * s.rise_rate * s.stroke_s;
    // The first flight is `2 u / g`; the train after it is `1 / (1 - e)` times
    // that, because each bounce returns `e` of the speed and so `e` of the
    // time.  Both from the ONE restitution this application has — `kRestitution`
    // — rather than from a second number that would be free to disagree with
    // the bounce the ball actually takes.
    s.flight_s = 2.0 * s.rise_rate / (gravity * (1.0 - kRestitution));
    s.max_ball_speed = (s.flight_s > 0.0) ? room_m / s.flight_s : 0.0;
    return s;
}

HopDemand stepHop(const HopSizing& size, HopCycle& cycle,
                  const HopFrame& frame, double dt) {
    if (!frame.enabled) {
        cycle = HopCycle{};
        return HopDemand{};
    }

    if (frame.airborne) {
        if (cycle.phase != HopPhase::Flight) {
            // A throw that let go is the only one that counts.  Counting the
            // window opening instead would report hops the plate never managed
            // — which is the thing this number exists to make visible.
            if (cycle.phase == HopPhase::Throw) ++cycle.thrown;
            cycle.phase = HopPhase::Flight;
            cycle.phase_s = 0.0;
        }
        cycle.phase_s += dt;
        return HopDemand{};   // the flight is #23's
    }

    switch (cycle.phase) {
        case HopPhase::Off:
        case HopPhase::Waiting:
        case HopPhase::Flight:
            // The ball is on the plate, so a throw is there to be armed: the
            // first frame after the hop is switched on, and every landing after
            // that.  Waiting is free — the tracking command is unchanged and
            // the ball is on the plate — and it is what keeps a hop from
            // coasting the ball off the rim or compounding a loop that is
            // already behind.  See the header.
            if (!frame.loop_has_room || frame.ball_speed > size.max_ball_speed) {
                cycle.phase = HopPhase::Waiting;
                cycle.phase_s = 0.0;
                return HopDemand{};
            }
            // Arming records the heave, and the whole hop happens between it
            // and `drop_m` below it.
            cycle.phase = HopPhase::Charge;
            cycle.phase_s = 0.0;
            cycle.z_arm = frame.plate_z;
            break;

        case HopPhase::Charge:
            // Down far enough, or out of window.  The window is the backstop
            // for a plate that cannot descend — a leg already on its lower
            // stop, or a command retreated back toward the level pose — so that
            // such a hop throws what it has rather than charging for ever.
            if (frame.plate_z <= cycle.z_arm - size.drop_m ||
                cycle.phase_s >= size.stroke_s) {
                cycle.phase = HopPhase::Throw;
                cycle.phase_s = 0.0;
            }
            break;

        case HopPhase::Throw:
            // Back at the heave it started from with the ball still on the
            // plate, or out of window getting there.  Start the cycle again
            // rather than keep pushing: whatever kept the ball down is still
            // there, and a plate parked high has spent the travel tilt wants
            // for nothing.  `z_arm` is deliberately NOT re-read — that is what
            // keeps a run of failed throws from marching the plate down.
            if (frame.plate_z >= cycle.z_arm || cycle.phase_s >= size.stroke_s) {
                cycle.phase = HopPhase::Charge;
                cycle.phase_s = 0.0;
            }
            break;
    }

    cycle.phase_s += dt;

    // **The stroke is a ramp, and that is the whole difference between the
    // plate carrying the ball off and the plate hitting it.**  See
    // `HopSizing::stroke_s` for the two measurements.
    const double ramp = std::min(1.0, cycle.phase_s / size.stroke_s);
    const double rate = ramp * size.rise_rate;
    return HopDemand{true, (cycle.phase == HopPhase::Throw) ? rate : -rate};
}

}  // namespace caliburn
