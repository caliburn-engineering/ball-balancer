// src/attract_mode.h
#pragma once

#include "setpoint_path.h"

#include <Eigen/Core>

namespace caliburn {

/// The demo driving itself for a visitor who has not arrived yet.
///
/// A motionless canvas is indistinguishable from a broken build, and a ball
/// sitting balanced reads as a still image of one.  So the application opens
/// with the loop already closed and the ball already tracing a circle: motion
/// that is continuous, obviously deliberate, and obviously being *produced*
/// rather than replayed.
///
/// **This used to be a disturbance schedule, and is not any more.**  The demo
/// opened with the ball displaced 72 mm and kicked it every four seconds, on
/// the argument that what distinguishes a control system from an animation is
/// *recovery*: something disturbs the ball and the loop puts it back.  The
/// argument is sound and the execution read as a fault.  Two reasons, and the
/// first is the one that matters:
///
/// A kick is legible only against a still baseline.  A ball parked at the
/// centre, nudged, returning — that is a story.  But the still moment between
/// kicks is also the moment the demo looks like a static image, so the demo
/// spent most of its time looking broken in order to make the other part
/// legible.  A circle needs no baseline: at every instant the ball is
/// somewhere it was told to be, and the plate is visibly working to keep it
/// there.  The loop is not demonstrating that it *can* respond, it is
/// demonstrating that it *is* responding, continuously.
///
/// The second reason is the opening displacement specifically.  Handing a
/// state feedback a 72 mm error at t = 0, with the legs exactly at home and no
/// servo history, is a step input: measured, the legs swung **20.7 degrees
/// apart within three frames** and the table dropped 4.5 mm before settling
/// inside 250 ms.  Correct, and it looks like the mechanism glitching.
/// Starting the ball ON the path at the path's own velocity makes both the
/// position and the velocity error zero at t = 0, and the same measurement
/// gives a 3.0 degree opening swing with the table height not moving at all.
///
/// What is given up is stated plainly, because it was a real acceptance
/// criterion: the demo no longer shows disturbance rejection unprompted.  The
/// visitor can still see it — the Nudge buttons are right there — and the
/// property is still pinned by `test_attract_mode`, which sweeps a 0.26 m/s
/// disturbance over every direction.  It is no longer *performed*.
///
/// See [#17](https://github.com/caliburn-engineering/caliburn/issues/17) for
/// the opening, and
/// [#24](https://github.com/caliburn-engineering/caliburn/issues/24) for the
/// paths it now opens on.

/// The path the demo opens tracing.
///
/// The circle rather than a cornered shape: it is the one whose tracking error
/// is smooth, so the opening reads as competence rather than as the ball
/// stumbling at every corner.  The corners are a better demonstration, and are
/// one dropdown away — but they are an argument the visitor should choose to
/// hear, not the first thing they see.
///
/// 120 mm on a 300 mm plate at a ten-second lap, which is 75 mm/s — a third of
/// the speed cap, and slow enough that the ball keeps up with 3.7 mm of mean
/// error.  Sized and paced to look effortless, because the opening's job is to
/// say "this is running", not to find the limits.
SetpointPath openingPath();

/// The ball state the demo opens on: **on** the path, moving **with** it.
///
/// Both halves matter and for the same reason.  On the path, so the position
/// error is zero; moving with it, so the velocity error is too — the reference
/// state carries the path's own velocity (see `BallReference` in
/// `auto_balance.h`), and a ball placed correctly but at rest would still hand
/// the loop a 75 mm/s error to answer on the first frame.
///
/// The result is an opening frame the controller has nothing to do about,
/// which is exactly what makes it look calm.
Eigen::Vector4d attractStart(const SetpointPath& path);

}  // namespace caliburn
