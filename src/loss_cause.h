// src/loss_cause.h
#pragma once

#include "sim_step.h"

#include <string>

namespace caliburn {

/// Why the demo lost the ball, named once, so that the most informative moment
/// in the demo is not a blink.
///
/// `resetBall()` used to teleport the ball back to centre with no banner, no
/// marker and no cause, and `ball_auto_reset_` defaulted true — so the loop
/// losing the ball rendered as the ball appearing in the middle again.  That is
/// #22's failure mode exactly: the demo looked fine and was not.  See
/// [#33](https://github.com/caliburn-engineering/caliburn/issues/33) and D13/D14
/// of `docs/plans/2026-09-09-bouncing-ball-control-decisions.md`.

/// The raw flags on the frame the ball crossed the rim.
///
/// **The frame itself, with no lookback, and that is measured rather than
/// assumed.**  The obvious worry is that the loop gives up long before the ball
/// is actually off the edge, so the flags at the crossing would have gone quiet
/// again and every loss would read as "rolled off".  It does not happen: over
/// 72 directions of a 1.30 m/s shove under Nominal — the softest tuning that
/// loses the ball at all — `clipped` and `saturated` are BOTH still raised on
/// the crossing frame in 72 of 72, and widening the window to 0.25, 0.5, 1.0 or
/// 2.0 s changes not one of the four counts.  The reason is mechanical: a ball
/// out at the rim is a large error, and a large error is a command the
/// mechanism cannot hold, right up to the frame it leaves on.
///
/// So there is no window, no latch and no decay constant to defend.  A rule
/// that reads one frame is one a test can put in a known state.
struct LossFlags {
    /// The plate's velocity Jacobian was too ill-conditioned for anything
    /// computed from its motion to be believed — `PlateMotion::rates_trustworthy`.
    bool rates_untrustworthy = false;

    /// The loop asked for a pose the mechanism will not hold, and was pulled
    /// back to one it will — `SimReport::clipped`.
    bool clipped = false;

    /// The ball was in the air: the plate had already let go of it —
    /// `SimReport::airborne`.
    bool separated = false;

    /// A leg command stood on a travel limit — `SimReport::saturated`.
    bool saturated = false;
};

/// The five causes, in the order `lossCause` resolves them.
enum class LossCause {
    RatesUntrustworthy,
    WorkspaceClipped,
    Separated,
    Saturated,
    RolledOff,
};

/// One cause, always, by a precedence derived from CAUSALITY rather than ranked
/// by severity: each entry is upstream of the next.
///
///   1. **rates untrustworthy** — the arithmetic cannot be believed at all, so
///      none of the three flags below is evidence of anything.
///   2. **workspace-clipped** — the mechanism refused the command, so what the
///      servos then did is a consequence rather than a cause.
///   3. **separated** — the plate left the ball, so where the ball went next was
///      not the plate's to decide.
///   4. **saturated** — the servo ran out of travel.
///   5. **rolled off** — nothing failed; the ball simply went too far.
///
/// **Why five and not three.**  Measured against the closed loop on `ce522fa`,
/// entries 1 and 3 could not be reached at any nudge the interface composes, and
/// collapsing the rule to the three that could was a live option — #33's own
/// comment puts it as one.  Re-measured on this tree, after #23 gave the landing
/// a coefficient of restitution, over 72 kick directions and 14.5 s per run:
///
/// | tuning | kick | untrusted | clipped | saturated | separated | lost |
/// |---|---|---|---|---|---|---|
/// | Detuned | 0.20 | 0 | 77 | 0 | 0 | 0 |
/// | Nominal | 0.20 | 0 | 502 | 92 | 0 | 0 |
/// | Aggressive | 0.20 | 0 | 606 | 72 | 0 | 0 |
/// | Nominal | 0.30 | 0 | 779 | 397 | **28** | 0 |
/// | Aggressive | 0.90 | **100** | 3287 | 2909 | 6 | 0 |
/// | Nominal | 1.30 | 0 | 1363 | 1147 | 0 | **72** |
///
/// `kMaxNudgeSpeed` is 0.20, so the first three rows are the whole envelope the
/// two Nudge buttons can compose — and in it the loop keeps the ball, clips, and
/// sometimes saturates.  Neither 1 nor 3 rises there.  `separated` first appears
/// at 0.30, half again the ceiling; before #23's bounce it was zero at every
/// speed measured to 0.90, so restitution moved it from unreachable to nearby.
/// `rates_untrustworthy` needs 0.90, and it is scarce even there: `legCommand`
/// and `stepServosOnPlate` both retreat to poses whose velocity Jacobian this
/// repository will vouch for (#29), so the flag survives only where the marched
/// pose and the analytic one have parted company.
///
/// **What keeps all five is that the closed loop is not the only driver.**  The
/// banner is now mostly for the hand-driven paths — the servo sliders and the
/// model panel — and those reach 3 and 5 easily: measured, a plate slammed to
/// two legs high and one low with the ball already running at 0.6 m/s loses it
/// on frame 7, in the air, which is `separated` at the crossing frame.  The same
/// slam from a ball at rest loses it on frame 49 with nothing raised at all,
/// which is `rolled off` and is the correct answer — nothing failed, the visitor
/// tilted the plate.  Entry 1 is the one no harness in this repo reaches through
/// the application; it is also the one that says *do not believe the other
/// four*, so a rule that dropped it would have no way to say the thing that
/// most needs saying.
///
/// The precedence is a claim about mechanism, so it is pinned by a test over all
/// sixteen flag combinations rather than sampled through the loop — see
/// `test_loss_cause`.
LossCause lossCause(const LossFlags& f);

/// The flags a finished step leaves behind.  Here rather than at the call site
/// because which field means which cause is part of the rule: `separated` is
/// `airborne`, and a panel deciding that for itself is a second opinion about
/// what the step reported.
///
/// `rates_untrustworthy` comes from the STATE rather than the report, because
/// that is where the plate's motion lives — the step does not copy it out, and
/// a copy is one more thing that can be taken from the wrong frame.
LossFlags lossFlags(const SimReport& r, const SimState& s);

/// Two or three words, for the line the visitor reads first.
const char* lossCauseLabel(LossCause c);

/// One sentence saying what that cause means, for the hover.  The label has to
/// fit a phone; the sentence does not have to.
const char* lossCauseSentence(LossCause c);

/// The raw flags, named and comma-separated, or "none".
///
/// **One sentence that can be wrong is better than five flags that cannot be
/// read** — but the flags are not thrown away for that.  They ride along in the
/// plot marker, so that the full picture survives for anyone who looks, and so
/// that a visitor who disagrees with the sentence can see what it was decided
/// from.
std::string lossFlagsLine(const LossFlags& f);

}  // namespace caliburn
