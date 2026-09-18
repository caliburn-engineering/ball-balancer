// tests/test_loss_cause.cpp
//
// The precedence that names one cause for a lost ball.  It is a claim about
// MECHANISM — rates untrustworthy is upstream of workspace-clipped, which is
// upstream of separated, which is upstream of saturated — and a claim about
// mechanism is exactly the kind that survives being stated in prose and quietly
// stops being true in the code.
//
// **Two harnesses, because the rule and the demo are two different claims.**
//
//   - The rule is a pure function of four booleans, so all sixteen inputs are
//     enumerable and every branch is reachable.  That matters: measured through
//     the closed loop, `rates_untrustworthy` does not rise at any nudge the
//     interface composes, so a test that drove only the demo would leave the
//     top of the precedence unpinned and nobody would notice it had broken.
//     See `lossCause`.
//   - The demo is what actually happens, and `lossFlags` is the seam between
//     the two.  So a second harness drives `stepSim` to a real loss and asks
//     what the visitor would be shown.  Without it the rule could be perfect
//     and read from the wrong fields.
#include "cascade_fixture.h"
#include "loss_cause.h"
#include "sim_step.h"
#include "test_helpers.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace caliburn;

namespace {

constexpr double kDeg = M_PI / 180.0;

LossFlags flags(bool untrusted, bool clipped, bool separated, bool saturated) {
    LossFlags f;
    f.rates_untrustworthy = untrusted;
    f.clipped = clipped;
    f.separated = separated;
    f.saturated = saturated;
    return f;
}

// --- The rule, over every input it has ---

// Sixteen combinations, each asserted against the causal order rather than
// against a table transcribed from it: the expected answer is recomputed here
// from "the first flag that is set, in order", which is what the precedence
// says.  A table would agree with `lossCause` by construction if both were
// edited together, which is the failure this is meant to catch.
void test_the_precedence_names_the_most_upstream_flag_that_is_set() {
    for (int bits = 0; bits < 16; ++bits) {
        const bool un = bits & 1, cl = bits & 2, se = bits & 4, sa = bits & 8;
        const LossFlags f = flags(un, cl, se, sa);

        const LossCause expected =
            un ? LossCause::RatesUntrustworthy :
            cl ? LossCause::WorkspaceClipped   :
            se ? LossCause::Separated          :
            sa ? LossCause::Saturated          :
                 LossCause::RolledOff;

        if (lossCause(f) != expected) {
            std::fprintf(stderr,
                         "FAIL: flags un=%d cl=%d se=%d sa=%d -> %s, expected %s\n",
                         un, cl, se, sa, lossCauseLabel(lossCause(f)),
                         lossCauseLabel(expected));
            std::exit(1);
        }
    }
}

// Each entry masks every entry below it, which is the whole content of "each is
// upstream of the next".  Stated separately from the sweep above because it is
// the property the ORDER exists for, and a sweep that agreed by accident would
// still pass this only if the order is right.
void test_each_cause_masks_the_ones_downstream_of_it() {
    ASSERT_TRUE(lossCause(flags(true, true, true, true)) ==
                LossCause::RatesUntrustworthy);
    ASSERT_TRUE(lossCause(flags(false, true, true, true)) ==
                LossCause::WorkspaceClipped);
    ASSERT_TRUE(lossCause(flags(false, false, true, true)) ==
                LossCause::Separated);
    ASSERT_TRUE(lossCause(flags(false, false, false, true)) ==
                LossCause::Saturated);
    ASSERT_TRUE(lossCause(flags(false, false, false, false)) ==
                LossCause::RolledOff);
}

// Nothing failed is a cause too.  The banner always says something, so
// "rolled off" has to be a real answer rather than the absence of one.
void test_no_flag_at_all_still_names_a_cause() {
    const LossCause c = lossCause(LossFlags{});
    ASSERT_TRUE(c == LossCause::RolledOff);
    ASSERT_TRUE(std::string(lossCauseLabel(c)) == "rolled off");
}

// Every cause is sayable.  A branch that returns an empty label is a banner
// that names nothing, and it would pass every assertion above.
void test_every_cause_has_a_label_and_a_sentence() {
    const LossCause all[] = {
        LossCause::RatesUntrustworthy, LossCause::WorkspaceClipped,
        LossCause::Separated, LossCause::Saturated, LossCause::RolledOff,
    };
    for (LossCause c : all) {
        ASSERT_TRUE(std::string(lossCauseLabel(c)).size() > 0);
        // Short enough to sit on one line of the control panel on a phone,
        // which is what #25 leaves this ticket responsible for.  The sentence
        // is the hover, and has room to be a sentence.
        ASSERT_TRUE(std::string(lossCauseLabel(c)).size() <= 20);
        ASSERT_TRUE(std::string(lossCauseSentence(c)).size() > 20);
    }
}

// The flags survive the summary.  One sentence that can be wrong is better than
// five flags that cannot be read, but only if the five are still there.
void test_the_raw_flags_are_all_reported_beside_the_cause() {
    const std::string all = lossFlagsLine(flags(true, true, true, true));
    ASSERT_TRUE(all.find("rates untrustworthy") != std::string::npos);
    ASSERT_TRUE(all.find("clipped") != std::string::npos);
    ASSERT_TRUE(all.find("separated") != std::string::npos);
    ASSERT_TRUE(all.find("saturated") != std::string::npos);

    // And a loss with nothing raised says so, rather than being blank — a blank
    // line in the marker reads as a marker that failed to record anything.
    ASSERT_TRUE(lossFlagsLine(LossFlags{}) == "none");
}

// --- The seam: which report field is which flag ---

// `lossFlags` is the only place the mapping lives, and every field of it can be
// read from the wrong place without the rule above noticing.  `separated` is
// `airborne` in particular: they are different words for the same fact, and
// `SimReport` has a `left_plate` right beside it that means something else
// entirely.
void test_the_flags_are_read_from_the_fields_they_name() {
    SimReport r;
    SimState s;
    r.clipped = true;
    r.saturated = false;
    r.airborne = true;
    r.left_plate = true;          // the loss itself, and NOT one of the flags
    s.motion.rates_trustworthy = false;

    const LossFlags f = lossFlags(r, s);
    ASSERT_TRUE(f.clipped);
    ASSERT_TRUE(!f.saturated);
    ASSERT_TRUE(f.separated);
    ASSERT_TRUE(f.rates_untrustworthy);

    // Trustworthy rates are the DEFAULT, so the negation has to be the right
    // way round: a sign slip here names every loss in the demo "rates
    // untrustworthy" and masks the other four causes permanently.
    s.motion.rates_trustworthy = true;
    ASSERT_TRUE(!lossFlags(r, s).rates_untrustworthy);
}

// --- The demo ---

// A real loss, through the same step the browser runs, and the cause the
// visitor would be shown.
//
// 1.30 m/s rather than something the Nudge buttons can compose, and that is the
// measurement rather than a shortcut: under the shipped tunings this loop does
// not lose the ball below about 1.20 m/s, and `kMaxNudgeSpeed` is 0.20.  The
// demo the visitor meets does not lose the ball — the banner exists for the
// hand-driven paths and for whatever a future tuning does — so reaching a loss
// at all means shoving harder than the interface will.
//
// What is pinned is the BRANCH: a loop that loses the ball while clipping
// against the holdable set says so, and does not say "rolled off".
void test_a_loss_under_the_loop_is_named_workspace_clipped() {
    const auto models = getBuiltinModels();
    const ModelEntry& e = cascadeModel(models);
    const SimPlate plate = cascadePlate(e.params);

    SimInput in;
    in.design = cascadeDesign(e.params);
    in.design.K = gainForPreset(e, presetNamed("Nominal"));

    int losses = 0;
    const int kDirections = 24;
    for (int d = 0; d < kDirections; ++d) {
        const double theta = 2.0 * M_PI * d / kDirections;
        SimState s = simStart(plate, in.design.home_leg_rad,
                              Eigen::Vector4d(0.06, -0.04, 0.0, 0.0));
        const int settle = static_cast<int>(2.5 / in.dt);
        const int total = settle + static_cast<int>(12.0 / in.dt);
        for (int k = 0; k < total; ++k) {
            if (k == settle && !s.ball.airborne) {
                s.ball.rolling(2) += 1.30 * std::cos(theta);
                s.ball.rolling(3) += 1.30 * std::sin(theta);
            }
            const SimReport f = stepSim(plate, in, s);
            if (!f.left_plate) continue;

            ++losses;
            const LossFlags fl = lossFlags(f, s);
            // Both are raised on the crossing frame, every direction — this is
            // the measurement `LossFlags` cites for having no lookback window.
            ASSERT_TRUE(fl.clipped);
            ASSERT_TRUE(fl.saturated);
            // And clipped is upstream, so that is the one word the banner gets.
            ASSERT_TRUE(lossCause(fl) == LossCause::WorkspaceClipped);
            break;
        }
    }
    ASSERT_EQ(losses, kDirections);
}

// The other half of the contract, and the half that is about the demo rather
// than about the banner: what the interface actually offers does NOT lose the
// ball.  `kMaxNudgeSpeed` is the hardest shove the two Nudge buttons compose,
// and Nominal keeps the ball through it from every direction — so no visitor
// meets this banner by pressing what is in front of them.
//
// Measured on the way past: inside that envelope the loop clips and sometimes
// saturates, and it never separates and never distrusts its own rates.  That is
// the fact that makes the top of the precedence unreachable through the loop,
// and the reason the hand-driven cases below exist.
void test_the_offered_envelope_does_not_lose_the_ball() {
    const auto models = getBuiltinModels();
    const ModelEntry& e = cascadeModel(models);
    const SimPlate plate = cascadePlate(e.params);

    SimInput in;
    in.design = cascadeDesign(e.params);
    in.design.K = gainForPreset(e, presetNamed("Nominal"));

    int clipped_frames = 0, separated_frames = 0, untrusted_frames = 0;
    const int kDirections = 72;
    for (int d = 0; d < kDirections; ++d) {
        const double theta = 2.0 * M_PI * d / kDirections;
        SimState s = simStart(plate, in.design.home_leg_rad,
                              Eigen::Vector4d(0.06, -0.04, 0.0, 0.0));
        const int settle = static_cast<int>(2.5 / in.dt);
        const int total = settle + static_cast<int>(12.0 / in.dt);
        for (int k = 0; k < total; ++k) {
            if (k == settle && !s.ball.airborne) {
                s.ball.rolling(2) += kMaxNudgeSpeed * std::cos(theta);
                s.ball.rolling(3) += kMaxNudgeSpeed * std::sin(theta);
            }
            const SimReport f = stepSim(plate, in, s);
            ASSERT_TRUE(!f.left_plate);
            const LossFlags fl = lossFlags(f, s);
            if (fl.clipped) ++clipped_frames;
            if (fl.separated) ++separated_frames;
            if (fl.rates_untrustworthy) ++untrusted_frames;
        }
    }
    ASSERT_TRUE(clipped_frames > 0);
    ASSERT_EQ(separated_frames, 0);
    ASSERT_EQ(untrusted_frames, 0);
}

// The two rows of `lossCause`'s table that the precedence's shape rests on, and
// they were prose until this test existed.
//
// D15 of the decision record is about exactly this failure — "the 3-of-24 number
// has gated two sessions and exists only in prose" — and the argument for
// keeping entries 1 and 3 is the same kind of claim: that they are reachable
// ABOVE the offered envelope even though nothing inside it raises them.  If they
// stop being reachable, the honest response is to re-argue the precedence, and
// that cannot happen if nothing notices.
//
// Both are asserted as "> 0" rather than at the measured counts.  The counts (28
// frames and 100 frames over 72 directions) are a property of the tuning, the
// contact model and the corner fillet all at once, and pinning them exactly
// would fail on every change to any of the three for no reason anybody could
// act on.  Reachability is the claim being made; reachability is what is pinned.
void test_the_two_scarce_causes_are_reachable_above_the_offered_envelope() {
    const auto models = getBuiltinModels();
    const ModelEntry& e = cascadeModel(models);
    const SimPlate plate = cascadePlate(e.params);

    struct Row { const char* preset; double kick; };
    // Nominal at 0.30 is where the plate still takes the ball off the surface —
    // #23's restitution gave that back, and the pre-bounce measurement on #33
    // recorded zero at every speed out to 0.90.  Aggressive at 0.90 is where the
    // marched pose and the analytic one part company far enough for the plate to
    // stop vouching for its own rates.
    const Row rows[] = {{"Nominal", 0.30}, {"Aggressive", 0.90}};

    int separated_frames = 0, untrusted_frames = 0;
    for (const Row& r : rows) {
        SimInput in;
        in.design = cascadeDesign(e.params);
        in.design.K = gainForPreset(e, presetNamed(r.preset));

        const int kDirections = 72;
        for (int d = 0; d < kDirections; ++d) {
            const double theta = 2.0 * M_PI * d / kDirections;
            SimState s = simStart(plate, in.design.home_leg_rad,
                                  Eigen::Vector4d(0.06, -0.04, 0.0, 0.0));
            const int settle = static_cast<int>(2.5 / in.dt);
            const int total = settle + static_cast<int>(12.0 / in.dt);
            for (int k = 0; k < total; ++k) {
                if (k == settle && !s.ball.airborne) {
                    s.ball.rolling(2) += r.kick * std::cos(theta);
                    s.ball.rolling(3) += r.kick * std::sin(theta);
                }
                const SimReport f = stepSim(plate, in, s);
                const LossFlags fl = lossFlags(f, s);
                if (fl.separated) ++separated_frames;
                if (fl.rates_untrustworthy) ++untrusted_frames;
                if (f.left_plate) break;
            }
        }
    }
    ASSERT_TRUE(separated_frames > 0);
    ASSERT_TRUE(untrusted_frames > 0);
}

// --- What this rule CANNOT do, pinned so that fixing it is noticed ---

// **Under the closed loop the precedence names the same cause every time, and
// this test exists to say so out loud.**
//
// Measured after the rule shipped: across three tunings and shove speeds from
// 1.20 to 5.00 m/s, 1,080 losses over 72 directions apiece, every single one
// carried the identical flag combination — `clipped` and `saturated`, with
// `rates_untrustworthy` and `separated` clear.  So `lossCause` returns
// `WorkspaceClipped` for all of them, and NO reordering of the precedence could
// do better, because the input does not vary.  Moving the anchor earlier does
// not help either: `clipped` is already raised by the time the ball is 75 mm
// out, which is a quarter of the way to the rim.
//
// The reason is mechanical rather than incidental.  `clipped` means the gain
// asked for a pose the mechanism will not hold, and a ball far off centre is a
// large error, so a large error is a command that cannot be held — continuously,
// for the whole excursion.  It is a "the loop is working hard" signal, not a
// fault signal, and sitting second in the precedence it masks `separated` and
// `saturated` permanently.
//
// So the banner's real information content under the loop is nil.  What it
// still discriminates is who was DRIVING — the hand-driven cases above reach
// `separated` and `rolled off` — which the visitor already knows.
//
// **This is a characterisation test, not a specification.**  It pins a
// limitation rather than a requirement, and the day it FAILS is the day
// somebody made the diagnosis discriminate, which is the wanted change.  When
// that happens: do not relax this test, delete it, and correct CONTEXT.md's
// "Losing the ball, and saying so" along with the banner's claim.
//
// See the discussion on
// [#33](https://github.com/caliburn-engineering/caliburn/issues/33).
void test_the_closed_loop_reports_one_constant_cause_today() {
    const auto models = getBuiltinModels();
    const ModelEntry& e = cascadeModel(models);
    const SimPlate plate = cascadePlate(e.params);

    int losses = 0;
    for (const char* preset : {"Detuned", "Nominal", "Aggressive"}) {
        for (double kick : {1.30, 2.00, 5.00}) {
            SimInput in;
            in.design = cascadeDesign(e.params);
            in.design.K = gainForPreset(e, presetNamed(preset));

            const int kDirections = 24;
            for (int d = 0; d < kDirections; ++d) {
                const double theta = 2.0 * M_PI * d / kDirections;
                SimState s = simStart(plate, in.design.home_leg_rad,
                                      Eigen::Vector4d(0.06, -0.04, 0.0, 0.0));
                const int settle = static_cast<int>(2.5 / in.dt);
                const int total = settle + static_cast<int>(12.0 / in.dt);
                for (int k = 0; k < total; ++k) {
                    if (k == settle && !s.ball.airborne) {
                        s.ball.rolling(2) += kick * std::cos(theta);
                        s.ball.rolling(3) += kick * std::sin(theta);
                    }
                    const SimReport f = stepSim(plate, in, s);
                    if (!f.left_plate) continue;

                    ++losses;
                    const LossFlags fl = lossFlags(f, s);
                    ASSERT_TRUE(fl.clipped);
                    ASSERT_TRUE(fl.saturated);
                    ASSERT_TRUE(!fl.rates_untrustworthy);
                    ASSERT_TRUE(!fl.separated);
                    ASSERT_TRUE(lossCause(fl) == LossCause::WorkspaceClipped);
                    break;
                }
            }
        }
    }
    // Every combination lost the ball, which is itself part of the claim: there
    // is no closed-loop loss in this range that reports anything else, because
    // there is no closed-loop loss in this range that got away.
    ASSERT_EQ(losses, 3 * 3 * 24);
}

// --- The hand-driven paths, which are what the banner is now mostly for ---

// The legs dragged apart with the ball already running: the plate drops away
// from underneath it, contact ends, and the ball is in the air when it crosses
// the rim.  `separated` at the crossing frame, and no loop is involved — with
// the loop open `clipped` and `saturated` are false by construction, because
// nobody asked the gain for anything.
//
// This is the only route in the repository that reaches branch 3 through the
// application rather than through the rule, which is why it is a test and not a
// note.  Branch 1 has no such route: `legCommand` and `stepServosOnPlate` both
// retreat to poses whose velocity Jacobian is believable (#29), so a plate the
// application will actually assemble almost never distrusts its own rates.  It
// stays in the precedence on the argument in `lossCause`, pinned by the sweep
// above rather than by a run.
void test_a_ball_thrown_off_by_hand_is_named_separated() {
    const auto models = getBuiltinModels();
    const ModelEntry& e = cascadeModel(models);
    const SimPlate plate = cascadePlate(e.params);

    SimInput in;
    in.design = cascadeDesign(e.params);
    in.closed_loop = false;   // the servo sliders own the legs
    // Two legs to their high stop and one to its low stop — the hardest tilt
    // the three sliders can be dragged to, which is what a visitor exploring
    // the mechanism does within about a second of finding them.
    in.open_loop_cmd_rad = {80.0 * kDeg, 80.0 * kDeg, 10.0 * kDeg};

    // The ball 200 mm out and rolling, so that it reaches the rim while the
    // plate is still dropping rather than after it has settled.
    SimState s = simStart(plate, in.design.home_leg_rad,
                          Eigen::Vector4d(0.20, 0.0, 0.6, 0.0));

    bool lost = false;
    for (int k = 0; k < static_cast<int>(6.0 / in.dt); ++k) {
        const SimReport f = stepSim(plate, in, s);
        if (!f.left_plate) continue;
        lost = true;
        const LossFlags fl = lossFlags(f, s);
        ASSERT_TRUE(fl.separated);
        ASSERT_TRUE(!fl.clipped);     // the loop is not driving
        ASSERT_TRUE(!fl.saturated);   // ...so it is not complaining either
        ASSERT_TRUE(lossCause(fl) == LossCause::Separated);
        break;
    }
    ASSERT_TRUE(lost);
}

// And the same slam from a ball at rest at the centre: it reaches the rim after
// the plate has finished moving, so nothing is raised and the cause is "rolled
// off".
//
// **That is the answer being correct, not the answer being empty.**  A visitor
// who tilts the plate over and watches the ball run off has not found a fault,
// and a banner that blamed the mechanism for it would be lying to them.  The
// case is here because it is the one the sentence is easiest to get wrong.
void test_a_ball_tilted_off_by_hand_is_named_rolled_off() {
    const auto models = getBuiltinModels();
    const ModelEntry& e = cascadeModel(models);
    const SimPlate plate = cascadePlate(e.params);

    SimInput in;
    in.design = cascadeDesign(e.params);
    in.closed_loop = false;
    in.open_loop_cmd_rad = {80.0 * kDeg, 10.0 * kDeg, 10.0 * kDeg};

    SimState s = simStart(plate, in.design.home_leg_rad);

    bool lost = false;
    for (int k = 0; k < static_cast<int>(6.0 / in.dt); ++k) {
        const SimReport f = stepSim(plate, in, s);
        if (!f.left_plate) continue;
        lost = true;
        const LossFlags fl = lossFlags(f, s);
        ASSERT_TRUE(lossCause(fl) == LossCause::RolledOff);
        ASSERT_TRUE(lossFlagsLine(fl) == "none");
        break;
    }
    ASSERT_TRUE(lost);
}

}  // namespace

int main() {
    test_the_precedence_names_the_most_upstream_flag_that_is_set();
    test_each_cause_masks_the_ones_downstream_of_it();
    test_no_flag_at_all_still_names_a_cause();
    test_every_cause_has_a_label_and_a_sentence();
    test_the_raw_flags_are_all_reported_beside_the_cause();
    test_the_flags_are_read_from_the_fields_they_name();
    test_a_loss_under_the_loop_is_named_workspace_clipped();
    test_the_offered_envelope_does_not_lose_the_ball();
    test_the_two_scarce_causes_are_reachable_above_the_offered_envelope();
    test_the_closed_loop_reports_one_constant_cause_today();
    test_a_ball_thrown_off_by_hand_is_named_separated();
    test_a_ball_tilted_off_by_hand_is_named_rolled_off();
    std::printf("test_loss_cause: all passed\n");
    return 0;
}
