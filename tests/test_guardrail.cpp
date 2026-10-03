// tests/test_guardrail.cpp
//
// End-to-end tests for the design-gating / engagement guardrail:
//   - `nextEngage` transitions (pure function tests)
//   - The invariant that a non-usable design is never handed to the execution
//     layer, exercised through the `stepSim` seam.
//
// Builds against `ball_dynamics` only.  `designToOffer` itself is not called
// here (it takes AppState which pulls imgui headers), but `Offer` is constructed
// directly to simulate what `designToOffer` returns for a failed solve.

#include "attract_mode.h"
#include "auto_balance.h"
#include "cascade_fixture.h"
#include "plate_design.h"
#include "sim_step.h"
#include "test_helpers.h"

#include <cstdio>
#include <cstdlib>

using namespace caliburn;

namespace {

constexpr double kDeg = M_PI / 180.0;
constexpr double kHome = 45.0 * kDeg;
constexpr double kDt = 1.0 / 60.0;

// ---------------------------------------------------------------------------
// nextEngage transitions — pure function tests
// ---------------------------------------------------------------------------

void test_auto_engage_fires_once() {
    // First usable frame: auto-engage fires, latch closes.
    Engage e{};
    e = nextEngage(e, true, true, LossKind::None);
    ASSERT_TRUE(e.engaged);
    ASSERT_TRUE(e.auto_engaged);

    // Visitor deliberately drops and then the design is still usable.
    // The latch is closed, so the loop must not re-fire.
    e.engaged = false;
    Engage e2 = nextEngage(e, true, true, LossKind::Deliberate);
    ASSERT_TRUE(!e2.engaged);
    ASSERT_TRUE(e2.auto_engaged);
}

void test_drop_on_unusable_design() {
    Engage e{true, true};
    Engage next = nextEngage(e, false, true, LossKind::None);
    ASSERT_TRUE(!next.engaged);
    ASSERT_TRUE(next.auto_engaged);  // latch survives a transient drop
}

void test_transient_reengages_when_design_recovers() {
    // Start engaged.
    Engage e{true, true};

    // Design goes away.
    e = nextEngage(e, false, true, LossKind::None);
    ASSERT_TRUE(!e.engaged);

    // Design recovers: Transient loss → re-engage.
    Engage recovered = nextEngage(e, true, true, LossKind::Transient);
    ASSERT_TRUE(recovered.engaged);
}

void test_deliberate_drop_stays_dropped_when_design_recovers() {
    Engage e{true, true};

    // Visitor drops (checkbox uncheck — caller records Deliberate).
    e.engaged = false;

    // Design usable, but loss is Deliberate: must stay dropped.
    Engage after = nextEngage(e, true, true, LossKind::Deliberate);
    ASSERT_TRUE(!after.engaged);
}

void test_no_engage_without_ball() {
    // Auto-engage requires ball_on.
    Engage e{};
    Engage next = nextEngage(e, true, false, LossKind::None);
    ASSERT_TRUE(!next.engaged);

    // Transient re-engage also requires ball_on.
    e.auto_engaged = true;
    Engage t = nextEngage(e, true, false, LossKind::Transient);
    ASSERT_TRUE(!t.engaged);
}

// ---------------------------------------------------------------------------
// End-to-end through stepSim
// ---------------------------------------------------------------------------

// Helper: advance engage state and record loss kind, mirroring what PlateView
// does every frame.
void advanceEngage(Engage& e, LossKind& loss,
                   bool usable, bool ball_on) {
    const Engage next = nextEngage(e, usable, ball_on, loss);
    if (e.engaged && !next.engaged)
        loss = LossKind::Transient;
    else if (next.engaged)
        loss = LossKind::None;
    e = next;
}

void test_closed_loop_keeps_ball_from_opening_position() {
    const auto& models = getBuiltinModels();
    const auto& cascade = cascadeModel(models);
    const auto& tk = cascadeKinematics();

    AutoBalanceDesign d = cascadeDesign(cascade.params);
    d.K = defaultGain(cascade);

    SimPlate plate = cascadePlate(cascade.params);
    SimState sim   = simStart(plate, kHome, attractStart(openingPath()));

    Engage engage{};
    LossKind loss = LossKind::None;

    const bool usable = gainFitsCascade(d) &&
                        samePlant(d.mechanism, d.gravity,
                                  plate.kinematics().params(), plate.gravity());

    for (int i = 0; i < 60; ++i) {
        advanceEngage(engage, loss, usable, true);

        SimInput in;
        in.dt          = kDt;
        in.design      = d;
        in.closed_loop = engage.engaged && usable;
        in.ball_enabled = true;

        SimReport rep = stepSim(plate, in, sim);
        ASSERT_TRUE(!rep.left_plate);
    }

    ASSERT_TRUE(engage.engaged);
}

void test_loop_drops_and_reengages_on_transient_design_loss() {
    const auto& models = getBuiltinModels();
    const auto& cascade = cascadeModel(models);

    AutoBalanceDesign good = cascadeDesign(cascade.params);
    good.K = defaultGain(cascade);

    // bad: empty K, gainFitsCascade returns false → usable = false
    AutoBalanceDesign bad = cascadeDesign(cascade.params);

    SimPlate plate = cascadePlate(cascade.params);
    SimState sim   = simStart(plate, kHome, attractStart(openingPath()));

    Engage engage{};
    LossKind loss = LossKind::None;

    // Phase 1: engage with good design (20 frames).
    for (int i = 0; i < 20; ++i) {
        const bool usable = gainFitsCascade(good) &&
                            samePlant(good.mechanism, good.gravity,
                                      plate.kinematics().params(), plate.gravity());
        advanceEngage(engage, loss, usable, true);
        SimInput in; in.dt = kDt; in.design = good;
        in.closed_loop = engage.engaged && usable;
        in.ball_enabled = true;
        stepSim(plate, in, sim);
    }
    ASSERT_TRUE(engage.engaged);

    // Phase 2: bad design (solve failed) — loop must drop.
    {
        const bool usable = gainFitsCascade(bad);
        advanceEngage(engage, loss, usable, true);
        SimInput in; in.dt = kDt; in.design = bad;
        in.closed_loop = false;
        in.ball_enabled = true;
        stepSim(plate, in, sim);
    }
    ASSERT_TRUE(!engage.engaged);
    ASSERT_EQ(loss, LossKind::Transient);

    // Phase 3: good design returns — loop must re-engage automatically.
    {
        const bool usable = gainFitsCascade(good) &&
                            samePlant(good.mechanism, good.gravity,
                                      plate.kinematics().params(), plate.gravity());
        advanceEngage(engage, loss, usable, true);
    }
    ASSERT_TRUE(engage.engaged);
    ASSERT_EQ(loss, LossKind::None);
}

void test_deliberate_drop_does_not_reengage_when_design_returns() {
    const auto& models = getBuiltinModels();
    const auto& cascade = cascadeModel(models);

    AutoBalanceDesign good = cascadeDesign(cascade.params);
    good.K = defaultGain(cascade);
    AutoBalanceDesign bad = cascadeDesign(cascade.params);  // empty K

    SimPlate plate = cascadePlate(cascade.params);
    SimState sim   = simStart(plate, kHome, attractStart(openingPath()));

    Engage engage{};
    LossKind loss = LossKind::None;

    // Engage.
    for (int i = 0; i < 5; ++i) {
        const bool usable = gainFitsCascade(good) &&
                            samePlant(good.mechanism, good.gravity,
                                      plate.kinematics().params(), plate.gravity());
        advanceEngage(engage, loss, usable, true);
        SimInput in; in.dt = kDt; in.design = good;
        in.closed_loop = engage.engaged && usable;
        in.ball_enabled = true;
        stepSim(plate, in, sim);
    }
    ASSERT_TRUE(engage.engaged);

    // Visitor deliberately drops (simulates checkbox uncheck).
    engage.engaged = false;
    loss = LossKind::Deliberate;

    // Design temporarily fails, then recovers.  Loop must stay dropped.
    {
        const bool usable_bad = gainFitsCascade(bad);
        engage = nextEngage(engage, usable_bad, true, loss);
        // loss stays Deliberate (caller doesn't change it — no transient here)
    }
    ASSERT_TRUE(!engage.engaged);

    {
        const bool usable_good = gainFitsCascade(good) &&
                                 samePlant(good.mechanism, good.gravity,
                                           plate.kinematics().params(),
                                           plate.gravity());
        engage = nextEngage(engage, usable_good, true, loss);
    }
    ASSERT_TRUE(!engage.engaged);  // deliberate: stays dropped
}

// ---------------------------------------------------------------------------
// AC3 end-to-end: non-offered design → engagement off → plate holds home pose
// ---------------------------------------------------------------------------

void test_failed_offer_plate_commands_home_not_stale_tilt() {
    // Asserts the full AC3 chain:
    //   (a) A non-offered design (with a reason, same as designToOffer returns
    //       on a failed LQR solve) drives engagement off.
    //   (b) Through the stepSim seam the plate is commanded to home — not
    //       frozen at whatever tilt the previous gain was asking for.
    //   (c) The ball does not leave the plate.
    //
    // The Offer is constructed directly to stand in for what designToOffer
    // produces when lqr_result.success is false: offered=false, reason set,
    // design.K empty (cascadeDesign without a K).
    const auto& models = getBuiltinModels();
    const auto& cascade = cascadeModel(models);

    AutoBalanceDesign good = cascadeDesign(cascade.params);
    good.K = defaultGain(cascade);

    SimPlate plate = cascadePlate(cascade.params);
    SimState sim   = simStart(plate, kHome, attractStart(openingPath()));

    const bool usable_good =
        gainFitsCascade(good) &&
        samePlant(good.mechanism, good.gravity,
                  plate.kinematics().params(), plate.gravity());

    Engage   engage{};
    LossKind loss = LossKind::None;

    // Phase 1: engage the loop with a good design.
    for (int i = 0; i < 20; ++i) {
        advanceEngage(engage, loss, usable_good, true);
        SimInput in; in.dt = kDt; in.design = good;
        in.closed_loop  = engage.engaged && usable_good;
        in.ball_enabled = true;
        stepSim(plate, in, sim);
    }
    ASSERT_TRUE(engage.engaged);

    // Simulate a failed solver result: offered=false, reason non-empty, K cleared.
    Offer fail_offer;
    fail_offer.offered = false;
    fail_offer.reason  = "the LQR solve failed";
    fail_offer.design  = cascadeDesign(cascade.params);  // K is empty

    ASSERT_TRUE(!fail_offer.offered);
    ASSERT_TRUE(!fail_offer.reason.empty());

    // Phase 2: failed design — loop must drop, plate must command home.
    for (int i = 0; i < 5; ++i) {
        const bool usable = gainFitsCascade(fail_offer.design);  // false: K empty
        advanceEngage(engage, loss, usable, true);

        SimInput in; in.dt = kDt; in.design = fail_offer.design;
        in.closed_loop  = engage.engaged && usable;  // false
        in.ball_enabled = true;
        const SimReport rep = stepSim(plate, in, sim);

        ASSERT_TRUE(!engage.engaged);
        // open_loop_cmd_rad defaults to home; closed_loop is false.
        for (int j = 0; j < 3; ++j)
            ASSERT_NEAR(rep.cmd_rad[j], kHome, 1e-9);
        ASSERT_TRUE(!rep.left_plate);
    }
}

}  // namespace

int main() {
    test_auto_engage_fires_once();
    test_drop_on_unusable_design();
    test_transient_reengages_when_design_recovers();
    test_deliberate_drop_stays_dropped_when_design_recovers();
    test_no_engage_without_ball();
    test_closed_loop_keeps_ball_from_opening_position();
    test_loop_drops_and_reengages_on_transient_design_loss();
    test_deliberate_drop_does_not_reengage_when_design_returns();
    test_failed_offer_plate_commands_home_not_stale_tilt();

    std::printf("test_guardrail: all tests passed\n");
    return 0;
}
