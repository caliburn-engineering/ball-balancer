// tests/cascade_fixture.h
#pragma once

// The three facts every closed-loop test on the Ball-Balancer Cascade needs:
// which preset it is, what gain the application opens on, and how big the ball
// is.  They live here rather than in each test file because a second
// transcription is how two tests come to disagree about the plant they are
// both claiming to measure — the same reason `cascadeMechanism` and friends
// are the single source for the application itself.
//
// **The simulation loop is no longer here either, and it used not to be
// anywhere.**  This comment used to argue that the loops were deliberately not
// shared — that `test_auto_balance` measures settling against a setpoint while
// `test_attract_mode` records a trace, and one harness serving both would take
// a parameter per difference.  That is still true of the MEASURING, which is
// why each file still has its own runner.  It was never true of the causal
// order, and the five copies of that drifted apart exactly as this file's own
// argument about the gain predicts they would: see `sim_step.h` and #30.
//
// So a runner here decides what to measure and when to shove the ball, and
// `stepSim` decides what a frame IS.

#include "analysis/lqr.h"
#include "analysis/model_library.h"
#include "auto_balance.h"
#include "ball_sim.h"
#include "rolling_dynamics.h"
#include "table_kinematics.h"
#include "test_helpers.h"

#include <Eigen/Core>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace caliburn {

/// The simulator's ball, not a copy of it: the friction dead band that
/// `test_auto_balance` asserts about is a claim on this exact rolling
/// resistance.
constexpr double kFixtureBallRadius = kPlateBall.radius;
constexpr double kFixtureBallFriction = kPlateBall.rolling_friction;

inline const ModelEntry& cascadeModel(const std::vector<ModelEntry>& models) {
    for (const auto& m : models)
        if (isCascadeModel(m)) return m;
    std::fprintf(stderr, "FAIL: no cascade model in the library\n");
    std::exit(1);
}

/// The plate itself.  `legCommand` needs it because the servo travel limits
/// are a box and the workspace is not — see issue #22 — so every test that
/// evaluates the loop needs a mechanism to check the command against.
inline const TableKinematics& cascadeKinematics() {
    static const TableKinematics tk(
        cascadeMechanism(cascadeModel(getBuiltinModels()).params));
    return tk;
}

inline Eigen::MatrixXd gainFor(const ModelEntry& e,
                               const Eigen::VectorXd& q,
                               const Eigen::VectorXd& r) {
    const LqrResult res = computeLQR(e.system,
                                     q.asDiagonal().toDenseMatrix(),
                                     r.asDiagonal().toDenseMatrix());
    ASSERT_TRUE(res.success);
    return res.K;
}

/// The preset named `name`, or a failed test.  Shared for the same reason
/// everything else in this file is: two test files that each look a preset up
/// their own way are two test files that can come to disagree about which
/// tuning they are measuring.  `test_attract_mode` used to reach for
/// `lqrPresets().back()`, which silently tracks whatever is listed last.
inline const LqrPreset& presetNamed(const char* name) {
    for (const auto& p : lqrPresets())
        if (std::string(p.name) == name) return p;
    std::fprintf(stderr, "FAIL: no preset named %s\n", name);
    std::exit(1);
}

/// A preset's gain against this plant.  The two weight vectors always travel
/// together — a preset IS the pair — so they are handed over as one thing.
inline Eigen::MatrixXd gainForPreset(const ModelEntry& e, const LqrPreset& p) {
    return gainFor(e, presetStateWeights(p, 7), presetInputWeights(p, 3));
}

/// The gain the application actually opens on — the weights from
/// `defaultLqrStateWeights` / `defaultLqrInputWeights`, solved against this
/// preset.  Every "does the demo work" claim is a claim about this matrix.
inline Eigen::MatrixXd defaultGain(const ModelEntry& e) {
    return gainFor(e, defaultLqrStateWeights(7), defaultLqrInputWeights(3));
}

}  // namespace caliburn
