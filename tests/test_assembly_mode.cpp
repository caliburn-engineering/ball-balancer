// tests/test_assembly_mode.cpp
//
// The 3-RRS constraint equations have TWO roots, and the residual norm cannot
// tell them apart.  Beside the assembly the machine is built in — table above
// its knees — there is a folded one with the table lying flat on the base, and
// for this mechanism's parameters that folded root is EXACT: `R_ground ==
// R_table` and `L1 == L2` make `z_c = 0` satisfy every leg's length
// constraint to machine precision, at every servo angle.
//
// Newton does not care which it lands on.  Seeded from a pose that a hard
// manoeuvre has carried down past the knee plane, it converges on the folded
// root in one iteration, reports success with a residual of 1e-11, and stays
// there for the rest of the session — the next frame is seeded from a pose
// that is already a root.  The plate then has no tilt authority at all, and
// every ball put on it rolls off.
//
// That is issue #22, and these are the tests that keep it fixed.
#include "table_kinematics.h"

#include "auto_balance.h"
#include "cascade_fixture.h"
#include "test_helpers.h"

#include <algorithm>
#include <cmath>
#include <random>

using namespace caliburn;

namespace {

constexpr double kDeg = M_PI / 180.0;

TableParams cascadeTable() {
    return cascadeMechanism(cascadeModel(getBuiltinModels()).params);
}

std::array<double, 3> all(double deg) {
    return {deg * kDeg, deg * kDeg, deg * kDeg};
}

// The defect itself, pinned as a property of the mechanism rather than as a
// memory.  If a future parameter change stops `z_c = 0` being a root, this
// test fails and the guard below can be reconsidered — which is the point.
void test_the_folded_pose_is_an_exact_root() {
    const TableKinematics tk(cascadeTable());
    for (double deg : {10.0, 30.0, 45.0, 60.0, 80.0}) {
        const Eigen::Vector3d f = tk.fk_residual(all(deg), TablePose{0.0, 0.0, 0.0});
        ASSERT_TRUE(f.norm() < 1e-12);
    }
}

// And it is a DIFFERENT root, not a rounding of the real one — nor a family of
// them.  The folded assembly is exactly one pose, the same at every servo
// angle, which is what lets the floor below be placed so generously.
void test_the_folded_assembly_is_one_pose_at_every_angle() {
    const TableKinematics tk(cascadeTable());
    for (double deg : {10.0, 30.0, 45.0, 60.0, 80.0}) {
        const FKResult folded = tk.forward_kinematics(all(deg), TablePose{0, 0, 0});
        ASSERT_TRUE(folded.converged);
        ASSERT_NEAR(folded.pose.z_c, 0.0, 1e-12);
        ASSERT_NEAR(folded.pose.phi, 0.0, 1e-12);
        ASSERT_NEAR(folded.pose.theta, 0.0, 1e-12);
        // While the built assembly moves with the servos, and is nowhere near.
        ASSERT_TRUE(tk.home_pose(deg * kDeg).z_c > 19.0 * tk.assembly_floor(all(deg)));
    }
}

// The separator, measured.  Over the servo box the built assembly clears the
// floor by 1.63x at worst and about 19.8x typically — 0.163 and 1.98 mean knee
// heights against a floor at 0.1 of one — while the folded root has to climb
// from zero to reach it at all.  The thin low tail is the point of the
// assertion: it is what stops the floor being raised on a hunch.
void test_the_floor_separates_the_two_assemblies() {
    const TableParams tp = cascadeTable();
    const TableKinematics tk(tp);
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> U(tp.alpha_min, tp.alpha_max);

    int built = 0, low_tail = 0;
    double worst_ratio = 1e9;
    for (int k = 0; k < 60000; ++k) {
        const std::array<double, 3> a = {U(rng), U(rng), U(rng)};
        const FKResult fk = tk.solve_pose(a, tk.home_pose((a[0]+a[1]+a[2]) / 3.0));
        if (!fk.converged) continue;   // no assembly exists for this triple
        ++built;
        const double ratio = fk.pose.z_c / tk.assembly_floor(a);
        worst_ratio = std::min(worst_ratio, ratio);
        if (ratio < 5.0) ++low_tail;   // below half a mean knee height
    }
    ASSERT_TRUE(built > 20000);

    // Above: nothing accepted comes near the floor.  The measured minimum sits
    // at 1.63 times it — 0.163 mean knee heights against a floor of 0.1.
    ASSERT_TRUE(worst_ratio > 1.3);

    // Below: the low tail is real, and small.  Both halves matter — if the
    // tail ever thickened, a floor chosen against a thin one would start
    // discarding poses the mechanism can genuinely reach.
    ASSERT_TRUE(low_tail > 0);
    ASSERT_TRUE(low_tail < built / 1000);
}

// The bug, reproduced and refused: seeded from the folded pose, `solve_pose`
// must climb back to the assembly the machine is built in.  Plain
// `forward_kinematics` does not, and cannot — it is a solver, not a mechanism.
void test_solve_pose_recovers_from_a_folded_seed() {
    const TableKinematics tk(cascadeTable());
    const std::array<double, 3> a = all(45.0);
    const TablePose folded{0.0, 0.0, 0.0};

    const FKResult naive = tk.forward_kinematics(a, folded);
    ASSERT_TRUE(naive.converged);              // it "succeeds" —
    ASSERT_NEAR(naive.pose.z_c, 0.0, 1e-9);    // — on the wrong root

    const FKResult fixed = tk.solve_pose(a, folded);
    ASSERT_TRUE(fixed.converged);
    ASSERT_NEAR(fixed.pose.z_c, tk.home_pose(45.0 * kDeg).z_c, 1e-9);
}

// Once folded, the old code stayed folded forever: every frame was seeded from
// a pose that was already a root, so Newton returned it in one iteration.
// Sixty seconds of that is what a visitor was seeing.
void test_a_folded_pose_does_not_persist_across_frames() {
    const TableKinematics tk(cascadeTable());
    TablePose seed{0.0, 0.0, 0.0};
    for (int frame = 0; frame < 600; ++frame) {
        const FKResult fk = tk.solve_pose(all(45.0), seed);
        ASSERT_TRUE(fk.converged);
        ASSERT_TRUE(fk.pose.z_c > tk.assembly_floor(all(45.0)));
        seed = fk.pose;
    }
}

// A leg triple with no assembly at all is reported as a failure, not as a
// folded plate.  The servo travel limits are a box; the workspace is not, so
// a per-leg clamp can ask for a configuration that does not exist.  The
// caller's answer is to keep the pose it had — a mechanism driven into a
// singularity binds and stops, it does not lie flat.
void test_an_unreachable_triple_fails_rather_than_folding() {
    const TableParams tp = cascadeTable();
    const TableKinematics tk(tp);
    // Found by sweeping: the servo box is roughly half reachable, and the
    // extremes of it are not.
    const std::array<double, 3> wide = {80.0 * kDeg, 10.0 * kDeg, 80.0 * kDeg};
    const FKResult fk = tk.solve_pose(wide, tk.home_pose(45.0 * kDeg));
    ASSERT_TRUE(!fk.converged);
    // And it did not quietly hand back the folded root as consolation.
    ASSERT_TRUE(!(fk.pose.z_c > 0.0 && fk.pose.z_c <= tk.assembly_floor(wide)));
}

// The second half of #22, at the seam nothing else covers.
//
// `legCommand`'s clip is tested in `test_auto_balance`, where the command
// lives.  This is the other one: `stepServos` lags each leg independently, so
// its result sits on the straight line from where the legs are to where they
// were told to go — and that line can leave the workspace even when both of
// its ends are inside it.  Clipping the command alone left two kick
// directions in every 360 with one frame, mid-flight, that had no assembly.
void test_the_servo_path_stays_inside_the_workspace() {
    const TableParams tp = cascadeTable();
    const TableKinematics tk(tp);

    // A command that IS assemblable, from legs that are, but far enough away
    // that the straight line between them leaves the workspace.
    const std::array<double, 3> from = all(45.0);
    const std::array<double, 3> to = {62.0 * kDeg, 32.0 * kDeg, 40.0 * kDeg};
    ASSERT_TRUE(tk.can_assemble(from));
    ASSERT_TRUE(tk.can_assemble(to));

    // Driven all the way there in steps, the plate always has a pose.
    std::array<double, 3> a = from;
    for (int k = 0; k < 600; ++k) {
        a = stepServosOnPlate(tk, a, to, 0.05, 1.0 / 60.0);
        ASSERT_TRUE(tk.can_assemble(a));
    }
}

// And the retreat reports a broken precondition rather than freezing quietly.
void test_a_retreat_from_an_unassemblable_safe_end_says_so() {
    const TableParams tp = cascadeTable();
    const TableKinematics tk(tp);

    const std::array<double, 3> nowhere = {80.0 * kDeg, 10.0 * kDeg, 80.0 * kDeg};
    ASSERT_TRUE(!tk.can_assemble(nowhere));

    const Retreat r = retreatToWorkspace(tk, nowhere, nowhere);
    ASSERT_TRUE(r.retreated);
    ASSERT_TRUE(r.safe_was_unassemblable);

    // Where the precondition holds, the flag stays down.
    const Retreat ok = retreatToWorkspace(tk, nowhere, all(45.0));
    ASSERT_TRUE(ok.retreated);
    ASSERT_TRUE(!ok.safe_was_unassemblable);
    ASSERT_TRUE(tk.can_assemble(ok.alpha_rad));
}

}  // namespace

int main() {
    test_the_folded_pose_is_an_exact_root();
    test_the_folded_assembly_is_one_pose_at_every_angle();
    test_the_floor_separates_the_two_assemblies();
    test_solve_pose_recovers_from_a_folded_seed();
    test_a_folded_pose_does_not_persist_across_frames();
    test_an_unreachable_triple_fails_rather_than_folding();
    test_the_servo_path_stays_inside_the_workspace();
    test_a_retreat_from_an_unassemblable_safe_end_says_so();
    std::printf("test_assembly_mode: all passed\n");
    return 0;
}
