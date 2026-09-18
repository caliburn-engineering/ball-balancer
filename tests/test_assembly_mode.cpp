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

    // Driven all the way there in steps, the plate always has a pose — and,
    // since #29, always one whose Jacobian it is allowed to believe.
    std::array<double, 3> a = from;
    for (int k = 0; k < 600; ++k) {
        a = stepServosOnPlate(tk, a, to, 0.05, 1.0 / 60.0);
        ASSERT_TRUE(tk.can_assemble(a));
        ASSERT_TRUE(tk.can_hold(a, kRatesUntrustworthyAbove));
    }
}

// The defect #29 turned out to be, pinned as a property of the mechanism.
//
// The folded root is not the only other assembly.  Approaching a
// direct-kinematics singularity the built assembly the plate is in meets a
// SECOND built one and the two swap, and `assembly_floor` cannot tell them
// apart because neither is anywhere near the floor: at the triple below they
// stand at 153 mm and 75 mm of heave against a floor of 8.5 mm — the nearer of
// them nine times clear of it — and both satisfy the constraint equations to
// better than 1e-9 (measured 1e-12 and 5e-16).
//
// This is the triple the aggressive tuning drove the plate to in #29, reached
// by the old retreat, which asked only whether an assembly EXISTED.  It does —
// two of them.
void test_a_second_built_assembly_waits_past_the_singularity() {
    const TableKinematics tk(cascadeTable());
    const std::array<double, 3> a = {22.93 * kDeg, 22.93 * kDeg, 67.07 * kDeg};

    const FKResult high = tk.solve_pose(a, tk.level_pose(a));
    const FKResult low = tk.solve_pose(
        a, TablePose{11.95 * kDeg, -6.89 * kDeg, 0.0755});

    ASSERT_TRUE(high.converged);
    ASSERT_TRUE(low.converged);
    // Measured at 1e-12 and 5e-16; the bar is 1e-9, loose enough that the
    // claim is "both are real roots" rather than a claim about the solver.
    ASSERT_TRUE(tk.fk_residual(a, high.pose).norm() < 1e-9);
    ASSERT_TRUE(tk.fk_residual(a, low.pose).norm() < 1e-9);

    // Two distinct poses, and the plate is tilted the OTHER WAY on the second:
    // the roll angles lean to opposite sides.  That is what the loop was
    // fighting — it asked for a tilt and got its mirror.
    ASSERT_TRUE(std::abs(high.pose.z_c - low.pose.z_c) > 0.05);
    ASSERT_TRUE(high.pose.phi * low.pose.phi < 0.0);

    // Neither is the folded root, so the guard that catches THAT sees nothing.
    ASSERT_TRUE(high.pose.z_c > tk.assembly_floor(a));
    ASSERT_TRUE(low.pose.z_c > tk.assembly_floor(a));

    // What does see it is the conditioning: both sit far past the line this
    // repository already draws at 20, which is why that line is what the
    // retreat is bounded by.
    ASSERT_TRUE(tk.condition_number(a, high.pose) > 100.0);
    ASSERT_TRUE(tk.condition_number(a, low.pose) > kRatesUntrustworthyAbove);
    ASSERT_TRUE(!tk.can_hold(a, kRatesUntrustworthyAbove));
    ASSERT_TRUE(tk.can_assemble(a));     // and `can_assemble` still says yes
}

// So the retreat stops short of it.  The same command that used to land on
// that triple now lands somewhere the plate can be believed, and the bound
// that stops it is `kRatesUntrustworthyAbove` rather than a number chosen to
// make this case go away.
void test_the_retreat_stops_at_the_conditioned_boundary() {
    const TableKinematics tk(cascadeTable());
    const std::array<double, 3> level = all(45.0);
    // Legs 1 and 2 at their floor and leg 3 at its ceiling: the hardest thing
    // the servo box can ask this plate for, and what the aggressive tuning's
    // command clamps to.
    const std::array<double, 3> corner = {10.0 * kDeg, 10.0 * kDeg, 80.0 * kDeg};
    const double span = 45.0 * kDeg - 10.0 * kDeg;

    const Retreat r = retreatToHoldable(tk, corner, level,
                                        kRatesUntrustworthyAbove);
    ASSERT_TRUE(r.retreated);
    ASSERT_TRUE(!r.safe_was_unholdable);
    ASSERT_TRUE(tk.can_hold(r.alpha_rad, kRatesUntrustworthyAbove));

    // It is a BOUND and not a refusal: the plate is still allowed nearly all of
    // that ray.  Measured, the retreat gives up at 0.6124 of it where an
    // unbounded one runs to 0.6307 — three per cent, and the singularity is
    // in it.
    const double scale = (45.0 * kDeg - r.alpha_rad[0]) / span;
    ASSERT_TRUE(scale > 0.55);
    ASSERT_TRUE(scale < 0.625);

    // And the bound is live: ask for a tighter one and the retreat gives up
    // sooner, ask for a looser one and it goes further.  A limit that moved
    // nothing would be a limit that was not being consulted.
    // Measured: 0.5210 at a limit of 10, 0.6272 at 40.
    const double tight = (45.0 * kDeg -
        retreatToHoldable(tk, corner, level, 10.0).alpha_rad[0]) / span;
    const double loose = (45.0 * kDeg -
        retreatToHoldable(tk, corner, level, 40.0).alpha_rad[0]) / span;
    ASSERT_TRUE(tight < scale);
    ASSERT_TRUE(loose > scale);
}

// And the retreat reports a broken precondition rather than freezing quietly.
void test_a_retreat_from_an_unassemblable_safe_end_says_so() {
    const TableParams tp = cascadeTable();
    const TableKinematics tk(tp);

    const std::array<double, 3> nowhere = {80.0 * kDeg, 10.0 * kDeg, 80.0 * kDeg};
    ASSERT_TRUE(!tk.can_assemble(nowhere));

    const Retreat r = retreatToHoldable(tk, nowhere, nowhere,
                                        kRatesUntrustworthyAbove);
    ASSERT_TRUE(r.retreated);
    ASSERT_TRUE(r.safe_was_unholdable);

    // Where the precondition holds, the flag stays down.
    const Retreat ok = retreatToHoldable(tk, nowhere, all(45.0),
                                         kRatesUntrustworthyAbove);
    ASSERT_TRUE(ok.retreated);
    ASSERT_TRUE(!ok.safe_was_unholdable);
    // What the retreat promises since #29 is holdable, which is the stronger
    // of the two — asserted as the promise rather than as its consequence.
    ASSERT_TRUE(tk.can_hold(ok.alpha_rad, kRatesUntrustworthyAbove));
}

// The home pose is the induction base case for both retreats: `legCommand`
// retreats toward it, and the legs start there.  Its condition number is the
// one number in that argument, so it is pinned rather than asserted in a
// comment.  Measured 4.7 — a quarter of the line at which the plate stops
// being believed, and inside the panel's own "Good" band.
void test_the_home_pose_is_somewhere_the_plate_can_be_held() {
    const TableKinematics tk(cascadeTable());
    const std::array<double, 3> home = all(45.0);
    const FKResult fk = tk.solve_pose(home, tk.level_pose(home));
    ASSERT_TRUE(fk.converged);
    ASSERT_NEAR(tk.condition_number(home, fk.pose), 4.7, 0.05);
    ASSERT_TRUE(tk.can_hold(home, kRatesUntrustworthyAbove));
}

// The march's own property, and the reason it replaced a plain bisection.
//
// Holdability is a reachability set intersected with a condition sublevel set,
// and the condition number rises toward a singularity and falls again past it
// — so a ray out of the home pose can run good, bad, good.  A bisection can
// converge into that third band and hand back a command the legs cannot travel
// to, because the servo step is clipped at the band.  Measured over 4000 random
// targets in the servo box, a plain bisection did exactly that on 1 of the 2147
// retreats it performed: rare, and the kind of rare that is a bug rather than a
// tolerance.
//
// The target below is that one.  Bisecting the whole ray lands on 0.9062 with
// four of two hundred samples behind it unholdable; marching first lands on
// 0.8837 with none.  So the claim is: whatever the retreat returns, every point
// between `safe` and it is holdable too — sampled far more finely than the
// march itself walks.
void test_the_retreat_leaves_no_unholdable_gap_behind_it() {
    const TableKinematics tk(cascadeTable());
    const std::array<double, 3> level = all(45.0);
    const std::array<double, 3> target = {13.7853 * kDeg, 35.0876 * kDeg,
                                          22.8554 * kDeg};

    const Retreat r = retreatToHoldable(tk, target, level,
                                        kRatesUntrustworthyAbove);
    ASSERT_TRUE(r.retreated);
    ASSERT_TRUE(tk.can_hold(r.alpha_rad, kRatesUntrustworthyAbove));

    for (int i = 0; i <= 400; ++i) {
        const double s = i / 400.0;
        std::array<double, 3> a{};
        for (int leg = 0; leg < 3; ++leg)
            a[leg] = level[leg] + s * (r.alpha_rad[leg] - level[leg]);
        ASSERT_TRUE(tk.can_hold(a, kRatesUntrustworthyAbove));
    }
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
    test_a_second_built_assembly_waits_past_the_singularity();
    test_the_retreat_stops_at_the_conditioned_boundary();
    test_a_retreat_from_an_unassemblable_safe_end_says_so();
    test_the_home_pose_is_somewhere_the_plate_can_be_held();
    test_the_retreat_leaves_no_unholdable_gap_behind_it();
    std::printf("test_assembly_mode: all passed\n");
    return 0;
}
