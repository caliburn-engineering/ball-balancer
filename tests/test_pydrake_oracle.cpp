// tests/test_pydrake_oracle.cpp
//
// Second live-oracle test: cross-checks oracle_fixture.h using pydrake's
// ContinuousAlgebraicRiccatiEquation (Schur-based, LAPACK) as an independent
// solver, giving three-way corroboration: the C++ sign-function solver,
// scipy's Schur solve, and pydrake's Schur solve.
//
// Gated by CALIBURN_LIVE_ORACLE cmake option; carries the "live-oracle" ctest
// label.  Exits with SKIP_RETURN_CODE (77) if pydrake is not importable, so
// the test reports "Skipped" rather than "Failed" in ctest's output.
//
// The comparison is tests/oracle_compare.h: code lines only, K and P as
// numbers to the default suite's relative tolerance, everything else exactly.
// The fresh Oracle line must name pydrake, so nothing but pydrake's solve can
// pass this test.
//
// The checked-in fixture is scipy's.  pydrake only checks it: a mismatch here
// with test_live_oracle passing means the two oracles disagree, which is a
// finding to investigate, not a fixture to overwrite with pydrake's output.
//
// To invoke:
//   cmake -S . -B build -DCALIBURN_LIVE_ORACLE=ON \
//         -DCALIBURN_ORACLE_PYTHON=.oracle_venv/bin/python3
//   cmake --build build -j2
//   ctest --test-dir build -L live-oracle --output-on-failure

#include "oracle_compare.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <unistd.h>
#include <vector>

// Paths injected at compile time by CMake.
#ifndef EXPORT_PLANTS_EXE
#  error "EXPORT_PLANTS_EXE must be defined by CMake"
#endif
#ifndef PYDRAKE_ORACLE_SCRIPT
#  error "PYDRAKE_ORACLE_SCRIPT must be defined by CMake"
#endif
#ifndef ORACLE_FIXTURE_FILE
#  error "ORACLE_FIXTURE_FILE must be defined by CMake"
#endif
#ifndef ORACLE_PYTHON
#  error "ORACLE_PYTHON must be defined by CMake"
#endif

static constexpr int SKIP_CODE = 77;

// Same relative tolerance as test_live_oracle and test_oracle_fixtures.
using namespace oracle_compare;

int main() {
    // Detect pydrake availability at runtime.
    if (std::system("\"" ORACLE_PYTHON "\" -c \"import pydrake\" 2>/dev/null") != 0) {
        std::printf("SKIP: pydrake does not import in %s; install pydrake there, or "
                    "point CALIBURN_ORACLE_PYTHON at a Python that has it\n",
                    ORACLE_PYTHON);
        return SKIP_CODE;
    }

    char plants_path[256];
    char fresh_path[256];
    std::snprintf(plants_path, sizeof(plants_path),
                  "/tmp/bb_pydrake_oracle_plants_%d.json", (int)getpid());
    std::snprintf(fresh_path, sizeof(fresh_path),
                  "/tmp/bb_pydrake_oracle_fresh_%d.h", (int)getpid());

    // Run export_plants to produce the plant JSON.
    {
        char cmd[1024];
        std::snprintf(cmd, sizeof(cmd), "\"%s\" > \"%s\"",
                      EXPORT_PLANTS_EXE, plants_path);
        if (std::system(cmd) != 0) {
            std::fprintf(stderr, "FAIL: export_plants exited non-zero\n");
            return 1;
        }
    }

    // Run the pydrake oracle generator to produce a fresh fixture header.
    {
        char cmd[1024];
        std::snprintf(cmd, sizeof(cmd), "\"%s\" \"%s\" \"%s\" > \"%s\"",
                      ORACLE_PYTHON, PYDRAKE_ORACLE_SCRIPT, plants_path, fresh_path);
        const int ret = std::system(cmd);
        std::remove(plants_path);
        if (ret != 0) {
            std::fprintf(stderr, "FAIL: gen_oracle_pydrake.py exited non-zero\n");
            return 1;
        }
    }

    const std::string fresh = read_file(fresh_path);
    std::remove(fresh_path);
    const std::string stored = read_file(ORACLE_FIXTURE_FILE);

    if (fresh.empty()) {
        std::fprintf(stderr, "FAIL: pydrake oracle generator produced no output\n");
        return 1;
    }
    if (stored.empty()) {
        std::fprintf(stderr,
            "FAIL: cannot read checked-in fixture at %s\n",
            ORACLE_FIXTURE_FILE);
        return 1;
    }

    const std::string oracle = oracle_line(fresh);
    if (oracle.rfind("// Oracle: pydrake", 0) != 0) {
        std::fprintf(stderr,
            "FAIL: the generator did not use pydrake (%s)\n"
            "  pydrake imported, so its path should have been taken.\n",
            oracle.empty() ? "no Oracle line" : oracle.c_str());
        return 1;
    }

    const int mismatches = compare(lines_of(fresh), lines_of(stored));
    if (mismatches > 0) {
        std::fprintf(stderr,
            "FAIL: %d line(s) of the pydrake-generated fixture differ from %s\n"
            "  (code lines only; K and P to relative tolerance %.0e)\n"
            "  The fixture is scipy's.  If test_live_oracle passes, the two\n"
            "  oracles disagree; do not regenerate the fixture with pydrake.\n",
            mismatches, ORACLE_FIXTURE_FILE, kRelTol);
        return 1;
    }

    std::printf("pydrake oracle check passed: pydrake agrees with oracle_fixture.h\n"
                "  (%s)\n", oracle.c_str());
    return 0;
}
