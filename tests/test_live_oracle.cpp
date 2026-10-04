// tests/test_live_oracle.cpp
//
// Live-oracle test: runs the oracle pipeline (export_plants →
// gen_oracle_fixtures.py) and cross-checks the freshly generated fixture
// against the checked-in oracle_fixture.h.
//
// Gated by CALIBURN_LIVE_ORACLE cmake option; carries the "live-oracle" ctest
// label.  Exits with SKIP_RETURN_CODE (77) if the oracle (scipy, imported by
// the Python that CALIBURN_ORACLE_PYTHON names) is not available, so the test
// reports "Skipped" rather than "Failed" in ctest's output.
//
// The comparison is line by line.  The K(i, j) and P(i, j) assignments are
// compared as numbers, to the same relative tolerance the default suite holds
// the solver to: a different scipy, LAPACK or CPU moves their last digits, and
// that is not a stale fixture.  Every other line — names, shapes, plant hashes
// — must match exactly.  The Date and Oracle metadata lines are skipped, but
// the fresh Oracle line must name scipy: the generator falls back to a
// pure-Python sign function when scipy will not import, and that fallback
// shares its algorithm with the solver under test, so it is no oracle.
//
// To invoke:
//   cmake -S . -B build -DCALIBURN_LIVE_ORACLE=ON \
//         -DCALIBURN_ORACLE_PYTHON=.oracle_venv/bin/python3
//   cmake --build build -j2
//   ctest --test-dir build -L live-oracle --output-on-failure

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

// Paths injected at compile time by CMake.
#ifndef EXPORT_PLANTS_EXE
#  error "EXPORT_PLANTS_EXE must be defined by CMake"
#endif
#ifndef GEN_ORACLE_SCRIPT
#  error "GEN_ORACLE_SCRIPT must be defined by CMake"
#endif
#ifndef ORACLE_FIXTURE_FILE
#  error "ORACLE_FIXTURE_FILE must be defined by CMake"
#endif
#ifndef ORACLE_PYTHON
#  error "ORACLE_PYTHON must be defined by CMake"
#endif

static constexpr int SKIP_CODE = 77;

// Matches tests/test_oracle_fixtures.cpp: the measured agreement between the
// solver and scipy is ~1e-14, so 1e-10 absorbs any platform's last digits.
static constexpr double kRelTol = 1e-10;

static std::string read_file(const char* path) {
    std::ifstream f(path);
    if (!f) return {};
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

static std::vector<std::string> lines_of(const std::string& text) {
    std::vector<std::string> out;
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("// Date  :", 0) == 0) continue;
        if (line.rfind("// Oracle:", 0) == 0) continue;
        out.push_back(line);
    }
    return out;
}

static std::string oracle_line(const std::string& text) {
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line))
        if (line.rfind("// Oracle:", 0) == 0) return line;
    return {};
}

// Splits "    K(0, 3) = -3.65;" into its "    K(0, 3) = " prefix and value.
// Returns false for any line that is not a K or P element assignment.
static bool split_assignment(const std::string& line, std::string& prefix,
                             double& value) {
    const size_t first = line.find_first_not_of(' ');
    if (first == std::string::npos) return false;
    if (line[first] != 'K' && line[first] != 'P') return false;
    if (line.compare(first + 1, 1, "(") != 0) return false;
    const size_t eq = line.find(" = ");
    if (eq == std::string::npos || line.empty() || line.back() != ';')
        return false;
    prefix = line.substr(0, eq + 3);
    const std::string number = line.substr(eq + 3, line.size() - eq - 4);
    char* end = nullptr;
    value = std::strtod(number.c_str(), &end);
    return end != number.c_str() && *end == '\0';
}

// Returns the number of mismatching lines, reporting each one.
static int compare(const std::vector<std::string>& fresh,
                   const std::vector<std::string>& stored) {
    if (fresh.size() != stored.size()) {
        std::fprintf(stderr,
            "FAIL: fresh fixture has %zu lines, checked-in has %zu\n",
            fresh.size(), stored.size());
        return 1;
    }
    int mismatches = 0;
    for (size_t i = 0; i < fresh.size(); ++i) {
        std::string fp, sp;
        double fv = 0.0, sv = 0.0;
        const bool fa = split_assignment(fresh[i], fp, fv);
        const bool sa = split_assignment(stored[i], sp, sv);
        bool same;
        if (fa && sa) {
            const double scale = std::max({std::abs(fv), std::abs(sv), 1.0});
            same = fp == sp && std::abs(fv - sv) <= kRelTol * scale;
        } else {
            same = fresh[i] == stored[i];
        }
        if (!same) {
            std::fprintf(stderr, "  fresh  : %s\n  stored : %s\n",
                         fresh[i].c_str(), stored[i].c_str());
            ++mismatches;
        }
    }
    return mismatches;
}

int main() {
    // Detect oracle availability at runtime.
    if (std::system("\"" ORACLE_PYTHON "\" -c \"import scipy\" 2>/dev/null") != 0) {
        std::printf("SKIP: scipy does not import in %s; install "
                    "tools/requirements.txt there, or point "
                    "CALIBURN_ORACLE_PYTHON at a Python that has it\n",
                    ORACLE_PYTHON);
        return SKIP_CODE;
    }

    char plants_path[256];
    char fresh_path[256];
    std::snprintf(plants_path, sizeof(plants_path),
                  "/tmp/bb_live_oracle_plants_%d.json", (int)getpid());
    std::snprintf(fresh_path, sizeof(fresh_path),
                  "/tmp/bb_live_oracle_fresh_%d.h", (int)getpid());

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

    // Run the oracle generator to produce a fresh fixture header.
    {
        char cmd[1024];
        std::snprintf(cmd, sizeof(cmd), "\"%s\" \"%s\" \"%s\" > \"%s\"",
                      ORACLE_PYTHON, GEN_ORACLE_SCRIPT, plants_path, fresh_path);
        const int ret = std::system(cmd);
        std::remove(plants_path);
        if (ret != 0) {
            std::fprintf(stderr, "FAIL: gen_oracle_fixtures.py exited non-zero\n");
            return 1;
        }
    }

    const std::string fresh = read_file(fresh_path);
    std::remove(fresh_path);
    const std::string stored = read_file(ORACLE_FIXTURE_FILE);

    if (fresh.empty()) {
        std::fprintf(stderr, "FAIL: oracle generator produced no output\n");
        return 1;
    }
    if (stored.empty()) {
        std::fprintf(stderr,
            "FAIL: cannot read checked-in fixture at %s\n",
            ORACLE_FIXTURE_FILE);
        return 1;
    }

    const std::string oracle = oracle_line(fresh);
    // A prefix check: the fallback's line reads "... (no scipy)".
    if (oracle.rfind("// Oracle: scipy", 0) != 0) {
        std::fprintf(stderr,
            "FAIL: the generator did not use scipy (%s)\n"
            "  scipy imported, so its fallback should not have run.\n",
            oracle.empty() ? "no Oracle line" : oracle.c_str());
        return 1;
    }

    const int mismatches = compare(lines_of(fresh), lines_of(stored));
    if (mismatches > 0) {
        std::fprintf(stderr,
            "FAIL: %d line(s) of the freshly generated fixture differ from %s\n"
            "  (K and P to relative tolerance %.0e; every other line exactly)\n"
            "  Regenerate with:\n"
            "    cmake --build build --target export_plants\n"
            "    ./build/export_plants | %s %s - > %s\n",
            mismatches, ORACLE_FIXTURE_FILE, kRelTol, ORACLE_PYTHON,
            GEN_ORACLE_SCRIPT, ORACLE_FIXTURE_FILE);
        return 1;
    }

    std::printf("Live oracle check passed: fresh fixture matches oracle_fixture.h\n"
                "  (%s)\n", oracle.c_str());
    return 0;
}
