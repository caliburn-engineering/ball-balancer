// tests/test_live_oracle.cpp
//
// Live-oracle test: runs the oracle pipeline (export_plants →
// gen_oracle_fixtures.py) and cross-checks the freshly generated fixture
// against the checked-in oracle_fixture.h.
//
// Gated by CALIBURN_LIVE_ORACLE cmake option; carries the "live-oracle" ctest
// label.  Exits with SKIP_RETURN_CODE (77) if the oracle (python3 + scipy) is
// not available on this machine, so the test reports "Skipped" rather than
// "Failed" in ctest's output.
//
// The two metadata lines that vary across regenerations (date and scipy
// version) are stripped before the text comparison, so a fixture that was
// generated on a different day or with a marginally different scipy still
// passes as long as the numerical constants are identical.
//
// To invoke:
//   cmake -S . -B build -DGLFW_BUILD_WAYLAND=OFF -DCALIBURN_LIVE_ORACLE=ON
//   cmake --build build -j
//   ctest --test-dir build -L live-oracle --output-on-failure

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>

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

static constexpr int SKIP_CODE = 77;

static std::string read_file(const char* path) {
    std::ifstream f(path);
    if (!f) return {};
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

// Strip the two metadata lines that differ between regenerations.
static std::string normalize(const std::string& text) {
    std::istringstream in(text);
    std::ostringstream out;
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("// Date  :", 0) == 0) continue;
        if (line.rfind("// Oracle:", 0) == 0) continue;
        out << line << '\n';
    }
    return out.str();
}

int main() {
    // Detect oracle availability at runtime.
    if (std::system("python3 -c \"import scipy\" 2>/dev/null") != 0) {
        std::printf("SKIP: python3 + scipy unavailable; "
                    "install scipy>=1.11 to run live-oracle tests\n");
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
        std::snprintf(cmd, sizeof(cmd), "python3 \"%s\" \"%s\" > \"%s\"",
                      GEN_ORACLE_SCRIPT, plants_path, fresh_path);
        const int ret = std::system(cmd);
        std::remove(plants_path);
        if (ret != 0) {
            std::fprintf(stderr, "FAIL: gen_oracle_fixtures.py exited non-zero\n");
            return 1;
        }
    }

    const std::string fresh  = normalize(read_file(fresh_path));
    std::remove(fresh_path);

    const std::string stored = normalize(read_file(ORACLE_FIXTURE_FILE));

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

    if (fresh != stored) {
        std::fprintf(stderr,
            "FAIL: freshly-generated fixture differs from %s\n"
            "  Regenerate with:\n"
            "    cmake --build build --target export_plants\n"
            "    ./build/export_plants | python3 %s - > %s\n",
            ORACLE_FIXTURE_FILE, GEN_ORACLE_SCRIPT, ORACLE_FIXTURE_FILE);
        return 1;
    }

    std::printf("Live oracle check passed: fresh fixture matches oracle_fixture.h\n");
    return 0;
}
