// tools/export_plants.cpp
//
// Walks the preset plants and emits each plant's A, B, Q, R in JSON so the
// Python oracle generator can solve CARE independently and produce fixture
// constants for the test suite.
//
// Usage:
//   ./export_plants > plants.json
//   python3 tools/gen_oracle_fixtures.py plants.json > tests/oracle_fixture.h
//
// The hash field in each record is a DJB2-64 digest of the concatenated
// matrix strings (full-precision %.17g rows).  A changed plant or weight set
// produces a different hash; comparing the header's stored hash against a
// fresh export is how downstream tools detect staleness without re-solving.
//
// Only the Ball-Balancer Cascade preset is exported; the other library entries
// are demonstration plants the product does not design LQR gains against, and
// exporting them would just add noise to the fixture.

#include "analysis/model_library.h"
#include "auto_balance.h"

#include <Eigen/Core>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using namespace caliburn;

namespace {

// Full-precision MATLAB-style matrix string (%.17g per element).
// The existing matrixToString uses %g (6 sig figs), which is not enough
// for an oracle: tiny rounding differences in the Jacobian coefficients
// propagate into the gain and make the fixture fail at tight tolerances.
std::string matrixFull(const Eigen::MatrixXd& m) {
    std::string s;
    char buf[32];
    for (int r = 0; r < m.rows(); ++r) {
        if (r > 0) s += "; ";
        for (int c = 0; c < m.cols(); ++c) {
            if (c > 0) s += " ";
            std::snprintf(buf, sizeof(buf), "%.17g", m(r, c));
            s += buf;
        }
    }
    return s;
}

// DJB2-64: a simple, portable hash for staleness detection.
// The Python generator and any downstream tool can reproduce it from the
// same string with no external library.
uint64_t djb2_64(const std::string& s) {
    uint64_t h = 5381;
    for (unsigned char c : s) h = h * 33u + c;
    return h;
}

// Escape a string for JSON output (handles the few characters that appear in
// matrix strings: digits, spaces, semicolons, minus signs, 'e', 'E', '.').
// No general-purpose escaping is needed for the values we emit.
void printJsonString(const std::string& s) {
    std::printf("\"");
    for (char c : s) {
        if (c == '"')       std::printf("\\\"");
        else if (c == '\\') std::printf("\\\\");
        else                std::putchar(c);
    }
    std::printf("\"");
}

struct PlantRecord {
    std::string name;
    uint64_t    hash;
    std::string A, B, Q, R;
    int n_states, n_inputs;
};

PlantRecord makeRecord(const std::string& name,
                       const LinearSystem& sys,
                       const Eigen::MatrixXd& Q,
                       const Eigen::MatrixXd& R) {
    PlantRecord rec;
    rec.name     = name;
    rec.n_states = sys.states();
    rec.n_inputs = sys.inputs();
    rec.A = matrixFull(sys.A);
    rec.B = matrixFull(sys.B);
    rec.Q = matrixFull(Q);
    rec.R = matrixFull(R);

    const std::string combined = rec.A + "|" + rec.B + "|" + rec.Q + "|" + rec.R;
    rec.hash = djb2_64(combined);
    return rec;
}

void printRecord(const PlantRecord& rec, bool last) {
    std::printf("    {\n");
    std::printf("      \"name\": ");         printJsonString(rec.name); std::printf(",\n");
    std::printf("      \"plant_hash\": \"%016llx\",\n",
                (unsigned long long)rec.hash);
    std::printf("      \"n_states\": %d,\n", rec.n_states);
    std::printf("      \"n_inputs\": %d,\n", rec.n_inputs);
    std::printf("      \"A\": ");  printJsonString(rec.A); std::printf(",\n");
    std::printf("      \"B\": ");  printJsonString(rec.B); std::printf(",\n");
    std::printf("      \"Q\": ");  printJsonString(rec.Q); std::printf(",\n");
    std::printf("      \"R\": ");  printJsonString(rec.R); std::printf("\n");
    std::printf("    }");
    if (!last) std::printf(",");
    std::printf("\n");
}

}  // namespace

int main() {
    const std::vector<ModelEntry> models = getBuiltinModels();

    // Find the cascade model -- the only preset the product designs LQR against.
    const ModelEntry* cascade = nullptr;
    for (const auto& m : models)
        if (isCascadeModel(m)) { cascade = &m; break; }
    if (!cascade) {
        std::fprintf(stderr, "export_plants: no cascade model found\n");
        return 1;
    }

    const auto& presets = lqrPresets();
    const int n = cascade->system.states();
    const int m = cascade->system.inputs();

    std::vector<PlantRecord> records;
    records.reserve(presets.size());
    for (const auto& p : presets) {
        const Eigen::MatrixXd Q = presetStateWeights(p, n).asDiagonal().toDenseMatrix();
        const Eigen::MatrixXd R = presetInputWeights(p, m).asDiagonal().toDenseMatrix();
        records.push_back(makeRecord(
            std::string("Ball-Balancer Cascade / ") + p.name,
            cascade->system, Q, R));
    }

    std::printf("{\n");
    std::printf("  \"generator\": \"export_plants\",\n");
    std::printf("  \"records\": [\n");
    for (std::size_t i = 0; i < records.size(); ++i)
        printRecord(records[i], i + 1 == records.size());
    std::printf("  ]\n");
    std::printf("}\n");
    return 0;
}
