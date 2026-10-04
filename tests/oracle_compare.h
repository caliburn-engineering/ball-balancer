// tests/oracle_compare.h
//
// The comparison both live-oracle tests make between a freshly generated
// fixture and the checked-in tests/oracle_fixture.h.
//
// Line by line, over the code lines only.  Comment lines are prose that each
// generator words its own way (which tool made it, when, how to rerun it);
// everything a test consumes is in code.  The K(i, j) and P(i, j) assignments
// compare as numbers, to the default suite's relative tolerance, because a
// different oracle, LAPACK or CPU moves their last digits.  Every other code
// line — function names, shapes, the kHash_* constants — must match exactly.
// The Oracle comment is read separately, so a test can insist on its oracle.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace oracle_compare {

// Matches tests/test_oracle_fixtures.cpp: the measured agreement between the
// solver and scipy is ~1e-14, so 1e-10 absorbs any platform's last digits.
constexpr double kRelTol = 1e-10;

inline std::string read_file(const char* path) {
    std::ifstream f(path);
    if (!f) return {};
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

inline std::vector<std::string> lines_of(const std::string& text) {
    std::vector<std::string> out;
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("//", 0) == 0) continue;
        out.push_back(line);
    }
    return out;
}

inline std::string oracle_line(const std::string& text) {
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line))
        if (line.rfind("// Oracle:", 0) == 0) return line;
    return {};
}

// Splits "    K(0, 3) = -3.65;" into its "    K(0, 3) = " prefix and value.
// Returns false for any line that is not a K or P element assignment.
inline bool split_assignment(const std::string& line, std::string& prefix,
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
inline int compare(const std::vector<std::string>& fresh,
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

}  // namespace oracle_compare
