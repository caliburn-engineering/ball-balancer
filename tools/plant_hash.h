// tools/plant_hash.h
//
// The plant hash baked into tests/oracle_fixture.h.  One definition, shared by
// the exporter that writes it and the staleness test that checks it, so a
// change to the scheme moves both sides together instead of reading as a
// stale fixture that regenerating cannot fix.
//
// The hash is a DJB2-64 digest of A, B, Q and R as full-precision (%.17g)
// MATLAB-style strings joined by "|".  DJB2 needs no library, so the Python
// generator and any other tool can reproduce it.

#pragma once

#include "linear_system.h"

#include <Eigen/Core>
#include <cstdint>
#include <cstdio>
#include <string>

namespace caliburn {

// Full-precision MATLAB-style matrix string (%.17g per element).
// matrixToString uses %g (6 sig figs), which is not enough for an oracle:
// tiny rounding differences in the Jacobian coefficients propagate into the
// gain and make the fixture fail at tight tolerances.
inline std::string matrixFull(const Eigen::MatrixXd& m) {
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

inline uint64_t djb2_64(const std::string& s) {
    uint64_t h = 5381;
    for (unsigned char c : s) h = h * 33u + c;
    return h;
}

inline uint64_t plantHash(const LinearSystem& sys,
                          const Eigen::MatrixXd& Q,
                          const Eigen::MatrixXd& R) {
    return djb2_64(matrixFull(sys.A) + "|" + matrixFull(sys.B) + "|" +
                   matrixFull(Q) + "|" + matrixFull(R));
}

}  // namespace caliburn
