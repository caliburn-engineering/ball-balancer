// From caliburn reference/test/assert_rel.h @ f60305f
#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>

#define ASSERT_REL_NEAR(actual, expected, tol)                                 \
    do {                                                                       \
        double a_ = (actual), e_ = (expected), t_ = (tol);                    \
        double scale_ = std::max({std::abs(a_), std::abs(e_), 1.0});          \
        if (std::abs(a_ - e_) > t_ * scale_) {                                \
            std::fprintf(stderr, "FAIL: %s:%d: %s = %g, expected %g (rel tol %g)\n", \
                         __FILE__, __LINE__, #actual, a_, e_, t_);             \
            std::exit(1);                                                      \
        }                                                                      \
    } while (0)

#ifdef EIGEN_CORE_H
#define ASSERT_MATRIX_REL_NEAR(A, B, tol)                                     \
    do {                                                                       \
        double na_ = (A).norm(), nb_ = (B).norm(), t_ = (tol);                \
        double scale_ = std::max({na_, nb_, 1.0});                             \
        if (((A) - (B)).norm() > t_ * scale_) {                               \
            std::fprintf(stderr, "FAIL: %s:%d: matrix rel error = %g (rel tol %g)\n", \
                         __FILE__, __LINE__, ((A) - (B)).norm() / scale_, t_); \
            std::exit(1);                                                      \
        }                                                                      \
    } while (0)
#endif
