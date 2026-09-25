#!/usr/bin/env python3
"""gen_oracle_fixtures.py -- Generate oracle fixture constants for CARE tests.

Usage (two-step):
    ./export_plants > plants.json
    python3 tools/gen_oracle_fixtures.py plants.json > tests/oracle_fixture.h

Or pipe directly:
    ./export_plants | python3 tools/gen_oracle_fixtures.py - > tests/oracle_fixture.h

The generator reads the JSON emitted by the C++ export_plants tool, solves the
continuous-time algebraic Riccati equation (CARE) for each record using an
independent solver, and emits a C++ header of fixture constants for the test
suite to compare against.

Solver selection (in order of preference):
  1. scipy.linalg.solve_continuous_are  -- Schur-based (LAPACK DGEES); the
     preferred oracle, independent of the C++ sign-function implementation.
  2. Pure-Python sign-function fallback  -- no external dependencies; converges
     on the same invariant subspace as scipy but via a different algorithm.
     Used automatically when scipy is not installed.

Pinned versions for reproducibility:
  scipy>=1.11, numpy>=1.24  (see requirements.txt)
"""

from __future__ import annotations

import json
import math
import sys
from datetime import date
from typing import Optional

# ── Matrix helpers (pure Python, no numpy) ──────────────────────────────────

def mat_rows(m):       return len(m)
def mat_cols(m):       return len(m[0])
def mat_el(m, r, c):  return m[r][c]

def mat_zero(n, k=None):
    k = k or n
    return [[0.0] * k for _ in range(n)]

def mat_eye(n):
    m = mat_zero(n)
    for i in range(n): m[i][i] = 1.0
    return m

def mat_mul(a, b):
    ra, ca = mat_rows(a), mat_cols(a)
    rb, cb = mat_rows(b), mat_cols(b)
    assert ca == rb
    c = mat_zero(ra, cb)
    for i in range(ra):
        for j in range(cb):
            s = 0.0
            for k in range(ca): s += a[i][k] * b[k][j]
            c[i][j] = s
    return c

def mat_add(a, b):
    ra, ca = mat_rows(a), mat_cols(a)
    return [[a[i][j] + b[i][j] for j in range(ca)] for i in range(ra)]

def mat_sub(a, b):
    ra, ca = mat_rows(a), mat_cols(a)
    return [[a[i][j] - b[i][j] for j in range(ca)] for i in range(ra)]

def mat_scale(a, s):
    return [[a[i][j] * s for j in range(mat_cols(a))] for i in range(mat_rows(a))]

def mat_transpose(a):
    ra, ca = mat_rows(a), mat_cols(a)
    return [[a[i][j] for i in range(ra)] for j in range(ca)]

def mat_norm(a):
    s = 0.0
    for row in a:
        for x in row: s += x * x
    return math.sqrt(s)

def mat_diag(d):
    n = len(d)
    m = mat_zero(n)
    for i in range(n): m[i][i] = d[i]
    return m

def mat_block(a, b, c, d):
    """[[a, b], [c, d]]"""
    ra, ca = mat_rows(a), mat_cols(a)
    rb, cb = mat_rows(b), mat_cols(b)
    rc, cc = mat_rows(c), mat_cols(c)
    rd, cd = mat_rows(d), mat_cols(d)
    assert ca == cc and rb == rd and ra == rb and ca == cc
    n = ra + rc
    k = ca + cb
    m = mat_zero(n, k)
    for i in range(ra):
        for j in range(ca): m[i][j]     = a[i][j]
        for j in range(cb): m[i][ca+j]  = b[i][j]
    for i in range(rc):
        for j in range(cc): m[ra+i][j]  = c[i][j]
        for j in range(cd): m[ra+i][ca+j] = d[i][j]
    return m

def mat_block_get(m, r0, c0, nrows, ncols):
    return [[m[r0+i][c0+j] for j in range(ncols)] for i in range(nrows)]

def mat_inv(m):
    """Gauss-Jordan inversion."""
    n = mat_rows(m)
    aug = [list(row) + [1.0 if i == j else 0.0 for j in range(n)]
           for i, row in enumerate(m)]
    for col in range(n):
        # Pivot
        best = col
        for row in range(col+1, n):
            if abs(aug[row][col]) > abs(aug[best][col]): best = row
        aug[col], aug[best] = aug[best], aug[col]
        piv = aug[col][col]
        if abs(piv) < 1e-300:
            raise ValueError("singular matrix in mat_inv")
        for j in range(2*n): aug[col][j] /= piv
        for row in range(n):
            if row == col: continue
            f = aug[row][col]
            for j in range(2*n): aug[row][j] -= f * aug[col][j]
    return [row[n:] for row in aug]

def mat_det(m):
    """LU-based determinant."""
    n = mat_rows(m)
    a = [list(row) for row in m]
    sign = 1.0
    for col in range(n):
        best = col
        for row in range(col+1, n):
            if abs(a[row][col]) > abs(a[best][col]): best = row
        if best != col:
            a[col], a[best] = a[best], a[col]
            sign = -sign
        if abs(a[col][col]) < 1e-300: return 0.0
        for row in range(col+1, n):
            f = a[row][col] / a[col][col]
            for j in range(col, n): a[row][j] -= f * a[col][j]
    prod = sign
    for i in range(n): prod *= a[i][i]
    return prod

# ── CARE via matrix sign function (pure Python) ─────────────────────────────

def mat_lstsq(lhs, rhs):
    """Least-squares solution X of lhs * X = rhs via normal equations.

    lhs is (2n x n), rhs is (2n x n).  Implements the normal-equation solve
    X = (lhs' lhs)^{-1} lhs' rhs.  When lhs has full column rank (which it
    does after sign-function convergence) this is the unique minimum-norm
    solution and agrees with numpy.linalg.lstsq.
    """
    lt = mat_transpose(lhs)
    gram = mat_mul(lt, lhs)      # n x n, positive definite when lhs has full rank
    rhs2 = mat_mul(lt, rhs)      # n x n
    return mat_mul(mat_inv(gram), rhs2)


def care_sign_function(A, B, Q, R):
    """Solve  A'P + PA - P B R^-1 B' P + Q = 0  for the stabilising P.

    Algorithm: form the Hamiltonian H = [[A, B R^-1 B'], [Q, -A']] and apply
    Newton's iteration for the matrix sign function with determinant scaling.
    Extract P from the converged sign matrix by solving the overdetermined
    system [W12; W22+I] P = [W11+I; W21] in the least-squares sense.

    This is the same algorithm and the same Hamiltonian convention as the C++
    lqr.cpp implementation, re-implemented independently in Python.  A bug
    in the Hamiltonian formation or the sign-function recovery manifests
    differently when the two implementations are compared in test_oracle_fixtures.
    """
    n = mat_rows(A)
    Rinv    = mat_inv(R)
    BRinvBt = mat_mul(mat_mul(B, Rinv), mat_transpose(B))  # B R^-1 B'
    negAt   = mat_scale(mat_transpose(A), -1.0)            # -A'

    # H = [[A, B R^-1 B'], [Q, -A']] -- same sign convention as C++ lqr.cpp
    H = mat_block(A, BRinvBt, Q, negAt)
    Z = [list(row) for row in H]

    TOL = 1e-9
    MAX_ITER = 150
    two_n = float(2 * n)

    for _ in range(MAX_ITER):
        Z_old = Z
        det = mat_det(Z)
        if not math.isfinite(det) or abs(det) < 1e-300:
            raise ValueError("sign function: singular iterate")
        scale = abs(det) ** (-1.0 / two_n)
        Z = mat_scale(Z, scale)
        Zinv = mat_inv(Z)
        # Z <- (Z + Z^-1) / 2  ≡  Z - (Z - Z^-1)/2
        diff = mat_sub(Z, Zinv)
        Z = mat_sub(Z, mat_scale(diff, 0.5))
        delta = mat_sub(Z, Z_old)
        if mat_norm(delta) < TOL * max(1.0, mat_norm(Z)):
            break

    # Extract P from sign(H).  The converged sign matrix has blocks:
    #   Z = [[W11, W12], [W21, W22]]
    # P solves the overdetermined 2n x n system (consistent in exact arithmetic):
    #   [W12; W22 + I] P = [W11 + I; W21]
    # Solved here via normal equations; the C++ uses SVD for the same system.
    W11 = mat_block_get(Z, 0, 0, n, n)
    W12 = mat_block_get(Z, 0, n, n, n)
    W21 = mat_block_get(Z, n, 0, n, n)
    W22 = mat_block_get(Z, n, n, n, n)

    # Build lhs = [W12; W22+I] and rhs = [W11+I; W21] as 2n x n matrices.
    lhs = [[W12[i][j] for j in range(n)] for i in range(n)] + \
          [[(W22[i][j] + (1.0 if i == j else 0.0)) for j in range(n)] for i in range(n)]
    rhs = [[(W11[i][j] + (1.0 if i == j else 0.0)) for j in range(n)] for i in range(n)] + \
          [[W21[i][j] for j in range(n)] for i in range(n)]

    P = mat_lstsq(lhs, rhs)

    # Symmetrise to remove the rounding asymmetry between triangles.
    for i in range(n):
        for j in range(n):
            avg = 0.5 * (P[i][j] + P[j][i])
            P[i][j] = avg
            P[j][i] = avg
    return P


def care_scipy(A_np, B_np, Q_np, R_np):
    """Solve CARE with scipy.linalg.solve_continuous_are (Schur-based, LAPACK)."""
    import scipy.linalg as la  # pylint: disable=import-outside-toplevel
    return la.solve_continuous_are(A_np, B_np, Q_np, R_np)


def solve_care(A, B, Q, R):
    """Solve CARE, preferring scipy's Schur solver; fall back to sign function.

    Returns P as a pure-Python list-of-lists.
    """
    try:
        import numpy as np                             # noqa: F401
        A_np = np.array(A)
        B_np = np.array(B)
        Q_np = np.array(Q)
        R_np = np.array(R)
        P_np = care_scipy(A_np, B_np, Q_np, R_np)
        return P_np.tolist()
    except ImportError:
        pass

    # numpy/scipy unavailable -- use pure-Python sign-function solver.
    return care_sign_function(A, B, Q, R)


def gain_from_P(P, B, R):
    """K = R^-1 B' P"""
    return mat_mul(mat_inv(R), mat_mul(mat_transpose(B), P))


# ── Matrix parsing ───────────────────────────────────────────────────────────

def parse_matlab_matrix(s: str):
    """Parse a MATLAB-style matrix string into a list-of-lists of floats."""
    rows = []
    for row_s in s.split(";"):
        row_s = row_s.strip()
        if not row_s:
            continue
        row = [float(x) for x in row_s.split()]
        rows.append(row)
    return rows


# ── C++ header emission ──────────────────────────────────────────────────────

def _matrix_cpp_array(m, name: str, indent: str = "    ") -> str:
    """Emit a C-style 2-D initialiser list for the matrix."""
    rows = mat_rows(m)
    cols = mat_cols(m)
    lines = [f"// {rows}x{cols}",
             f"static const double {name}[{rows}][{cols}] = {{"]
    for i, row in enumerate(m):
        vals = ", ".join(f"{v:.17g}" for v in row)
        comma = "," if i + 1 < rows else ""
        lines.append(f"    {{{vals}}}{comma}")
    lines.append("};")
    return ("\n" + indent).join(lines)


def _ident(name: str) -> str:
    """Turn 'Ball-Balancer Cascade / Detuned' into 'cascade_detuned'."""
    # Strip the common prefix and produce a snake_case ident for the remainder.
    s = name.lower()
    # "ball-balancer cascade / detuned" → "cascade / detuned" → "cascade_detuned"
    s = s.replace("ball-balancer ", "")   # → "cascade / detuned"
    s = s.replace("/ ", "")              # → "cascade detuned"
    s = s.replace(" ", "_")             # → "cascade_detuned"
    s = s.replace("-", "_")
    return s.strip("_")


def emit_header(records: list[dict]) -> str:
    today = str(date.today())

    try:
        import scipy
        oracle_id = f"scipy {scipy.__version__} solve_continuous_are (Schur / LAPACK)"
    except ImportError:
        oracle_id = "pure-Python matrix sign function (no scipy)"

    lines = [
        "// tests/oracle_fixture.h",
        "// Generated by tools/gen_oracle_fixtures.py -- DO NOT EDIT.",
        "// Re-run the pipeline to regenerate:",
        "//   cmake --build build --target export_plants",
        "//   ./build/export_plants | python3 tools/gen_oracle_fixtures.py - > tests/oracle_fixture.h",
        "//",
        f"// Oracle: {oracle_id}",
        f"// Date  : {today}",
        "//",
        "// Each record: K (gain, m x n) and P (Riccati solution, n x n).",
        "// The plant_hash identifies the A, B, Q, R matrices the gain was",
        "// solved against; re-exporting and comparing hashes detects staleness.",
        "#pragma once",
        "",
        "#include <Eigen/Core>",
        "",
        "namespace caliburn::oracle {",
        "",
    ]

    for rec in records:
        ident = _ident(rec["name"])
        A  = parse_matlab_matrix(rec["A"])
        B  = parse_matlab_matrix(rec["B"])
        Q  = parse_matlab_matrix(rec["Q"])
        R  = parse_matlab_matrix(rec["R"])

        P = solve_care(A, B, Q, R)
        K = gain_from_P(P, B, R)

        n = mat_rows(A)
        m = mat_rows(R)
        lines += [
            f"// ── {rec['name']} ─────────────────────────────────────────",
            f"// plant_hash : {rec['plant_hash']}",
            f"// K : {m}x{n} gain (u = -K x), P : {n}x{n} Riccati solution",
            "",
            f"inline Eigen::MatrixXd oracle_K_{ident}() {{",
            f"    Eigen::MatrixXd K({m}, {n});",
        ]
        for i in range(m):
            for j in range(n):
                lines.append(f"    K({i}, {j}) = {K[i][j]:.17g};")
        lines += [
            "    return K;",
            "}",
            "",
            f"inline Eigen::MatrixXd oracle_P_{ident}() {{",
            f"    Eigen::MatrixXd P({n}, {n});",
        ]
        for i in range(n):
            for j in range(n):
                lines.append(f"    P({i}, {j}) = {P[i][j]:.17g};")
        lines += [
            "    return P;",
            "}",
            "",
            f"constexpr const char* kHash_{ident} = \"{rec['plant_hash']}\";",
            "",
        ]

    lines += [
        "}  // namespace caliburn::oracle",
        "",
    ]
    return "\n".join(lines)


# ── Entry point ──────────────────────────────────────────────────────────────

def main() -> None:
    if len(sys.argv) < 2 or sys.argv[1] == "-":
        data = json.load(sys.stdin)
    else:
        with open(sys.argv[1]) as fh:
            data = json.load(fh)

    records = data["records"]
    print(emit_header(records), end="")


if __name__ == "__main__":
    main()
