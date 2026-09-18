// src/analysis/controller_builders.h
#pragma once

#include "linear_system.h"
#include <Eigen/Core>
#include <vector>

namespace caliburn {

// Parameter blocks live here, not in app_state.h: app_state.h includes imgui.h,
// which analysis_lib and every test target cannot see.  app_state.h already
// pulls in six analysis headers, so it takes a seventh.  Fields stay `float`
// so ImGui sliders bind them directly; the builders widen to double.

struct PIDParams {
    float Kp = 1.0f;
    float Ki = 0.0f;
    float Kd = 0.0f;
    float tau_f = 0.01f;  // derivative filter time constant.  Precondition: > 0.
};

enum class CompensatorMode { Lead, Lag, LeadLag };

struct LeadLagParams {
    CompensatorMode mode = CompensatorMode::Lead;
    float Kc = 1.0f;
    float alpha_lead = 0.1f;  // < 1.  Precondition: > 0.
    float T_lead = 1.0f;      // Precondition: > 0.
    float alpha_lag = 10.0f;  // > 1.  Precondition: > 0.
    float T_lag = 1.0f;       // Precondition: > 0.
};

// One control loop: reads e[out] = r[out] - y[out], drives u[in].
// Both parameter blocks are carried inline and LoopKind selects which is read,
// so switching kind loses neither tuning nor the pairing — the pairing being
// the expensive decision.  See issue #2 for the rejected alternatives.
struct Loop {
    int out = 0;
    int in = 0;
    PIDParams pid{};
    LeadLagParams leadlag{};
};

// Not ControllerType: that enum lives in app_state.h and carries StateSpace and
// GainMatrix, which mean nothing to the assembler.
enum class LoopKind { PID, LeadLag };

// C(s) = Kp + Ki/s + Kd*s/(tau_f*s + 1), minimal realization.
LinearSystem buildPID(const PIDParams& p);

// C(s) = Kc * (s + 1/T) / (s + 1/(alpha*T)), minimal realization; the cascade
// mode is the two sections in series.
//
// Note the constant: the design spec writes Kc*(Ts+1)/(alpha*T*s+1) in prose but
// gives a realization computing the pole-zero form above — the two differ by a
// factor of alpha, which Kc absorbs.  The realization is what was verified, so
// it is what ships, and the pole-zero form is the contract the tests assert.
// Both forms agree at alpha = 1 (pure gain Kc).
LinearSystem buildLeadLag(const LeadLagParams& p);

// Assemble the m-out / p-in controller for a loop list against `plant`.
// The plant is passed rather than (int p, int m): two same-typed ints transpose
// silently and build a wrongly-shaped controller, and this matches the existing
// stateFeedbackClose(const LinearSystem&, const MatrixXd&).
LinearSystem buildLoopController(const std::vector<Loop>& loops,
                                 LoopKind kind,
                                 const LinearSystem& plant);

}  // namespace caliburn
