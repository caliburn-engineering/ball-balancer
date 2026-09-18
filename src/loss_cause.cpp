// src/loss_cause.cpp
#include "loss_cause.h"

namespace caliburn {

LossCause lossCause(const LossFlags& f) {
    // Written as a chain rather than as a table, because a chain is what the
    // rule IS: each test is reached only when everything upstream of it has
    // been ruled out, so the order on the page is the order of the claim.
    if (f.rates_untrustworthy) return LossCause::RatesUntrustworthy;
    if (f.clipped)             return LossCause::WorkspaceClipped;
    if (f.separated)           return LossCause::Separated;
    if (f.saturated)           return LossCause::Saturated;
    return LossCause::RolledOff;
}

LossFlags lossFlags(const SimReport& r, const SimState& s) {
    LossFlags f;
    // The plate's own credibility, from the motion it was last assembled with.
    f.rates_untrustworthy = !s.motion.rates_trustworthy;
    f.clipped = r.clipped;
    // `airborne`, and NOT `left_plate`: the two sit beside each other in
    // `SimReport` and mean opposite halves of this question.  `left_plate` is
    // the loss being explained; `airborne` is one of the things that might
    // explain it.
    f.separated = r.airborne;
    f.saturated = r.saturated;
    return f;
}

const char* lossCauseLabel(LossCause c) {
    switch (c) {
        case LossCause::RatesUntrustworthy: return "rates untrusted";
        case LossCause::WorkspaceClipped:   return "workspace-clipped";
        case LossCause::Separated:          return "separated";
        case LossCause::Saturated:          return "saturated";
        case LossCause::RolledOff:          return "rolled off";
    }
    return "rolled off";
}

const char* lossCauseSentence(LossCause c) {
    switch (c) {
        case LossCause::RatesUntrustworthy:
            return "The plate was too near a kinematic singularity for its own "
                   "velocity Jacobian to be believed, so nothing computed from "
                   "how it was moving — including whether it was still holding "
                   "the ball — is evidence of anything.";
        case LossCause::WorkspaceClipped:
            return "The loop asked the mechanism for a pose it cannot hold and "
                   "was pulled back to one it can, so the tilt that reached the "
                   "ball was not the tilt the gain wanted.  The plate ran out "
                   "of workspace before the gain ran out of authority.";
        case LossCause::Separated:
            return "The plate accelerated away from the ball faster than "
                   "gravity, so contact ended and the ball was flying when it "
                   "crossed the rim.  A plate that has let go cannot steer what "
                   "it let go of.";
        case LossCause::Saturated:
            return "A leg command stood on a travel limit, so the loop was "
                   "asking for more tilt than the servos have left to give.";
        case LossCause::RolledOff:
            return "Nothing failed.  The mechanism held every command, the legs "
                   "had travel in hand and the plate never let go — the ball "
                   "simply went further than the plate is wide.";
    }
    return "";
}

std::string lossFlagsLine(const LossFlags& f) {
    std::string s;
    auto add = [&s](const char* name) {
        if (!s.empty()) s += ", ";
        s += name;
    };
    // In the precedence's own order, so that the line reads as the working the
    // cause was chosen from rather than as an unordered set.
    if (f.rates_untrustworthy) add("rates untrustworthy");
    if (f.clipped)             add("clipped");
    if (f.separated)           add("separated");
    if (f.saturated)           add("saturated");
    return s.empty() ? std::string("none") : s;
}

}  // namespace caliburn
