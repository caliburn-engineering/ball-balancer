// src/plate_design.h
#pragma once

#include "auto_balance.h"

#include <string>
#include <vector>

namespace caliburn {

struct AppState;
struct ModelEntry;

/// The design surface's current answer, as seen from the plate.
///
/// `offered` is the model panel's claim that LQR is selected, the solve
/// succeeded, and the plant is the cascade.  `reason` names the first
/// objection when it is not, so the panel can say why without the plate having
/// to guess.  `design` carries servo parameters regardless of `offered`, so
/// `servo_tau` and `home_leg_rad` reach the simulation every frame.
///
/// The plate adds two further checks in `setDesign`: the gain's shape, and
/// whether the mechanism in `design` is the one being simulated.  Those belong
/// at the seam only the plate can see — no caller outside the plate knows its
/// geometry or gravity.
struct Offer {
    AutoBalanceDesign design;
    bool offered = false;
    std::string reason;
};

/// The model panel's offer decision as a pure free function.
///
/// Assembles `design` from the current cascade parameter list, checks the
/// controller type and the LQR solve result, and returns an `Offer` with
/// `offered = true` only when all three pass.
///
/// No GL, ImGui or GLFW state is read or written: this is a computation over
/// the application's logical state, and the same answer would come out the same
/// way in a headless harness.  The two checks only the plate can make — gain
/// shape and mechanism match — are left to `PlateView::setDesign`.
Offer designToOffer(const AppState& state,
                    const std::vector<ModelEntry>& presets);

}  // namespace caliburn
