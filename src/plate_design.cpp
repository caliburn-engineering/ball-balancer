// src/plate_design.cpp
#include "plate_design.h"
#include "app_state.h"
#include "analysis/model_library.h"
#include "sim_step.h"

namespace caliburn {

Offer designToOffer(const AppState& state,
                    const std::vector<ModelEntry>& presets) {
    Offer offer;

    const bool is_cascade =
        state.preset_index >= 0 &&
        state.preset_index < static_cast<int>(presets.size()) &&
        isCascadeModel(presets[state.preset_index]);

    // The operating point — servo tau, home angle, mechanism, gravity — is
    // always assembled from the current params so those values reach the plate
    // every frame, whether or not a gain is on offer.
    if (is_cascade)
        offer.design = cascadeDesign(state.current_params);

    if (!is_cascade) {
        offer.reason = "plant is not the Ball-Balancer Cascade";
    } else if (state.ctrl_type != ControllerType::LQR) {
        offer.reason = "select LQR as the controller type";
    } else if (!state.lqr_result.success) {
        offer.reason = "the LQR solve failed";
    } else {
        offer.design.K = state.lqr_result.K;
        offer.offered = true;
    }

    return offer;
}

}  // namespace caliburn
