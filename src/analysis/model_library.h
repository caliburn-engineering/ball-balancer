// src/analysis/model_library.h
#pragma once

#include "linear_system.h"
#include "table_kinematics.h"

#include <Eigen/Core>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace caliburn {

struct PhysicalParam {
    std::string name;      // Display name ("Mass")
    std::string symbol;    // Short symbol ("m")
    std::string unit;      // Unit string ("kg")
    float value;           // Default/current value
    float min_val;         // Slider minimum
    float max_val;         // Slider maximum
    bool logarithmic = false;
};

struct DerivationStep {
    std::string title;     // "Equations of Motion (ODE)"
    std::string content;   // Multi-line explanation with equations
};

struct ModelEntry {
    std::string name;
    std::string description;
    LinearSystem system;

    // Physical modeling (optional — empty means no physical params)
    std::vector<PhysicalParam> params;
    std::vector<DerivationStep> derivation;
    std::function<LinearSystem(const std::vector<PhysicalParam>&)> builder;
};

// Built-in preset models.
std::vector<ModelEntry> getBuiltinModels();

// True for the Ball-Balancer Cascade, the one preset whose states, inputs and
// parameters the plate view knows how to drive.  Matched by name because that
// is the only thing that identifies a plant: `7 states, 3 inputs` is a shape,
// and shapes coincide.
bool isCascadeModel(const ModelEntry& entry);

// The mechanism a cascade parameter set describes.
//
// One function, three callers: the model builder linearises about it, the
// application hands it to the plate to check the gain was designed for the
// plate being simulated, and the tests rebuild it.  Transcribing it a second
// time is how the geometry check comes to reject two identical plates, which
// costs the balance loop entirely and says nothing on screen.
//
// `params` must be a cascade parameter list; anything missing falls back to
// that parameter's default rather than to zero, so a partial list yields the
// stock 3-RRS table and not a degenerate one.
TableParams cascadeMechanism(const std::vector<PhysicalParam>& params);

// Gravity, from the same list.  Not part of the mechanism — the plate's
// geometry and the acceleration acting on the ball are independent, and only
// one of them is a shape — but it is equally part of the plant a gain was
// designed against.
double cascadeGravity(const std::vector<PhysicalParam>& params);

// The first-order leg lag [s], floored where the model would divide by it.
// The simulated servo reads the same function, so the two cannot differ at the
// slider's bottom end.
double cascadeServoTau(const std::vector<PhysicalParam>& params);

// The leg angle the plant is linearised about [rad].  States 0..2 are
// deviations from this, so the plate has to regulate to the same number.
double cascadeHomeLegAngle(const std::vector<PhysicalParam>& params);

// The preset the application opens on.  Resolved by name, not by position:
// this is a ball-balancer product, and a plain index would silently select a
// different plant the next time the library gains an entry.  Falls back to 0.
int defaultModelIndex(const std::vector<ModelEntry>& models);

// Parse a matrix from a MATLAB-style string: "0 1; -2 -3"
// Rows separated by ';', columns by whitespace.
// Returns nullopt on parse failure.
std::optional<Eigen::MatrixXd> parseMatrix(const std::string& text);

// Format a matrix as MATLAB-style string: "0 1; -2 -3"
std::string matrixToString(const Eigen::MatrixXd& mat);

}  // namespace caliburn
