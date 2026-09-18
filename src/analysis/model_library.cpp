// src/analysis/model_library.cpp
#include "model_library.h"

// The cascade model reads its mechanism block from the same kinematics the
// 3D view draws, rather than carrying a transcribed copy of the Jacobian.
#include "table_kinematics.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>

namespace caliburn {

namespace {
constexpr double kPi = 3.14159265358979323846;
}

namespace {

// The cascade's parameter list, read by symbol.  By symbol and not by index:
// the list is edited by hand and a reordered entry would otherwise silently
// swap two physical quantities.
double cascadeParam(const std::vector<PhysicalParam>& params,
                    const char* symbol, double fallback) {
    for (const auto& p : params)
        if (p.symbol == symbol) return static_cast<double>(p.value);
    return fallback;
}

// The servo travel the mechanism is built with.  Fixed rather than exposed:
// they are the hardware's limits, not a design choice, and the plate view and
// the model have to agree on them for a saturating command to mean anything.
constexpr double kAlphaMin = 10.0 * kPi / 180.0;
constexpr double kAlphaMax = 80.0 * kPi / 180.0;

}  // namespace

bool isCascadeModel(const ModelEntry& entry) {
    return entry.name.rfind("Ball-Balancer Cascade", 0) == 0;
}

TableParams cascadeMechanism(const std::vector<PhysicalParam>& params) {
    TableParams tp;
    tp.R_ground  = cascadeParam(params, "Rg", 0.300);
    tp.R_table   = cascadeParam(params, "Rt", 0.300);
    tp.L1        = cascadeParam(params, "L1", 0.150);
    tp.L2        = cascadeParam(params, "L2", 0.150);
    tp.alpha_min = kAlphaMin;
    tp.alpha_max = kAlphaMax;
    return tp;
}

double cascadeGravity(const std::vector<PhysicalParam>& params) {
    return cascadeParam(params, "g", 9.81);
}

double cascadeServoTau(const std::vector<PhysicalParam>& params) {
    return std::max(1e-3, cascadeParam(params, "\xcf\x84", 0.05));
}

double cascadeHomeLegAngle(const std::vector<PhysicalParam>& params) {
    return cascadeParam(params, "a0", 45.0) * kPi / 180.0;
}

int defaultModelIndex(const std::vector<ModelEntry>& models) {
    for (std::size_t i = 0; i < models.size(); ++i)
        if (models[i].name == "Ball-Balancer Cascade") return static_cast<int>(i);
    return 0;
}

std::optional<Eigen::MatrixXd> parseMatrix(const std::string& text) {
    std::vector<std::vector<double>> rows;
    std::istringstream row_stream(text);
    std::string row_str;

    while (std::getline(row_stream, row_str, ';')) {
        std::istringstream col_stream(row_str);
        std::vector<double> row;
        double val;
        while (col_stream >> val) {
            row.push_back(val);
        }
        // Check for non-numeric content remaining
        if (col_stream.fail() && !col_stream.eof()) {
            return std::nullopt;
        }
        if (!row.empty()) {
            rows.push_back(row);
        }
    }

    if (rows.empty()) return std::nullopt;

    size_t cols = rows[0].size();
    for (const auto& r : rows) {
        if (r.size() != cols) return std::nullopt;
    }

    Eigen::MatrixXd mat(rows.size(), cols);
    for (size_t i = 0; i < rows.size(); ++i) {
        for (size_t j = 0; j < cols; ++j) {
            mat(static_cast<int>(i), static_cast<int>(j)) = rows[i][j];
        }
    }
    return mat;
}

std::string matrixToString(const Eigen::MatrixXd& mat) {
    std::string text;
    for (int i = 0; i < mat.rows(); ++i) {
        if (i > 0) text += "; ";
        for (int j = 0; j < mat.cols(); ++j) {
            if (j > 0) text += " ";
            char num[32];
            std::snprintf(num, sizeof(num), "%g", mat(i, j));
            text += num;
        }
    }
    return text;
}

std::vector<ModelEntry> getBuiltinModels() {
    std::vector<ModelEntry> models;

    // =====================================================================
    // 1. First-Order (K=1, tau=1)
    // =====================================================================
    {
        ModelEntry entry;
        entry.name = "First-Order";
        entry.description = "1/(s+1), \xcf\x84=1s";

        entry.params = {
            {"Gain", "K", "", 1.0f, 0.01f, 10.0f, false},
            {"Time constant", "\xcf\x84", "s", 1.0f, 0.01f, 100.0f, true},
        };

        entry.builder = [](const std::vector<PhysicalParam>& p) -> LinearSystem {
            double K = p[0].value, tau = p[1].value;
            LinearSystem sys;
            sys.A = (Eigen::MatrixXd(1, 1) << -1.0 / tau).finished();
            sys.B = (Eigen::MatrixXd(1, 1) << 1.0 / tau).finished();
            sys.C = (Eigen::MatrixXd(1, 1) << K).finished();
            sys.D = Eigen::MatrixXd::Zero(1, 1);
            return sys;
        };
        entry.system = entry.builder(entry.params);

        entry.derivation = {
            {"Physical System",
             "A first-order system has one energy storage element.\n"
             "Examples:\n"
             "  - RC low-pass filter (capacitor stores charge)\n"
             "  - Thermal system (mass stores heat)\n"
             "  - Hydraulic tank (volume stores fluid)\n\n"
             "A single time constant \xcf\x84 governs the speed of response.\n"
             "DC gain K sets the steady-state output per unit input."},

            {"Governing ODE",
             "The first-order ODE relating input u(t) to output y(t):\n\n"
             "  \xcf\x84 \xc2\xb7 dy/dt + y = K \xc2\xb7 u\n\n"
             "where:\n"
             "  y(t) = output (voltage, temperature, level, ...)\n"
             "  u(t) = input (source, heat flux, flow, ...)\n"
             "  \xcf\x84 = time constant [s]\n"
             "  K = DC gain [-]"},

            {"Laplace Transform",
             "Apply Laplace transform (zero initial conditions):\n\n"
             "  \xcf\x84\xc2\xb7s\xc2\xb7Y(s) + Y(s) = K\xc2\xb7U(s)\n"
             "  Y(s)\xc2\xb7(\xcf\x84s + 1) = K\xc2\xb7U(s)\n\n"
             "Transfer function:\n"
             "  G(s) = Y(s)/U(s) = K / (\xcf\x84s + 1)\n\n"
             "Properties:\n"
             "  Pole: s = -1/\xcf\x84 (stable for \xcf\x84 > 0)\n"
             "  DC gain: G(0) = K\n"
             "  Bandwidth: \xcf\x89_bw = 1/\xcf\x84 rad/s\n"
             "  Step response settles in ~4\xcf\x84 seconds"},

            {"State-Space Form",
             "Define state x such that y = K\xc2\xb7x:\n\n"
             "  dx/dt = -(1/\xcf\x84)\xc2\xb7x + (1/\xcf\x84)\xc2\xb7u\n"
             "  y = K\xc2\xb7x\n\n"
             "Matrix form:\n"
             "  A = [-1/\xcf\x84]    B = [1/\xcf\x84]\n"
             "  C = [K]       D = [0]\n\n"
             "The eigenvalue of A is -1/\xcf\x84 = the pole of G(s)."},
        };

        models.push_back(std::move(entry));
    }

    // =====================================================================
    // 2. Second-Order Mass-Spring-Damper (m=1, b=1.4, k=1)
    // =====================================================================
    {
        ModelEntry entry;
        entry.name = "Second-Order";
        entry.description = "Mass-spring-damper (\xcf\x89n=1, \xce\xb6=0.7)";

        entry.params = {
            {"Mass", "m", "kg", 1.0f, 0.1f, 100.0f, true},
            {"Damping", "b", "N\xc2\xb7s/m", 1.4f, 0.0f, 50.0f, false},
            {"Stiffness", "k", "N/m", 1.0f, 0.01f, 100.0f, true},
        };

        entry.builder = [](const std::vector<PhysicalParam>& p) -> LinearSystem {
            double m = p[0].value, b = p[1].value, k = p[2].value;
            LinearSystem sys;
            sys.A = (Eigen::MatrixXd(2, 2) << 0, 1, -k / m, -b / m).finished();
            sys.B = (Eigen::MatrixXd(2, 1) << 0, 1.0 / m).finished();
            sys.C = (Eigen::MatrixXd(1, 2) << 1, 0).finished();
            sys.D = Eigen::MatrixXd::Zero(1, 1);
            return sys;
        };
        entry.system = entry.builder(entry.params);

        entry.derivation = {
            {"Physical System",
             "A mass m connected to a fixed wall by a spring (stiffness k)\n"
             "and a viscous damper (coefficient b). An external force F\n"
             "is applied to the mass.\n\n"
             "This models many real systems:\n"
             "  - Vehicle suspension element\n"
             "  - Mechanical actuator with compliance\n"
             "  - Loudspeaker cone dynamics\n"
             "  - Structural vibration mode"},

            {"Equations of Motion (ODE)",
             "Newton's second law (sum of forces = m\xc2\xb7" "a):\n\n"
             "  m\xc2\xb7" "d\xc2\xb2x/dt\xc2\xb2 = -k\xc2\xb7x - b\xc2\xb7" "dx/dt + F\n\n"
             "Rearranged:\n"
             "  m\xc2\xb7x'' + b\xc2\xb7x' + k\xc2\xb7x = F\n\n"
             "where:\n"
             "  x = displacement from equilibrium [m]\n"
             "  F = applied force [N]\n"
             "  m = mass [kg]\n"
             "  b = viscous damping coefficient [N\xc2\xb7s/m]\n"
             "  k = spring stiffness [N/m]"},

            {"Laplace Transform",
             "Applying Laplace transform:\n\n"
             "  (m\xc2\xb7s\xc2\xb2 + b\xc2\xb7s + k)\xc2\xb7X(s) = F(s)\n\n"
             "Transfer function:\n"
             "  G(s) = X(s)/F(s) = 1 / (m\xc2\xb7s\xc2\xb2 + b\xc2\xb7s + k)\n\n"
             "Standard form (divide by m):\n"
             "  G(s) = (1/m) / (s\xc2\xb2 + 2\xce\xb6\xcf\x89n\xc2\xb7s + \xcf\x89n\xc2\xb2)\n\n"
             "Natural frequency and damping ratio:\n"
             "  \xcf\x89n = sqrt(k/m)       [rad/s]\n"
             "  \xce\xb6 = b / (2\xc2\xb7sqrt(m\xc2\xb7k))  [-]\n\n"
             "Poles: s = -\xce\xb6\xcf\x89n \xc2\xb1 \xcf\x89n\xc2\xb7sqrt(\xce\xb6\xc2\xb2 - 1)\n"
             "  \xce\xb6 < 1: underdamped (oscillates)\n"
             "  \xce\xb6 = 1: critically damped\n"
             "  \xce\xb6 > 1: overdamped (no oscillation)"},

            {"State-Space Form",
             "Choose states: x1 = position, x2 = velocity\n\n"
             "  x1' = x2\n"
             "  x2' = -(k/m)\xc2\xb7x1 - (b/m)\xc2\xb7x2 + (1/m)\xc2\xb7" "F\n\n"
             "Matrix form:\n"
             "  A = [  0       1   ]    B = [  0  ]\n"
             "      [-k/m   -b/m   ]        [1/m  ]\n\n"
             "  C = [1  0]              D = [0]\n\n"
             "Eigenvalues of A are the poles of G(s)."},
        };

        models.push_back(std::move(entry));
    }

    // =====================================================================
    // 3. Ball-Balancer (ball rolling on tilting plate)
    // =====================================================================
    {
        ModelEntry entry;
        entry.name = "Ball-Balancer";
        entry.description = "Ball on plate, 4 states, 2 in, 2 out";

        entry.params = {
            {"Gravity", "g", "m/s\xc2\xb2", 9.81f, 1.0f, 20.0f, false},
        };

        entry.builder = [](const std::vector<PhysicalParam>& p) -> LinearSystem {
            double g = p[0].value;
            constexpr double k = 5.0 / 7.0;  // solid sphere rolling factor
            LinearSystem sys;
            sys.A = Eigen::MatrixXd::Zero(4, 4);
            sys.A(0, 2) = 1.0;
            sys.A(1, 3) = 1.0;
            sys.B = Eigen::MatrixXd::Zero(4, 2);
            sys.B(2, 1) = k * g;
            sys.B(3, 0) = k * g;
            sys.C = Eigen::MatrixXd::Zero(2, 4);
            sys.C(0, 0) = 1.0;
            sys.C(1, 1) = 1.0;
            sys.D = Eigen::MatrixXd::Zero(2, 2);
            return sys;
        };
        entry.system = entry.builder(entry.params);

        entry.derivation = {
            {"Physical System",
             "A ball rolling freely on a flat plate that can tilt\n"
             "in two axes. The plate tilt angles \xce\xb8x and \xce\xb8y are\n"
             "the control inputs. Ball position (x, y) is measured.\n\n"
             "The ball rolls without slipping, so rotational inertia\n"
             "adds to the effective translational mass."},

            {"Equations of Motion (ODE)",
             "For a sphere rolling without slip on an incline:\n\n"
             "  (m + I/r\xc2\xb2)\xc2\xb7" "a = m\xc2\xb7g\xc2\xb7sin(\xce\xb8)\n\n"
             "For a solid sphere, I = (2/5)\xc2\xb7m\xc2\xb7r\xc2\xb2:\n"
             "  (m + 2m/5)\xc2\xb7" "a = m\xc2\xb7g\xc2\xb7sin(\xce\xb8)\n"
             "  (7/5)\xc2\xb7" "a = g\xc2\xb7sin(\xce\xb8)\n\n"
             "Linearize for small angles (sin\xce\xb8 ~= \xce\xb8):\n"
             "  x'' = (5g/7)\xc2\xb7\xce\xb8y     (tilt about y moves ball in x)\n"
             "  y'' = (5g/7)\xc2\xb7\xce\xb8x     (tilt about x moves ball in y)"},

            {"Transfer Function",
             "Each axis is a double integrator with gain:\n\n"
             "  G(s) = X(s)/\xce\x98(s) = (5g/7) / s\xc2\xb2\n\n"
             "Two poles at s = 0 (marginally stable).\n"
             "No damping \xe2\x80\x94 the ball never stops on its own.\n"
             "Active feedback control is essential."},

            {"State-Space Form",
             "States: [x, y, x', y']\n"
             "Inputs: [\xce\xb8x, \xce\xb8y]\n"
             "Outputs: [x, y]\n\n"
             "  A = [0 0 1 0]    B = [  0      0   ]\n"
             "      [0 0 0 1]        [  0      0   ]\n"
             "      [0 0 0 0]        [  0    5g/7  ]\n"
             "      [0 0 0 0]        [5g/7    0    ]\n\n"
             "  C = [1 0 0 0]    D = [0 0]\n"
             "      [0 1 0 0]        [0 0]\n\n"
             "Double integrator: needs at least PD control."},
        };

        models.push_back(std::move(entry));
    }

    // =====================================================================
    // 4. Ball-Balancer Cascade (legs -> plate -> ball)
    // =====================================================================
    //
    // The plant the product is actually about.  Model 3 takes plate tilt as
    // its input, which is not a thing anyone can command: the real actuators
    // are three servo-driven legs, and the tilt they produce is the *middle*
    // of the chain.  This model spans the whole chain, so a controller
    // designed against it commands what the hardware exposes.
    //
    // The mechanism block is not transcribed here.  It is read from the same
    // TableKinematics the 3D view draws, at the same home pose, so the model
    // and the picture cannot disagree — and the geometry sliders propagate
    // into the controller design rather than only into the rendering.
    {
        ModelEntry entry;
        entry.name = "Ball-Balancer Cascade";
        entry.description = "Legs -> plate -> ball, 7 states, 3 in, 5 out";

        entry.params = {
            {"Gravity",          "g",             "m/s\xc2\xb2", 9.81f,  1.0f,  20.0f, false},
            {"Servo lag",        "\xcf\x84",      "s",           0.05f,  0.01f,  0.50f, true},
            {"Ground radius",    "Rg",            "m",           0.300f, 0.10f,  0.60f, false},
            {"Table radius",     "Rt",            "m",           0.300f, 0.10f,  0.60f, false},
            {"Lower leg",        "L1",            "m",           0.150f, 0.05f,  0.40f, false},
            {"Upper leg",        "L2",            "m",           0.150f, 0.05f,  0.40f, false},
            {"Home leg angle",   "a0",            "deg",        45.0f,  20.0f,  70.0f, false},
        };

        entry.builder = [](const std::vector<PhysicalParam>& p) -> LinearSystem {
            const double g     = cascadeGravity(p);
            const double tau   = cascadeServoTau(p);
            const double a0    = cascadeHomeLegAngle(p);
            constexpr double k = 5.0 / 7.0;  // solid sphere rolling factor

            // Not transcribed here: the same function the application hands to
            // the plate, so the mechanism the gain is designed against and the
            // one it is checked against cannot be two different tables.
            const TableKinematics tk(cascadeMechanism(p));
            const TablePose home = tk.home_pose(a0);
            Eigen::Matrix3d Jv = tk.velocity_jacobian({a0, a0, a0}, home);

            // A geometry the sliders can reach but the mechanism cannot close
            // makes J_pose singular, and the solve returns infinities.  Fall
            // back to a decoupled unit map rather than handing NaN to the
            // Riccati solver, which would fail with a mystifying message.
            if (!Jv.allFinite()) Jv = Eigen::Matrix3d::Identity();

            LinearSystem sys;
            sys.A = Eigen::MatrixXd::Zero(7, 7);
            sys.B = Eigen::MatrixXd::Zero(7, 3);
            sys.C = Eigen::MatrixXd::Zero(5, 7);
            sys.D = Eigen::MatrixXd::Zero(5, 3);

            // States 0..2: leg angle deviation from home, first-order servo.
            for (int i = 0; i < 3; ++i) {
                sys.A(i, i) = -1.0 / tau;
                sys.B(i, i) =  1.0 / tau;
            }

            // States 3..6: ball [x, y, vx, vy] in the plate frame.
            sys.A(3, 5) = 1.0;
            sys.A(4, 6) = 1.0;

            // Tilt reaches the ball through the mechanism.  Pitch (row 1 of
            // Jv) drives x; roll (row 0) drives y with a sign flip, because
            // R = Ry(theta)*Rx(phi) lifts the +Y edge as phi grows.  Same
            // convention as ballTiltFromPose — see CONTEXT.md.
            for (int j = 0; j < 3; ++j) {
                sys.A(5, j) =  k * g * Jv(1, j);
                sys.A(6, j) = -k * g * Jv(0, j);
            }

            // Outputs: roll, pitch, heave, then ball x and y.
            sys.C.block(0, 0, 3, 3) = Jv;
            sys.C(3, 3) = 1.0;
            sys.C(4, 4) = 1.0;
            return sys;
        };
        entry.system = entry.builder(entry.params);

        entry.derivation = {
            {"Physical System",
             "A 3-RRS parallel mechanism: three servo-driven legs support a\n"
             "circular table with three degrees of freedom \xe2\x80\x94 roll (\xcf\x86), pitch\n"
             "(\xce\xb8) and heave (z). A ball rolls freely on the table.\n\n"
             "The legs are the only actuators. Plate tilt is not commanded\n"
             "directly; it is what the leg angles produce. Model 3 starts at\n"
             "the tilt and so describes a plant nobody can drive \xe2\x80\x94 this one\n"
             "starts at the leg commands, which is what the hardware takes."},

            {"Mechanism Kinematics",
             "Leg angles and table pose are tied by three loop-closure\n"
             "constraints, f(\xce\xb1, p) = 0, with p = [\xcf\x86, \xce\xb8, z]. Implicit\n"
             "differentiation gives the map from leg motion to pose motion:\n\n"
             "  J_pose\xc2\xb7\xce\xb4p + J_\xce\xb1\xc2\xb7\xce\xb4\xce\xb1 = 0\n"
             "  \xce\xb4p = -J_pose\xe2\x81\xbb\xc2\xb9\xc2\xb7J_\xce\xb1\xc2\xb7\xce\xb4\xce\xb1 = Jv\xc2\xb7\xce\xb4\xce\xb1\n\n"
             "Jv is evaluated once at the home pose (all legs at \xce\xb1\xe2\x82\x80) and is\n"
             "fully coupled: every leg moves all three pose axes. Its\n"
             "condition number is the Jacobian Analysis readout in the Plate\n"
             "Control panel \xe2\x80\x94 near a singularity the plant loses authority\n"
             "over one pose direction and the design degrades with it."},

            {"Servo Dynamics",
             "Each leg is a position-commanded servo with first-order lag:\n\n"
             "  \xce\xb4\xce\xb1\xe1\xb5\xa2' = (1/\xcf\x84)\xc2\xb7(\xce\xb4\xce\xb1\xe1\xb5\xa2,cmd - \xce\xb4\xce\xb1\xe1\xb5\xa2)\n\n"
             "Without this the leg angle would be an algebraic input, the\n"
             "mechanism would contribute no states, and there would be nothing\n"
             "for a state-feedback design to act on. The lag is also what\n"
             "makes an aggressive tuning saturate rather than track."},

            {"Ball Dynamics",
             "A solid sphere rolling without slip on an incline:\n\n"
             "  (m + I/r\xc2\xb2)\xc2\xb7" "a = m\xc2\xb7g\xc2\xb7sin(\xce\xb8),  I = (2/5)m\xc2\xb7r\xc2\xb2\n"
             "  (7/5)\xc2\xb7" "a = g\xc2\xb7sin(\xce\xb8)\n\n"
             "Linearised for small angles:\n\n"
             "  x'' = +(5g/7)\xc2\xb7\xce\xb8      (pitch runs the ball along +x)\n"
             "  y'' = -(5g/7)\xc2\xb7\xcf\x86      (roll runs the ball along -y)\n\n"
             "The sign on \xcf\x86 is not cosmetic. The table rotation is\n"
             "R = Ry(\xce\xb8)\xc2\xb7Rx(\xcf\x86), so a positive \xcf\x86 LIFTS the +y edge and the\n"
             "ball runs downhill the other way. Getting it wrong yields a\n"
             "model that compiles, simulates, and rolls uphill."},

            {"Cascade / State-Space Form",
             "Substituting \xce\xb4p = Jv\xc2\xb7\xce\xb4\xce\xb1 into the ball equations couples the\n"
             "two halves without adding pose states \xe2\x80\x94 pose is algebraic in \xce\xb4\xce\xb1.\n\n"
             "States:  [\xce\xb4\xce\xb1\xe2\x82\x81, \xce\xb4\xce\xb1\xe2\x82\x82, \xce\xb4\xce\xb1\xe2\x82\x83, x, y, x', y']\n"
             "Inputs:  [\xce\xb4\xce\xb1\xe2\x82\x81,cmd, \xce\xb4\xce\xb1\xe2\x82\x82,cmd, \xce\xb4\xce\xb1\xe2\x82\x83,cmd]\n"
             "Outputs: [\xcf\x86, \xce\xb8, z, x, y]\n\n"
             "  A = [ -I/\xcf\x84        0    0 ]   B = [ I/\xcf\x84 ]\n"
             "      [    0         0    I ]       [  0  ]\n"
             "      [ M\xc2\xb7Jv         0    0 ]       [  0  ]\n\n"
             "  M = (5g/7)\xc2\xb7[ +row\xe2\x82\x81(Jv) ; -row\xe2\x82\x80(Jv) ]\n\n"
             "  C = [ Jv   0   0 ]           D = 0\n"
             "      [ 0    I   0 ]\n\n"
             "Three servo poles at -1/\xcf\x84 and four at the origin: the ball is a\n"
             "double integrator in each axis and will never settle on its own.\n"
             "Controllable from three legs, so LQR has a stabilising solution."},
        };

        models.push_back(std::move(entry));
    }

    // =====================================================================
    // 5. Inverted Pendulum on Cart
    // =====================================================================
    {
        ModelEntry entry;
        entry.name = "Inverted Pendulum";
        entry.description = "Cart-pendulum (unstable)";

        entry.params = {
            {"Cart mass", "M", "kg", 1.0f, 0.1f, 50.0f, true},
            {"Pendulum mass", "m", "kg", 0.1f, 0.01f, 10.0f, true},
            {"Pendulum length", "l", "m", 0.5f, 0.05f, 5.0f, true},
            {"Gravity", "g", "m/s\xc2\xb2", 9.81f, 1.0f, 20.0f, false},
        };

        entry.builder = [](const std::vector<PhysicalParam>& p) -> LinearSystem {
            double M = p[0].value, m = p[1].value;
            double l = p[2].value, g = p[3].value;
            LinearSystem sys;
            sys.A = Eigen::MatrixXd::Zero(4, 4);
            sys.A(0, 1) = 1.0;
            sys.A(1, 2) = -m * g / M;
            sys.A(2, 3) = 1.0;
            sys.A(3, 2) = (M + m) * g / (M * l);
            sys.B = (Eigen::MatrixXd(4, 1) << 0, 1.0 / M, 0,
                     -1.0 / (M * l))
                        .finished();
            sys.C = Eigen::MatrixXd::Zero(2, 4);
            sys.C(0, 0) = 1.0;
            sys.C(1, 2) = 1.0;
            sys.D = Eigen::MatrixXd::Zero(2, 1);
            return sys;
        };
        entry.system = entry.builder(entry.params);

        entry.derivation = {
            {"Physical System",
             "A rigid pendulum (mass m, length l) mounted on a cart\n"
             "(mass M). The cart moves on a frictionless track driven\n"
             "by horizontal force F. The pendulum is balanced in the\n"
             "inverted (upright) position.\n\n"
             "Classic benchmark for control design \xe2\x80\x94 inherently unstable."},

            {"Equations of Motion (ODE)",
             "Nonlinear Newton-Euler equations:\n\n"
             "  (M+m)\xc2\xb7x'' + m\xc2\xb7l\xc2\xb7\xce\xb8''\xc2\xb7" "cos\xce\xb8 - m\xc2\xb7l\xc2\xb7\xce\xb8'\xc2\xb2\xc2\xb7sin\xce\xb8 = F\n"
             "  m\xc2\xb7l\xc2\xb7x''\xc2\xb7" "cos\xce\xb8 + m\xc2\xb7l\xc2\xb2\xc2\xb7\xce\xb8'' - m\xc2\xb7g\xc2\xb7l\xc2\xb7sin\xce\xb8 = 0\n\n"
             "Linearize about \xce\xb8=0 (upright):\n"
             "  cos\xce\xb8 ~= 1, sin\xce\xb8 ~= \xce\xb8, \xce\xb8'\xc2\xb2 ~= 0\n\n"
             "  (M+m)\xc2\xb7x'' + m\xc2\xb7l\xc2\xb7\xce\xb8'' = F\n"
             "  x'' + l\xc2\xb7\xce\xb8'' = g\xc2\xb7\xce\xb8\n\n"
             "Solve for accelerations:\n"
             "  x'' = -(m\xc2\xb7g/M)\xc2\xb7\xce\xb8 + (1/M)\xc2\xb7" "F\n"
             "  \xce\xb8'' = ((M+m)\xc2\xb7g/(M\xc2\xb7l))\xc2\xb7\xce\xb8 - (1/(M\xc2\xb7l))\xc2\xb7" "F"},

            {"Stability Analysis",
             "The linearized system has a positive eigenvalue:\n"
             "  s = +sqrt((M+m)\xc2\xb7g/(M\xc2\xb7l))\n\n"
             "This pole in the right half-plane makes the system\n"
             "open-loop unstable \xe2\x80\x94 the pendulum falls without\n"
             "active control.\n\n"
             "Increasing m/M ratio or l makes the unstable pole\n"
             "slower (easier to control). Decreasing l makes it\n"
             "faster (harder to balance)."},

            {"State-Space Form",
             "States: [x, x', \xce\xb8, \xce\xb8']\n"
             "Input: F (force on cart)\n"
             "Outputs: [x, \xce\xb8]\n\n"
             "  A = [0      1         0            0   ]\n"
             "      [0      0      -m\xc2\xb7g/M         0   ]\n"
             "      [0      0         0            1   ]\n"
             "      [0      0   (M+m)g/(Ml)        0   ]\n\n"
             "  B = [0; 1/M; 0; -1/(Ml)]\n"
             "  C = [1 0 0 0; 0 0 1 0]\n"
             "  D = [0; 0]"},
        };

        models.push_back(std::move(entry));
    }

    // =====================================================================
    // 6. Quarter-Car Suspension
    // =====================================================================
    {
        ModelEntry entry;
        entry.name = "Quarter-Car";
        entry.description = "Suspension (sprung + unsprung)";

        entry.params = {
            {"Sprung mass", "ms", "kg", 300.0f, 50.0f, 2000.0f, true},
            {"Unsprung mass", "mu", "kg", 50.0f, 10.0f, 200.0f, true},
            {"Spring stiffness", "ks", "N/m", 20000.0f, 1000.0f, 200000.0f, true},
            {"Tire stiffness", "kt", "N/m", 200000.0f, 50000.0f, 500000.0f, true},
            {"Damping", "bs", "N\xc2\xb7s/m", 1000.0f, 100.0f, 10000.0f, true},
        };

        entry.builder = [](const std::vector<PhysicalParam>& p) -> LinearSystem {
            double ms = p[0].value, mu = p[1].value;
            double ks = p[2].value, kt = p[3].value, bs = p[4].value;
            LinearSystem sys;
            sys.A = Eigen::MatrixXd::Zero(4, 4);
            sys.A(0, 1) = 1.0;
            sys.A(1, 0) = -ks / ms;
            sys.A(1, 1) = -bs / ms;
            sys.A(1, 2) = ks / ms;
            sys.A(1, 3) = bs / ms;
            sys.A(2, 3) = 1.0;
            sys.A(3, 0) = ks / mu;
            sys.A(3, 1) = bs / mu;
            sys.A(3, 2) = -(ks + kt) / mu;
            sys.A(3, 3) = -bs / mu;
            sys.B = (Eigen::MatrixXd(4, 1) << 0, 0, 0, kt / mu).finished();
            sys.C = Eigen::MatrixXd::Zero(2, 4);
            sys.C(0, 0) = 1.0;
            sys.C(1, 2) = 1.0;
            sys.D = Eigen::MatrixXd::Zero(2, 1);
            return sys;
        };
        entry.system = entry.builder(entry.params);

        entry.derivation = {
            {"Physical System",
             "Two-mass model of a vehicle suspension:\n"
             "  - Sprung mass ms: car body, passengers, cargo\n"
             "  - Unsprung mass mu: wheel, tire, axle, brakes\n"
             "  - Suspension spring ks and damper bs connect the masses\n"
             "  - Tire stiffness kt connects wheel to road\n\n"
             "Road surface profile zr is the disturbance input."},

            {"Equations of Motion (ODE)",
             "Newton's second law for each mass:\n\n"
             "Sprung mass (car body):\n"
             "  ms\xc2\xb7zs'' = -ks\xc2\xb7(zs - zu) - bs\xc2\xb7(zs' - zu')\n\n"
             "Unsprung mass (wheel):\n"
             "  mu\xc2\xb7zu'' = ks\xc2\xb7(zs - zu) + bs\xc2\xb7(zs' - zu') - kt\xc2\xb7(zu - zr)\n\n"
             "where:\n"
             "  zs = sprung mass displacement [m]\n"
             "  zu = unsprung mass displacement [m]\n"
             "  zr = road profile input [m]"},

            {"Design Trade-offs",
             "Increasing suspension stiffness (ks):\n"
             "  + Reduces body roll, better handling\n"
             "  - Harsher ride, less comfort\n\n"
             "Increasing damping (bs):\n"
             "  + Controls resonance peaks\n"
             "  - Transmits more high-frequency vibration\n\n"
             "Key frequencies:\n"
             "  Body bounce: ~(1/2\xcf\x80)\xc2\xb7sqrt(ks/ms) ~ 1 Hz\n"
             "  Wheel hop:   ~(1/2\xcf\x80)\xc2\xb7sqrt(kt/mu) ~ 10 Hz"},

            {"State-Space Form",
             "States: [zs, zs', zu, zu']\n"
             "Input: zr (road profile)\n"
             "Outputs: [zs, zu]\n\n"
             "  A = [  0        1          0          0     ]\n"
             "      [-ks/ms  -bs/ms     ks/ms      bs/ms   ]\n"
             "      [  0        0          0          1     ]\n"
             "      [ ks/mu   bs/mu  -(ks+kt)/mu  -bs/mu   ]\n\n"
             "  B = [0; 0; 0; kt/mu]\n"
             "  C = [1 0 0 0; 0 0 1 0]\n"
             "  D = [0; 0]"},
        };

        models.push_back(std::move(entry));
    }

    // =====================================================================
    // 7. Double Mass-Spring-Damper
    // =====================================================================
    {
        ModelEntry entry;
        entry.name = "Double Mass-Spring";
        entry.description = "Two coupled oscillators";

        entry.params = {
            {"Mass 1", "m1", "kg", 1.0f, 0.1f, 50.0f, true},
            {"Mass 2", "m2", "kg", 1.0f, 0.1f, 50.0f, true},
            {"Spring 1", "k1", "N/m", 1.0f, 0.01f, 100.0f, true},
            {"Spring 2", "k2", "N/m", 1.0f, 0.01f, 100.0f, true},
            {"Damper 1", "b1", "N\xc2\xb7s/m", 0.1f, 0.0f, 10.0f, false},
            {"Damper 2", "b2", "N\xc2\xb7s/m", 0.1f, 0.0f, 10.0f, false},
        };

        entry.builder = [](const std::vector<PhysicalParam>& p) -> LinearSystem {
            double m1 = p[0].value, m2 = p[1].value;
            double k1 = p[2].value, k2 = p[3].value;
            double b1 = p[4].value, b2 = p[5].value;
            LinearSystem sys;
            sys.A = Eigen::MatrixXd::Zero(4, 4);
            sys.A(0, 1) = 1.0;
            sys.A(1, 0) = -(k1 + k2) / m1;
            sys.A(1, 1) = -(b1 + b2) / m1;
            sys.A(1, 2) = k2 / m1;
            sys.A(1, 3) = b2 / m1;
            sys.A(2, 3) = 1.0;
            sys.A(3, 0) = k2 / m2;
            sys.A(3, 1) = b2 / m2;
            sys.A(3, 2) = -k2 / m2;
            sys.A(3, 3) = -b2 / m2;
            sys.B = (Eigen::MatrixXd(4, 1) << 0, 1.0 / m1, 0, 0).finished();
            sys.C = Eigen::MatrixXd::Zero(2, 4);
            sys.C(0, 0) = 1.0;
            sys.C(1, 2) = 1.0;
            sys.D = Eigen::MatrixXd::Zero(2, 1);
            return sys;
        };
        entry.system = entry.builder(entry.params);

        entry.derivation = {
            {"Physical System",
             "Two masses connected by springs and dampers in series.\n"
             "Mass 1 is attached to a wall by spring k1 and damper b1,\n"
             "and coupled to mass 2 by spring k2 and damper b2.\n"
             "Force F is applied to mass 1.\n\n"
             "Models: coupled oscillators, multi-story buildings,\n"
             "drivetrain compliance, multi-DOF vibration."},

            {"Equations of Motion (ODE)",
             "Newton's second law for each mass:\n\n"
             "Mass 1:\n"
             "  m1\xc2\xb7x1'' = -k1\xc2\xb7x1 - b1\xc2\xb7x1'\n"
             "             + k2\xc2\xb7(x2 - x1) + b2\xc2\xb7(x2' - x1') + F\n\n"
             "Expand:\n"
             "  m1\xc2\xb7x1'' = -(k1+k2)\xc2\xb7x1 - (b1+b2)\xc2\xb7x1'\n"
             "             + k2\xc2\xb7x2 + b2\xc2\xb7x2' + F\n\n"
             "Mass 2:\n"
             "  m2\xc2\xb7x2'' = -k2\xc2\xb7(x2 - x1) - b2\xc2\xb7(x2' - x1')\n"
             "  m2\xc2\xb7x2'' = k2\xc2\xb7x1 + b2\xc2\xb7x1' - k2\xc2\xb7x2 - b2\xc2\xb7x2'"},

            {"Modal Analysis",
             "The system has two natural modes:\n"
             "  Mode 1 (in-phase): both masses move together\n"
             "  Mode 2 (out-of-phase): masses move opposite\n\n"
             "Natural frequencies depend on mass and stiffness ratios.\n"
             "With equal masses and springs, the modes are well-separated\n"
             "and visible as two resonance peaks in the Bode plot.\n\n"
             "Changing the mass ratio shifts peak separation.\n"
             "Increasing damping reduces peak height and broadens them."},

            {"State-Space Form",
             "States: [x1, x1', x2, x2']\n"
             "Input: F (force on mass 1)\n"
             "Outputs: [x1, x2]\n\n"
             "  A = [    0          1        0       0    ]\n"
             "      [-(k1+k2)/m1 -(b1+b2)/m1 k2/m1 b2/m1]\n"
             "      [    0          0        0       1    ]\n"
             "      [  k2/m2      b2/m2   -k2/m2 -b2/m2  ]\n\n"
             "  B = [0; 1/m1; 0; 0]\n"
             "  C = [1 0 0 0; 0 0 1 0]\n"
             "  D = [0; 0]"},
        };

        models.push_back(std::move(entry));
    }

    return models;
}

}  // namespace caliburn
