#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <cstddef>
#include <optional>

#include "obstacles.hpp"   // MPCObsSet + compute_lateral_bounds

// ---------------------------------
// MPC configuration (kept as-is)
// ---------------------------------
struct MPCParams {
    int    N   = 20;
    double dt  = 0.1;
    double L   = 2.7;

    // weights
    double wy    = 0.20;
    double wpsi  = 0.02;
    double wv    = 0.10;
    double wa    = 0.05;
    double wdd   = 0.10;   // effort on ddelta
    double wda   = 0.00;   // slew a     (u_k - u_{k-1})
    double wddd  = 0.00;   // slew ddelta
    double wyf   = 8.0;    // terminal ey
    double wpsif = 6.0;    // terminal epsi

    // bounds
    double a_min = -5.0, a_max = 3.0;
    double ddelta_max = 0.25;     // |steer rate| [rad/s]
    double delta_max  = 0.40;     // steering angle cap [rad]
    double v_min = 0.0,  v_max = 40.0;

    // lateral hard band (also used when no obstacles)
    double ey_max = 10.0;
};

// Preview point used by your sim
struct PreviewPoint {
    double kappa = 0.0;   // road curvature at s_k
    double v_ref = 0.0;   // desired speed at s_k
};

// Horizon preview (your sim fills ref.hp[i].{kappa,v_ref})
struct MPCRef {
    std::vector<PreviewPoint> hp;  // length N (or >= N)
    std::vector<double> ey_ref;      // size N+1 (state at k=0..N)
    double ey_ref_N{0.0};            // terminal, if you prefer separate
};

// State / Input / Output
struct MPCState {
    double ey    = 0.0;   // lateral error [m]
    double epsi  = 0.0;   // heading error [rad]
    double v     = 0.0;   // speed [m/s]
    double delta = 0.0;   // steering angle [rad]
};

struct MPCControl {
    double a      = 0.0;  // accel [m/s^2]
    double ddelta = 0.0;  // steering rate [rad/s]
    bool   ok     = false;
};

// ---------------------------------
// LTV MPC Controller
// ---------------------------------
class LTV_MPC {
public:
    explicit LTV_MPC(const MPCParams& p): P(p) {
        nom_.x.resize(P.N + 1);
        nom_.u.resize(P.N);
    }

    void setObstacleConstraints(std::optional<MPCObsSet> s) { obs_ = std::move(s); }

    // Optional warm-start
    void setNominal(const std::vector<MPCState>& x_nom,
                    const std::vector<Eigen::Vector2d>& u_nom) {
        nom_.x = x_nom; nom_.u = u_nom;
    }

    // Build linearized discrete model x_{k+1} = A_k x_k + B_k u_k + c_k
    void buildLinearization(const MPCRef& ref);

    // Assemble and solve QP; returns first input to apply
    MPCControl solveQP(const MPCState& x0, const MPCRef& ref);

    // Backward-compatible wrapper (your sim calls mpc.solve(...))
    MPCControl solve(const MPCState& x0, const MPCRef& ref) {
        buildLinearization(ref);
        return solveQP(x0, ref);
    }

    // helpers
    static void angleWrap(double& a);

    void setCorridorBounds(const std::vector<double>& lo,
        const std::vector<double>& up);

private:
    MPCParams P;

    struct LinModel {
        std::vector<Eigen::Matrix<double,4,4>> A;  // ey, epsi, v, delta
        std::vector<Eigen::Matrix<double,4,2>> B;  // a, ddelta
        std::vector<Eigen::Matrix<double,4,1>> c;  // affine
    } lm_;

    struct Nominal {
        std::vector<MPCState>        x;   // N+1
        std::vector<Eigen::Vector2d> u;   // N
    } nom_;

    std::vector<double> ey_lo_provided_;
    std::vector<double> ey_up_provided_;
    bool have_corridor_bounds_{false};

    std::optional<MPCObsSet> obs_;        // per-step ey half-spaces
};
