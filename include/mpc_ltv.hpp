#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <cstddef>

struct MPCParams {
    int    N   = 20;
    double dt  = 0.1;
    double L   = 2.7;

    // weights
    double wy    = 4.0;
    double wpsi  = 3.0;
    double wv    = 1.0;
    double wa    = 0.1;
    double wdd   = 0.05;
    double wda   = 0.0;   // slew a
    double wddd  = 0.0;   // slew ddelta
    double wyf   = 8.0;
    double wpsif = 6.0;

    // bounds
    double a_min = -3.0, a_max = 3.0;
    double ddelta_max = 0.5;        // rad/s
    double v_min = 0.0,  v_max = 40.0;
    double delta_max = 0.4;         // rad

    // NEW: soft-constraint params
    double ey_max = 1.5;            // lateral band (m)
    double w_sigma_ey = 1e4;        // heavy penalty on ey slack
    double w_sigma_delta = 1e4;     // heavy penalty on delta slack
};

struct PreviewPoint {
    double kappa;     // road curvature at s_k
    double v_ref;     // desired speed at s_k
};

struct MPCRef {
    // length N sequences
    std::vector<PreviewPoint> hp;
};

struct MPCState {
    double ey;      // lateral error [m]
    double epsi;    // heading error [rad]
    double v;       // speed [m/s]
    double delta;   // steering angle [rad]
};

struct MPCControl {
    double a;        // accel [m/s^2]
    double ddelta;   // steering rate [rad/s]
    bool   ok = false;
};

class LTV_MPC {
public:
    explicit LTV_MPC(const MPCParams& p);

    // Solve one step given current state x0 and horizon preview (kappa, v_ref)
    MPCControl solve(const MPCState& x0, const MPCRef& ref);

    // warm-start hooks (optional)
    void resetWarmStart();

private:
    MPCParams P;

    // Nominal rollout for linearization (kept shallow for simplicity)
    struct Nominal {
        std::vector<MPCState> x;   // size N+1
        std::vector<Eigen::Vector2d> u; // [a, ddelta], size N
    } nom_;

    // Build linearized discrete model x_{k+1} = A x_k + B u_k + c
    void buildLinearization(const MPCRef& ref);

    // Assemble QP matrices (H, g, Aeq, beq, Ain, lbA, ubA, lb, ub)
    // and call OSQP; returns first control.
    MPCControl solveQP(const MPCState& x0, const MPCRef& ref);

    // helpers
    static void angleWrap(double& a);
};
