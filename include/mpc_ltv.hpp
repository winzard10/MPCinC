#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <cstddef>

struct MPCParams {
    int    N = 20;        // horizon steps
    double dt = 0.1;      // [s]
    double L  = 2.7;      // wheelbase [m]

    // Weights
    double wy = 8.0, wpsi = 4.0, wv = 0.5;             // state tracking
    double wa = 0.2, wdd = 1.0;                        // input effort: a, ddelta
    double wda = 0.1, wddd = 0.5;                      // input slew (Δu)
    double wyf = 20.0, wpsif = 10.0;                   // terminal

    // Limits
    double a_min = -6.0, a_max = 2.5;                  // [m/s^2]
    double ddelta_max = 0.7;                           // [rad/s]
    double delta_max  = 0.5;                           // [rad]
    double v_min = 0.0, v_max = 40.0;                  // [m/s]
    double ey_max = 1.2;                               // soft lane bound (unused in hard form below)
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
