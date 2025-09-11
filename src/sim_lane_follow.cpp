#include "centerline_map.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <filesystem>
#include "mpc_ltv.hpp"


namespace fs = std::filesystem;

struct VehicleParams {
    double L{2.7};     // wheelbase [m]
    double dt{0.1};    // step [s]
};

struct Limits {
    double delta_max{0.5};     // ~28.6 deg
    double ddelta_max{0.7};    // rad/s
    double a_min{-6.0}, a_max{2.5};
};

struct State {
    double s{0.0};     // along-road (we’ll record s_proj from project())
    double x{0.0}, y{0.0}, psi{0.0};
    double v{25.0};
    double delta{0.0};
};

struct Control {
    double a{0.0};       // accel [m/s^2]
    double ddelta{0.0};  // steering rate [rad/s]
};

static inline double clamp(double u, double lo, double hi) {
    return std::min(std::max(u, lo), hi);
}

// Kinematic bicycle (pose integrated for logging)
static State step_vehicle(const State& s, const Control& u, const VehicleParams& vp, const Limits& lim)
{
    State n = s;
    const double dt = vp.dt;

    const double a = clamp(u.a, lim.a_min, lim.a_max);
    const double ddelta = clamp(u.ddelta, -lim.ddelta_max, lim.ddelta_max);

    n.v     = std::max(0.0, s.v + a * dt);
    n.s     = s.s + n.v * dt;

    n.delta = clamp(s.delta + ddelta * dt, -lim.delta_max, lim.delta_max);

    n.x   = s.x + s.v * std::cos(s.psi) * dt;
    n.y   = s.y + s.v * std::sin(s.psi) * dt;
    n.psi = s.psi + (s.v / vp.L) * std::tan(s.delta) * dt;

    return n;
}

int main(int argc, char** argv)
{
    // ----- Load map -----
    fs::path csv_path = (argc > 1) ? fs::path(argv[1]) : fs::path("data/lane_centerlines.csv");
    CenterlineMap map;
    if (!map.load_csv(csv_path.string())) {
        std::cerr << "Failed to load centerlines: " << csv_path << "\n";
        return 1;
    }
    std::cout << "Loaded map: " << csv_path << "  (s ∈ [" << map.s_min() << ", " << map.s_max() << "])\n";

    VehicleParams vp; Limits lim;
    State st;

    // initialize on RIGHT lane at start
    st.s = map.s_min();
    auto initR = map.right_lane_at(st.s);
    st.x = initR.x; st.y = initR.y; st.psi = initR.psi;
    st.v = std::min(28.0, initR.v_ref); st.delta = 0.0;

    std::ofstream log("sim_log.csv");
    log << std::fixed << std::setprecision(6);
    log << "t,s,x,y,psi,v,delta,a_cmd,ddelta_cmd,ey,epsi,dv,"
           "v_ref,x_ref,y_ref,psi_ref\n";

    const double T = 60.0;
    const int    N = static_cast<int>(T / vp.dt);

    MPCParams mpcp;
    mpcp.N = 20; mpcp.dt = vp.dt; mpcp.L = vp.L;
    LTV_MPC mpc(mpcp);

    for (int k = 0; k <= N; ++k) {
        double t = k * vp.dt;

        // --- Frenet projection onto RIGHT lane ---
        auto proj = map.project(st.x, st.y, CenterlineMap::LaneRef::Right);
        auto cref = map.center_at(proj.s_proj);  // tangent & v_ref

        // Errors for logging (and for building MPC state)
        double ey   = proj.ey;
        double epsi = st.psi - proj.psi_ref;  // actual - reference (wrap is optional for small errors)
        while (epsi >  M_PI) epsi -= 2*M_PI;
        while (epsi <= -M_PI) epsi += 2*M_PI;
        double dv   = st.v - cref.v_ref;

        // ---- Build MPC preview (kappa, v_ref) along horizon ----
        MPCRef pref;
        pref.hp.resize(mpcp.N);
        double s = proj.s_proj;
        for (int i = 0; i < mpcp.N; ++i) {
            // simple forward preview using current speed; you can use cref.v_ref if you prefer
            s = std::min(s + std::max(1e-3, st.v) * mpcp.dt, map.s_max());
            auto c = map.center_at(s);
            pref.hp[i] = { c.kappa, c.v_ref };
        }

        // ---- Current state in Frenet error coordinates ----
        MPCState xk{ ey, epsi, st.v, st.delta };

        // ---- Solve MPC ----
        MPCControl u_mpc = mpc.solve(xk, pref);

        // Fallback if solver failed
        if (!u_mpc.ok) { u_mpc.a = 0.0; u_mpc.ddelta = 0.0; }

        // ---- Log (use s_proj so the CSV s column is the along-path coordinate)
        log << t << "," << proj.s_proj << "," << st.x << "," << st.y << ","
            << st.psi << "," << st.v << "," << st.delta << ","
            << u_mpc.a << "," << u_mpc.ddelta << ","
            << ey << "," << epsi << "," << dv << ","
            << cref.v_ref << "," << proj.x_ref << "," << proj.y_ref << "," << proj.psi_ref << "\n";

        // ---- Advance vehicle dynamics (step_vehicle expects Control)
        Control u_cmd;
        u_cmd.a = u_mpc.a;
        u_cmd.ddelta = u_mpc.ddelta;
        st = step_vehicle(st, u_cmd, vp, lim);

        // optional: stop at end of map
        if (proj.s_proj > map.s_max() - 1.0) break;

    }
    log.close();

    std::cout << "Wrote sim_log.csv (Frenet errors). Plot again.\n";
    return 0;
}
