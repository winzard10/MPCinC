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

static CenterlineMap::LaneRef parse_lane(const std::string& s) {
    std::string t = s;
    std::transform(t.begin(), t.end(), t.begin(), ::tolower);
    if (t == "right") return CenterlineMap::LaneRef::Right;
    if (t == "left")  return CenterlineMap::LaneRef::Left;
    return CenterlineMap::LaneRef::Center;
}

static double smoothstep01(double u) {
    // 3u^2 - 2u^3, clamped to [0,1]
    if (u <= 0.0) return 0.0;
    if (u >= 1.0) return 1.0;
    return u*u*(3.0 - 2.0*u);
}

// sample y of a lane at s
static double lane_y_at(const CenterlineMap& map, double s, CenterlineMap::LaneRef which) {
    if (which == CenterlineMap::LaneRef::Right) return map.right_lane_at(s).y;
    if (which == CenterlineMap::LaneRef::Left)  return map.left_lane_at(s).y;
    return map.center_at(s).y;
}


int main(int argc, char** argv)
{
    // ---- CLI options ----
    CenterlineMap::LaneRef lane_from = CenterlineMap::LaneRef::Right;
    CenterlineMap::LaneRef lane_to   = CenterlineMap::LaneRef::Right;
    double t_change = 1e18;   // default: no lane change
    double T_change = 3.0;
    std::string map_file = "data/lane_centerlines.csv";

    for (int i=1; i<argc; ++i) {
        std::string a = argv[i];
        if (a == "--lane-from" && i+1<argc)      lane_from = parse_lane(argv[++i]);
        else if (a == "--lane-to" && i+1<argc)   lane_to   = parse_lane(argv[++i]);
        else if (a == "--t-change" && i+1<argc)  t_change  = std::stod(argv[++i]);
        else if (a == "--T-change" && i+1<argc)  T_change  = std::stod(argv[++i]);
        else if (a == "--map" && i+1<argc)       map_file  = argv[++i];
    }

    // ---- Load map (declare 'map' before using it) ----
    CenterlineMap map;
    if (!map.load_csv(map_file)) {
        std::cerr << "Failed to load centerlines: " << map_file << "\n";
        return 1;
    }
    std::cout << "Loaded map: " << map_file
            << "  (s ∈ [" << map.s_min() << ", " << map.s_max() << "])\n";


    VehicleParams vp; Limits lim;
    State st;

    // initialize on RIGHT lane at start
    st.s = map.s_min();
    auto initR = map.right_lane_at(st.s);
    st.x = initR.x; st.y = initR.y; st.psi = initR.psi;
    st.v = std::min(28.0, initR.v_ref); st.delta = 0.0;

    std::ofstream log("sim_log.csv");
    log << std::fixed << std::setprecision(6);
    log << "t,s,x,y,psi,v,delta,a_cmd,ddelta_cmd,ey,epsi,dv,v_ref,x_ref,y_ref,psi_ref,alpha\n";

    const double T = 60.0;
    const int    N = static_cast<int>(T / vp.dt);

    MPCParams mpcp;
    mpcp.N = 20; mpcp.dt = vp.dt; mpcp.L = vp.L;
    LTV_MPC mpc(mpcp);

    for (int k = 0; k <= N; ++k) {
        double t = k * vp.dt;

        // --- Project to centerline first to get s, tangent, psi_ref ---
        auto projC = map.project(st.x, st.y);  // centerline projection (no lane arg)
        auto cref  = map.center_at(projC.s_proj);

        // Smooth blend factor alpha(t) for the lane change
        double alpha = smoothstep01( (t - t_change) / T_change );

        // Blend y reference between lanes at the same s
        double y_from = lane_y_at(map, projC.s_proj, lane_from);
        double y_to   = lane_y_at(map, projC.s_proj, lane_to);
        double y_ref  = (1.0 - alpha) * y_from + alpha * y_to;

        // We keep x_ref on centerline projection for numerical stability
        double x_ref  = projC.x_ref;
        double psi_ref = cref.psi;

        // Signed lateral error to blended ref using centerline normal
        double nx = -std::sin(psi_ref);
        double ny =  std::cos(psi_ref);
        double ey  = (st.x - x_ref) * nx + (st.y - y_ref) * ny;

        // Heading error relative to road tangent
        double epsi = st.psi - psi_ref;
        while (epsi >  M_PI) epsi -= 2*M_PI;
        while (epsi <= -M_PI) epsi += 2*M_PI;

        // Speed error for logging
        double dv = st.v - cref.v_ref;

        // ---- Build MPC preview (kappa, v_ref) along horizon ----
        MPCRef pref;
        pref.hp.resize(mpcp.N);
        double s_h = projC.s_proj;
        for (int i = 0; i < mpcp.N; ++i) {
            s_h = std::min(s_h + std::max(1e-3, st.v) * mpcp.dt, map.s_max());
            auto c = map.center_at(s_h);
            pref.hp[i].kappa = c.kappa;
            pref.hp[i].v_ref = c.v_ref;
        }

        // ---- Current state in Frenet error coordinates ----
        MPCState xk{ ey, epsi, st.v, st.delta };

        // ---- Solve MPC ----
        MPCControl u_mpc = mpc.solve(xk, pref);

        // Fallback if solver failed
        if (!u_mpc.ok) { u_mpc.a = 0.0; u_mpc.ddelta = 0.0; }

        // ---- Log (use s_proj so the CSV s column is the along-path coordinate)
        log << t << "," << projC.s_proj << "," << st.x << "," << st.y << ","
            << st.psi << "," << st.v << "," << st.delta << ","
            << u_mpc.a << "," << u_mpc.ddelta << ","
            << ey << "," << epsi << "," << dv << ","
            << cref.v_ref << "," << x_ref << "," << y_ref << "," << psi_ref << ","
            << alpha << "\n";   // NEW

        // ---- Advance vehicle dynamics (step_vehicle expects Control)
        Control u_cmd;
        u_cmd.a = u_mpc.a;
        u_cmd.ddelta = u_mpc.ddelta;
        st = step_vehicle(st, u_cmd, vp, lim);

        // optional: stop at end of map
        if (projC.s_proj > map.s_max() - 1.0) break;

    }
    log.close();

    std::cout << "Wrote sim_log.csv (Frenet errors). Plot again.\n";
    return 0;
}
