#include "centerline_map.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <filesystem>
#include <Eigen/Dense>
#include "mpc_ltv.hpp"
#include "obstacles.hpp"


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
    log << "t,s,x,y,psi,v,delta,a_cmd,ddelta_cmd,ey,epsi,dv,v_ref,x_ref,y_ref,psi_ref,alpha,dmin\n";

    // ---- Load obstacles (optional) ----
    Obstacles obstacles;
    std::string obs_file = "data/obstacles.csv";
    if (!std::filesystem::exists(obs_file))
        obs_file = "data/obstacles_example.csv";
    if (std::filesystem::exists(obs_file)) {
        if (obstacles.load_csv(obs_file)) {
            std::cout << "Loaded obstacles: " << obs_file
                    << " (" << obstacles.items.size() << " items)\n";
        } else {
            std::cout << "Failed to parse obstacles file: " << obs_file << "\n";
        }
    } else {
        std::cout << "No obstacles file found; running without obstacles.\n";
    }

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

        // === BEGIN: obstacle corridor constraints (stable) ===
        MPCObsSet oc;
        oc.obs.resize(mpcp.N);

        const double t0     = t;        // current sim time
        const double margin = 0.8;      // clearance added to obstacle radius [m]
        const double L_look = 200.0;     // only obstacles 0..L_look m ahead
        const double slope  = 0.30;     // max corridor change per step [m/step]

        // Precompute ref position and lane frames along the horizon
        std::vector<Eigen::Vector2d> p_ref_h(mpcp.N);
        std::vector<double>          psi_ref_h(mpcp.N);
        std::vector<Eigen::Vector2d> n_lane_h(mpcp.N);
        std::vector<Eigen::Vector2d> t_hat_h(mpcp.N);

        double s_for_h = projC.s_proj;
        for (int k = 0; k < mpcp.N; ++k) {
            s_for_h = std::min(s_for_h + std::max(0.1, st.v) * mpcp.dt, map.s_max());
            const auto c = map.center_at(s_for_h);            // needs {x,y,psi}
            p_ref_h[k]   = Eigen::Vector2d(c.x, c.y);
            psi_ref_h[k] = c.psi;
            n_lane_h[k]  = Eigen::Vector2d(-std::sin(psi_ref_h[k]),  std::cos(psi_ref_h[k]));
            t_hat_h[k]   = Eigen::Vector2d( std::cos(psi_ref_h[k]),  std::sin(psi_ref_h[k]));
        }

        // Start with a wide default corridor (your nominal ey bounds)
        std::vector<double> ey_min(mpcp.N, -mpcp.ey_max);
        std::vector<double> ey_max(mpcp.N,  +mpcp.ey_max);

        // Tighten the corridor with obstacles
        for (int k = 0; k < mpcp.N; ++k) {
            const double tk = t0 + (k + 1) * mpcp.dt;
            auto act = obstacles.active_at(tk);
            for (const auto& a : act) {
                const Eigen::Vector2d p_ob(a.x, a.y);
                const Eigen::Vector2d d = p_ob - p_ref_h[k];
                const double x_lon = t_hat_h[k].dot(d);  // along-lane distance
                if (x_lon < 0.0 || x_lon > L_look) continue;  // ignore behind/too far

                const double y_lat = n_lane_h[k].dot(d); // +left / -right
                const double Rplus = a.radius + margin;

                if (y_lat >= 0.0) {
                    // obstacle on LEFT → upper bound
                    ey_max[k] = std::min(ey_max[k], y_lat - Rplus);
                } else {
                    // obstacle on RIGHT → lower bound
                    ey_min[k] = std::max(ey_min[k], y_lat + Rplus);
                }
            }
        }

        // Smooth the corridor to prevent step-to-step flipping (forward + backward pass)
        for (int k = 1; k < mpcp.N; ++k) {
            ey_max[k] = std::min(ey_max[k], ey_max[k-1] + slope);
            ey_min[k] = std::max(ey_min[k], ey_min[k-1] - slope);
        }
        for (int k = mpcp.N - 2; k >= 0; --k) {
            ey_max[k] = std::min(ey_max[k], ey_max[k+1] + slope);
            ey_min[k] = std::max(ey_min[k], ey_min[k+1] - slope);
        }

        // Ensure corridor never inverts; if too tight, open a tiny gap (epsilon)
        const double eps = 0.05;
        for (int k = 0; k < mpcp.N; ++k) {
            if (ey_min[k] > ey_max[k] - eps) {
                const double mid = 0.5*(ey_min[k] + ey_max[k]);
                ey_min[k] = mid - 0.5*eps;
                ey_max[k] = mid + 0.5*eps;
            }
            // Convert corridor to inequalities using existing sigma_ey(k):
            // ey >= ey_min  →  (+1)*ey + σ >=  ey_min
            // ey ≤ ey_max   →  (-1)*ey + σ ≥ -ey_max
            oc.obs[k].push_back(ObsIneq{ +1.0,  ey_min[k]  + 1.2});
            oc.obs[k].push_back(ObsIneq{ -1.0, -ey_max[k]  + 1.2});
        }

        // Install constraints for this solve
        mpc.setObstacleConstraints(std::move(oc));
        // === END: obstacle corridor constraints ===

        // ---- Solve MPC ----
        MPCControl u_mpc = mpc.solve(xk, pref);

        // Fallback if solver failed
        if (!u_mpc.ok) { u_mpc.a = 0.0; u_mpc.ddelta = 0.0; }

        // ---- Minimum distance to obstacles (at current sim time t) ----
        double dmin = std::numeric_limits<double>::infinity();
        {
            auto active = obstacles.active_at(t);  // use 't' (your current sim time)
            for (const auto& a : active) {
                double dx = st.x - a.x;            // use st.x / st.y (state), not x/y
                double dy = st.y - a.y;
                double d  = std::sqrt(dx*dx + dy*dy) - a.radius;
                if (d < dmin) dmin = d;
            }
            if (!std::isfinite(dmin)) dmin = std::numeric_limits<double>::quiet_NaN();
        }

        // ---- Log (use s_proj so the CSV 's' is the along-path coordinate)
        log << t << "," << projC.s_proj << "," << st.x << "," << st.y << ","
            << st.psi << "," << st.v << "," << st.delta << ","
            << u_mpc.a << "," << u_mpc.ddelta << ","
            << ey << "," << epsi << "," << dv << ","
            << cref.v_ref << "," << x_ref << "," << y_ref << "," << psi_ref << ","
            << alpha << "," << dmin << "\n";

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
