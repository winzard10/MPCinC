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

// ---------- Corridor graph planner (Algorithm 2 & 3 inspired) ----------
struct CorrOut {
    std::vector<double> lo, up;    // corridor bounds per step
    std::vector<double> ey_ref;    // chosen center path (same length as horizon)
  };
  
  // make raw corridor [lo, up] by sweeping obstacles (disk inflated by margin)
  static inline void build_raw_corridor(
      const std::vector<Eigen::Vector2d>& p_ref_h,
      const std::vector<double>&          psi_ref_h,
      double t0, double dt,
      const Obstacles& obstacles,
      double margin, double L_look, double ey_cap,
      std::vector<double>& lo, std::vector<double>& up)
  {
    const int N = (int)p_ref_h.size();
    lo.assign(N, -ey_cap);
    up.assign(N, +ey_cap);
  
    for (int k=0; k<N; ++k) {
      const double tk = t0 + (k+1)*dt;
      const double psi = psi_ref_h[k];
      const double tx = std::cos(psi), ty = std::sin(psi);
      const double nx = -ty,           ny =  ty;
  
      for (const auto& a : obstacles.active_at(tk)) {
        const double dx = a.x - p_ref_h[k].x();
        const double dy = a.y - p_ref_h[k].y();
        const double ds = dx*tx + dy*ty;                   // along-lane
        if (ds < -2.0 || ds > L_look) continue;
  
        const double ey = dx*nx + dy*ny;                   // lateral
        const double inflate = a.radius + margin;
        if (ey >= 0.0)   up[k] = std::min(up[k], ey - inflate);  // tighten left
        else             lo[k] = std::max(lo[k], ey + inflate);  // tighten right
      }
      if (lo[k] > up[k]) { const double m=0.5*(lo[k]+up[k]); lo[k]=up[k]=m; }
    }
  }
  
  // simple temporal smoothing of bounds (slew-rate limit)
  static inline void smooth_bounds(std::vector<double>& lo,
                                   std::vector<double>& up,
                                   double slope_per_step)
  {
    const int N = (int)lo.size();
    for (int k=1; k<N; ++k) {
      lo[k] = std::max(lo[k], lo[k-1]-slope_per_step);
      up[k] = std::min(up[k], up[k-1]+slope_per_step);
      if (lo[k] > up[k]) { const double m=0.5*(lo[k]+up[k]); lo[k]=up[k]=m; }
    }
    for (int k=N-2; k>=0; --k) {
      lo[k] = std::max(lo[k], lo[k+1]-slope_per_step);
      up[k] = std::min(up[k], up[k+1]+slope_per_step);
      if (lo[k] > up[k]) { const double m=0.5*(lo[k]+up[k]); lo[k]=up[k]=m; }
    }
  }
  
  // check that the straight segment between (i, ey_i) and (j, ey_j)
  // stays inside [lo,up] for all integer m in [i..j]
  static inline bool segment_stays_inside(const std::vector<double>& lo,
                                          const std::vector<double>& up,
                                          int i, int j, double ey_i, double ey_j)
  {
    if (j <= i) return false;
    const int len = j - i;
    for (int m=i; m<=j; ++m) {
      const double tau = (double)(m - i) / (double)len;
      const double ey  = (1.0 - tau)*ey_i + tau*ey_j;
      if (ey < lo[m] - 1e-6 || ey > up[m] + 1e-6) return false;
    }
    return true;
  }
  
  // graph search with heading-change cost; nodes = {lo, mid, up} at each step
  static inline CorrOut plan_corridor_graph(
      const std::vector<double>& s_grid,     // centerline s for each step (monotone)
      const std::vector<double>& lo_raw,
      const std::vector<double>& up_raw,
      double slope_per_step)
  {
    const int N = (int)s_grid.size();
    std::vector<double> lo = lo_raw, up = up_raw;
    smooth_bounds(lo, up, slope_per_step);
  
    // candidate ey values per step
    auto mid = [&](int k){ return 0.5*(lo[k]+up[k]); };
    const int C = 3; // {lo, mid, up}
  
    // Dijkstra on layered graph (k,c) → (k’,c’), k’>k
    struct Node { double cost; int prev_k, prev_c; };
    std::vector<std::array<Node,3>> best(N);
    for (int c=0;c<C;++c) best[0][c] = {0.0, -1, -1};
  
    // precompute s-spacing
    auto ey_of = [&](int k,int c)->double{
      return c==0? lo[k] : (c==1? mid(k) : up[k]);
    };
  
    for (int k=0; k<N-1; ++k) {
      for (int c=0; c<C; ++c) {
        const double ey_k = ey_of(k,c);
        const double base = best[k][c].cost;
        if (!std::isfinite(base)) continue;
  
        // try to jump to any future step k’ in [k+1 .. N-1]
        for (int kp=k+1; kp<N; ++kp) {
          const double ds = std::max(1e-3, s_grid[kp] - s_grid[k]);
          for (int cp=0; cp<C; ++cp) {
            const double ey_kp = ey_of(kp,cp);
            if (!segment_stays_inside(lo, up, k, kp, ey_k, ey_kp)) continue;
  
            const double dtheta = std::atan2(ey_kp - ey_k, ds); // ≈ eq (10) with c=0
            const double cost   = base + std::abs(dtheta);
  
            if (best[kp][cp].prev_k < 0 || cost < best[kp][cp].cost) {
              best[kp][cp] = {cost, k, c};
            }
          }
        }
      }
    }
  
    // pick best terminal
    int kf = N-1, cf = 0;
    for (int c=1; c<C; ++c)
      if (best[kf][c].prev_k >= 0 && best[kf][c].cost < best[kf][cf].cost) cf = c;
  
    // backtrack ey_ref
    std::vector<double> ey_ref(N);
    int k = kf, c = cf;
    while (k >= 0) {
      ey_ref[k] = ey_of(k,c);
      const int pk = best[k][c].prev_k, pc = best[k][c].prev_c;
      if (pk < 0) break;
      k = pk; c = pc;
    }
    // fill any leading holes (shouldn’t happen, but safe)
    for (int i=0; i<k; ++i) ey_ref[i] = ey_ref[k];
  
    return { std::move(lo), std::move(up), std::move(ey_ref) };
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

    const double T = 180.0;
    const int    N = static_cast<int>(T / vp.dt);

    MPCParams mpcp;
    mpcp.N = 20; mpcp.dt = vp.dt; mpcp.L = vp.L;
    LTV_MPC mpc(mpcp);

    for (int k = 0; k <= N; ++k) {
    const double t = k * vp.dt;

    // --- Project current pose to centerline & reference at s ---
    auto projC  = map.project(st.x, st.y);      // expects .s_proj and .x_ref
    auto cref   = map.center_at(projC.s_proj);
    const double psi_ref = cref.psi;

    // Optional lane-change blend (kept for your CLI options)
    const double alpha  = smoothstep01((t - t_change)/T_change);
    const double y_from = lane_y_at(map, projC.s_proj, lane_from);
    const double y_to   = lane_y_at(map, projC.s_proj, lane_to);
    const double y_ref  = (1.0 - alpha) * y_from + alpha * y_to;

    // Frenet errors (left-positive)
    const double nx = -std::sin(psi_ref);
    const double ny =  std::cos(psi_ref);
    const double ey    = (st.x - projC.x_ref) * nx + (st.y - y_ref) * ny;
    double epsi = st.psi - psi_ref;
    while (epsi >  M_PI) epsi -= 2*M_PI;
    while (epsi <=-M_PI) epsi += 2*M_PI;

    const double dv = st.v - cref.v_ref;

    // ---- Build MPC preview and frames in ONE pass (no uninitialized vectors) ----
    MPCRef pref; pref.hp.resize(mpcp.N);
    std::vector<Eigen::Vector2d> p_ref_h(mpcp.N);
    std::vector<double>          psi_ref_h(mpcp.N);
    std::vector<double>          s_grid(mpcp.N);

    double s_h = projC.s_proj;
    for (int i = 0; i < mpcp.N; ++i) {
        s_grid[i] = s_h;
        const auto c = map.center_at(s_h);

        // preview signals for MPC
        pref.hp[i].kappa = c.kappa;
        pref.hp[i].v_ref = c.v_ref;

        // frames for corridor
        p_ref_h[i]   = Eigen::Vector2d(c.x, c.y);
        psi_ref_h[i] = c.psi;

        // advance station using v_ref (stable preview spacing)
        const double v_step = std::max(1e-3, c.v_ref);
        s_h = std::min(s_h + v_step * mpcp.dt, map.s_max());
    }

    // ---- Current state in Frenet error coordinates ----
    MPCState xk{ ey, epsi, st.v, st.delta };

    // === BEGIN: Corridor planning (graph-based) ===
    const double t0     = t;
    const double margin = 1.2;   // inflate obstacle radius by margin [m]
    const double L_look = 70.0;  // consider obstacles up to this far ahead [m]
    const double slope  = 0.25;  // corridor slew per step [m/step]

    // 1) raw bounds from obstacles (Algorithm-2 tightened road bounds)
    std::vector<double> lo_raw, up_raw;
    build_raw_corridor(p_ref_h, psi_ref_h, t0, mpcp.dt,
                       obstacles, margin, L_look, mpcp.ey_max,
                       lo_raw, up_raw);

    // 2) corridor graph (Algorithm-3 style) with heading-change cost
    CorrOut cor = plan_corridor_graph(s_grid, lo_raw, up_raw, slope);

    // center of corridor works well; if you have a graph path, use that instead
    auto mid = [&](int i){ return 0.5*(cor.lo[i] + cor.up[i]); };

    pref.ey_ref.resize(mpcp.N + 1);
    for (int i = 0; i < mpcp.N; ++i) pref.ey_ref[i] = mid(i);
    pref.ey_ref[mpcp.N] = pref.ey_ref[mpcp.N-1];  // terminal
    pref.ey_ref_N = pref.ey_ref.back();

    // 3) install as hard inequalities ey ∈ [lo, up]
    MPCObsSet oc; oc.obs.assign(mpcp.N, {});
    for (int i = 0; i < mpcp.N; ++i) {
        oc.obs[i].push_back({ +1.0,  cor.lo[i]  });  // ey >= lo
        oc.obs[i].push_back({ -1.0, -cor.up[i]  });  // ey <= up
    }
    mpc.setObstacleConstraints(std::move(oc));
    // === END: Corridor planning ===

    // ---- Solve MPC ----
    MPCControl u_mpc = mpc.solve(xk, pref);
    if (!u_mpc.ok) { u_mpc.a = 0.0; u_mpc.ddelta = 0.0; }  // safe fallback

    // ---- Minimum distance to currently active obstacles (for logging) ----
    double dmin = std::numeric_limits<double>::infinity();
    for (const auto& a : obstacles.active_at(t)) {
        const double dx = st.x - a.x;
        const double dy = st.y - a.y;
        const double d  = std::sqrt(dx*dx + dy*dy) - a.radius;
        if (d < dmin) dmin = d;
    }
    if (!std::isfinite(dmin))
        dmin = std::numeric_limits<double>::quiet_NaN();

    // ---- Log one row per step ----
    log << t << "," << projC.s_proj << "," << st.x << "," << st.y << ","
        << st.psi << "," << st.v << "," << st.delta << ","
        << u_mpc.a << "," << u_mpc.ddelta << ","
        << ey << "," << epsi << "," << dv << ","
        << cref.v_ref << "," << projC.x_ref << "," << y_ref << "," << psi_ref << ","
        << alpha << "," << dmin << "\n";

    // ---- Advance vehicle dynamics ----
    Control u_cmd{ u_mpc.a, u_mpc.ddelta };
    st = step_vehicle(st, u_cmd, vp, lim);

    // stop if we reach end of map
    if (projC.s_proj > map.s_max() - 1.0) break;
}

}