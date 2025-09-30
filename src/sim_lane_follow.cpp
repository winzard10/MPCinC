#include "centerline_map.hpp"
#include "mpc_ltv.hpp"
#include "obstacles.hpp"
#include "corridor_planner.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

// ---------- Basic types ----------
struct VehicleParams {
  double L{3.0};   // wheelbase [m]
  double track_w{1.7}; // track width [m]
  double dt{0.1};  // step [s]
};

struct Limits {
  double delta_max{0.5};   // ~28.6 deg
  double ddelta_max{0.7};  // rad/s
  double a_min{-6.0}, a_max{2.5};
};

struct State {
  double s{0.0};   // along-road (we’ll log s_proj from project())
  double x{0.0}, y{0.0}, psi{0.0};
  double v{25.0};
  double delta{0.0};
};

struct Control {
  double a{0.0};       // accel [m/s^2]
  double ddelta{0.0};  // steering rate [rad/s]
};

struct CLI {
  CenterlineMap::LaneRef lane_from{CenterlineMap::LaneRef::Right};
  CenterlineMap::LaneRef lane_to  {CenterlineMap::LaneRef::Right};
  double t_change{1e18};     // default: no lane change
  double T_change{3.0};
  std::string map_file{"data/lane_centerlines.csv"};
  std::string obs_file_primary{"data/obstacles.csv"};
  std::string obs_file_fallback{"data/obstacles_example.csv"};
  std::string log_file{"sim_log.csv"};
  double T{180.0};
};

// ---------- Small math & utils ----------
static inline double clamp(double u, double lo, double hi) {
  return std::min(std::max(u, lo), hi);
}

static inline double smoothstep01(double u) {
  // 3u^2 - 2u^3, clamped to [0,1]
  if (u <= 0.0) return 0.0;
  if (u >= 1.0) return 1.0;
  return u * u * (3.0 - 2.0 * u);
}

static inline double wrapAngle(double a) {
  while (a >  M_PI) a -= 2.0 * M_PI;
  while (a <=-M_PI) a += 2.0 * M_PI;
  return a;
}

static inline CenterlineMap::LaneRef parseLane(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), ::tolower);
  if (s == "right") return CenterlineMap::LaneRef::Right;
  if (s == "left")  return CenterlineMap::LaneRef::Left;
  return CenterlineMap::LaneRef::Center;
}

// ---------- Vehicle model ----------
static State stepVehicle(const State& s, const Control& u,
                         const VehicleParams& vp, const Limits& lim) {
  State n = s;
  const double dt = vp.dt;

  const double a      = clamp(u.a,     lim.a_min, lim.a_max);
  const double ddelta = clamp(u.ddelta, -lim.ddelta_max, lim.ddelta_max);

  n.v     = std::max(0.0, s.v + a * dt);
  n.s     = s.s + n.v * dt;
  n.delta = clamp(s.delta + ddelta * dt, -lim.delta_max, lim.delta_max);

  // Kinematic bicycle forward Euler for pose (integrate for logging)
  n.x   = s.x + s.v * std::cos(s.psi) * dt;
  n.y   = s.y + s.v * std::sin(s.psi) * dt;
  n.psi = s.psi + (s.v / vp.L) * std::tan(s.delta) * dt;

  return n;
}

// ---------- Lane pose helpers ----------
struct LanePose {
  double x{0}, y{0}, psi{0}, kappa{0}, v_ref{0}, lane_width{3.7};
};

static inline LanePose lanePoseAt(const CenterlineMap& map, double s,
                                  CenterlineMap::LaneRef which) {
  if (which == CenterlineMap::LaneRef::Right) {
    auto r = map.right_lane_at(s); return {r.x, r.y, r.psi, r.kappa, r.v_ref, r.lane_width};
  } else if (which == CenterlineMap::LaneRef::Left) {
    auto l = map.left_lane_at(s);  return {l.x, l.y, l.psi, l.kappa, l.v_ref, l.lane_width};
  }
  auto c = map.center_at(s);       return {c.x, c.y, c.psi, c.kappa, c.v_ref, c.lane_width};
}

// Compute lateral offset (ey) of a point (x,y) from the centerline
double lateralOffsetFromCenterline(const CenterlineMap& map, double x, double y, CenterlineMap::LaneRef which) {
  // Project obstacle position onto centerline
  auto proj = map.project(x, y);    // expects .s_proj, .x_ref, .y_ref, .psi, etc.

  LanePose cref = lanePoseAt(map, proj.s_proj, which);

  // Normal vector to centerline heading (left-positive)
  const double nx = -std::sin(cref.psi);
  const double ny =  std::cos(cref.psi);

  // Lateral offset (dot product with normal)
  const double ey = (x - cref.x) * nx + (y - cref.y) * ny;
  return ey;
}

// ---------- CLI parsing ----------
static void parseCLI(int argc, char** argv, CLI& cli) {
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--lane-from"   && i + 1 < argc) cli.lane_from = parseLane(argv[++i]);
    else if (a == "--lane-to"   && i + 1 < argc) cli.lane_to   = parseLane(argv[++i]);
    else if (a == "--t-change"  && i + 1 < argc) cli.t_change  = std::stod(argv[++i]);
    else if (a == "--T-change"  && i + 1 < argc) cli.T_change  = std::stod(argv[++i]);
    else if (a == "--map"       && i + 1 < argc) cli.map_file  = argv[++i];
    else if (a == "--log"       && i + 1 < argc) cli.log_file  = argv[++i];
    else if (a == "--T"         && i + 1 < argc) cli.T         = std::stod(argv[++i]);
    else if (a == "--help" || a == "-h") {
      std::cout <<
        "Usage: sim_main [options]\n"
        "  --lane-from {left|center|right}\n"
        "  --lane-to   {left|center|right}\n"
        "  --t-change  <sec>     (default: no change)\n"
        "  --T-change  <sec>     (default: 3.0)\n"
        "  --map       <file>    (default: data/lane_centerlines.csv)\n"
        "  --log       <file>    (default: sim_log.csv)\n"
        "  --T         <sec>     (total sim time, default: 180)\n";
      std::exit(0);
    }
  }
}

// ---------- Main ----------
int main(int argc, char** argv) {
  // --- CLI ---
  CLI cli;
  parseCLI(argc, argv, cli);

  // --- Map ---
  CenterlineMap map;
  if (!map.load_csv(cli.map_file)) {
    std::cerr << "Failed to load centerlines: " << cli.map_file << "\n";
    return 1;
  }
  std::cout << "Loaded map: " << cli.map_file
            << "  (s ∈ [" << map.s_min() << ", " << map.s_max() << "])\n";

  // --- Obstacles (optional) ---
  Obstacles obstacles;
  std::string obs_file = fs::exists(cli.obs_file_primary)
                         ? cli.obs_file_primary
                         : (fs::exists(cli.obs_file_fallback) ? cli.obs_file_fallback : "");
  if (!obs_file.empty()) {
    if (obstacles.load_csv(obs_file)) {
      std::cout << "Loaded obstacles: " << obs_file
                << " (" << obstacles.items.size() << " items)\n";
    } else {
      std::cout << "Failed to parse obstacles file: " << obs_file << "\n";
    }
  } else {
    std::cout << "No obstacles file found; running without obstacles.\n";
  }

  // --- Vehicle & MPC ---
  VehicleParams vp; Limits lim; State st;
  MPCParams mpcp; mpcp.N = 20; mpcp.dt = vp.dt; mpcp.L = vp.L;
  LTV_MPC mpc(mpcp);

  // Initialize vehicle pose on chosen starting lane at s_min
  st.s = map.s_min();
  {
    LanePose p0 = lanePoseAt(map, st.s, cli.lane_from);
    st = State{st.s, p0.x, p0.y, p0.psi, std::min(28.0, p0.v_ref), 0.0};
  }

  // --- Logging ---
  std::ofstream log(cli.log_file);
  if (!log) {
    std::cerr << "Cannot open log file: " << cli.log_file << "\n";
    return 1;
  }
  log << std::fixed << std::setprecision(6);
  log << "t,s,x,y,psi,v,delta,a_cmd,ddelta_cmd,ey,epsi,dv,"
         "v_ref,x_ref,y_ref,psi_ref,alpha,dmin\n";

  // --- Simulation loop ---
  const int steps = static_cast<int>(cli.T / vp.dt);
  for (int k = 0; k <= steps; ++k) {
    const double t = k * vp.dt;

    // Project current pose to centerline for reference heading (almost not relevant anymore)
    auto projC = map.project(st.x, st.y);     // expects .s_proj and .x_ref etc.
    auto cref  = map.center_at(projC.s_proj);
    // const double psi_ref = cref.psi;

    // Lane-change blend for target XY (interpolate between lanes)
    const double alpha = smoothstep01((t - cli.t_change) / cli.T_change);
    LanePose fromP = lanePoseAt(map, projC.s_proj, cli.lane_from);
    LanePose toP   = lanePoseAt(map, projC.s_proj, cli.lane_to);

    auto blendYaw = [](double psi1, double psi2, double a) {
      // Wrap difference to (-pi, pi]
      const double d = std::atan2(std::sin(psi2 - psi1), std::cos(psi2 - psi1));
      return psi1 + a * d;
  };
  
    // Blended reference heading for the lane-change path
    const double psi_ref = blendYaw(fromP.psi, toP.psi, alpha);
    const double x_ref = (1.0 - alpha) * fromP.x + alpha * toP.x;
    const double y_ref = (1.0 - alpha) * fromP.y + alpha * toP.y;

    // Frenet errors (left-positive w.r.t centerline heading)
    const double nx = -std::sin(psi_ref);
    const double ny =  std::cos(psi_ref);
    const double ey   = (st.x - x_ref) * nx + (st.y - y_ref) * ny;
    const double epsi = wrapAngle(st.psi - psi_ref);
    const double dv   = st.v - cref.v_ref;

    // --- MPC preview (reference horizon & corridor frames) ---
    MPCRef pref; pref.hp.resize(mpcp.N);
    std::vector<Eigen::Vector2d> p_ref_h(mpcp.N);
    std::vector<double>          psi_ref_h(mpcp.N);
    std::vector<double>          v_ref_h(mpcp.N);
    std::vector<double>          s_grid(mpcp.N);
    std::vector<CenterlineMap::LaneRef> lane_ref_h(mpcp.N);
    std::vector<double>          lane_w_h(mpcp.N);

    double s_h = projC.s_proj;
    double t_h = t;

    for (int i = 0; i < mpcp.N; ++i) {
      s_grid[i] = s_h;

      CenterlineMap::LaneRef which;
      if (t_h < cli.t_change) which = CenterlineMap::LaneRef::Right;
      else if (t_h > cli.t_change + cli.T_change) which = CenterlineMap::LaneRef::Left;
      else which = CenterlineMap::LaneRef::Center;

      LanePose c = lanePoseAt(map, s_h, which);

      // MPC preview
      pref.hp[i].kappa = c.kappa;
      pref.hp[i].v_ref = c.v_ref;

      // Corridor frame + metadata
      p_ref_h[i]   = {c.x, c.y};
      psi_ref_h[i] = c.psi;
      lane_ref_h[i]= which;
      lane_w_h[i]  = std::max(0.1, c.lane_width); // guard
      
      for (const auto &a : obstacles.active_at(t)) {
        double ey = lateralOffsetFromCenterline(map, a.x, a.y, which);
        obstacles.items[a.idx].ey_obs = ey;  // write back to the real object
      }

      // advance
      const double v_step = std::max(1e-3, c.v_ref);
      s_h = std::min(s_h + v_step * mpcp.dt, map.s_max());
      t_h += mpcp.dt;
    }

    std::vector<double> lo(mpcp.N), up(mpcp.N);
    for (int i = 0; i < mpcp.N; ++i) {
      const double w = lane_w_h[i];
      switch (lane_ref_h[i]) {
        case CenterlineMap::LaneRef::Right:
          lo[i] = - (0.5) * w;
          up[i] = + (1.5) * w;
          break;
        case CenterlineMap::LaneRef::Left:
          lo[i] = - (1.5) * w;
          up[i] = + (0.5) * w;
          break;
        default: // Center
          lo[i] = - w;
          up[i] = + w;
          break;
      }
    }

    // Current state in Frenet error coordinates
    MPCState xk{ey, epsi, st.v, st.delta};

    // --- Corridor planning (Algorithm-2 & 3 inspired) ---
    const double t0     = t;
    const double margin = 0.1;   // inflate obstacle radius by margin [m]
    const double L_look = 70.0;  // look-ahead distance for obstacles [m]
    const double slope  = 0.25;  // bound slew per step [m/step]

    corridor::buildRawBounds(p_ref_h, psi_ref_h, lane_w_h, lane_ref_h, t0, mpcp.dt,
                             obstacles, vp.track_w, margin, L_look,
                             lo, up);

    corridor::Output cor = corridor::planGraph(s_grid, lo, up, slope);

    pref.ey_ref.resize(mpcp.N + 1);
    for (int i = 0; i < mpcp.N; ++i)
        pref.ey_ref[i] = cor.ey_ref[i];   // <-- use planned path
    pref.ey_ref.back() = pref.ey_ref[mpcp.N - 1];
    pref.ey_ref_N = pref.ey_ref.back();

    mpc.setCorridorBounds(cor.lo, cor.up);

    // --- Solve MPC ---
    MPCControl u_mpc = mpc.solve(xk, pref);
    if (!u_mpc.ok) { u_mpc.a = 0.0; u_mpc.ddelta = 0.0; }  // safe fallback

    // --- Measure min distance to currently active obstacles (for logging) ---
    double dmin = std::numeric_limits<double>::infinity();
    for (const auto& a : obstacles.active_at(t)) {
      const double dx = st.x - a.x;
      const double dy = st.y - a.y;
      const double d  = std::sqrt(dx * dx + dy * dy) - a.radius;
      if (d < dmin) dmin = d;
    }
    if (!std::isfinite(dmin))
      dmin = std::numeric_limits<double>::quiet_NaN();

    // --- Log ---
    log << t << "," << projC.s_proj << "," << st.x << "," << st.y << ","
        << st.psi << "," << st.v << "," << st.delta << ","
        << u_mpc.a << "," << u_mpc.ddelta << ","
        << ey << "," << epsi << "," << dv << ","
        << cref.v_ref << "," << x_ref << "," << y_ref << "," << psi_ref << ","
        << alpha << "," << dmin << "\n";

    // --- Advance vehicle dynamics ---
    const Control u_cmd{u_mpc.a, u_mpc.ddelta};
    st = stepVehicle(st, u_cmd, vp, lim);

    // Stop if we reach end of map
    if (projC.s_proj > map.s_max() - 1.0) break;
  }

  std::cout << "Log written to: " << cli.log_file << "\n";
  return 0;
}
