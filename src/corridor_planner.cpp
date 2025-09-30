#include "corridor_planner.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace corridor {

// Build raw lateral bounds [lo, up] by sweeping obstacles (disks inflated by margin)
void buildRawBounds(
    const std::vector<Eigen::Vector2d>& p_ref_h,
    const std::vector<double>&          psi_ref_h,
    const std::vector<double>&          lane_w_h,
    const std::vector<CenterlineMap::LaneRef>& lane_ref_h,
    double t0, double dt,
    const Obstacles& obstacles,
    double track_w,
    double margin, double L_look,
    std::vector<double>& lo, std::vector<double>& up) {

  const int N = static_cast<int>(p_ref_h.size());

  for (int k = 0; k < N; ++k) {
    const double tk  = t0 + (k + 1) * dt;
    const double psi = psi_ref_h[k];
    const double tx  = std::cos(psi), ty = std::sin(psi);


    const double lo_k = lo[k];
    const double up_k = up[k];

    for (const auto& obs : obstacles.active_at(tk)) {
      const double dx = obs.x - p_ref_h[k].x();
      const double dy = obs.y - p_ref_h[k].y();
      const double ds = dx * tx + dy * ty;               // along-lane
      if (ds < -L_look/3.0 || ds > L_look) continue;

      double ey_obs = obstacles.items[obs.idx].ey_obs;

      if (lane_ref_h[k] == CenterlineMap::LaneRef::Right) {
        if (ey_obs <= 0.0) lo[k] = std::max(lo[k], lo_k + lane_w_h[k]/2.0 + (ey_obs + obs.radius + track_w/2.0 + margin));
        else {
          if (ey_obs >= lane_w_h[k]/2.0) {up[k] = std::min(up[k], lo_k+ lane_w_h[k]/2.0 +(ey_obs-obs.radius-track_w/2.0 - margin));}
          else if (lane_w_h[k]/2.0 > -(ey_obs - obs.radius - track_w - margin)) {up[k] = std::min(up[k], up_k - 3.0*lane_w_h[k]/2.0 + (ey_obs - obs.radius - track_w/2.0 - margin));}
          else {lo[k] = std::max(lo[k], lo_k + lane_w_h[k]/2.0 + (ey_obs + obs.radius + track_w/2.0 + margin));}
        }
            
      } else if (lane_ref_h[k] == CenterlineMap::LaneRef::Left) {
          if (ey_obs >= 0.0)  up[k] = std::min(up[k], up_k - lane_w_h[k]/2 + (ey_obs - obs.radius - track_w/2.0 - margin));
          else { 
            if (ey_obs <= -lane_w_h[k]/2.0) {lo[k] = std::max(lo[k], up_k - lane_w_h[k]/2.0 + (ey_obs + obs.radius + track_w/2.0 + margin));  }
            else if (lane_w_h[k]/2.0 > ey_obs + obs.radius + track_w + margin) {lo[k] = std::max(lo[k], lo_k + 3.0*lane_w_h[k]/2.0 + (ey_obs + obs.radius + track_w/2.0 + margin)); }
            else { up[k] = std::min(up[k], up_k - lane_w_h[k]/2.0 + (ey_obs - obs.radius - track_w/2.0 - margin)); }
          }
      } else { 
          if (ey_obs >= 0.0) up[k] = std::min(up[k], up_k + (ey_obs - obs.radius - track_w/2.0 - margin));
          else               lo[k] = std::max(lo[k], lo_k + (ey_obs + obs.radius + track_w/2.0 + margin));
      }
    }

    if (lo[k] > up[k]) {  // safety clamp
      const double mid = 0.5 * (lo[k] + up[k]);
      lo[k] = up[k] = mid;
    }
  }
}

// simple temporal smoothing (slew-rate limit on bounds)
static inline void smoothBounds(std::vector<double>& lo,
                                std::vector<double>& up,
                                double slope_per_step) {
  const int N = static_cast<int>(lo.size());
  for (int k = 1; k < N; ++k) {
    lo[k] = std::max(lo[k], lo[k - 1] - slope_per_step);
    up[k] = std::min(up[k], up[k - 1] + slope_per_step);
    if (lo[k] > up[k]) { const double m = 0.5 * (lo[k] + up[k]); lo[k] = up[k] = m; }
  }
  for (int k = N - 2; k >= 0; --k) {
    lo[k] = std::max(lo[k], lo[k + 1] - slope_per_step);
    up[k] = std::min(up[k], up[k + 1] + slope_per_step);
    if (lo[k] > up[k]) { const double m = 0.5 * (lo[k] + up[k]); lo[k] = up[k] = m; }
  }
}

// check straight segment (i, ey_i) → (j, ey_j) stays inside [lo,up]
static inline bool segmentInside(const std::vector<double>& lo,
                                 const std::vector<double>& up,
                                 int i, int j, double ey_i, double ey_j) {
  if (j <= i) return false;
  const int len = j - i;
  for (int m = i; m <= j; ++m) {
    const double tau = static_cast<double>(m - i) / static_cast<double>(len);
    const double ey  = (1.0 - tau) * ey_i + tau * ey_j;
    if (ey < lo[m] - 1e-6 || ey > up[m] + 1e-6) return false;
  }
  return true;
}

Output planGraph(
    const std::vector<double>& s_grid,     // centerline s per step (monotone)
    const std::vector<double>& lo_raw,
    const std::vector<double>& up_raw,
    double slope_per_step) {

  const int N = static_cast<int>(s_grid.size());

  // 1) smooth the bounds
  std::vector<double> lo = lo_raw, up = up_raw;
  smoothBounds(lo, up, slope_per_step);

  // 2) layered graph over {lo, mid, up}
  auto mid = [&](int k) { return 0.5 * (lo[k] + up[k]); };
  constexpr int C = 3; // candidates: 0=lo, 1=mid, 2=up

  auto eyOf = [&](int k, int c) -> double {
    return (c == 0) ? lo[k] : ((c == 1) ? mid(k) : up[k]);
  };

  struct Node { double cost; int prev_k; int prev_c; };
  std::vector<std::array<Node, C>> best(N);
  for (int c = 0; c < C; ++c) best[0][c] = {0.0, -1, -1};

  for (int k = 0; k < N - 1; ++k) {
    for (int c = 0; c < C; ++c) {
      const double ey_k = eyOf(k, c);
      const double base = best[k][c].cost;
      if (!std::isfinite(base)) continue;

      for (int kp = k + 1; kp < N; ++kp) {
        const double ds = std::max(1e-3, s_grid[kp] - s_grid[k]);
        for (int cp = 0; cp < C; ++cp) {
          const double ey_kp = eyOf(kp, cp);
          if (!segmentInside(lo, up, k, kp, ey_k, ey_kp)) continue;

          const double dtheta = std::atan2(ey_kp - ey_k, ds);  // heading change proxy
          const double cost   = base + std::abs(dtheta);

          if (best[kp][cp].prev_k < 0 || cost < best[kp][cp].cost) {
            best[kp][cp] = {cost, k, c};
          }
        }
      }
    }
  }

  // 3) pick best terminal & backtrack ey_ref
  int kf = N - 1, cf = 0;
  for (int c = 1; c < C; ++c)
    if (best[kf][c].prev_k >= 0 && best[kf][c].cost < best[kf][cf].cost) cf = c;

  std::vector<double> ey_ref(N);
  int k = kf, c = cf;
  while (k >= 0) {
    ey_ref[k] = eyOf(k, c);
    const int pk = best[k][c].prev_k, pc = best[k][c].prev_c;
    if (pk < 0) break;
    k = pk; c = pc;
  }
  for (int i = 0; i < k; ++i) ey_ref[i] = ey_ref[k];  // fill leading holes (safety)

  return {std::move(lo), std::move(up), std::move(ey_ref)};
}

} // namespace corridor

