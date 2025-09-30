#pragma once
#include <vector>
#include <Eigen/Dense>
#include "centerline_map.hpp"
#include "obstacles.hpp"

namespace corridor {

struct Output {
  std::vector<double> lo, up;    // corridor bounds per step
  std::vector<double> ey_ref;    // chosen center/path per step
};

// Build raw lateral bounds [lo, up] by sweeping obstacles (disks inflated by margin)
void buildRawBounds(
  const std::vector<Eigen::Vector2d>& p_ref_h,
  const std::vector<double>&          psi_ref_h,
  const std::vector<double>&          lane_w_h,
  const std::vector<CenterlineMap::LaneRef>& lane_ref_h,
  double t0, double dt,
  const Obstacles& obstacles,
  double track_w, double margin, double L_look, // by value
  std::vector<double>& lo, std::vector<double>& up);

// Plan a smooth centerline/path within [lo, up] using a layered graph DP
Output planGraph(const std::vector<double>& s_grid,
  const std::vector<double>& lo_raw,
  const std::vector<double>& up_raw,
  double slope_per_step);

} // namespace corridor
