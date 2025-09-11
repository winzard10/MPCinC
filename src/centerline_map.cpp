#include "centerline_map.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>

static inline double clamp(double x, double a, double b) {
    return std::min(std::max(x, a), b);
}

bool CenterlineMap::load_csv(const std::string& path) {
    s_.clear(); xr_.clear(); yr_.clear(); xl_.clear(); yl_.clear();
    psi_.clear(); kappa_.clear(); vref_.clear(); yc_.clear();

    std::ifstream fin(path);
    if (!fin) return false;

    std::string line;
    if (!std::getline(fin, line)) return false; // header

    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string cell; std::vector<double> v;
        while (std::getline(ss, cell, ',')) v.push_back(cell.empty()? 0.0 : std::stod(cell));
        if (v.size() < 8) continue;

        s_.push_back(v[0]); xr_.push_back(v[1]); yr_.push_back(v[2]);
        xl_.push_back(v[3]); yl_.push_back(v[4]);
        psi_.push_back(v[5]); kappa_.push_back(v[6]); vref_.push_back(v[7]);
    }
    if (s_.size() < 2) return false;

    // road midline y = (yr + yl)/2
    yc_.resize(yr_.size());
    for (std::size_t i = 0; i < yr_.size(); ++i) yc_[i] = 0.5*(yr_[i] + yl_[i]);

    // s must be strictly increasing
    for (std::size_t i = 1; i < s_.size(); ++i)
        if (!(s_[i] > s_[i-1])) return false;

    return true;
}

CenterlineMap::Row CenterlineMap::row(std::size_t i) const {
    return {s_[i], xr_[i], yr_[i], xl_[i], yl_[i], psi_[i], kappa_[i], vref_[i]};
}

std::size_t CenterlineMap::upper_index(double s) const {
    if (s <= s_.front()) return 1;
    if (s >= s_.back())  return s_.size() - 1;
    auto it = std::upper_bound(s_.begin(), s_.end(), s);
    return std::distance(s_.begin(), it);
}

CenterlineMap::Row CenterlineMap::sample(double s) const {
    std::size_t j = upper_index(s), i = j - 1;
    double s0 = s_[i], s1 = s_[j];
    double t = (clamp(s, s0, s1) - s0) / (s1 - s0);
    auto L = [t](double a, double b){ return a + t*(b-a); };
    return {clamp(s, s_.front(), s_.back()),
            L(xr_[i],xr_[j]), L(yr_[i],yr_[j]),
            L(xl_[i],xl_[j]), L(yl_[i],yl_[j]),
            L(psi_[i],psi_[j]), L(kappa_[i],kappa_[j]), L(vref_[i],vref_[j])};
}

CenterlineMap::LanePose CenterlineMap::right_lane_at(double s) const {
    auto r = sample(s);
    return {r.x_right, r.y_right, r.psi_center, r.kappa_center, r.v_ref};
}

CenterlineMap::LanePose CenterlineMap::left_lane_at(double s) const {
    auto r = sample(s);
    return {r.x_left, r.y_left, r.psi_center, r.kappa_center, r.v_ref};
}

CenterlineMap::CenterPose CenterlineMap::center_at(double s) const {
    std::size_t j = upper_index(s), i = j - 1;
    double s0 = s_[i], s1 = s_[j];
    double t = (clamp(s, s0, s1) - s0) / (s1 - s0);
    auto L = [t](double a, double b){ return a + t*(b-a); };

    CenterPose c;
    c.x     = L(s_[i], s_[j]);      // x ≈ s in this dataset
    c.y     = L(yc_[i], yc_[j]);    // midline y
    c.psi   = L(psi_[i], psi_[j]);
    c.kappa = L(kappa_[i], kappa_[j]);
    c.v_ref = L(vref_[i], vref_[j]);
    return c;
}

// pick which y-curve to project onto
static inline const std::vector<double>& pick_y_curve(
    const std::vector<double>& yr,
    const std::vector<double>& yl,
    const std::vector<double>& yc,
    CenterlineMap::LaneRef which)
{
    switch (which) {
        case CenterlineMap::LaneRef::Right:  return yr;  // right lane centerline
        case CenterlineMap::LaneRef::Left:   return yl;  // left lane centerline
        default:                             return yc;  // road midline
    }
}

// Project (x,y) to nearest segment on the chosen reference curve in a local window
CenterlineMap::ProjectResult CenterlineMap::project(double x, double y, LaneRef which) const {
    ProjectResult out;

    const auto& ycurve = pick_y_curve(yr_, yl_, yc_, which);

    // initial window around x (~s)
    std::size_t j  = upper_index(x);
    std::size_t n  = s_.size();
    std::size_t i0 = (j > 50) ? j - 50 : 1;
    std::size_t i1 = std::min(n - 1, j + 50);

    double best_d2 = 1e300;
    double best_s  = s_.front();
    double best_x  = s_.front();      // x ≈ s
    double best_y  = ycurve.front();

    for (std::size_t i = i0; i < i1; ++i) {
        // segment endpoints (x is s; y from selected curve)
        double x0 = s_[i],   y0 = ycurve[i];
        double x1 = s_[i+1], y1 = ycurve[i+1];

        double dx = x1 - x0, dy = y1 - y0;
        double seg2 = dx*dx + dy*dy;
        if (seg2 <= 1e-12) continue;

        // projection parameter in [0,1]
        double t = ((x - x0)*dx + (y - y0)*dy) / seg2;
        t = clamp(t, 0.0, 1.0);

        double px = x0 + t*dx;
        double py = y0 + t*dy;
        double d2 = (x - px)*(x - px) + (y - py)*(y - py);

        if (d2 < best_d2) {
            best_d2 = d2;
            best_s  = s_[i] + t*(s_[i+1] - s_[i]);
            best_x  = px;
            best_y  = py;
        }
    }

    // tangent & v_ref from road midline (smooth shared tangent)
    auto c = center_at(best_s);

    out.s_proj  = best_s;
    out.x_ref   = best_x;
    out.y_ref   = best_y;
    out.psi_ref = c.psi;

    // signed lateral error (+ left of tangent)
    double nx = -std::sin(out.psi_ref);
    double ny =  std::cos(out.psi_ref);
    out.ey = (x - out.x_ref)*nx + (y - out.y_ref)*ny;

    return out;
}
