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
    psi_.clear(); kappa_.clear(); vref_.clear(); xc_.clear(); yc_.clear();
    xl_border_.clear(); yl_border_.clear(); xr_border_.clear(); yr_border_.clear();

    std::ifstream fin(path);
    if (!fin) return false;

    std::string line;
    if (!std::getline(fin, line)) return false; // header

    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string cell; std::vector<double> v;
        while (std::getline(ss, cell, ',')) v.push_back(cell.empty()? 0.0 : std::stod(cell));
        if (v.size() < 14) continue;

        s_.push_back(v[0]); xr_.push_back(v[1]); yr_.push_back(v[2]);
        xl_.push_back(v[3]); yl_.push_back(v[4]);
        psi_.push_back(v[5]); kappa_.push_back(v[6]); vref_.push_back(v[7]);
        xc_.push_back(v[8]); yc_.push_back(v[9]);
        xl_border_.push_back(v[10]); yl_border_.push_back(v[11]);
        xr_border_.push_back(v[12]); yr_border_.push_back(v[13]);
    }
    if (s_.size() < 2) return false;

    // s must be strictly increasing
    for (std::size_t i = 1; i < s_.size(); ++i)
        if (!(s_[i] > s_[i-1])) return false;

    return true;
}

CenterlineMap::Row CenterlineMap::row(std::size_t i) const {
    return {s_[i], xr_[i], yr_[i], xl_[i], yl_[i], psi_[i], kappa_[i], vref_[i], xc_[i], yc_[i],
            xl_border_[i], yl_border_[i], xr_border_[i], yr_border_[i]};
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
    auto r = sample(s);
    return {r.x_center, r.y_center, r.psi_center, r.kappa_center, r.v_ref};
}

// pick which left, right, center curve to project onto
using XYRef = std::pair<const std::vector<double>&, const std::vector<double>&>;

static inline XYRef pick_xy_curve(
    const std::vector<double>& xr, const std::vector<double>& xl, const std::vector<double>& xc,
    const std::vector<double>& yr, const std::vector<double>& yl, const std::vector<double>& yc,
    CenterlineMap::LaneRef which)
{
    switch (which) {
        case CenterlineMap::LaneRef::Right: return {xr, yr};
        case CenterlineMap::LaneRef::Left:  return {xl, yl};
        default:                            return {xc, yc};
    }
}


CenterlineMap::ProjectResult CenterlineMap::project(double x, double y, LaneRef which) const {
    ProjectResult out;

    const auto& curve  = pick_xy_curve(xr_, xl_, xc_, yr_, yl_, yc_, which);
    const auto& xcurve = curve.first;
    const auto& ycurve = curve.second;

    const std::size_t n = s_.size();
    if (n < 2) {
        // Fallback: no segments to project onto; snap to closest vertex if present.
        out.s_proj  = n ? s_.front() : 0.0;
        out.x_ref   = n ? xcurve.front() : x;
        out.y_ref   = n ? ycurve.front() : y;
        out.psi_ref = center_at(out.s_proj).psi;
        double nx = -std::sin(out.psi_ref), ny = std::cos(out.psi_ref);
        out.ey = (x - out.x_ref)*nx + (y - out.y_ref)*ny;
        return out;
    }

    // initial window around x (using index j into s_)
    std::size_t j  = upper_index(x);                       // assumes this is consistent with s_/curve ordering
    std::size_t i0 = (j > 50) ? (j - 50) : 0;              // include segment [0,1]
    std::size_t i1 = std::min(n - 1, j + 50);              // last vertex index; last segment is [i1-1, i1]

    double best_d2 = 1e300;
    double best_s  = s_.front();
    double best_x  = xcurve.front();
    double best_y  = ycurve.front();

    for (std::size_t i = i0; i < i1; ++i) {
        // segment endpoints (x from selected curve, y from selected curve)
        const double x0 = xcurve[i],     y0 = ycurve[i];
        const double x1 = xcurve[i + 1], y1 = ycurve[i + 1];

        const double dx = x1 - x0, dy = y1 - y0;
        const double seg2 = dx*dx + dy*dy;
        if (seg2 <= 1e-12) continue; // skip degenerate segment

        // projection parameter clamped to [0,1]
        double t = ((x - x0)*dx + (y - y0)*dy) / seg2;
        t = clamp(t, 0.0, 1.0);

        const double px = x0 + t*dx;
        const double py = y0 + t*dy;
        const double d2 = (x - px)*(x - px) + (y - py)*(y - py);

        if (d2 < best_d2) {
            best_d2 = d2;
            best_s  = s_[i] + t*(s_[i + 1] - s_[i]);
            best_x  = px;
            best_y  = py;
        }
    }

    // tangent & v_ref from road midline (shared smooth tangent)
    const auto c = center_at(best_s);

    out.s_proj  = best_s;
    out.x_ref   = best_x;
    out.y_ref   = best_y;
    out.psi_ref = c.psi;

    // signed lateral error (+ left of tangent)
    const double nx = -std::sin(out.psi_ref);
    const double ny =  std::cos(out.psi_ref);
    out.ey = (x - out.x_ref)*nx + (y - out.y_ref)*ny;

    return out;
}
