#include "obstacles.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <limits>

#ifndef OSQP_INFTY
#define OSQP_INFTY 1e20
#endif

static inline std::string lower(std::string s){
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return char(std::tolower(c)); });
    return s;
}

bool Obstacles::load_csv(const std::string& path) {
    items.clear();
    std::ifstream fin(path);
    if(!fin) return false;

    std::string line;
    // try to read header if present
    if (!std::getline(fin, line)) return false;

    while(std::getline(fin, line)){
        if(line.empty()) continue;
        std::stringstream ss(line);
        std::string cell; std::vector<std::string> cols;
        while(std::getline(ss, cell, ',')) cols.push_back(cell);
        if(cols.size() < 9) continue;

        Obstacle ob;
        try {
            ob.id = std::stoi(cols[0]);
            const auto kind = lower(cols[1]);
            ob.kind = (kind=="moving") ? Obstacle::Kind::Moving : Obstacle::Kind::Static;
            ob.x0 = std::stod(cols[2]); ob.y0 = std::stod(cols[3]);
            ob.vx = std::stod(cols[4]); ob.vy = std::stod(cols[5]);
            ob.radius  = std::stod(cols[6]);
            ob.t_start = std::stod(cols[7]);
            ob.t_end   = std::stod(cols[8]);
        } catch (...) {
            continue;
        }
        items.push_back(ob);
    }
    return true;
}

std::optional<std::pair<double,double>>
Obstacles::position_of(const Obstacle& ob, double t) const {
    if (t < ob.t_start || t > ob.t_end) return std::nullopt;
    if (ob.kind == Obstacle::Kind::Moving) {
        const double dt = t - ob.t_start;
        return std::make_pair(ob.x0 + ob.vx*dt, ob.y0 + ob.vy*dt);
    }
    return std::make_pair(ob.x0, ob.y0);
}

std::vector<Obstacles::Active> Obstacles::active_at(double t) const {
    std::vector<Active> out;
    out.reserve(items.size());
    for (const auto& ob : items) {
        auto p = position_of(ob, t);
        if (!p) continue;
        out.push_back(Active{ob.id, p->first, p->second, ob.radius});
    }
    return out;
}

// ---------------------------------------------------------
// MPC-facing utility: half-spaces → per-step ey bounds
// ---------------------------------------------------------
void compute_lateral_bounds(const MPCObsSet& set,
                            int N,
                            double ey_cap,
                            std::vector<double>& ey_upper,
                            std::vector<double>& ey_lower) {
    if (N < 0) N = 0;
    const int K = std::min<int>(N, static_cast<int>(set.obs.size()));

    ey_upper.assign(N, +OSQP_INFTY);
    ey_lower.assign(N, -OSQP_INFTY);

    for (int k = 0; k < K; ++k) {
        for (const auto& ineq : set.obs[k]) {
            const double a = ineq.a, b = ineq.b;
            if (a > 0.0) { // ey >= b/a
                ey_lower[k] = std::max(ey_lower[k], b / a);
            } else if (a < 0.0) { // ey <= b/a
                ey_upper[k] = std::min(ey_upper[k], b / a);
            }
        }
        ey_upper[k] = std::min(ey_upper[k], +ey_cap);
        ey_lower[k] = std::max(ey_lower[k], -ey_cap);
        if (ey_lower[k] > ey_upper[k]) {
            const double mid = 0.5 * (ey_lower[k] + ey_upper[k]);
            ey_lower[k] = mid; ey_upper[k] = mid;
        }
    }
}
