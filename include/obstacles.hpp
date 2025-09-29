#pragma once
#include <vector>
#include <string>
#include <optional>

// -----------------------------
// Obstacle CSV model & queries
// -----------------------------
struct Obstacle {
    int    id = -1;
    enum class Kind { Static, Moving } kind = Kind::Static;
    double x0 = 0.0, y0 = 0.0;          // initial position
    double vx = 0.0, vy = 0.0;          // velocity if moving [m/s]
    double radius = 0.5;                // safety radius (already includes vehicle inflation)
    double t_start = 0.0, t_end = 1e9;  // active time window [s]
};

struct Obstacles {
    // Load from CSV with header:
    // id,kind,x0,y0,vx,vy,radius,t_start,t_end   (kind ∈ {static,moving})
    bool load_csv(const std::string& path);

    // Position at time t (for moving obstacles).
    // Returns std::nullopt if t is outside [t_start, t_end].
    std::optional<std::pair<double,double>> position_of(const Obstacle& ob, double t) const;

    // Get list of active obstacles at time t with their positions.
    struct Active { int id; double x; double y; double radius; };
    std::vector<Active> active_at(double t) const;

    std::vector<Obstacle> items;
};

// ---------------------------------------------
// MPC-facing half-spaces and ey bound utility
// ---------------------------------------------
// Each inequality encodes a half-space in lateral error: a * ey >= b.
struct ObsIneq { double a; double b; };

struct MPCObsSet {
    // obs[k] = list of half-spaces to apply at step k (k = 0..N-1)
    std::vector<std::vector<ObsIneq>> obs;
};
