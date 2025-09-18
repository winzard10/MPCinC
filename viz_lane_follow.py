#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ---------- inputs ----------
sim_path  = sys.argv[1] if len(sys.argv) > 1 else "sim_log.csv"
map_path  = "data/lane_centerlines.csv"
# Prefer enhanced file if present
if Path("data/lane_centerlines_enhanced.csv").exists():
    map_path = "data/lane_centerlines_enhanced.csv"

# ---------- load sim log ----------
data = np.genfromtxt(sim_path, delimiter=",", names=True)

t      = data["t"]
x      = data["x"]
y      = data["y"]
v      = data["v"]
delta  = data["delta"]
a_cmd  = data["a_cmd"]
ddcmd  = data["ddelta_cmd"]
ey     = data["ey"]
epsi   = data["epsi"]
dv     = data["dv"]
v_ref  = data["v_ref"]
x_ref  = data["x_ref"]    # (possibly blended) reference centerline
y_ref  = data["y_ref"]

print(f"Loaded: {sim_path}")
print(f"Duration: {t[-1]:.2f}s, samples: {len(t)}")
print(f"RMS ey: {np.sqrt(np.mean(ey**2)):.4f} m,  RMS epsi: {np.sqrt(np.mean(epsi**2)):.5f} rad")

# ---------- try to load lane map ----------
have_map = Path(map_path).exists()
lane = None
if have_map:
    lane = np.genfromtxt(map_path, delimiter=",", names=True)
    cols = set(lane.dtype.names or [])

    # Base lane center curves
    def get(name_candidates):
        for nm in name_candidates:
            if nm in cols: return lane[nm]
        raise KeyError(name_candidates[0])

    # Support both legacy and enhanced CSVs
    xr = get(["x_right", "xr"]) if have_map else None
    yr = get(["y_right", "yr"]) if have_map else None
    xl = get(["x_left",  "xl"]) if have_map else None
    yl = get(["y_left",  "yl"]) if have_map else None

    # Optional enhanced columns
    xcenter = lane["x_centerline"] if "x_centerline" in cols else 0.5*(xl + xr)
    ycenter = lane["y_centerline"] if "y_centerline" in cols else 0.5*(yl + yr)

    xlb = lane["x_left_border"]  if "x_left_border"  in cols else None
    ylb = lane["y_left_border"]  if "y_left_border"  in cols else None
    xrb = lane["x_right_border"] if "x_right_border" in cols else None
    yrb = lane["y_right_border"] if "y_right_border" in cols else None
else:
    xr = yr = xl = yl = xcenter = ycenter = xlb = ylb = xrb = yrb = None

# ---------- figure ----------
fig, axs = plt.subplots(3, 2, figsize=(12, 10))
axs = axs.ravel()

# 1) XY trajectory + borders + divider
if have_map:
    # Draw borders if available; otherwise draw lane centers and shade corridor
    if xlb is not None and xrb is not None:
        axs[0].plot(xlb, ylb, label="Left lane outer border", linewidth=3)
        axs[0].plot(xrb, yrb, label="Right lane outer border", linewidth=3)
        # Light corridor shading between lane centers (fallback if borders far apart)
        try:
            # Sort by x for clean fill (use left centerline as base)
            idxL = np.argsort(xl); idxR = np.argsort(xr)
            xL, yL = xl[idxL], yl[idxL]
            xR, yR = xr[idxR], yr[idxR]
            yR_on_L = np.interp(xL, xR, yR)
            axs[0].fill_between(xL, yL, yR_on_L, alpha=0.06, zorder=0)
        except Exception:
            pass
    else:
        # Legacy: plot left/right lane centers as boundaries and shade
        axs[0].plot(xl, yl, label="Left lane center", linewidth=1)
        axs[0].plot(xr, yr, label="Right lane center", linewidth=1)
        try:
            idxL = np.argsort(xl); idxR = np.argsort(xr)
            xL, yL = xl[idxL], yl[idxL]
            xR, yR = xr[idxR], yr[idxR]
            yR_on_L = np.interp(xL, xR, yR)
            axs[0].fill_between(xL, yL, yR_on_L, alpha=0.06, zorder=0)
        except Exception:
            pass

    # Lane divider / centerline between lanes
    axs[0].plot(xcenter, ycenter, linestyle='--', linewidth=1.5, label="Lane divider (centerline)" )

# Vehicle reference & actual
axs[0].plot(x_ref, y_ref, linewidth=2, label="Reference path")
axs[0].plot(x, y, "--", label="Vehicle path")
axs[0].set_title("Trajectory (XY)")
axs[0].set_xlabel("X [m]"); axs[0].set_ylabel("Y [m]")
axs[0].legend(); axs[0].grid(True)

# 2) Speed tracking
axs[1].plot(t, v_ref, label="v_ref")
axs[1].plot(t, v, "--", label="v")
axs[1].set_title("Speed Tracking")
axs[1].set_xlabel("Time [s]"); axs[1].set_ylabel("Speed [m/s]")
axs[1].legend(); axs[1].grid(True)

# 3) Errors
axs[2].plot(t, ey, label="e_y [m]")
axs[2].plot(t, epsi, "--", label="e_psi [rad]")
axs[2].set_title("Tracking Errors")
axs[2].set_xlabel("Time [s]"); axs[2].set_ylabel("Error")
axs[2].legend(); axs[2].grid(True)

# 4) Acceleration
axs[3].plot(t, a_cmd)
axs[3].set_title("Acceleration Command")
axs[3].set_xlabel("Time [s]"); axs[3].set_ylabel("a [m/s²]")
axs[3].grid(True)

# 5) Steering rate
axs[4].plot(t, ddcmd)
axs[4].set_title("Steering Rate Command")
axs[4].set_xlabel("Time [s]"); axs[4].set_ylabel("d(delta)/dt [rad/s]")
axs[4].grid(True)

# 6) Steering angle
axs[5].plot(t, delta)
axs[5].set_title("Steering Angle")
axs[5].set_xlabel("Time [s]"); axs[5].set_ylabel("delta [rad]")
axs[5].grid(True)

plt.tight_layout()
plt.show()
