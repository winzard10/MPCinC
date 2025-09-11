#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

path = sys.argv[1] if len(sys.argv) > 1 else "sim_log.csv"
data = np.genfromtxt(path, delimiter=",", names=True)

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
x_ref  = data["x_ref"]
y_ref  = data["y_ref"]

# ---- quick stats ----
print("Loaded:", path)
print(f"Duration: {t[-1]:.2f}s, samples: {len(t)}")
print(f"RMS ey: {np.sqrt(np.mean(ey**2)):.3f} m,  RMS epsi: {np.sqrt(np.mean(epsi**2)):.4f} rad")

# ---- subplot layout ----
fig, axs = plt.subplots(3, 2, figsize=(12, 10))
axs = axs.ravel()

# 1) XY trajectory
axs[0].plot(x_ref, y_ref, label="Lane center", linewidth=2)
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
