# Highway MPC Simulation

This project implements a **lane-following highway controller** using **Linear Time-Varying MPC (LTV-MPC)** in C++, with visualization in Python.  
The workflow is: **build → run simulations → visualize results**.

---

## Quick Commands

```bash
rm -rf build
cmake -S . -B build
cmake --build build -j

# Run simulations
./build/sim_lane_follow
./build/centerline_demo data/lane_centerlines.csv

# Visualize results
python3 viz_lane_follow.py
