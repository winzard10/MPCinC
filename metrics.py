import pandas as pd
import numpy as np
from pathlib import Path

def summarize(sim_csv: str = "sim_log.csv", save_csv: str | None = "logs/summary_metrics.csv"):
    df = pd.read_csv(sim_csv)
    # Expected columns from your sim:
    # t,s,x,y,psi,v,delta,a_cmd,ddelta_cmd,ey,epsi,dv,v_ref,x_ref,y_ref,psi_ref
    
    metrics = {}
    def pct95(x): return np.percentile(np.abs(x), 95)
    def meanabs(x): return np.mean(np.abs(x))
    def rmse(err): return np.sqrt(np.mean(np.square(err)))
    
    metrics["|ey|_mean"]   = meanabs(df["ey"])
    metrics["|ey|_p95"]    = pct95(df["ey"])
    metrics["|epsi|_mean"] = meanabs(df["epsi"])
    metrics["|epsi|_p95"]  = pct95(df["epsi"])
    # speed tracking error
    if "v_ref" in df.columns:
        verr = df["v"] - df["v_ref"]
        metrics["speed_RMSE"] = rmse(verr)
        metrics["speed_mean_err"] = verr.mean()
    # control effort
    if "a_cmd" in df.columns:
        metrics["|a|_mean"] = meanabs(df["a_cmd"])
        metrics["|a|_p95"]  = pct95(df["a_cmd"])
    if "ddelta_cmd" in df.columns:
        metrics["|ddelta|_mean"] = meanabs(df["ddelta_cmd"])
        metrics["|ddelta|_p95"]  = pct95(df["ddelta_cmd"])
    # Steering angle usage
    metrics["|delta|_mean"] = meanabs(df["delta"])
    metrics["|delta|_p95"]  = pct95(df["delta"])
    
    out = pd.DataFrame([metrics])
    if save_csv:
        Path(save_csv).parent.mkdir(parents=True, exist_ok=True)
        # append or create
        if Path(save_csv).exists():
            prev = pd.read_csv(save_csv)
            out_all = pd.concat([prev, out], ignore_index=True)
            out_all.to_csv(save_csv, index=False)
        else:
            out.to_csv(save_csv, index=False)
    return out

def lane_change_metrics(sim_csv="sim_log.csv", map_csv="data/lane_centerlines.csv",
                        margin=0.2, epsilon=0.05, L=2.7):
    df = pd.read_csv(sim_csv)
    # require alpha column (0→1 blend)
    alpha = df.get("alpha", None)
    if alpha is None:
        raise ValueError("alpha column not found; run the lane-change sim that logs alpha.")

    # window of the change: where alpha in (0,1)
    mask = (alpha > 1e-6) & (alpha < 1-1e-6)
    if mask.sum() == 0:
        return {"note": "no lane change detected (alpha constant)"}

    t = df["t"].values
    w = np.where(mask)[0]
    t_start, t_end = t[w[0]], t[w[-1]]

    # KPIs in the change window
    ey = df["ey"].values
    epsi = df["epsi"].values
    peak_ey = float(np.max(np.abs(ey[w])))
    peak_epsi = float(np.max(np.abs(epsi[w])))

    # settling after change: next 3 s
    post = (t >= t_end) & (t <= t_end + 3.0)
    settled = np.all(np.abs(ey[post]) < epsilon) if post.any() else False
    settle_time = float(t_end - t[w[0]])  # duration of the blend itself

    # comfort proxies
    v = df["v"].values
    ddelta = df["ddelta_cmd"].values
    delta = df["delta"].values
    ay = (v**2 / L) * delta
    peak_ay = float(np.max(np.abs(ay[w])))
    jerk_proxy = (v**2 / L) * np.abs(ddelta)
    peak_jerk = float(np.max(jerk_proxy[w]))

    # constraint usage
    # (optional: write a_max, delta_max into CSV once so we can read them)
    frac_sat_delta = float(np.mean(np.isclose(np.abs(delta[w]), np.max(np.abs(delta[w])), atol=1e-6)))

    # corridor violations: use map borders if available
    viol_count = 0
    try:
        lane = pd.read_csv(map_csv)
        # Interp left/right y on reference x_ref (from sim log)
        x_ref = df["x_ref"].values
        xl, yl = lane["xl"].values, lane["yl"].values
        xr, yr = lane["xr"].values, lane["yr"].values
        # sort by x for interpolation
        idxL = np.argsort(xl); idxR = np.argsort(xr)
        yL_on_ref = np.interp(x_ref, xl[idxL], yl[idxL])
        yR_on_ref = np.interp(x_ref, xr[idxR], yr[idxR])
        y = df["y"].values
        # count outside band ± margin
        low = np.minimum(yL_on_ref, yR_on_ref) - margin
        high = np.maximum(yL_on_ref, yR_on_ref) + margin
        viol_count = int(np.sum((y < low) | (y > high)))
    except Exception:
        pass

    return {
        "t_change_start": float(t_start),
        "t_change_end": float(t_end),
        "change_duration": float(t_end - t_start),
        "peak_ey": peak_ey,
        "peak_epsi": peak_epsi,
        "settled_within_3s": bool(settled),
        "peak_ay": peak_ay,
        "peak_jerk_proxy": peak_jerk,
        "violations_count": viol_count,
        "frac_sat_delta": frac_sat_delta,
    }

if __name__ == "__main__":
    summary = summarize()
    print(summary.to_string(index=False))