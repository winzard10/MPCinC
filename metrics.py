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

if __name__ == "__main__":
    summary = summarize()
    print(summary.to_string(index=False))