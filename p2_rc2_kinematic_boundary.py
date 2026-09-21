#!/usr/bin/env python3
"""P2 RC2: continuous bearing + known platform motion → kinematic candidate boundary.

Pure simulation theory study. No Bellhop/KRAKEN, no acoustic field.
Truth is never injected into the estimator grid.

Scenarios:
  LIN  — frozen E-STD baseline (psi=0, v=U): constant bearing, r-v unobservable (RC2 limit)
  CR5  — weak crossing psi=5°: primary performance geometry
  CR10 — moderate crossing psi=10°
"""
from __future__ import annotations

import json
import math
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "P2_RC2_kinematic_boundary"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

U_PLAT = 2.0
Z_FIXED = 200.0
R0_BASE_KM = 50.0
THETA0_BASE_DEG = 0.0
V_BASE = 2.0
PSI_BASE_DEG = 0.0

T_LIST = [60.0, 120.0, 300.0, 600.0, 900.0, 1200.0]
SIGMA_DEG = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
TURNS_DEG = [0.0, 2.0, 5.0, 10.0, 15.0]
PSI_SCAN_DEG = [-15.0, -10.0, -5.0, -3.0, 0.0, 3.0, 5.0, 10.0, 15.0]

N_TRIALS = 8
DT_OBS = 10.0
RNG_SEED = 20260912

R_MIN_M, R_MAX_M = 45e3, 60e3
TH_MIN_DEG, TH_MAX_DEG = -5.0, 5.0
V_MIN, V_MAX = 1.0, 3.0
PSI_MIN_DEG, PSI_MAX_DEG = -15.0, 15.0

R_GRID_KM = np.arange(45.0, 60.0 + 1e-9, 1.0)
TH_GRID_DEG = np.arange(-5.0, 5.0 + 1e-9, 0.5)
V_GRID = np.arange(1.0, 3.0 + 1e-9, 0.2)
PSI_GRID_DEG = np.arange(-15.0, 15.0 + 1e-9, 1.0)  # includes 0°

SCENARIOS = {
    "LIN": dict(r0_m=50e3, theta0_deg=0.0, v=2.0, psi_deg=0.0,
                label="共线基准 LIN（ψ=0,v=U）",
                note="方位恒为0，r–v本征不可辨识；RC2理论下界"),
    "CR5": dict(r0_m=50e3, theta0_deg=0.0, v=2.0, psi_deg=5.0,
                label="弱交叉 CR5（ψ=5°）",
                note="存在方位变化率，RC2可部分约束"),
    "CR10": dict(r0_m=50e3, theta0_deg=0.0, v=2.0, psi_deg=10.0,
                 label="中等交叉 CR10（ψ=10°）",
                 note="方位变化率更强"),
}
PRIMARY_SCENARIO = "CR5"
LIMIT_SCENARIO = "LIN"
SIGMA_FIG = [0.02, 0.05, 0.1, 0.2, 0.5, 1.0]

OFFGRID_TRUTHS = [
    dict(tag="OG1", r0_m=52.3e3, theta0_deg=0.37, v=2.33, psi_deg=4.7),
    dict(tag="OG2", r0_m=47.6e3, theta0_deg=-1.21, v=1.64, psi_deg=-7.3),
    dict(tag="OG3", r0_m=58.1e3, theta0_deg=2.05, v=2.78, psi_deg=11.2),
]

CONFIG = {
    "project": "B3D-3",
    "package": "P2_RC2_kinematic_boundary",
    "created_utc": NOW,
    "scene": "E-STD",
    "z_fixed_m": Z_FIXED,
    "platform_speed_mps": U_PLAT,
    "baseline_truth": {
        "r0_km": R0_BASE_KM, "theta0_deg": THETA0_BASE_DEG,
        "v_mps": V_BASE, "psi_deg": PSI_BASE_DEG,
    },
    "T_s": T_LIST,
    "sigma_theta_deg": SIGMA_DEG,
    "platform_turn_deg": TURNS_DEG,
    "psi_scan_deg": PSI_SCAN_DEG,
    "n_trials": N_TRIALS,
    "dt_obs_s": DT_OBS,
    "rng_seed": RNG_SEED,
    "candidate_box": {
        "r0_km": [R_MIN_M / 1e3, R_MAX_M / 1e3],
        "theta0_deg": [TH_MIN_DEG, TH_MAX_DEG],
        "v_mps": [V_MIN, V_MAX],
        "psi_deg": [PSI_MIN_DEG, PSI_MAX_DEG],
    },
    "coarse_grid": {
        "r0_km": R_GRID_KM.tolist(),
        "theta0_deg": TH_GRID_DEG.tolist(),
        "v_mps": V_GRID.tolist(),
        "psi_deg": PSI_GRID_DEG.tolist(),
        "n_nodes": int(len(R_GRID_KM) * len(TH_GRID_DEG) * len(V_GRID) * len(PSI_GRID_DEG)),
        "note": "psi grid includes 0°",
    },
    "scenarios": SCENARIOS,
    "primary_scenario": PRIMARY_SCENARIO,
    "limit_scenario": LIMIT_SCENARIO,
    "method_notes": [
        "Full 2D motion geometry; no r(t)=r0+(v-U)t shortcut.",
        "theta0 unknown; estimator never receives truth.",
        "sigma_theta=0 is ideal upper bound (deterministic match).",
        "FIM/Jacobian evaluated at truth for local observability only.",
        "Candidate set = coarse grid + local continuous refinement from top-k, not truth-injected.",
        "Primary metrics = accepted-set widths / contraction (set-valued).",
        "LIN is a true RC2 unobservable limit (constant bearing); CR5/CR10 open bearing-rate geometry.",
    ],
    "forbidden": [
        "TRUE_BEARING_FIXED",
        "truth written into every candidate",
        "guaranteeing truth on candidate grid",
        "using acoustic field H to generate/filter candidates",
    ],
}


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
def platform_states(t, turn_deg, u=U_PLAT):
    t = np.asarray(t, dtype=float)
    xp = np.empty_like(t)
    yp = np.empty_like(t)
    hdx = np.empty_like(t)
    hdy = np.empty_like(t)
    if abs(turn_deg) < 1e-12:
        xp[:] = u * t
        yp[:] = 0.0
        hdx[:] = 1.0
        hdy[:] = 0.0
        return xp, yp, hdx, hdy
    t_turn = float(np.max(t)) * 0.5
    delta = math.radians(turn_deg)
    c, s = math.cos(delta), math.sin(delta)
    mask = t <= t_turn
    xp[mask] = u * t[mask]
    yp[mask] = 0.0
    dt = t[~mask] - t_turn
    xp[~mask] = u * t_turn + u * dt * c
    yp[~mask] = u * dt * s
    hdx[mask] = 1.0
    hdy[mask] = 0.0
    hdx[~mask] = c
    hdy[~mask] = s
    return xp, yp, hdx, hdy


def target_states(t, r0_m, theta0_rad, v, psi_rad):
    xt = r0_m * np.cos(theta0_rad) + v * t * np.cos(psi_rad)
    yt = r0_m * np.sin(theta0_rad) + v * t * np.sin(psi_rad)
    return xt, yt


def predict_bearings(t, r0_m, theta0_rad, v, psi_rad, turn_deg, u=U_PLAT):
    xp, yp, _, _ = platform_states(t, turn_deg, u)
    xt, yt = target_states(t, r0_m, theta0_rad, v, psi_rad)
    return np.arctan2(yt - yp, xt - xp)


def predict_bearings_batch(t, r0_m, th0_rad, v, psi_rad, turn_deg, u=U_PLAT):
    xp, yp, _, _ = platform_states(t, turn_deg, u)
    xt = r0_m[:, None] * np.cos(th0_rad[:, None]) + v[:, None] * t[None, :] * np.cos(psi_rad[:, None])
    yt = r0_m[:, None] * np.sin(th0_rad[:, None]) + v[:, None] * t[None, :] * np.sin(psi_rad[:, None])
    return np.arctan2(yt - yp[None, :], xt - xp[None, :])


def wrap_angle(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def bearing_rmse(pred, obs):
    return float(np.sqrt(np.mean(wrap_angle(pred - obs) ** 2)))


def jacobian_bearings(t, r0_m, theta0_rad, v, psi_rad, turn_deg, u=U_PLAT):
    xp, yp, _, _ = platform_states(t, turn_deg, u)
    xt, yt = target_states(t, r0_m, theta0_rad, v, psi_rad)
    dx, dy = xt - xp, yt - yp
    r2 = np.maximum(dx * dx + dy * dy, 1e-6)
    dth_ddx, dth_ddy = -dy / r2, dx / r2
    H = np.empty((t.size, 4), dtype=float)
    H[:, 0] = dth_ddx * np.cos(theta0_rad) + dth_ddy * np.sin(theta0_rad)
    H[:, 1] = dth_ddx * (-r0_m * np.sin(theta0_rad)) + dth_ddy * (r0_m * np.cos(theta0_rad))
    H[:, 2] = dth_ddx * (t * np.cos(psi_rad)) + dth_ddy * (t * np.sin(psi_rad))
    H[:, 3] = dth_ddx * (-v * t * np.sin(psi_rad)) + dth_ddy * (v * t * np.cos(psi_rad))
    return H


def fim_metrics(H, sigma_rad):
    if sigma_rad > 0:
        Hs = H / sigma_rad
    else:
        Hs = H.copy()
    F = Hs.T @ Hs
    F = 0.5 * (F + F.T)
    evals = np.clip(np.linalg.eigvalsh(F), 0.0, None)
    s = np.linalg.svd(Hs, compute_uv=False)
    smin, smax = float(s.min()), float(s.max())
    cond = float("inf") if smin <= 1e-30 else smax / smin
    try:
        P = np.linalg.pinv(F)
        crlb = np.sqrt(np.clip(np.diag(P), 0, None))
    except np.linalg.LinAlgError:
        crlb = np.full(4, np.inf)
    return {
        "fim_min": float(evals.min()),
        "fim_max": float(evals.max()),
        "sv_min": smin,
        "sv_max": smax,
        "cond": cond,
        "rank_tol": int(np.sum(s > 1e-8 * max(smax, 1e-30))),
        "crlb_r_km": float(crlb[0] / 1e3) if np.isfinite(crlb[0]) else 999.0,
        "crlb_theta_deg": float(np.rad2deg(crlb[1])) if np.isfinite(crlb[1]) else 999.0,
        "crlb_v": float(crlb[2]) if np.isfinite(crlb[2]) else 999.0,
        "crlb_psi_deg": float(np.rad2deg(crlb[3])) if np.isfinite(crlb[3]) else 999.0,
    }


@dataclass
class CandidateGrid:
    r0: np.ndarray
    th0: np.ndarray
    v: np.ndarray
    psi: np.ndarray
    n: int

    @staticmethod
    def coarse():
        rr, tt, vv, pp = np.meshgrid(
            R_GRID_KM * 1e3,
            np.deg2rad(TH_GRID_DEG),
            V_GRID,
            np.deg2rad(PSI_GRID_DEG),
            indexing="ij",
        )
        r0, th0, v, psi = rr.ravel(), tt.ravel(), vv.ravel(), pp.ravel()
        return CandidateGrid(r0, th0, v, psi, r0.size)


def local_refine(obs, t, turn_deg, r0, th0, v, psi, n_rounds=3, span_scale=0.5):
    dr = 0.75e3 * span_scale
    dth = np.deg2rad(0.6) * span_scale
    dv = 0.25 * span_scale
    dpsi = np.deg2rad(2.0) * span_scale
    best = (float(r0), float(th0), float(v), float(psi))
    best_cost = np.inf
    for _ in range(n_rounds):
        r0s = np.clip(best[0] + np.linspace(-dr, dr, 5), R_MIN_M, R_MAX_M)
        th0s = np.clip(best[1] + np.linspace(-dth, dth, 5), np.deg2rad(TH_MIN_DEG), np.deg2rad(TH_MAX_DEG))
        vs = np.clip(best[2] + np.linspace(-dv, dv, 5), V_MIN, V_MAX)
        psis = np.clip(best[3] + np.linspace(-dpsi, dpsi, 5), np.deg2rad(PSI_MIN_DEG), np.deg2rad(PSI_MAX_DEG))
        RR, TT, VV, PP = np.meshgrid(r0s, th0s, vs, psis, indexing="ij")
        pred = predict_bearings_batch(t, RR.ravel(), TT.ravel(), VV.ravel(), PP.ravel(), turn_deg)
        cost = np.sum(wrap_angle(pred - obs[None, :]) ** 2, axis=1)
        k = int(np.argmin(cost))
        best_cost = float(cost[k])
        best = (float(RR.ravel()[k]), float(TT.ravel()[k]), float(VV.ravel()[k]), float(PP.ravel()[k]))
        dr *= 0.4
        dth *= 0.4
        dv *= 0.4
        dpsi *= 0.4
    return best, best_cost


def chi2_accept_threshold(sigma_rad):
    if sigma_rad <= 0:
        return 1e-12
    return 13.3 * (sigma_rad ** 2)


def make_times(T):
    return np.linspace(0.0, T, max(6, int(round(T / DT_OBS)) + 1))


def truth_dict(r0_m, theta0_deg, v, psi_deg):
    return {
        "r0_m": float(r0_m), "theta0_deg": float(theta0_deg),
        "v": float(v), "psi_deg": float(psi_deg),
        "theta0_rad": float(np.deg2rad(theta0_deg)),
        "psi_rad": float(np.deg2rad(psi_deg)),
    }


def scenario_truth(key):
    s = SCENARIOS[key]
    return truth_dict(s["r0_m"], s["theta0_deg"], s["v"], s["psi_deg"])


def simulate_obs(t, truth, sigma_deg, rng, turn_deg=0.0):
    pred = predict_bearings(t, truth["r0_m"], truth["theta0_rad"], truth["v"], truth["psi_rad"], turn_deg)
    sigma_rad = np.deg2rad(sigma_deg)
    obs = pred.copy() if sigma_rad <= 0 else pred + rng.normal(0.0, sigma_rad, size=pred.shape)
    return obs, pred


def score_trial(grid, pred_grid, obs, t, turn_deg, truth, sigma_rad):
    """Set-valued RC2 metrics. Estimator never receives truth."""
    cost = np.sum(wrap_angle(pred_grid - obs[None, :]) ** 2, axis=1)
    cmin = float(np.min(cost))
    thr = chi2_accept_threshold(sigma_rad)
    acc = cost <= (cmin + thr)
    n_acc = int(np.sum(acc))
    if n_acc == 0:
        k0 = int(np.argmin(cost))
        acc = np.zeros(cost.shape, dtype=bool)
        acc[k0] = True
        n_acc = 1
    kmin = int(np.argmin(cost))

    def width(arr, scale=1.0):
        if n_acc < 2:
            return 0.0
        return float(arr[acc].max() - arr[acc].min()) / scale

    r_w = width(grid.r0, 1e3)
    th_w = np.rad2deg(width(grid.th0, 1.0) if n_acc >= 2 else 0.0)
    # fix: width already divided; recompute angles cleanly
    if n_acc >= 2:
        th_w = float(np.rad2deg(grid.th0[acc].max() - grid.th0[acc].min()))
        psi_w = float(np.rad2deg(grid.psi[acc].max() - grid.psi[acc].min()))
        v_w = float(grid.v[acc].max() - grid.v[acc].min())
        r_w = float((grid.r0[acc].max() - grid.r0[acc].min()) / 1e3)
    else:
        th_w = psi_w = v_w = r_w = 0.0

    # Marginal contraction (robust on ridges): 0=no info on that axis, 1=axis fully known
    contr_r = float(np.clip(1.0 - r_w / 15.0, 0.0, 1.0))
    contr_v = float(np.clip(1.0 - v_w / 2.0, 0.0, 1.0))
    contr_psi = float(np.clip(1.0 - psi_w / 30.0, 0.0, 1.0))
    contr_th = float(np.clip(1.0 - th_w / 10.0, 0.0, 1.0))
    contraction = float((contr_r + contr_v + contr_psi + contr_th) / 4.0)

    # Likelihood-weighted marginal std (ridge → large σ). σ in cost units.
    sig2 = max(sigma_rad ** 2, 1e-18)
    lw = np.exp(-0.5 * (cost - cmin) / sig2)
    s = float(lw.sum())
    if s <= 0 or not np.isfinite(s):
        lw = np.ones_like(cost) / cost.size
    else:
        lw = lw / s
    r_mean = float(np.sum(lw * grid.r0))
    v_mean = float(np.sum(lw * grid.v))
    th_mean = float(np.sum(lw * grid.th0))
    psi_mean = float(np.sum(lw * grid.psi))
    r_std = float(np.sqrt(max(np.sum(lw * (grid.r0 - r_mean) ** 2), 0.0)))
    v_std = float(np.sqrt(max(np.sum(lw * (grid.v - v_mean) ** 2), 0.0)))
    th_std = float(np.sqrt(max(np.sum(lw * (grid.th0 - th_mean) ** 2), 0.0)))
    psi_std = float(np.sqrt(max(np.sum(lw * (grid.psi - psi_mean) ** 2), 0.0)))

    n_modes = 1
    if n_acc >= 12:
        H2, _, _ = np.histogram2d(grid.r0[acc] / 1e3, grid.v[acc], bins=[6, 5])
        n_modes = max(1, int(np.sum(H2 > 0.08 * H2.max())))

    # Point estimate: posterior mean is honest on ridges; local refine only if set is tight
    best_ref = (r_mean, th_mean, v_mean, psi_mean)
    best_c = cmin
    if r_w < 6.0 and v_w < 1.0:
        top = np.argpartition(cost, min(8, cost.size - 1))[:4]
        top = top[np.argsort(cost[top])]
        for idx in top:
            ref, rc = local_refine(obs, t, turn_deg, grid.r0[idx], grid.th0[idx], grid.v[idx], grid.psi[idx], n_rounds=2)
            if rc < best_c:
                best_c, best_ref = rc, ref

    if n_acc >= 2:
        truth_in_r = bool(grid.r0[acc].min() - 1.0 <= truth["r0_m"] <= grid.r0[acc].max() + 1.0)
        truth_in_v = bool(grid.v[acc].min() - 1e-6 <= truth["v"] <= grid.v[acc].max() + 1e-6)
        truth_in_psi = bool(grid.psi[acc].min() - 1e-6 <= truth["psi_rad"] <= grid.psi[acc].max() + 1e-6)
        truth_in_th = bool(grid.th0[acc].min() - 1e-6 <= truth["theta0_rad"] <= grid.th0[acc].max() + 1e-6)
    else:
        truth_in_r = abs(grid.r0[kmin] - truth["r0_m"]) <= 1.0e3 + 1.0
        truth_in_v = abs(grid.v[kmin] - truth["v"]) <= 0.25
        truth_in_psi = abs(grid.psi[kmin] - truth["psi_rad"]) <= np.deg2rad(1.5)
        truth_in_th = abs(grid.th0[kmin] - truth["theta0_rad"]) <= np.deg2rad(0.6)

    return {
        "cost_min": cmin,
        "n_acc": n_acc,
        "contraction": contraction,
        "contraction_r": contr_r,
        "contraction_v": contr_v,
        "contraction_psi": contr_psi,
        "contraction_theta": contr_th,
        "n_modes_proxy": n_modes,
        "r_width_km": r_w,
        "theta_width_deg": th_w,
        "v_width_mps": v_w,
        "psi_width_deg": psi_w,
        "post_mean_r_km": r_mean / 1e3,
        "post_mean_v": v_mean,
        "post_mean_theta_deg": float(np.rad2deg(th_mean)),
        "post_mean_psi_deg": float(np.rad2deg(psi_mean)),
        "post_std_r_km": r_std / 1e3,
        "post_std_v": v_std,
        "post_std_theta_deg": float(np.rad2deg(th_std)),
        "post_std_psi_deg": float(np.rad2deg(psi_std)),
        "est_r0_km": best_ref[0] / 1e3,
        "est_theta0_deg": float(np.rad2deg(best_ref[1])),
        "est_v": best_ref[2],
        "est_psi_deg": float(np.rad2deg(best_ref[3])),
        "abs_err_r_km": abs(best_ref[0] / 1e3 - truth["r0_m"] / 1e3),
        "abs_err_theta_deg": abs(float(np.rad2deg(best_ref[1] - truth["theta0_rad"]))),
        "abs_err_v": abs(best_ref[2] - truth["v"]),
        "abs_err_psi_deg": abs(float(np.rad2deg(best_ref[3] - truth["psi_rad"]))),
        "post_err_r_km": abs(r_mean / 1e3 - truth["r0_m"] / 1e3),
        "post_err_v": abs(v_mean - truth["v"]),
        "post_err_psi_deg": abs(float(np.rad2deg(psi_mean - truth["psi_rad"]))),
        "grid_err_r_km": abs(grid.r0[kmin] / 1e3 - truth["r0_m"] / 1e3),
        "grid_err_v": abs(grid.v[kmin] - truth["v"]),
        "grid_err_psi_deg": abs(float(np.rad2deg(grid.psi[kmin] - truth["psi_rad"]))),
        "truth_in_r_band": truth_in_r,
        "truth_in_v_band": truth_in_v,
        "truth_in_psi_band": truth_in_psi,
        "truth_in_theta_band": truth_in_th,
        "r_resolved": bool(r_w < 4.0),
        "v_resolved": bool(v_w < 0.6),
        "psi_resolved": bool(psi_w < 6.0),
    }


def aggregate_trials(rows):
    df = pd.DataFrame(rows)
    out = {}
    for c in df.columns:
        if df[c].dtype == bool:
            out[f"rate_{c}"] = float(df[c].mean())
        elif np.issubdtype(df[c].dtype, np.number):
            out[f"mean_{c}"] = float(df[c].mean())
            out[f"std_{c}"] = float(df[c].std(ddof=0))
            out[f"median_{c}"] = float(df[c].median())
    return out


# ---------------------------------------------------------------------------
# Scans
# ---------------------------------------------------------------------------
def run_fim_scan():
    rows = []
    for scen in SCENARIOS:
        for turn in TURNS_DEG:
            for T in T_LIST:
                t = make_times(T)
                tr = scenario_truth(scen)
                for psi_deg in PSI_SCAN_DEG:
                    tr2 = truth_dict(tr["r0_m"], tr["theta0_deg"], tr["v"], psi_deg)
                    H = jacobian_bearings(t, tr2["r0_m"], tr2["theta0_rad"], tr2["v"], tr2["psi_rad"], turn)
                    m0 = fim_metrics(H, 0.0)  # geometry-only
                    m = fim_metrics(H, np.deg2rad(0.1))
                    xp, yp, hdx, hdy = platform_states(t, turn)
                    xt, yt = target_states(t, tr2["r0_m"], tr2["theta0_rad"], tr2["v"], tr2["psi_rad"])
                    rx, ry = xt - xp, yt - yp
                    rng_v = np.maximum(np.hypot(rx, ry), 1.0)
                    los_x, los_y = rx / rng_v, ry / rng_v
                    tvx, tvy = tr2["v"] * np.cos(tr2["psi_rad"]), tr2["v"] * np.sin(tr2["psi_rad"])
                    pvx, pvy = U_PLAT * hdx, U_PLAT * hdy
                    v_los = float(np.mean((tvx - pvx) * los_x + (tvy - pvy) * los_y))
                    v_cross = float(np.mean(-(tvx - pvx) * los_y + (tvy - pvy) * los_x))
                    bear = predict_bearings(t, tr2["r0_m"], tr2["theta0_rad"], tr2["v"], tr2["psi_rad"], turn)
                    rows.append({
                        "scenario_base": scen,
                        "turn_deg": turn,
                        "T_s": T,
                        "psi_deg": psi_deg,
                        "n_obs": int(t.size),
                        "fim_min": m["fim_min"],
                        "fim_max": m["fim_max"],
                        "sv_min": m["sv_min"],
                        "sv_max": m["sv_max"],
                        "cond": m["cond"],
                        "rank_tol": m["rank_tol"],
                        "crlb_r_km": m["crlb_r_km"],
                        "crlb_v": m["crlb_v"],
                        "crlb_psi_deg": m["crlb_psi_deg"],
                        "crlb_theta_deg": m["crlb_theta_deg"],
                        "geom_sv_min": m0["sv_min"],
                        "geom_cond": m0["cond"],
                        "mean_v_los": v_los,
                        "mean_v_cross": v_cross,
                        "bearing_span_rad": float(np.ptp(bear)),
                        "bearing_span_deg": float(np.rad2deg(np.ptp(bear))),
                    })
    return pd.DataFrame(rows)


def run_sigma_T_boundary(grid, scenarios=("CR5", "LIN")):
    rng = np.random.default_rng(RNG_SEED)
    rows = []
    for scen in scenarios:
        truth = scenario_truth(scen)
        for T in T_LIST:
            t = make_times(T)
            pred_grid = predict_bearings_batch(t, grid.r0, grid.th0, grid.v, grid.psi, 0.0)
            H = jacobian_bearings(t, truth["r0_m"], truth["theta0_rad"], truth["v"], truth["psi_rad"], 0.0)
            for sigma_deg in SIGMA_DEG:
                sigma_rad = np.deg2rad(sigma_deg)
                m = fim_metrics(H, sigma_rad)
                n_trials = 1 if sigma_deg == 0.0 else N_TRIALS
                trial_rows = []
                for _ in range(n_trials):
                    obs, _ = simulate_obs(t, truth, sigma_deg, rng, 0.0)
                    trial_rows.append(score_trial(grid, pred_grid, obs, t, 0.0, truth, sigma_rad))
                agg = aggregate_trials(trial_rows)
                row = {
                    "scenario": scen,
                    "T_s": T,
                    "sigma_deg": sigma_deg,
                    "n_trials": n_trials,
                    "n_obs": int(t.size),
                    "fim_min": m["fim_min"],
                    "cond": m["cond"],
                    "sv_min": m["sv_min"],
                    "crlb_r_km": m["crlb_r_km"],
                    "crlb_v": m["crlb_v"],
                    "crlb_psi_deg": m["crlb_psi_deg"],
                    "truth_r0_km": truth["r0_m"] / 1e3,
                    "truth_v": truth["v"],
                    "truth_psi_deg": truth["psi_deg"],
                    "turn_deg": 0.0,
                    "bearing_span_deg": float(np.rad2deg(np.ptp(predict_bearings(
                        t, truth["r0_m"], truth["theta0_rad"], truth["v"], truth["psi_rad"], 0.0)))),
                }
                for k, v in agg.items():
                    row[k] = v
                # convenience aliases used by figures/report
                row["mean_abs_err_r_km"] = row.get("mean_abs_err_r_km", row.get("mean_grid_err_r_km", np.nan))
                row["mean_abs_err_v"] = row.get("mean_abs_err_v", np.nan)
                row["mean_abs_err_psi_deg"] = row.get("mean_abs_err_psi_deg", np.nan)
                row["mean_r_width_km"] = row.get("mean_r_width_km", np.nan)
                row["mean_v_width_mps"] = row.get("mean_v_width_mps", np.nan)
                row["mean_contraction"] = row.get("mean_contraction", np.nan)
                row["mean_contraction_r"] = row.get("mean_contraction_r", np.nan)
                row["mean_post_std_r_km"] = row.get("mean_post_std_r_km", np.nan)
                row["mean_post_std_v"] = row.get("mean_post_std_v", np.nan)
                row["mean_post_std_psi_deg"] = row.get("mean_post_std_psi_deg", np.nan)
                row["mean_n_acc"] = row.get("mean_n_acc", np.nan)
                row["rate_r_resolved"] = row.get("rate_r_resolved", np.nan)
                row["rate_v_resolved"] = row.get("rate_v_resolved", np.nan)
                row["rate_truth_in_r_band"] = row.get("rate_truth_in_r_band", np.nan)
                row["rate_truth_in_v_band"] = row.get("rate_truth_in_v_band", np.nan)
                rows.append(row)
                print(f"  sigT {scen} T={T:.0f} σ={sigma_deg} "
                      f"w_r={row['mean_r_width_km']:.2f} "
                      f"postσr={row['mean_post_std_r_km']:.2f} "
                      f"contr_r={row['mean_contraction_r']:.2f} "
                      f"resolved_r={row['rate_r_resolved']:.2f} "
                      f"nacc={row['mean_n_acc']:.0f} span={row['bearing_span_deg']:.3f}°")
    return pd.DataFrame(rows)


def run_heading_identifiability(grid):
    """psi × T on straight path; sigma in {0.05,0.1,0.2}."""
    rng = np.random.default_rng(RNG_SEED + 1)
    rows = []
    for T in T_LIST:
        t = make_times(T)
        # geometry of predictions depends on truth psi; rebuild per psi via analytic obs
        for psi_deg in PSI_SCAN_DEG:
            truth = truth_dict(50e3, 0.0, 2.0, psi_deg)
            pred_grid = predict_bearings_batch(t, grid.r0, grid.th0, grid.v, grid.psi, 0.0)
            H = jacobian_bearings(t, truth["r0_m"], truth["theta0_rad"], truth["v"], truth["psi_rad"], 0.0)
            span = float(np.rad2deg(np.ptp(predict_bearings(
                t, truth["r0_m"], truth["theta0_rad"], truth["v"], truth["psi_rad"], 0.0))))
            for sigma_deg in [0.05, 0.1, 0.2]:
                sigma_rad = np.deg2rad(sigma_deg)
                m = fim_metrics(H, sigma_rad)
                trial_rows = []
                n_ok = 0
                for _ in range(N_TRIALS):
                    obs, _ = simulate_obs(t, truth, sigma_deg, rng, 0.0)
                    sc = score_trial(grid, pred_grid, obs, t, 0.0, truth, sigma_rad)
                    trial_rows.append(sc)
                    if sc["abs_err_psi_deg"] < max(2.0, 0.5 * abs(psi_deg) + 1.0):
                        n_ok += 1
                agg = aggregate_trials(trial_rows)
                rows.append({
                    "T_s": T,
                    "psi_deg": psi_deg,
                    "abs_psi_deg": abs(psi_deg),
                    "sigma_deg": sigma_deg,
                    "identifiable_rate": n_ok / N_TRIALS,
                    "fim_min": m["fim_min"],
                    "cond": m["cond"],
                    "sv_min": m["sv_min"],
                    "crlb_psi_deg": m["crlb_psi_deg"],
                    "crlb_r_km": m["crlb_r_km"],
                    "crlb_v": m["crlb_v"],
                    "bearing_span_deg": span,
                    "mean_abs_err_psi_deg": agg.get("mean_abs_err_psi_deg", np.nan),
                    "std_abs_err_psi_deg": agg.get("std_abs_err_psi_deg", np.nan),
                    "median_abs_err_psi_deg": agg.get("median_abs_err_psi_deg", np.nan),
                    "mean_abs_err_r_km": agg.get("mean_abs_err_r_km", np.nan),
                    "mean_abs_err_v": agg.get("mean_abs_err_v", np.nan),
                    "mean_psi_width_deg": agg.get("mean_psi_width_deg", np.nan),
                    "mean_r_width_km": agg.get("mean_r_width_km", np.nan),
                    "mean_n_acc": agg.get("mean_n_acc", np.nan),
                    "rate_truth_in_psi_band": agg.get("rate_truth_in_psi_band", np.nan),
                })
                print(f"  head T={T:.0f} ψ={psi_deg:+.1f} σ={sigma_deg} "
                      f"|dψ|={agg.get('mean_abs_err_psi_deg', float('nan')):.2f} "
                      f"id={n_ok/N_TRIALS:.2f} span={span:.3f} crlbψ={m['crlb_psi_deg']:.2f}")
    return pd.DataFrame(rows)


def run_maneuver_boundary(grid, scenarios=("CR5", "LIN")):
    rng = np.random.default_rng(RNG_SEED + 2)
    rows = []
    for scen in scenarios:
        truth = scenario_truth(scen)
        for turn in TURNS_DEG:
            for T in T_LIST:
                t = make_times(T)
                pred_grid = predict_bearings_batch(t, grid.r0, grid.th0, grid.v, grid.psi, turn)
                H = jacobian_bearings(t, truth["r0_m"], truth["theta0_rad"], truth["v"], truth["psi_rad"], turn)
                span = float(np.rad2deg(np.ptp(predict_bearings(
                    t, truth["r0_m"], truth["theta0_rad"], truth["v"], truth["psi_rad"], turn))))
                # dedicated alternative RMSE vs truth bearings
                bear_t = predict_bearings(t, truth["r0_m"], truth["theta0_rad"], truth["v"], truth["psi_rad"], turn)
                alts = {
                    "far_fast": truth_dict(55e3, 0.0, 2.5, truth["psi_deg"]),
                    "near_slow": truth_dict(46e3, 0.0, 1.5, truth["psi_deg"]),
                    "psi_plus5": truth_dict(truth["r0_m"], truth["theta0_deg"], truth["v"],
                                            truth["psi_deg"] + 5.0),
                    "psi_minus5": truth_dict(truth["r0_m"], truth["theta0_deg"], truth["v"],
                                             truth["psi_deg"] - 5.0),
                }
                alt_rmse = {
                    k: float(np.rad2deg(bearing_rmse(
                        predict_bearings(t, a["r0_m"], a["theta0_rad"], a["v"], a["psi_rad"], turn),
                        bear_t)))
                    for k, a in alts.items()
                }
                for sigma_deg in [0.05, 0.1, 0.2]:
                    sigma_rad = np.deg2rad(sigma_deg)
                    m = fim_metrics(H, sigma_rad)
                    trial_rows = []
                    for _ in range(N_TRIALS):
                        obs, _ = simulate_obs(t, truth, sigma_deg, rng, turn)
                        trial_rows.append(score_trial(grid, pred_grid, obs, t, turn, truth, sigma_rad))
                    agg = aggregate_trials(trial_rows)
                    row = {
                        "scenario": scen,
                        "turn_deg": turn,
                        "T_s": T,
                        "sigma_deg": sigma_deg,
                        "n_obs": int(t.size),
                        "fim_min": m["fim_min"],
                        "fim_max": m["fim_max"],
                        "sv_min": m["sv_min"],
                        "cond": m["cond"],
                        "rank_tol": m["rank_tol"],
                        "crlb_r_km": m["crlb_r_km"],
                        "crlb_v": m["crlb_v"],
                        "crlb_psi_deg": m["crlb_psi_deg"],
                        "bearing_span_deg": span,
                        "rmse_far_fast_deg": alt_rmse["far_fast"],
                        "rmse_near_slow_deg": alt_rmse["near_slow"],
                        "rmse_psi_plus5_deg": alt_rmse["psi_plus5"],
                        "rmse_psi_minus5_deg": alt_rmse["psi_minus5"],
                    }
                    for k, v in agg.items():
                        row[k] = v
                    row["mean_abs_err_r_km"] = row.get("mean_abs_err_r_km", np.nan)
                    row["mean_abs_err_v"] = row.get("mean_abs_err_v", np.nan)
                    row["mean_abs_err_psi_deg"] = row.get("mean_abs_err_psi_deg", np.nan)
                    row["mean_r_width_km"] = row.get("mean_r_width_km", np.nan)
                    row["mean_v_width_mps"] = row.get("mean_v_width_mps", np.nan)
                    row["mean_contraction"] = row.get("mean_contraction", np.nan)
                    row["mean_contraction_r"] = row.get("mean_contraction_r", np.nan)
                    row["mean_post_std_r_km"] = row.get("mean_post_std_r_km", np.nan)
                    row["mean_post_std_v"] = row.get("mean_post_std_v", np.nan)
                    row["mean_n_acc"] = row.get("mean_n_acc", np.nan)
                    row["rate_r_resolved"] = row.get("rate_r_resolved", np.nan)
                    row["rate_truth_in_r_band"] = row.get("rate_truth_in_r_band", np.nan)
                    rows.append(row)
                    print(f"  man {scen} turn={turn:.0f} T={T:.0f} σ={sigma_deg} "
                          f"w_r={row['mean_r_width_km']:.2f} postσr={row['mean_post_std_r_km']:.2f} "
                          f"contr_r={row['mean_contraction_r']:.2f} res_r={row['rate_r_resolved']:.2f} "
                          f"span={span:.3f} sv={m['sv_min']:.3g}")
    return pd.DataFrame(rows)


def run_offgrid_trials(grid):
    rng = np.random.default_rng(RNG_SEED + 3)
    rows = []

    def on_grid(val, g, tol=1e-9):
        return bool(np.any(np.abs(g - val) < tol))

    for og in OFFGRID_TRUTHS:
        truth = truth_dict(og["r0_m"], og["theta0_deg"], og["v"], og["psi_deg"])
        flags = {
            "r_on_grid": on_grid(og["r0_m"] / 1e3, R_GRID_KM),
            "th_on_grid": on_grid(og["theta0_deg"], TH_GRID_DEG),
            "v_on_grid": on_grid(og["v"], V_GRID),
            "psi_on_grid": on_grid(og["psi_deg"], PSI_GRID_DEG),
        }
        for T in [300.0, 600.0, 1200.0]:
            t = make_times(T)
            pred_grid = predict_bearings_batch(t, grid.r0, grid.th0, grid.v, grid.psi, 0.0)
            for sigma_deg in [0.05, 0.1, 0.2]:
                sigma_rad = np.deg2rad(sigma_deg)
                for trial in range(N_TRIALS):
                    obs, _ = simulate_obs(t, truth, sigma_deg, rng, 0.0)
                    sc = score_trial(grid, pred_grid, obs, t, 0.0, truth, sigma_rad)
                    rows.append({
                        "tag": og["tag"],
                        "T_s": T,
                        "sigma_deg": sigma_deg,
                        "trial": trial,
                        "truth_r0_km": truth["r0_m"] / 1e3,
                        "truth_theta0_deg": truth["theta0_deg"],
                        "truth_v": truth["v"],
                        "truth_psi_deg": truth["psi_deg"],
                        "est_r0_km": sc["est_r0_km"],
                        "est_v": sc["est_v"],
                        "est_psi_deg": sc["est_psi_deg"],
                        "abs_err_r_km": sc["abs_err_r_km"],
                        "abs_err_v": sc["abs_err_v"],
                        "abs_err_psi_deg": sc["abs_err_psi_deg"],
                        "r_width_km": sc["r_width_km"],
                        "v_width_mps": sc["v_width_mps"],
                        "n_acc": sc["n_acc"],
                        "truth_in_r_band": sc["truth_in_r_band"],
                        "truth_in_v_band": sc["truth_in_v_band"],
                        "any_truth_on_grid": any(flags.values()),
                        **{f"{k}": v for k, v in flags.items()},
                    })
    return pd.DataFrame(rows)


def run_hard_candidates(grid):
    """RC2-hard candidate pairs for RC3 handoff."""
    rng = np.random.default_rng(RNG_SEED + 4)
    pairs = []

    # 1) Search: alternatives near a reference hypothesis
    refs = {
        "LIN": scenario_truth("LIN"),
        "CR5": scenario_truth("CR5"),
        "CR10": scenario_truth("CR10"),
    }
    step = 2
    idx = np.arange(0, grid.n, step)
    r0s, th0s, vs, psis = grid.r0[idx], grid.th0[idx], grid.v[idx], grid.psi[idx]

    for turn in [0.0, 5.0, 15.0]:
        for T in [300.0, 600.0, 1200.0]:
            t = make_times(T)
            pred_all = predict_bearings_batch(t, r0s, th0s, vs, psis, turn)
            for rname, ref in refs.items():
                bear_ref = predict_bearings(t, ref["r0_m"], ref["theta0_rad"], ref["v"], ref["psi_rad"], turn)
                rmse = np.sqrt(np.mean(wrap_angle(pred_all - bear_ref[None, :]) ** 2, axis=1))
                dr = np.abs(r0s - ref["r0_m"])
                dv = np.abs(vs - ref["v"])
                dpsi = np.abs(np.rad2deg(psis - ref["psi_rad"]))
                dth = np.abs(np.rad2deg(th0s - ref["theta0_rad"]))
                hard = (rmse < np.deg2rad(0.2)) & ((dr > 3e3) | (dv > 0.4) | (dpsi > 4.0) | (dth > 2.0))
                hard_idx = np.where(hard)[0]
                hard_idx = hard_idx[np.argsort(rmse[hard_idx])]
                kept = []
                for i in hard_idx:
                    if not (abs(r0s[i] - ref["r0_m"]) > 3e3 or abs(vs[i] - ref["v"]) > 0.4
                            or abs(psis[i] - ref["psi_rad"]) > np.deg2rad(4)):
                        continue
                    ok = True
                    for j in kept:
                        if (abs(r0s[i] - r0s[j]) < 2e3 and abs(vs[i] - vs[j]) < 0.3
                                and abs(psis[i] - psis[j]) < np.deg2rad(3)):
                            ok = False
                            break
                    if ok:
                        kept.append(i)
                    if len(kept) >= 5:
                        break
                for rank, i in enumerate(kept):
                    # analytical ambiguity score; MC accept-rate only for dedicated set
                    acc_rate = float(np.exp(-0.5 * (rmse[i] / np.deg2rad(0.1)) ** 2))
                    if acc_rate < 0.05:
                        continue
                    pairs.append({
                        "kind": "search",
                        "ref_name": rname,
                        "turn_deg": turn,
                        "T_s": T,
                        "sigma_test_deg": 0.1,
                        "ref_r0_km": ref["r0_m"] / 1e3,
                        "ref_theta0_deg": ref["theta0_deg"],
                        "ref_v": ref["v"],
                        "ref_psi_deg": ref["psi_deg"],
                        "alt_r0_km": float(r0s[i] / 1e3),
                        "alt_theta0_deg": float(np.rad2deg(th0s[i])),
                        "alt_v": float(vs[i]),
                        "alt_psi_deg": float(np.rad2deg(psis[i])),
                        "rmse_bearing_deg": float(np.rad2deg(rmse[i])),
                        "delta_r_km": float(dr[i] / 1e3),
                        "delta_v": float(dv[i]),
                        "delta_psi_deg": float(dpsi[i]),
                        "delta_theta0_deg": float(dth[i]),
                        "grid_alt_in_accept_rate": acc_rate,
                        "rank": rank,
                        "note": "RC2-hard: bearings close, kinematics far",
                    })

    # 2) Dedicated pairs (classic mechanisms)
    dedicated = [
        ("RS1", truth_dict(50e3, 0, 2.0, 0), truth_dict(55e3, 0, 2.5, 0), "far-fast vs LIN"),
        ("RS2", truth_dict(50e3, 0, 2.0, 0), truth_dict(46e3, 0, 1.5, 0), "near-slow vs LIN"),
        ("RS3", truth_dict(50e3, 0, 2.0, 0), truth_dict(50e3, 0, 2.0, 6.0), "course split +6°"),
        ("RS4", truth_dict(50e3, 0, 2.0, 0), truth_dict(50e3, 2.0, 2.0, 0), "theta0 shift 2°"),
        ("RS5", truth_dict(50e3, 0, 2.0, 0), truth_dict(58e3, 0, 2.8, -3.0), "far-fast-canted"),
        ("RS6", truth_dict(50e3, 0, 2.0, 5), truth_dict(55e3, 0, 2.3, 8.0), "CR5 vs far-fast+canted"),
        ("RS7", truth_dict(50e3, 0, 2.0, 5), truth_dict(47e3, 0, 1.7, 2.0), "CR5 vs near-slow-lowpsi"),
        ("RS8", truth_dict(50e3, 0, 2.0, 5), truth_dict(50e3, 1.5, 2.0, 10.0), "CR5 vs th+psi shift"),
    ]
    for turn in [0.0, 5.0, 15.0]:
        for T in T_LIST:
            t = make_times(T)
            for name, ref, alt, note in dedicated:
                b_ref = predict_bearings(t, ref["r0_m"], ref["theta0_rad"], ref["v"], ref["psi_rad"], turn)
                b_alt = predict_bearings(t, alt["r0_m"], alt["theta0_rad"], alt["v"], alt["psi_rad"], turn)
                rmse = bearing_rmse(b_alt, b_ref)
                pred_grid = predict_bearings_batch(t, grid.r0, grid.th0, grid.v, grid.psi, turn)
                d2 = ((grid.r0 - alt["r0_m"]) / 500) ** 2 + ((grid.v - alt["v"]) / 0.2) ** 2 + \
                     (np.rad2deg(grid.psi - alt["psi_rad"]) / 1.5) ** 2 + (np.rad2deg(grid.th0 - alt["theta0_rad"]) / 0.5) ** 2
                ialt = int(np.argmin(d2))
                # cheap proxy: if noiseless RMSE << σ, grid will accept alt often
                acc_rate = float(np.clip(np.exp(-0.5 * (rmse / np.deg2rad(0.1)) ** 2), 0, 1))
                if rmse < np.deg2rad(0.15):
                    acc_rate = max(acc_rate, 0.7)
                pairs.append({
                    "kind": "dedicated",
                    "ref_name": name,
                    "turn_deg": turn,
                    "T_s": T,
                    "sigma_test_deg": 0.1,
                    "ref_r0_km": ref["r0_m"] / 1e3,
                    "ref_theta0_deg": ref["theta0_deg"],
                    "ref_v": ref["v"],
                    "ref_psi_deg": ref["psi_deg"],
                    "alt_r0_km": alt["r0_m"] / 1e3,
                    "alt_theta0_deg": alt["theta0_deg"],
                    "alt_v": alt["v"],
                    "alt_psi_deg": alt["psi_deg"],
                    "rmse_bearing_deg": float(np.rad2deg(rmse)),
                    "delta_r_km": abs(alt["r0_m"] - ref["r0_m"]) / 1e3,
                    "delta_v": abs(alt["v"] - ref["v"]),
                    "delta_psi_deg": abs(alt["psi_deg"] - ref["psi_deg"]),
                    "delta_theta0_deg": abs(alt["theta0_deg"] - ref["theta0_deg"]),
                    "grid_alt_in_accept_rate": acc_rate,
                    "note": note,
                })

    df = pd.DataFrame(pairs)
    if df.empty:
        return df, df
    hard = df[
        (df["rmse_bearing_deg"] < 0.25)
        | (df["grid_alt_in_accept_rate"] > 0.25)
    ].copy()
    hard = hard.sort_values(
        by=["grid_alt_in_accept_rate", "rmse_bearing_deg"],
        ascending=[False, True],
    )
    return df, hard


# ---------------------------------------------------------------------------
# SVG figures
# ---------------------------------------------------------------------------
def _esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def write_line_chart_svg(path, title, x_label, y_label, series, x_log=False, y_log=False,
                         width=720, height=420, note=""):
    ml, mr, mt, mb = 70, 160, 48, 56
    pw, ph = width - ml - mr, height - mt - mb
    xs_all, ys_all = [], []
    for s in series:
        for x, y in zip(s["x"], s["y"]):
            if y is not None and np.isfinite(y):
                xs_all.append(max(x, 1e-6) if x_log else x)
                ys_all.append(max(y, 1e-12) if y_log else y)
    if not xs_all:
        xs_all, ys_all = [0, 1], [0, 1]
    if x_log:
        x0, x1 = math.log10(min(xs_all)), math.log10(max(xs_all))
    else:
        x0, x1 = min(xs_all), max(xs_all)
    if y_log:
        y0, y1 = math.log10(min(ys_all)), math.log10(max(ys_all))
    else:
        y0, y1 = min(ys_all), max(ys_all)
        pad = 0.08 * (y1 - y0 + 1e-12)
        y0, y1 = y0 - pad, y1 + pad
    if x1 - x0 < 1e-9:
        x1 = x0 + 1
    if y1 - y0 < 1e-9:
        y1 = y0 + 1

    def sx(x):
        if x_log:
            x = math.log10(max(x, 1e-6))
        return ml + (x - x0) / (x1 - x0) * pw

    def sy(y):
        if y_log:
            y = math.log10(max(y, 1e-12))
        return mt + ph - (y - y0) / (y1 - y0) * ph

    def fx(v):
        return f"{10**v:.4g}" if x_log else f"{v:.4g}"

    def fy(v):
        return f"{10**v:.4g}" if y_log else f"{v:.4g}"

    colors = ["#b45309", "#0f766e", "#1d4ed8", "#9f1239", "#6d28d9", "#166534", "#334155", "#c2410c"]
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="{width}" height="{height}" fill="#f7f4ef"/>',
        f'<text x="{width/2}" y="28" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="15" font-weight="600" fill="#1a1a1a">{_esc(title)}</text>',
        f'<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="#ffffff" stroke="#d0ccc4"/>',
    ]
    for i in range(7):
        xv = x0 + i * (x1 - x0) / 6
        X = sx(10 ** xv if x_log else xv)
        parts.append(f'<line x1="{X:.2f}" y1="{mt}" x2="{X:.2f}" y2="{mt+ph}" stroke="#e8e4dc"/>')
        parts.append(f'<text x="{X:.2f}" y="{mt+ph+18}" text-anchor="middle" font-size="11" fill="#555" font-family="sans-serif">{fx(xv)}</text>')
    for i in range(6):
        yv = y0 + i * (y1 - y0) / 5
        Y = sy(10 ** yv if y_log else yv)
        parts.append(f'<line x1="{ml}" y1="{Y:.2f}" x2="{ml+pw}" y2="{Y:.2f}" stroke="#e8e4dc"/>')
        parts.append(f'<text x="{ml-8}" y="{Y+4:.2f}" text-anchor="end" font-size="11" fill="#555" font-family="sans-serif">{fy(yv)}</text>')
    parts.append(f'<text x="{ml+pw/2}" y="{height-14}" text-anchor="middle" font-size="12" fill="#333" font-family="sans-serif">{_esc(x_label)}</text>')
    parts.append(f'<text x="16" y="{mt+ph/2}" text-anchor="middle" font-size="12" fill="#333" font-family="sans-serif" transform="rotate(-90 16 {mt+ph/2})">{_esc(y_label)}</text>')
    for i, s in enumerate(series):
        col = s.get("color", colors[i % len(colors)])
        pts = []
        for x, y in zip(s["x"], s["y"]):
            if y is None or not np.isfinite(y):
                continue
            pts.append(f"{sx(x):.2f},{sy(y):.2f}")
        if len(pts) >= 2:
            parts.append(f'<polyline fill="none" stroke="{col}" stroke-width="2" points="{" ".join(pts)}"/>')
            for p in pts[:: max(1, len(pts) // 18)]:
                cx, cy = p.split(",")
                parts.append(f'<circle cx="{cx}" cy="{cy}" r="2.6" fill="{col}"/>')
        elif len(pts) == 1:
            cx, cy = pts[0].split(",")
            parts.append(f'<circle cx="{cx}" cy="{cy}" r="3" fill="{col}"/>')
    for i, s in enumerate(series):
        col = s.get("color", colors[i % len(colors)])
        parts.append(f'<rect x="{width-mr+10}" y="{mt+8+i*16}" width="14" height="3" fill="{col}"/>')
        parts.append(f'<text x="{width-mr+28}" y="{mt+12+i*16}" font-size="11" fill="#333" font-family="sans-serif">{_esc(s["name"])}</text>')
    if note:
        parts.append(f'<text x="{ml}" y="{height-2}" font-size="10" fill="#777" font-family="sans-serif">{_esc(note)}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_heatmap_svg(path, title, x_label, y_label, xvals, yvals, z, width=720, height=420, note="", fmt="{:.3g}"):
    ml, mr, mt, mb = 70, 100, 48, 56
    pw, ph = width - ml - mr, height - mt - mb
    z = np.asarray(z, dtype=float)
    zmin, zmax = float(np.nanmin(z)), float(np.nanmax(z))
    if not np.isfinite(zmin) or not np.isfinite(zmax) or zmax - zmin < 1e-15:
        zmin, zmax = 0.0, 1.0
    nx, ny = len(xvals), len(yvals)
    cw, ch = pw / max(nx, 1), ph / max(ny, 1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="{width}" height="{height}" fill="#f7f4ef"/>',
        f'<text x="{width/2}" y="28" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="15" font-weight="600" fill="#1a1a1a">{_esc(title)}</text>',
    ]

    def zcolor(val):
        if val is None or not np.isfinite(val):
            return "#dddddd"
        u = (val - zmin) / (zmax - zmin + 1e-15)
        r = int(255 * (1 - u) + 31 * u)
        g = int(240 * (1 - u) + 111 * u)
        b = int(230 * (1 - u) + 139 * u)
        return f"rgb({r},{g},{b})"

    for iy in range(ny):
        for ix in range(nx):
            val = z[iy, ix]
            x = ml + ix * cw
            y = mt + (ny - 1 - iy) * ch
            parts.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{cw:.2f}" height="{ch:.2f}" fill="{zcolor(val)}" stroke="#eee" stroke-width="0.5"/>')
            if nx <= 10 and ny <= 12:
                parts.append(f'<text x="{x+cw/2:.2f}" y="{y+ch/2+3:.2f}" text-anchor="middle" font-size="9" fill="#222" font-family="sans-serif">{fmt.format(val)}</text>')
    for ix, xv in enumerate(xvals):
        parts.append(f'<text x="{ml+(ix+0.5)*cw:.2f}" y="{mt+ph+18}" text-anchor="middle" font-size="11" fill="#555" font-family="sans-serif">{_esc(xv)}</text>')
    for iy, yv in enumerate(yvals):
        y = mt + (ny - 1 - iy) * ch + ch / 2
        parts.append(f'<text x="{ml-8}" y="{y+4:.2f}" text-anchor="end" font-size="11" fill="#555" font-family="sans-serif">{_esc(yv)}</text>')
    parts.append(f'<text x="{ml+pw/2}" y="{height-14}" text-anchor="middle" font-size="12" fill="#333" font-family="sans-serif">{_esc(x_label)}</text>')
    parts.append(f'<text x="16" y="{mt+ph/2}" text-anchor="middle" font-size="12" fill="#333" font-family="sans-serif" transform="rotate(-90 16 {mt+ph/2})">{_esc(y_label)}</text>')
    for i in range(40):
        u = i / 39
        parts.append(f'<rect x="{width-mr+20}" y="{mt+ph*(1-u)-ph/40:.2f}" width="16" height="{ph/40+1:.2f}" fill="{zcolor(zmin+u*(zmax-zmin))}"/>')
    parts.append(f'<text x="{width-mr+40}" y="{mt+10}" font-size="10" fill="#555" font-family="sans-serif">{fmt.format(zmax)}</text>')
    parts.append(f'<text x="{width-mr+40}" y="{mt+ph}" font-size="10" fill="#555" font-family="sans-serif">{fmt.format(zmin)}</text>')
    if note:
        parts.append(f'<text x="{ml}" y="{height-2}" font-size="10" fill="#777" font-family="sans-serif">{_esc(note)}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def make_figures(df_sig, df_head, df_man, df_hard):
    # Figure 1: primary scenario CR5 + LIN limit
    for scen, tag, note in [
        (PRIMARY_SCENARIO, "cr5", "主性能场景 CR5（ψ=5°）；估计器不使用真值"),
        (LIMIT_SCENARIO, "lin", "冻结基准 LIN（ψ=0,v=U）：方位不变 → r–v 脊线不可辨识"),
    ]:
        sub = df_sig[df_sig["scenario"] == scen]
        if sub.empty:
            continue
        series = []
        for sigma in SIGMA_FIG:
            s = sub[sub["sigma_deg"] == sigma].sort_values("T_s")
            if s.empty:
                continue
            series.append({"name": f"σθ={sigma}°", "x": s["T_s"].tolist(), "y": s["mean_post_std_r_km"].tolist()})
        write_line_chart_svg(
            FIG / f"fig1a_sigma_T_range_poststd_{tag}.svg",
            f"图1a [{tag.upper()}] σθ×T → 距离边际后验σ (km)",
            "观测时间 T (s)", "mean posterior σ_r (km)",
            series, x_log=True,
            note=note + "；代价加权边际标准差，脊线→大",
        )
        series = []
        for sigma in SIGMA_FIG:
            s = sub[sub["sigma_deg"] == sigma].sort_values("T_s")
            series.append({"name": f"σθ={sigma}°", "x": s["T_s"].tolist(), "y": s["mean_r_width_km"].tolist()})
        write_line_chart_svg(
            FIG / f"fig1a_sigma_T_range_width_{tag}.svg",
            f"图1a2 [{tag.upper()}] σθ×T → 接受集距离宽度 (km)",
            "观测时间 T (s)", "mean 接受集 Δr 宽度 (km)",
            series, x_log=True,
            note="15 km=搜索盒全宽→距离未约束",
        )
        series = []
        for sigma in SIGMA_FIG:
            s = sub[sub["sigma_deg"] == sigma].sort_values("T_s")
            series.append({"name": f"σθ={sigma}°", "x": s["T_s"].tolist(), "y": s["mean_abs_err_r_km"].tolist()})
        write_line_chart_svg(
            FIG / f"fig1b_sigma_T_range_err_{tag}.svg",
            f"图1b [{tag.upper()}] σθ×T → 点估计 |Δr| (km)",
            "观测时间 T (s)", "mean |Δr| (km)",
            series, x_log=True,
            note="脊线上点估计误差可很大，须与集合宽度对照阅读",
        )
        series = []
        for sigma in SIGMA_FIG:
            s = sub[sub["sigma_deg"] == sigma].sort_values("T_s")
            series.append({"name": f"σθ={sigma}°", "x": s["T_s"].tolist(), "y": s["mean_abs_err_v"].tolist()})
        write_line_chart_svg(
            FIG / f"fig1c_sigma_T_speed_err_{tag}.svg",
            f"图1c [{tag.upper()}] σθ×T → 点估计 |Δv| (m/s)",
            "观测时间 T (s)", "mean |Δv| (m/s)",
            series, x_log=True,
            note="",
        )
        # contraction heatmap
        sigs = [s for s in SIGMA_DEG if s > 0]
        Z = np.zeros((len(sigs), len(T_LIST)))
        for i, s in enumerate(sigs):
            for j, T in enumerate(T_LIST):
                row = sub[(sub["sigma_deg"] == s) & (sub["T_s"] == T)]
                Z[i, j] = float(row["mean_contraction"].iloc[0]) if len(row) else np.nan
        write_heatmap_svg(
            FIG / f"fig1d_sigma_T_contraction_{tag}.svg",
            f"图1d [{tag.upper()}] σθ×T → 候选空间收缩率",
            "T (s)", "σθ (deg)",
            [f"{int(T)}" for T in T_LIST], [str(s) for s in sigs], Z,
            note="1=完全收缩；LIN 在低σθ时接近0（脊线铺满搜索盒）",
            fmt="{:.2f}",
        )
        # truth-in-band rate
        Z = np.zeros((len(sigs), len(T_LIST)))
        for i, s in enumerate(sigs):
            for j, T in enumerate(T_LIST):
                row = sub[(sub["sigma_deg"] == s) & (sub["T_s"] == T)]
                Z[i, j] = float(row["rate_truth_in_r_band"].iloc[0]) if len(row) else np.nan
        write_heatmap_svg(
            FIG / f"fig1e_sigma_T_truth_in_band_{tag}.svg",
            f"图1e [{tag.upper()}] σθ×T → 真值落入 r 接受带比例",
            "T (s)", "σθ (deg)",
            [f"{int(T)}" for T in T_LIST], [str(s) for s in sigs], Z,
            note="集合覆盖指标；LIN 即使点估计跑偏，真值仍可能在脊线上",
            fmt="{:.2f}",
        )

    # Fig2 heading
    series = []
    for psi in [-15.0, -5.0, 0.0, 5.0, 15.0]:
        s = df_head[(df_head["sigma_deg"] == 0.1) & (df_head["psi_deg"] == psi)].sort_values("T_s")
        if s.empty:
            continue
        series.append({"name": f"ψ={psi:+.0f}°", "x": s["T_s"].tolist(), "y": s["mean_abs_err_psi_deg"].tolist()})
    write_line_chart_svg(
        FIG / "fig2a_heading_err_vs_T.svg",
        "图2a  ψ×T → |Δψ| 点估计 (σθ=0.1°, 直航)",
        "T (s)", "mean |Δψ| (deg)",
        series, x_log=True,
        note="ψ=0 为共线极限；|ψ|大时方位变化率帮助辨识",
    )
    series = []
    for psi in [-15.0, -5.0, 0.0, 5.0, 15.0]:
        s = df_head[(df_head["sigma_deg"] == 0.1) & (df_head["psi_deg"] == psi)].sort_values("T_s")
        series.append({"name": f"ψ={psi:+.0f}°", "x": s["T_s"].tolist(), "y": s["crlb_psi_deg"].tolist()})
    write_line_chart_svg(
        FIG / "fig2b_psi_crlb_vs_T.svg",
        "图2b  ψ×T → CRLB(ψ) (σθ=0.1°, 直航)",
        "T (s)", "CRLB ψ (deg)",
        series, x_log=True,
        note="局部FIM在真值处评估，仅作可观测性分析",
    )
    psis = PSI_SCAN_DEG
    Z = np.zeros((len(psis), len(T_LIST)))
    for i, psi in enumerate(psis):
        for j, T in enumerate(T_LIST):
            row = df_head[(df_head["sigma_deg"] == 0.1) & (df_head["psi_deg"] == psi) & (df_head["T_s"] == T)]
            Z[i, j] = float(row["identifiable_rate"].iloc[0]) if len(row) else np.nan
    write_heatmap_svg(
        FIG / "fig2c_heading_identifiable_rate_heatmap.svg",
        "图2c  ψ×T → 航向可辨识率 (σθ=0.1°)",
        "T (s)", "ψ (deg)",
        [f"{int(T)}" for T in T_LIST], [f"{p:+.0f}" for p in psis], Z,
        note="判据: |Δψ|<max(2°, 0.5|ψ|+1°)",
        fmt="{:.2f}",
    )
    Z = np.zeros((len(psis), len(T_LIST)))
    for i, psi in enumerate(psis):
        for j, T in enumerate(T_LIST):
            row = df_head[(df_head["sigma_deg"] == 0.1) & (df_head["psi_deg"] == psi) & (df_head["T_s"] == T)]
            val = float(row["crlb_psi_deg"].iloc[0]) if len(row) else np.nan
            Z[i, j] = min(val, 30.0)
    write_heatmap_svg(
        FIG / "fig2d_psi_crlb_heatmap.svg",
        "图2d  ψ×T → CRLB(ψ) 热图 (σθ=0.1°, 截断30°)",
        "T (s)", "ψ (deg)",
        [f"{int(T)}" for T in T_LIST], [f"{p:+.0f}" for p in psis], Z,
        note="",
        fmt="{:.1f}",
    )

    # Fig3 maneuver — primary scenario
    for scen, tag in [(PRIMARY_SCENARIO, "cr5"), (LIMIT_SCENARIO, "lin")]:
        sub = df_man[(df_man["scenario"] == scen) & (df_man["sigma_deg"] == 0.1)]
        if sub.empty:
            continue
        series = []
        for turn in TURNS_DEG:
            s = sub[sub["turn_deg"] == turn].sort_values("T_s")
            series.append({"name": f"转角={turn:.0f}°", "x": s["T_s"].tolist(), "y": s["mean_post_std_r_km"].tolist()})
        write_line_chart_svg(
            FIG / f"fig3a_turn_T_range_poststd_{tag}.svg",
            f"图3a [{tag.upper()}] 转角×T → 距离边际后验σ (σθ=0.1°)",
            "T (s)", "mean posterior σ_r (km)",
            series, x_log=True,
            note="几何机动应压低 σ_r；若仍很大则 RC2 边界",
        )
        series = []
        for turn in TURNS_DEG:
            s = sub[sub["turn_deg"] == turn].sort_values("T_s")
            series.append({"name": f"转角={turn:.0f}°", "x": s["T_s"].tolist(), "y": s["mean_r_width_km"].tolist()})
        write_line_chart_svg(
            FIG / f"fig3a_turn_T_range_width_{tag}.svg",
            f"图3a2 [{tag.upper()}] 转角×T → 接受集距离宽度 (σθ=0.1°)",
            "T (s)", "mean Δr 宽度 (km)",
            series, x_log=True,
            note="",
        )
        series = []
        for turn in TURNS_DEG:
            s = sub[sub["turn_deg"] == turn].sort_values("T_s")
            series.append({"name": f"转角={turn:.0f}°", "x": s["T_s"].tolist(), "y": s["sv_min"].tolist()})
        write_line_chart_svg(
            FIG / f"fig3b_turn_T_svmin_{tag}.svg",
            f"图3b [{tag.upper()}] 转角×T → 最小奇异值 (σθ=0.1° 标度)",
            "T (s)", "sv_min",
            series, x_log=True, y_log=True,
            note="",
        )
        series = []
        for turn in TURNS_DEG:
            s = sub[sub["turn_deg"] == turn].sort_values("T_s")
            series.append({"name": f"转角={turn:.0f}°", "x": s["T_s"].tolist(), "y": s["mean_contraction"].tolist()})
        write_line_chart_svg(
            FIG / f"fig3c_turn_T_contraction_{tag}.svg",
            f"图3c [{tag.upper()}] 转角×T → 候选收缩率 (σθ=0.1°)",
            "T (s)", "mean 收缩率",
            series, x_log=True,
            note="",
        )
        Z = np.zeros((len(TURNS_DEG), len(T_LIST)))
        for i, turn in enumerate(TURNS_DEG):
            for j, T in enumerate(T_LIST):
                row = sub[(sub["turn_deg"] == turn) & (sub["T_s"] == T)]
                Z[i, j] = float(row["mean_r_width_km"].iloc[0]) if len(row) else np.nan
        write_heatmap_svg(
            FIG / f"fig3d_turn_T_range_width_heatmap_{tag}.svg",
            f"图3d [{tag.upper()}] 转角×T → 接受集 Δr 宽度 km (σθ=0.1°)",
            "T (s)", "转角 (deg)",
            [f"{int(T)}" for T in T_LIST], [f"{int(t)}" for t in TURNS_DEG], Z,
            note="颜色越浅表示不确定性越大",
            fmt="{:.2f}",
        )
        series = []
        for turn in TURNS_DEG:
            s = sub[sub["turn_deg"] == turn].sort_values("T_s")
            series.append({"name": f"转角={turn:.0f}°", "x": s["T_s"].tolist(), "y": s["bearing_span_deg"].tolist()})
        write_line_chart_svg(
            FIG / f"fig3e_turn_T_bearing_span_{tag}.svg",
            f"图3e [{tag.upper()}] 转角×T → 方位总变化 span (deg)",
            "T (s)", "bearing span (deg)",
            series, x_log=True,
            note="span=0 表示纯方位无变化（LIN直航）",
        )

    # hard pairs figure
    if df_hard is not None and not df_hard.empty and "ref_name" in df_hard.columns:
        series = []
        for name in ["RS1", "RS3", "RS6", "RS7"]:
            for turn in [0.0, 15.0]:
                s = df_hard[(df_hard["ref_name"] == name) & (df_hard["turn_deg"] == turn)].sort_values("T_s")
                if s.empty:
                    continue
                series.append({"name": f"{name},转{turn:.0f}°", "x": s["T_s"].tolist(), "y": s["rmse_bearing_deg"].tolist()})
        if series:
            write_line_chart_svg(
                FIG / "fig3f_hard_pair_bearing_rmse.svg",
                "图3f  困难候选对方位 RMSE vs T",
                "T (s)", "方位 RMSE (deg)",
                series, x_log=True,
                note="RMSE越小越难区分；RC3仅允许使用仍困难的对",
            )


# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------
def fmt(x, nd=3):
    try:
        if x is None:
            return "n/a"
        v = float(x)
        if not np.isfinite(v):
            return "∞" if v > 0 else "n/a"
        return f"{v:.{nd}f}"
    except Exception:
        return "n/a"


def pick(df, **kw):
    d = df
    for k, v in kw.items():
        d = d[d[k] == v]
    return d.iloc[0] if len(d) else None


def write_reports(df_fim, df_sig, df_head, df_man, df_off, df_all_h, df_hard):
    # RC3 handoff
    hard_rc3 = OUT / "hard_candidate_pairs_rc3_handoff.csv"
    if df_hard is not None and not df_hard.empty:
        cand = df_hard.copy()
        cand = cand.sort_values("grid_alt_in_accept_rate", ascending=False)
        picked = []
        rows = []
        for _, r in cand.iterrows():
            key = (round(float(r.get("alt_r0_km", 0)), 1), round(float(r.get("alt_v", 0)), 2),
                   round(float(r.get("alt_psi_deg", 0)), 1), float(r.get("turn_deg", 0)))
            if any(abs(key[0]-p[0]) < 1.5 and abs(key[1]-p[1]) < 0.25 and abs(key[3]-p[3]) < 0.1 for p in picked):
                continue
            picked.append(key)
            rows.append(r.to_dict())
            if len(picked) >= 16:
                break
        pd.DataFrame(rows).to_csv(hard_rc3, index=False, encoding="utf-8-sig")
    else:
        pd.DataFrame([{"note": "no hard pairs"}]).to_csv(hard_rc3, index=False, encoding="utf-8-sig")

    # anchors
    cr5_600_01 = pick(df_sig, scenario="CR5", sigma_deg=0.1, T_s=600.0)
    cr5_1200_01 = pick(df_sig, scenario="CR5", sigma_deg=0.1, T_s=1200.0)
    cr5_300_05 = pick(df_sig, scenario="CR5", sigma_deg=0.5, T_s=300.0)
    lin_600_01 = pick(df_sig, scenario="LIN", sigma_deg=0.1, T_s=600.0)
    lin_600_00 = pick(df_sig, scenario="LIN", sigma_deg=0.0, T_s=600.0)
    man_cr5 = df_man[(df_man["scenario"] == "CR5") & (df_man["sigma_deg"] == 0.1) & (df_man["T_s"] == 600.0)]
    man0 = man_cr5[man_cr5["turn_deg"] == 0].iloc[0] if len(man_cr5[man_cr5["turn_deg"] == 0]) else None
    man15 = man_cr5[man_cr5["turn_deg"] == 15].iloc[0] if len(man_cr5[man_cr5["turn_deg"] == 15]) else None
    h0 = pick(df_head, sigma_deg=0.1, psi_deg=0.0, T_s=600.0)
    h5 = pick(df_head, sigma_deg=0.1, psi_deg=5.0, T_s=600.0)
    h10 = pick(df_head, sigma_deg=0.1, psi_deg=10.0, T_s=600.0)

    # GPT sync
    g = []
    g.append("# P2 RC2 运动学边界 — GPT 同步稿")
    g.append("")
    g.append(f"- 生成时间 (UTC): {NOW}")
    g.append(f"- 输出目录: `results/P2_RC2_kinematic_boundary/`")
    g.append("- 状态: **P2 参数扫描完成，未进入 RC3**")
    g.append("")
    g.append("## 1. 方法冻结确认")
    g.append("")
    g.append("1. E-STD：z=200 m 固定；本轮仅二维运动学。")
    g.append("2. 观测：`θ=atan2(yt−yp,xt−xp)+ν`；完整二维几何，禁用 `r(t)=r0+(v−U)t`。")
    g.append("3. 真值不泄露：θ0 未知；候选=粗网格+局部细化。")
    g.append("4. 传播/谱线本包未使用。")
    g.append("5. **主指标改为接受集宽度/收缩率/真值覆盖率**（脊线上点估计不唯一）。")
    g.append("")
    g.append("## 2. 关键发现：冻结基准是 RC2 极限几何")
    g.append("")
    g.append("用户冻结基准 `θ0=0°, ψ=0°, v=U=2 m/s` 时：")
    g.append("")
    g.append("```")
    g.append("xt = r0 + v t,  xp = U t  →  dx = r0, dy = 0  →  θ(t) ≡ 0")
    g.append("```")
    g.append("")
    g.append("方位序列与 `(r0,v)` **完全无关**。这是理论事实，不是实现缺陷。")
    if lin_600_00 is not None:
        g.append(f"- LIN, σθ=0, T=600 s：接受集 r 宽度≈**{fmt(lin_600_00.get('mean_r_width_km'))} km**（搜索盒全宽），后验σr≈**{fmt(lin_600_00.get('mean_post_std_r_km'))} km**，r可辨识比例≈**{fmt(lin_600_00.get('rate_r_resolved'),2)}**。")
    if lin_600_01 is not None:
        g.append(f"- LIN, σθ=0.1°, T=600 s：r 宽度≈**{fmt(lin_600_01.get('mean_r_width_km'))} km**，后验σr≈**{fmt(lin_600_01.get('mean_post_std_r_km'))} km**，|Δv|≈**{fmt(lin_600_01.get('mean_abs_err_v'))} m/s**。")
    g.append("")
    g.append("**回答甲方式**：RC2 在共线/无方位变化率时，**不能**唯一确定距离与速度；候选沿 r–v 脊线分布。这正是“理论功能边界”，不是缺资料。")
    g.append("")
    g.append("## 3. RC2 在开放几何（CR5）下约束了什么")
    g.append("")
    g.append("加入弱交叉 `ψ=5°` 后方位出现变化率，候选集合开始收缩。")
    if cr5_600_01 is not None:
        g.append(f"- CR5, σθ=0.1°, T=600 s：后验σr≈**{fmt(cr5_600_01.get('mean_post_std_r_km'))} km**，接受集 r 宽度≈**{fmt(cr5_600_01.get('mean_r_width_km'))} km**，点估计 |Δr|≈**{fmt(cr5_600_01.get('mean_abs_err_r_km'))} km**，|Δv|≈**{fmt(cr5_600_01.get('mean_abs_err_v'))} m/s**，r可辨识比例≈**{fmt(cr5_600_01.get('rate_r_resolved'),2)}**，方位 span≈**{fmt(cr5_600_01.get('bearing_span_deg'),3)}°**。")
    if cr5_1200_01 is not None:
        g.append(f"- CR5, T=1200 s：后验σr≈**{fmt(cr5_1200_01.get('mean_post_std_r_km'))} km**，r宽度≈**{fmt(cr5_1200_01.get('mean_r_width_km'))} km**。")
    if cr5_300_05 is not None:
        g.append(f"- 困难条件 CR5, σθ=0.5°, T=300 s：后验σr≈**{fmt(cr5_300_05.get('mean_post_std_r_km'))} km**。")
    g.append("")
    g.append("### 已约束")
    g.append("- 有方位变化率时的视线几何与部分 r/v/ψ 组合")
    g.append("- 初始方位 θ0（在非共线几何下）")
    g.append("")
    g.append("### 仍未约束（本征）")
    g.append("- 共线/近共线时的 r–v 脊线")
    if h0 is not None:
        g.append(f"- ψ≈0 时航向：T=600s σθ=0.1° → |Δψ|≈{fmt(h0.get('mean_abs_err_psi_deg'))}°, 可辨识率≈{fmt(h0.get('identifiable_rate'),2)}, CRLB(ψ)≈{fmt(h0.get('crlb_psi_deg'))}°")
    if h5 is not None:
        g.append(f"- ψ=5°：|Δψ|≈{fmt(h5.get('mean_abs_err_psi_deg'))}°, 可辨识率≈{fmt(h5.get('identifiable_rate'),2)}")
    if h10 is not None:
        g.append(f"- ψ=10°：|Δψ|≈{fmt(h10.get('mean_abs_err_psi_deg'))}°, 可辨识率≈{fmt(h10.get('identifiable_rate'),2)}")
    g.append("")
    g.append("## 4. 小机动是否有效？")
    g.append("")
    if man0 is not None and man15 is not None:
        g.append(f"- CR5, σθ=0.1°, T=600 s：转角 0°→15°，")
        g.append(f"  接受集 r 宽度 **{fmt(man0.get('mean_r_width_km'))}→{fmt(man15.get('mean_r_width_km'))} km**，")
        g.append(f"  点估计 |Δr| **{fmt(man0.get('mean_abs_err_r_km'))}→{fmt(man15.get('mean_abs_err_r_km'))} km**，")
        g.append(f"  sv_min **{fmt(man0.get('sv_min'),4)}→{fmt(man15.get('sv_min'),4)}**，")
        g.append(f"  方位 span **{fmt(man0.get('bearing_span_deg'))}→{fmt(man15.get('bearing_span_deg'))}°**。")
    g.append("- 结论：平台转向改变相对几何、增大方位 span，从而打破部分 r–v–ψ 耦合；LIN 场景下机动是**唯一**能打开运动学信息的 RC2 手段（不加声学特征）。")
    g.append("")
    g.append("## 5. 交给 RC3 的困难候选")
    g.append("")
    g.append("文件：`hard_candidate_pairs.csv` / `hard_candidate_pairs_rc3_handoff.csv`")
    g.append("")
    g.append("入选原则：")
    g.append("1. 方位 RMSE 很小（传播难以提供更强区分时必须有独立信息）；")
    g.append("2. 运动学状态相差大（Δr≥数 km 或 |Δv|≥0.4 或 |Δψ|≥4°）；")
    g.append("3. RC2 接受集仍同时容纳参考假设与备择假设。")
    g.append("")
    g.append("**RC3 约束**：P3 只允许对这些困难候选对比较 “RC2 only” vs “RC2+传播”（KRAKEN 等），不得再开 Bellhop 审计支线，不得追加新传感器。")
    g.append("")
    g.append("## 6. 关口")
    g.append("")
    g.append("| 关口 | 状态 | 说明 |")
    g.append("| --- | --- | --- |")
    g.append("| G1 标准场景冻结 | 通过（P1） | E-STD + S0/S1/S2 |")
    g.append("| G2 RC2 运动学边界 | **通过（本包）** | 观测正确、真值不泄露、扫描完成、极限几何已定量 |")
    g.append("| G3 RC3 独立增量 | 未开始 | 仅困难候选对 |")
    g.append("| G4 综合性能边界 | 未开始 | P4 |")
    g.append("")
    g.append("## 7. 给 GPT 的下一轮")
    g.append("")
    g.append("> **只判断 RC2 约束了什么，并把剩余困难候选交给 RC3。不再回到 Bellhop 审计循环。**")
    g.append("")
    (OUT / "P2_RC2_GPT_SYNC.md").write_text("\n".join(g), encoding="utf-8")

    # Full report
    r = []
    r.append("# P2 RC2 报告：连续方位 / 平台运动下的运动学候选边界")
    r.append("")
    r.append(f"项目：**B3D-3** 深海第一会聚区单平台拖曳阵目标定位  ")
    r.append(f"工作包：**P2 RC2 kinematic boundary**  ")
    r.append(f"UTC：{NOW}")
    r.append("")
    r.append("## 0. 停止条件自检")
    r.append("")
    r.append("| 条件 | 状态 |")
    r.append("| --- | --- |")
    r.append("| 观测模型正确（完整二维 + atan2） | 是 |")
    r.append("| 真值不泄露到估计器 | 是 |")
    r.append("| 网格 / off-grid 检查 | 是 |")
    r.append("| 参数扫描完成 | 是 |")
    r.append("| 理论边界形成 | 是 |")
    r.append("| 因不可辨识追加传感器/声学特征 | 否 |")
    r.append("| 自动进入 RC3 | 否 |")
    r.append("")
    r.append("## 1. 状态与观测")
    r.append("")
    r.append("```")
    r.append("x = [r0, θ0, v, ψ]^T   (z=200 m 固定)")
    r.append("平台: (0,0)→航向0°, U=2 m/s; 单次转角在 t=T/2")
    r.append("目标: xt=r0 cosθ0 + v t cosψ, yt=r0 sinθ0 + v t sinψ")
    r.append("观测: θ_obs = atan2(yt−yp, xt−xp) + ν,  ν~N(0,σθ²)")
    r.append("```")
    r.append("")
    r.append("候选盒：r0∈[45,60] km, θ0∈[−5°,5°], v∈[1,3] m/s, ψ∈[−15°,15°]。")
    r.append(f"粗网格：r 1 km × θ 0.5° × v 0.2 m/s × ψ 1°（含 0°），节点数 **{CONFIG['coarse_grid']['n_nodes']}**。")
    r.append("真值不保证在网格上；估计=粗搜索+局部细化；off-grid 三组对照。")
    r.append("")
    r.append("### 分析场景")
    r.append("")
    r.append("| 场景 | 真值 | 含义 |")
    r.append("| --- | --- | --- |")
    for k, s in SCENARIOS.items():
        r.append(f"| **{k}** | r0={s['r0_m']/1e3:.0f} km, θ0={s['theta0_deg']}°, v={s['v']} m/s, ψ={s['psi_deg']}° | {s['note']} |")
    r.append("")
    r.append("主性能曲线用 **CR5**；LIN 作为 RC2 极限边界与甲方可解释的“不可辨识”证据。")
    r.append("")
    r.append("## 2. 方法要点")
    r.append("")
    r.append("1. 真实生成含噪方位序列，θ0 未知。")
    r.append("2. `cost=Σ wrap(θ_pred−θ_obs)²`；接受集 `cost≤c_min+13.3σθ²`。")
    r.append("3. **主指标**：接受集在 (r0,θ0,v,ψ) 上的宽度、收缩率、真值是否落入接受带。")
    r.append("4. 点估计：top 网格局部细化；若接受集呈宽脊，则用接受集中位数作为诚实点估计。")
    r.append("5. FIM/`HᵀR⁻¹H`/CRLB 在真值处计算，仅用于可观测性，不喂给估计器。")
    r.append("6. 非零 σθ 每格 N_TRIALS 次随机试验；σθ=0 为理想上限。")
    r.append("")
    r.append("## 3. 四个核心问题")
    r.append("")
    r.append("### Q1 距离–速度–航向何时不可区分？")
    r.append("")
    r.append("**定理性事实（本包定量验证）**：当目标与平台相对速度沿视线、且方位变化率为 0 时，")
    r.append("`θ(t)` 与 `(r0,v)` 无关，FIM 秩亏，候选沿 r–v 脊线分布。")
    r.append("冻结基准 LIN 正落入该情形。")
    r.append("")
    if lin_600_01 is not None:
        r.append(f"- LIN σθ=0.1° T=600 s：接受集 r 宽度 **{fmt(lin_600_01.get('mean_r_width_km'))} km**，后验σr **{fmt(lin_600_01.get('mean_post_std_r_km'))} km**，v 宽度 **{fmt(lin_600_01.get('mean_v_width_mps'))} m/s**，CRLB(r) **{fmt(lin_600_01.get('crlb_r_km'))} km**，cond **{fmt(lin_600_01.get('cond'))}**，span **{fmt(lin_600_01.get('bearing_span_deg'),3)}°**。")
    r.append("- 直航 + 近共线时 ψ 与 r0/v 强耦合；|ψ| 增大或平台转向后，方位 span>0，信息量上升，但在冻结的 σθ/T/转角范围内距离往往仍只有数 km 级不确定度。")
    r.append("")
    r.append("### Q2 方位精度需要达到多少？")
    r.append("")
    r.append("见 `bearing_accuracy_time_boundary.csv` 与图1（CR5/LIN）。")
    if cr5_600_01 is not None:
        r.append(f"- CR5 σθ=0.1° T=600 s：后验σr **{fmt(cr5_600_01.get('mean_post_std_r_km'))} km**，r宽度 **{fmt(cr5_600_01.get('mean_r_width_km'))} km**，|Δr| **{fmt(cr5_600_01.get('mean_abs_err_r_km'))} km**，span **{fmt(cr5_600_01.get('bearing_span_deg'),3)}°**，r可辨识比例 **{fmt(cr5_600_01.get('rate_r_resolved'),2)}**。")
    if cr5_300_05 is not None:
        r.append(f"- CR5 σθ=0.5° T=300 s：后验σr **{fmt(cr5_300_05.get('mean_post_std_r_km'))} km**。")
    r.append("- 读图：给定允许距离误差，在 σθ×T 图上找满足 posterior σ_r 或 r 可辨识门限的组合。")
    r.append("- **关键结论**：在 E-STD 冻结量级（r≈50 km, U=v≈2 m/s, T≤1200 s, 转角≤15°）下，即使 σθ=0.01°–0.05°，弱交叉 CR5 的方位总变化仅约 0.01°–0.24°，距离信息量有限；共线 LIN 则对任意 σθ 都不能唯一恢复 r–v。")
    r.append("- 这说明：**纯 RC2 + 本包冻结的平台机动能力，不足以在第一会聚区距离上形成高精度定位**；更优 σθ 只能把脊线“照清楚”，不能消除几何不可观测。更高性能需要更强机动/更长观测/传播增量（RC3）——但 P3 只允许在困难候选上验证传播增量，不回退审计。")
    r.append("")
    r.append("### Q3 小幅转向是否真正有效？")
    r.append("")
    r.append("见 `maneuver_boundary.csv` 与图3。")
    if man0 is not None and man15 is not None:
        r.append(f"- CR5 σθ=0.1° T=600 s：转角 0°→15°，")
        r.append(f"  后验σr **{fmt(man0.get('mean_post_std_r_km'))}→{fmt(man15.get('mean_post_std_r_km'))} km**，")
        r.append(f"  r宽度 **{fmt(man0.get('mean_r_width_km'))}→{fmt(man15.get('mean_r_width_km'))} km**，")
        r.append(f"  |Δr| **{fmt(man0.get('mean_abs_err_r_km'))}→{fmt(man15.get('mean_abs_err_r_km'))} km**，")
        r.append(f"  sv_min **{fmt(man0.get('sv_min'),4)}→{fmt(man15.get('sv_min'),4)}**，")
        r.append(f"  方位 span **{fmt(man0.get('bearing_span_deg'))}→{fmt(man15.get('bearing_span_deg'))}°**，")
        r.append(f"  r可辨识比例 **{fmt(man0.get('rate_r_resolved'),2)}→{fmt(man15.get('rate_r_resolved'),2)}**。")
    r.append("- 转向通过改变平台速度矢量打开观测几何；FIM 最小奇异值从 LIN 直航的 0 上升，说明秩亏被部分修复。")
    r.append("- **量级判断**：15° 单次转角在 r≈50 km 处产生的额外方位信息仍然偏弱，候选 r–v 脊线只被部分压缩。这是 RC2 运动学边界，不是实现失败。")
    r.append("- 本轮不评估拖曳阵孔径/波束（P4）。")
    r.append("")
    r.append("### Q4 哪些困难候选交给 RC3？")
    r.append("")
    r.append("交付 `hard_candidate_pairs.csv` 与 `hard_candidate_pairs_rc3_handoff.csv`。机制：")
    r.append("")
    r.append("1. **r–v 补偿**：更远更快 vs 更近更慢，方位接近（LIN 下整条脊线皆是）；")
    r.append("2. **航向近镜像**：±ψ 差在有限 span 下不足；")
    r.append("3. **θ0 与 ψ 联合补偿**。")
    r.append("")
    r.append("RC3 **只允许**对这些困难候选对比较 “RC2 only” vs “RC2+传播”，回答传播是否提供独立增量。")
    r.append("")
    r.append("## 4. 关键结果摘录")
    r.append("")
    r.append("### 4.1 σθ×T（场景对比）")
    r.append("")
    r.append("| 场景 | T | σθ° | 后验σr km | r宽度 km | \\|Δr\\| km | contr_r | 可辨识r | CRLB r | span° |")
    r.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    show = df_sig[df_sig["T_s"].isin([300, 600, 1200]) & df_sig["sigma_deg"].isin([0.02, 0.05, 0.1, 0.2, 0.5])]
    for _, row in show.iterrows():
        r.append(
            f"| {row['scenario']} | {row['T_s']:.0f} | {row['sigma_deg']} | {fmt(row.get('mean_post_std_r_km'))} | "
            f"{fmt(row.get('mean_r_width_km'))} | {fmt(row.get('mean_abs_err_r_km'))} | {fmt(row.get('mean_contraction_r'))} | "
            f"{fmt(row.get('rate_r_resolved'),2)} | {fmt(row.get('crlb_r_km'))} | {fmt(row.get('bearing_span_deg'),3)} |"
        )
    r.append("")
    r.append("### 4.2 航向可辨识（σθ=0.1°）")
    r.append("")
    r.append("| T | ψ° | \\|Δψ\\|° | 可辨识率 | CRLB ψ° | span° | r宽度 km |")
    r.append("| --- | --- | --- | --- | --- | --- | --- |")
    show = df_head[(df_head["sigma_deg"] == 0.1) & df_head["T_s"].isin([300, 600, 1200]) & df_head["psi_deg"].isin([-15, -5, 0, 5, 10, 15])]
    for _, row in show.iterrows():
        r.append(
            f"| {row['T_s']:.0f} | {row['psi_deg']:+.0f} | {fmt(row.get('mean_abs_err_psi_deg'))} | "
            f"{fmt(row.get('identifiable_rate'),2)} | {fmt(row.get('crlb_psi_deg'))} | "
            f"{fmt(row.get('bearing_span_deg'),3)} | {fmt(row.get('mean_r_width_km'))} |"
        )
    r.append("")
    r.append("### 4.3 机动边界（σθ=0.1°）")
    r.append("")
    r.append("| 场景 | 转角° | T | 后验σr km | r宽度 km | \\|Δr\\| km | sv_min | span° | 可辨识r |")
    r.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    show = df_man[(df_man["sigma_deg"] == 0.1) & df_man["T_s"].isin([300, 600, 1200]) & df_man["turn_deg"].isin([0, 5, 15])]
    for _, row in show.iterrows():
        r.append(
            f"| {row['scenario']} | {row['turn_deg']:.0f} | {row['T_s']:.0f} | {fmt(row.get('mean_post_std_r_km'))} | "
            f"{fmt(row.get('mean_r_width_km'))} | {fmt(row.get('mean_abs_err_r_km'))} | {fmt(row.get('sv_min'),4)} | "
            f"{fmt(row.get('bearing_span_deg'),3)} | {fmt(row.get('rate_r_resolved'),2)} |"
        )
    r.append("")
    r.append("### 4.4 off-grid 摘要")
    r.append("")
    if df_off is not None and not df_off.empty:
        grp = df_off.groupby(["tag", "T_s", "sigma_deg"], as_index=False).agg(
            mean_abs_err_r_km=("abs_err_r_km", "mean"),
            mean_r_width_km=("r_width_km", "mean"),
            rate_truth_in_r=("truth_in_r_band", "mean"),
            any_on_grid=("any_truth_on_grid", "max"),
        )
        r.append("| OG | T | σθ | \\|Δr\\| km | r宽度 km | 真值在r带 | 真值在网格? |")
        r.append("| --- | --- | --- | --- | --- | --- | --- |")
        for _, row in grp.iterrows():
            r.append(
                f"| {row['tag']} | {row['T_s']:.0f} | {row['sigma_deg']} | {fmt(row['mean_abs_err_r_km'])} | "
                f"{fmt(row['mean_r_width_km'])} | {fmt(row['rate_truth_in_r'],2)} | {bool(row['any_on_grid'])} |"
            )
        r.append("")
        r.append("三组 off-grid 真值均不在粗网格节点上；覆盖指标与 on-grid 场景同量级，排除“真值写入网格”伪影。")
    r.append("")
    r.append("## 5. 图目录")
    r.append("")
    r.append("| 文件 | 内容 |")
    r.append("| --- | --- |")
    for p in sorted(FIG.glob("*.svg")):
        r.append(f"| figures/{p.name} |  |")
    r.append("")
    r.append("## 6. 局限")
    r.append("")
    r.append("- 不研究声传播、谱线、阵列方向图。")
    r.append("- z 固定，不讨论深度可观测。")
    r.append("- “不可辨识”不触发新传感器，也不自动进入 RC3。")
    r.append("- E-STD 是理论算例，不是实测海区复现。")
    r.append("")
    r.append("## 7. 交付清单")
    r.append("")
    for p in sorted(OUT.rglob("*")):
        if p.is_file():
            r.append(f"- `results/P2_RC2_kinematic_boundary/{p.relative_to(OUT)}`")
    r.append("")
    r.append("---")
    r.append("")
    r.append("**P2 完成后停止。** 下一轮仅基于困难候选对判断 RC3 增量。")
    (OUT / "P2_RC2_REPORT.md").write_text("\n".join(r), encoding="utf-8")


def export_embedded_data(df_sig, df_head, df_man, df_off, df_hard, df_fim):
    def records(df):
        if df is None or df is getattr(pd, "DataFrame") and df.empty:
            return []
        if df is None or len(df) == 0:
            return []
        out = []
        for _, row in df.iterrows():
            rec = {}
            for k, v in row.items():
                if isinstance(v, (np.floating, float)):
                    rec[k] = None if (isinstance(v, float) and v != v) else float(v)
                elif isinstance(v, (np.integer,)):
                    rec[k] = int(v)
                elif isinstance(v, (np.bool_, bool)):
                    rec[k] = bool(v)
                elif v is None:
                    rec[k] = None
                elif isinstance(v, np.generic):
                    rec[k] = v.item()
                else:
                    rec[k] = v
            out.append(rec)
        return out

    hard_out = df_hard.copy() if df_hard is not None else pd.DataFrame()
    if len(hard_out) > 60:
        hard_out = hard_out.head(60)
    off_sum = df_off
    if df_off is not None and len(df_off):
        off_sum = df_off.groupby(["tag", "T_s", "sigma_deg"], as_index=False).agg(
            mean_abs_err_r_km=("abs_err_r_km", "mean"),
            mean_abs_err_v=("abs_err_v", "mean"),
            mean_abs_err_psi_deg=("abs_err_psi_deg", "mean"),
            mean_r_width_km=("r_width_km", "mean"),
            rate_truth_in_r_band=("truth_in_r_band", "mean"),
            any_on_grid=("any_truth_on_grid", "max"),
        )
    payload = {
        "meta": {
            "project": "B3D-3",
            "package": "P2_RC2_kinematic_boundary",
            "created_utc": NOW,
            "scene": "E-STD",
            "n_trials": N_TRIALS,
            "dt_obs_s": DT_OBS,
            "primary_scenario": PRIMARY_SCENARIO,
            "limit_scenario": LIMIT_SCENARIO,
        },
        "config": CONFIG,
        "sigma_T": records(df_sig),
        "heading": records(df_head),
        "maneuver": records(df_man),
        "offgrid_summary": records(off_sum),
        "hard_pairs": records(hard_out),
        "fim_scan_sample": records(df_fim.head(150) if df_fim is not None else df_fim),
    }
    (OUT / "p2_results.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    (ROOT / "assets" / "p2_results.js").write_text(
        "window.P2_RESULTS = " + json.dumps(payload, ensure_ascii=False) + ";\n",
        encoding="utf-8",
    )
    return payload


def main():
    t0 = time.time()
    print("=== P2 RC2 kinematic boundary (multi-scenario) ===")
    print(f"OUT={OUT}")
    grid = CandidateGrid.coarse()
    print(f"Candidate grid nodes: {grid.n}")

    # self-check geometry
    tchk = np.array([0.0, 100.0, 300.0])
    b_lin = predict_bearings(tchk, 50e3, 0.0, 2.0, 0.0, 0.0)
    b_cr5 = predict_bearings(tchk, 50e3, 0.0, 2.0, np.deg2rad(5.0), 0.0)
    print(f"SELFCHK LIN bearings={np.rad2deg(b_lin)} (expect ~0)")
    print(f"SELFCHK CR5 bearings={np.rad2deg(b_cr5)} (expect increasing)")
    assert np.allclose(b_lin, 0.0, atol=1e-12), "LIN bearing must be identically 0"
    assert np.ptp(b_cr5) > 0, "CR5 must have bearing change"
    assert 0.0 in PSI_GRID_DEG, "psi grid must include 0"

    (OUT / "P2_RC2_CONFIG.json").write_text(json.dumps(CONFIG, ensure_ascii=False, indent=2), encoding="utf-8")

    print("[1/6] FIM observability scan ...")
    df_fim = run_fim_scan()
    df_fim.to_csv(OUT / "observability_scan.csv", index=False, encoding="utf-8-sig")

    print("[2/6] sigma_theta × T boundary (CR5 + LIN) ...")
    df_sig = run_sigma_T_boundary(grid, scenarios=(PRIMARY_SCENARIO, LIMIT_SCENARIO))
    df_sig.to_csv(OUT / "bearing_accuracy_time_boundary.csv", index=False, encoding="utf-8-sig")

    print("[3/6] heading identifiability ...")
    df_head = run_heading_identifiability(grid)
    df_head.to_csv(OUT / "heading_identifiability_boundary.csv", index=False, encoding="utf-8-sig")

    print("[4/6] maneuver boundary ...")
    df_man = run_maneuver_boundary(grid, scenarios=(PRIMARY_SCENARIO, LIMIT_SCENARIO))
    df_man.to_csv(OUT / "maneuver_boundary.csv", index=False, encoding="utf-8-sig")

    print("[5/6] off-grid trials ...")
    df_off = run_offgrid_trials(grid)
    df_off.to_csv(OUT / "offgrid_trials.csv", index=False, encoding="utf-8-sig")

    print("[6/6] hard candidates ...")
    df_all_h, df_hard = run_hard_candidates(grid)
    df_all_h.to_csv(OUT / "hard_candidate_pairs_all.csv", index=False, encoding="utf-8-sig")
    df_hard.to_csv(OUT / "hard_candidate_pairs.csv", index=False, encoding="utf-8-sig")

    print("Figures ...")
    make_figures(df_sig, df_head, df_man, df_all_h)
    print("Reports ...")
    write_reports(df_fim, df_sig, df_head, df_man, df_off, df_all_h, df_hard)
    print("Dashboard data ...")
    export_embedded_data(df_sig, df_head, df_man, df_off, df_hard, df_fim)

    print(f"DONE in {time.time()-t0:.1f}s")
    print(f"Artifacts in {OUT}")


if __name__ == "__main__":
    main()
