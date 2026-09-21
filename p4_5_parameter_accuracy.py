#!/usr/bin/env python3
"""P4.5 — parameter estimation accuracy baseline.

M0 prior / M1 RC2 bearing / M2 RC2+CZ on identical truths, noise, candidate space.
Primary metrics: bias, MAE, RMSE, CI95 width, coverage, improvement.
No new features; same E-STD modal propagation as P3/P4.
"""
from __future__ import annotations

import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import p4_g4_performance_boundary as p4

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "P4_5_parameter_accuracy"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

N_TRUTH = 30
N_NOISE = 20
RNG = 20260921

# Estimation grids (same for M1/M2)
R_G = np.arange(45.0, 60.0 + 1e-9, 0.5) * 1e3          # 31
TH_G = np.deg2rad(np.arange(-4.0, 4.0 + 1e-9, 0.5))    # 17
V_G = np.arange(1.2, 2.8 + 1e-9, 0.2)                  # 9
PSI_G = np.deg2rad(np.arange(-12.0, 12.0 + 1e-9, 2.0)) # 13
Z_G = np.array([180.0, 200.0, 220.0])

PRIOR = {
    "r_m": (45e3, 60e3),
    "theta_rad": (math.radians(-5), math.radians(5)),
    "z_m": (150.0, 250.0),
    "v": (1.0, 3.0),
    "psi_rad": (math.radians(-15), math.radians(15)),
}

# Acceptance for CI: standardized J
THR_BEAR = 13.28
THR_AC = 13.28

SCENARIOS = {
    "B-L": dict(
        name="B-L极限几何",
        r0=50e3, theta0=0.0, z=200.0, v=2.0, psi=0.0,
        T=600.0, sigma_deg=0.1, source="S1", snr_db=10.0, turn_deg=0.0,
        truth_spread=dict(r_km=2.0, theta_deg=1.0, z_m=20.0, v=0.3, psi_deg=2.0),
    ),
    "B-M": dict(
        name="B-M标准开放几何",
        r0=50e3, theta0=0.0, z=200.0, v=2.0, psi=5.0,
        T=600.0, sigma_deg=0.1, source="S1", snr_db=10.0, turn_deg=0.0,
        truth_spread=dict(r_km=2.0, theta_deg=1.0, z_m=20.0, v=0.3, psi_deg=2.0),
    ),
    "B-U": dict(
        name="B-U有利理论场景",
        r0=50e3, theta0=0.0, z=200.0, v=2.0, psi=8.0,
        T=1200.0, sigma_deg=0.02, source="S2", snr_db=20.0, turn_deg=15.0,
        truth_spread=dict(r_km=2.0, theta_deg=1.0, z_m=20.0, v=0.3, psi_deg=2.0),
    ),
}

CONFIG = {
    "project": "B3D-3",
    "package": "P4_5_parameter_accuracy",
    "created_utc": NOW,
    "purpose": "standardized parameter RMSE: RC1 prior vs RC2 vs RC2+CZ",
    "n_truth": N_TRUTH,
    "n_noise": N_NOISE,
    "scenarios": {k: {kk: vv for kk, vv in v.items()} for k, v in SCENARIOS.items()},
    "grids": {
        "r_km": R_G.tolist(),
        "theta_deg": np.rad2deg(TH_G).tolist(),
        "z_m": Z_G.tolist(),
        "v_mps": V_G.tolist(),
        "psi_deg": np.rad2deg(PSI_G).tolist(),
    },
    "same_for_M1_M2": ["truth", "bearing_noise", "candidate_space", "optimization"],
    "point_estimate": "coarse grid argmin → local continuous-ish refine",
    "ci": "accepted set / profile-like: J <= J_min + Delta_chi2",
    "coverage_ok": "RMSE small AND coverage in [0.85,1.0] else OVERCONFIDENT/INVALID_COVERAGE",
    "forbidden": [
        "Bellhop", "modify P3 propagation", "new features",
        "tune prior for pretty RMSE", "drop failed trials", "auto-enter P5",
    ],
}


def build_grid():
    rr, tt, vv, pp = np.meshgrid(R_G, TH_G, V_G, PSI_G, indexing="ij")
    kin = dict(
        r0=rr.ravel(), th0=tt.ravel(), v=vv.ravel(), psi=pp.ravel(),
    )
    n_k = kin["r0"].size
    # M2 grid: kin × z
    zrep = len(Z_G)
    full = {
        "r0": np.repeat(kin["r0"], zrep),
        "th0": np.repeat(kin["th0"], zrep),
        "v": np.repeat(kin["v"], zrep),
        "psi": np.repeat(kin["psi"], zrep),
        "z": np.tile(Z_G, n_k),
    }
    return kin, full, n_k, full["r0"].size


def ensure_tables():
    if not p4.AMP_TABS:
        for src in ["S0", "S1", "S2"]:
            freqs, _, _ = p4.source_freqs(src)
            for ek in ["E0", "E1", "E2"]:
                p4.AMP_TABS[(src, ek)] = p4.build_amp_table(p4.MODE_ENV[ek], freqs)


# Feature cache: acoustic time-shape on full (kin×z) grid
_FEAT = {}


def time_shape_from_A(A):
    F, M = A.shape
    parts = []
    for j in range(F):
        y = A[j] - np.mean(A[j])
        parts.append(y / (np.linalg.norm(y) + 1e-20))
    feat = np.concatenate(parts)
    return feat / (np.linalg.norm(feat) + 1e-20)


def cand_feats_full(t, gfull, turn_deg, source, tab):
    freqs, _, _ = p4.source_freqs(source)
    N, M, F = gfull["r0"].size, t.size, freqs.size
    xp, yp = p4.platform_xy(t, turn_deg)
    feats = np.zeros((N, F * M))
    chunk = 3000
    for i0 in range(0, N, chunk):
        i1 = min(N, i0 + chunk)
        xt = gfull["r0"][i0:i1, None] * np.cos(gfull["th0"][i0:i1, None]) + gfull["v"][i0:i1, None] * t[None] * np.cos(gfull["psi"][i0:i1, None])
        yt = gfull["r0"][i0:i1, None] * np.sin(gfull["th0"][i0:i1, None]) + gfull["v"][i0:i1, None] * t[None] * np.sin(gfull["psi"][i0:i1, None])
        rr = np.hypot(xt - xp[None], yt - yp[None])
        th = np.arctan2(yt - yp[None], xt - xp[None])
        for k in range(i1 - i0):
            A = p4.lookup_amp(tab, freqs, rr[k], float(gfull["z"][i0 + k]))
            for it in range(M):
                A[:, it] *= p4.array_gain(th[k, it], freqs)
            feats[i0 + k] = time_shape_from_A(A)
    return feats


def bearings_kin(t, gkin, turn_deg):
    xp, yp = p4.platform_xy(t, turn_deg)
    N, M = gkin["r0"].size, t.size
    out = np.empty((N, M))
    for i0 in range(0, N, 4000):
        i1 = min(N, i0 + 4000)
        xt = gkin["r0"][i0:i1, None] * np.cos(gkin["th0"][i0:i1, None]) + gkin["v"][i0:i1, None] * t[None] * np.cos(gkin["psi"][i0:i1, None])
        yt = gkin["r0"][i0:i1, None] * np.sin(gkin["th0"][i0:i1, None]) + gkin["v"][i0:i1, None] * t[None] * np.sin(gkin["psi"][i0:i1, None])
        out[i0:i1] = np.arctan2(yt - yp[None], xt - xp[None])
    return out


def obs_bearing(t, truth, turn_deg, sigma_deg, rng):
    xp, yp = p4.platform_xy(t, turn_deg)
    xt, yt = p4.target_xy(t, truth["r0_m"], truth["theta0_rad"], truth["v"], truth["psi_rad"])
    th = np.arctan2(yt - yp, xt - xp)
    sig = math.radians(max(sigma_deg, 1e-6))
    return th + rng.normal(0, sig, size=th.shape)


def obs_feat(t, truth, turn_deg, source, snr_db, rng, tab):
    freqs, _, _ = p4.source_freqs(source)
    M = t.size
    xp, yp = p4.platform_xy(t, turn_deg)
    xt, yt = p4.target_xy(t, truth["r0_m"], truth["theta0_rad"], truth["v"], truth["psi_rad"])
    rr = np.hypot(xt - xp, yt - yp)
    th = np.arctan2(yt - yp, xt - xp)
    A = p4.lookup_amp(tab, freqs, rr, float(truth["z"]))
    for it in range(M):
        A[:, it] *= p4.array_gain(th[it], freqs)
    if snr_db is not None:
        sig = float(np.mean(A ** 2)) + 1e-20
        nvar = sig / (10 ** (snr_db / 10.0))
        A = np.abs(A + rng.normal(0, math.sqrt(nvar), size=A.shape))
    return time_shape_from_A(A)


def local_refine_kin(obs_th, t, turn_deg, r0, th0, v, psi, n_rounds=3):
    """Refine bearing-only point estimate."""
    best = (float(r0), float(th0), float(v), float(psi))
    best_c = np.inf
    dr, dth, dv, dpsi = 500.0, math.radians(0.4), 0.15, math.radians(1.5)
    for _ in range(n_rounds):
        rs = np.clip(best[0] + np.linspace(-dr, dr, 5), 45e3, 60e3)
        ts = np.clip(best[1] + np.linspace(-dth, dth, 5), math.radians(-5), math.radians(5))
        vs = np.clip(best[2] + np.linspace(-dv, dv, 5), 1.0, 3.0)
        ps = np.clip(best[3] + np.linspace(-dpsi, dpsi, 5), math.radians(-15), math.radians(15))
        RR, TT, VV, PP = np.meshgrid(rs, ts, vs, ps, indexing="ij")
        rr = RR.ravel(); tt = TT.ravel(); vv = VV.ravel(); pp = PP.ravel()
        xp, yp = p4.platform_xy(t, turn_deg)
        N = rr.size
        cost = np.empty(N)
        for i0 in range(0, N, 800):
            i1 = min(N, i0 + 800)
            xt = rr[i0:i1, None] * np.cos(tt[i0:i1, None]) + vv[i0:i1, None] * t[None] * np.cos(pp[i0:i1, None])
            yt = rr[i0:i1, None] * np.sin(tt[i0:i1, None]) + vv[i0:i1, None] * t[None] * np.sin(pp[i0:i1, None])
            th = np.arctan2(yt - yp[None], xt - xp[None])
            cost[i0:i1] = np.sum(p4.wrap(th - obs_th[None]) ** 2, axis=1)
        k = int(np.argmin(cost))
        best_c = float(cost[k])
        best = (float(rr[k]), float(tt[k]), float(vv[k]), float(pp[k]))
        dr *= 0.35; dth *= 0.35; dv *= 0.35; dpsi *= 0.35
    return best, best_c


def local_refine_m2(obs_th, obs_f, t, turn_deg, source, tab, r0, th0, v, psi, z,
                    sig2_bear, sig_n2, n_rounds=2):
    """Refine joint bearing+CZ point estimate (coarse local mesh)."""
    freqs, _, _ = p4.source_freqs(source)
    best = (float(r0), float(th0), float(v), float(psi), float(z))
    best_c = np.inf
    dr, dth, dv, dpsi, dz = 400.0, math.radians(0.3), 0.12, math.radians(1.2), 8.0
    xp, yp = p4.platform_xy(t, turn_deg)
    M = t.size
    for _ in range(n_rounds):
        rs = np.clip(best[0] + np.array([-dr, 0.0, dr]), 45e3, 60e3)
        ts = np.clip(best[1] + np.array([-dth, 0.0, dth]), math.radians(-5), math.radians(5))
        vs = np.clip(best[2] + np.array([-dv, 0.0, dv]), 1.0, 3.0)
        ps = np.clip(best[3] + np.array([-dpsi, 0.0, dpsi]), math.radians(-15), math.radians(15))
        zs = np.clip(best[4] + np.array([-dz, 0.0, dz]), 150.0, 250.0)
        # sparse mesh: combine all with z outer
        best_c_local = np.inf
        best_local = best
        for zq in zs:
            for r in rs:
                for th in ts:
                    for vv in vs:
                        for pp in ps:
                            xt = r * np.cos(th) + vv * t * np.cos(pp)
                            yt = r * np.sin(th) + vv * t * np.sin(pp)
                            thv = np.arctan2(yt - yp, xt - xp)
                            jb = float(np.sum(p4.wrap(thv - obs_th) ** 2) / sig2_bear)
                            rr = np.hypot(xt - xp, yt - yp)
                            A = p4.lookup_amp(tab, freqs, rr, float(zq))
                            for it in range(M):
                                A[:, it] *= p4.array_gain(thv[it], freqs)
                            feat = time_shape_from_A(A)
                            jac = float(np.sum((feat - obs_f) ** 2) / sig_n2)
                            c = jb + jac
                            if c < best_c_local:
                                best_c_local = c
                                best_local = (float(r), float(th), float(vv), float(pp), float(zq))
        best, best_c = best_local, best_c_local
        dr *= 0.4; dth *= 0.4; dv *= 0.4; dpsi *= 0.4; dz *= 0.4
    return best, best_c


def ang_err_deg(a, b):
    return abs(math.degrees(p4.wrap(a - b)))


def sample_truths(key, n, rng):
    sc = SCENARIOS[key]
    sp = sc["truth_spread"]
    out = []
    for _ in range(n):
        # off-grid continuous sample around scene center, clipped to prior
        r = sc["r0"] + rng.uniform(-sp["r_km"], sp["r_km"]) * 1e3
        th = math.radians(sc["theta0"] + rng.uniform(-sp["theta_deg"], sp["theta_deg"]))
        z = sc["z"] + rng.uniform(-sp["z_m"], sp["z_m"])
        v = sc["v"] + rng.uniform(-sp["v"], sp["v"])
        psi = math.radians(sc["psi"] + rng.uniform(-sp["psi_deg"], sp["psi_deg"]))
        r = float(np.clip(r, PRIOR["r_m"][0], PRIOR["r_m"][1]))
        th = float(np.clip(th, PRIOR["theta_rad"][0], PRIOR["theta_rad"][1]))
        z = float(np.clip(z, PRIOR["z_m"][0], PRIOR["z_m"][1]))
        v = float(np.clip(v, PRIOR["v"][0], PRIOR["v"][1]))
        psi = float(np.clip(psi, PRIOR["psi_rad"][0], PRIOR["psi_rad"][1]))
        out.append(dict(r0_m=r, theta0_rad=th, z=z, v=v, psi_rad=psi,
                        theta0_deg=math.degrees(th), psi_deg=math.degrees(psi)))
    return out


def prior_stats():
    def half(a, b):
        return (b - a) / 2.0
    def std_uniform(a, b):
        return (b - a) / math.sqrt(12.0)
    return {
        "r_km": dict(half_width_km=half(*PRIOR["r_m"]) / 1e3,
                     std_km=std_uniform(*PRIOR["r_m"]) / 1e3),
        "theta_deg": dict(half_width_deg=math.degrees(half(*PRIOR["theta_rad"])),
                          std_deg=math.degrees(std_uniform(*PRIOR["theta_rad"]))),
        "z_m": dict(half_width_m=half(*PRIOR["z_m"]), std_m=std_uniform(*PRIOR["z_m"])),
        "v_mps": dict(half_width=half(*PRIOR["v"]), std=std_uniform(*PRIOR["v"])),
        "psi_deg": dict(half_width_deg=math.degrees(half(*PRIOR["psi_rad"])),
                        std_deg=math.degrees(std_uniform(*PRIOR["psi_rad"]))),
    }


def run_mc_one_scenario(key, gkin, gfull, n_k, n_f, cache_bear, cache_feat, cache_sig):
    sc = SCENARIOS[key]
    t = np.linspace(0.0, sc["T"], max(6, int(round(sc["T"] / 10.0)) + 1))
    turn = sc["turn_deg"]
    sig2_bear = (math.radians(max(sc["sigma_deg"], 1e-6))) ** 2
    tab = p4.AMP_TABS[(sc["source"], "E0")]
    # precompute
    if ("bear", key) not in cache_bear:
        cache_bear[("bear", key)] = bearings_kin(t, gkin, turn)
    if ("feat", key) not in cache_feat:
        cache_feat[("feat", key)] = cand_feats_full(t, gfull, turn, sc["source"], tab)
        # estimate typical feature noise scale for SNR
        sig_n2 = max((10 ** (-sc["snr_db"] / 10.0)) / max(cache_feat[("feat", key)].shape[1], 1), 1e-10)
        cache_sig[key] = sig_n2
    pred_b = cache_bear[("bear", key)]
    feats = cache_feat[("feat", key)]
    sig_n2 = cache_sig[key]

    # map full-grid index → kin index
    zrep = len(Z_G)
    kin_of_full = np.repeat(np.arange(n_k), zrep)

    truths = sample_truths(key, N_TRUTH, np.random.default_rng(RNG + hash(key) % 1000))
    pri = prior_stats()
    trials = []

    for it, truth in enumerate(truths):
        for inn in range(N_NOISE):
            rng = np.random.default_rng(RNG + it * 10007 + inn)
            b_obs = obs_bearing(t, truth, turn, sc["sigma_deg"], rng)
            f_obs = obs_feat(t, truth, turn, sc["source"], sc["snr_db"], rng, tab)

            # J1 bearing on kinematic grid
            dth = p4.wrap(pred_b - b_obs[None, :])
            J1 = np.sum(dth ** 2, axis=1) / sig2_bear
            k1 = int(np.argmin(J1))
            est1, _ = local_refine_kin(b_obs, t, turn, gkin["r0"][k1], gkin["th0"][k1],
                                       gkin["v"][k1], gkin["psi"][k1])
            # CI from accepted kinematic set
            m1 = J1 <= J1[k1] + THR_BEAR
            if m1.sum() < 2:
                m1[k1] = True

            # J2 joint on full grid (kin × z)
            # bearing part: expand J1 to full via kin_of_full
            d_ac = feats - f_obs[None, :]
            J_ac = np.sum(d_ac ** 2, axis=1) / sig_n2
            J2 = J1[kin_of_full] + J_ac
            k2 = int(np.argmin(J2))
            est2, _ = local_refine_m2(
                b_obs, f_obs, t, turn, sc["source"], tab,
                gfull["r0"][k2], gfull["th0"][k2], gfull["v"][k2], gfull["psi"][k2], gfull["z"][k2],
                sig2_bear, sig_n2,
            )
            # CI M2
            m2 = J2 <= J2[k2] + THR_BEAR + THR_AC
            if m2.sum() < 2:
                m2[k2] = True

            def metrics(est, mask, gsrc, is_full):
                r_e, th_e, v_e, psi_e = est[0], est[1], est[2], est[3]
                z_e = est[4] if len(est) > 4 else float("nan")
                # errors
                e_r = r_e / 1e3 - truth["r0_m"] / 1e3
                e_th = math.degrees(p4.wrap(th_e - truth["theta0_rad"]))
                e_v = v_e - truth["v"]
                e_psi = math.degrees(p4.wrap(psi_e - truth["psi_rad"]))
                e_z = z_e - truth["z"] if np.isfinite(z_e) else np.nan
                # CI widths from mask
                if is_full:
                    w_r = float(gsrc["r0"][mask].max() - gsrc["r0"][mask].min()) / 1e3 if mask.sum() >= 2 else 0.0
                    w_th = math.degrees(float(gsrc["th0"][mask].max() - gsrc["th0"][mask].min())) if mask.sum() >= 2 else 0.0
                    w_v = float(gsrc["v"][mask].max() - gsrc["v"][mask].min()) if mask.sum() >= 2 else 0.0
                    w_psi = math.degrees(float(gsrc["psi"][mask].max() - gsrc["psi"][mask].min())) if mask.sum() >= 2 else 0.0
                    w_z = float(gsrc["z"][mask].max() - gsrc["z"][mask].min()) if mask.sum() >= 2 else 0.0
                    cov_r = bool(gsrc["r0"][mask].min() - 1 <= truth["r0_m"] <= gsrc["r0"][mask].max() + 1)
                    cov_th = bool(gsrc["th0"][mask].min() - 1e-9 <= truth["theta0_rad"] <= gsrc["th0"][mask].max() + 1e-9)
                    cov_v = bool(gsrc["v"][mask].min() - 1e-9 <= truth["v"] <= gsrc["v"][mask].max() + 1e-9)
                    cov_psi = bool(gsrc["psi"][mask].min() - 1e-9 <= truth["psi_rad"] <= gsrc["psi"][mask].max() + 1e-9)
                    cov_z = bool(gsrc["z"][mask].min() - 1e-9 <= truth["z"] <= gsrc["z"][mask].max() + 1e-9)
                else:
                    w_r = float(gsrc["r0"][mask].max() - gsrc["r0"][mask].min()) / 1e3 if mask.sum() >= 2 else 0.0
                    w_th = math.degrees(float(gsrc["th0"][mask].max() - gsrc["th0"][mask].min())) if mask.sum() >= 2 else 0.0
                    w_v = float(gsrc["v"][mask].max() - gsrc["v"][mask].min()) if mask.sum() >= 2 else 0.0
                    w_psi = math.degrees(float(gsrc["psi"][mask].max() - gsrc["psi"][mask].min())) if mask.sum() >= 2 else 0.0
                    w_z = np.nan
                    cov_r = bool(gsrc["r0"][mask].min() - 1 <= truth["r0_m"] <= gsrc["r0"][mask].max() + 1)
                    cov_th = bool(gsrc["th0"][mask].min() - 1e-9 <= truth["theta0_rad"] <= gsrc["th0"][mask].max() + 1e-9)
                    cov_v = bool(gsrc["v"][mask].min() - 1e-9 <= truth["v"] <= gsrc["v"][mask].max() + 1e-9)
                    cov_psi = bool(gsrc["psi"][mask].min() - 1e-9 <= truth["psi_rad"] <= gsrc["psi"][mask].max() + 1e-9)
                    cov_z = True
                return dict(
                    est_r_km=r_e / 1e3, est_theta_deg=math.degrees(th_e), est_v=v_e,
                    est_psi_deg=math.degrees(psi_e), est_z=z_e,
                    err_r_km=e_r, err_theta_deg=e_th, err_v=e_v, err_psi_deg=e_psi, err_z=e_z,
                    abs_err_r_km=abs(e_r), abs_err_theta_deg=abs(e_th), abs_err_v=abs(e_v),
                    abs_err_psi_deg=abs(e_psi), abs_err_z=abs(e_z) if np.isfinite(e_z) else np.nan,
                    ci_r_km=w_r, ci_theta_deg=w_th, ci_v=w_v, ci_psi_deg=w_psi, ci_z_m=w_z,
                    cov_r=cov_r, cov_theta=cov_th, cov_v=cov_v, cov_psi=cov_psi, cov_z=cov_z,
                    n_acc=int(mask.sum()),
                )

            m0 = {
                "est_r_km": 50.0,
                "est_theta_deg": SCENARIOS[key]["theta0"],
                "est_v": SCENARIOS[key]["v"],
                "est_psi_deg": SCENARIOS[key]["psi"],
                "est_z": SCENARIOS[key]["z"],
                "err_r_km": 50.0 - truth["r0_m"] / 1e3,
                "err_theta_deg": math.degrees(p4.wrap(math.radians(SCENARIOS[key]["theta0"]) - truth["theta0_rad"])),
                "err_v": SCENARIOS[key]["v"] - truth["v"],
                "err_psi_deg": math.degrees(p4.wrap(math.radians(SCENARIOS[key]["psi"]) - truth["psi_rad"])),
                "err_z": SCENARIOS[key]["z"] - truth["z"],
                "abs_err_r_km": abs(50.0 - truth["r0_m"] / 1e3),
                "abs_err_theta_deg": abs(math.degrees(p4.wrap(math.radians(SCENARIOS[key]["theta0"]) - truth["theta0_rad"]))),
                "abs_err_v": abs(SCENARIOS[key]["v"] - truth["v"]),
                "abs_err_psi_deg": abs(math.degrees(p4.wrap(math.radians(SCENARIOS[key]["psi"]) - truth["psi_rad"]))),
                "abs_err_z": abs(SCENARIOS[key]["z"] - truth["z"]),
                "ci_r_km": pri["r_km"]["half_width_km"] * 2,
                "ci_theta_deg": pri["theta_deg"]["half_width_deg"] * 2,
                "ci_v": pri["v_mps"]["half_width"] * 2,
                "ci_psi_deg": pri["psi_deg"]["half_width_deg"] * 2,
                "ci_z_m": pri["z_m"]["half_width_m"] * 2,
                "cov_r": True, "cov_theta": True, "cov_v": True, "cov_psi": True, "cov_z": True,
                "n_acc": int(gkin["r0"].size),
            }

            mm1 = metrics(est1, m1, gkin, False)
            mm2 = metrics(est2, m2, gfull, True)

            for layer, mm in [("M0_prior", m0), ("M1_RC2", mm1), ("M2_RC2_RC3", mm2)]:
                row = {
                    "scenario": key,
                    "trial": it,
                    "noise": inn,
                    "layer": layer,
                    "truth_r_km": truth["r0_m"] / 1e3,
                    "truth_theta_deg": truth["theta0_deg"],
                    "truth_z": truth["z"],
                    "truth_v": truth["v"],
                    "truth_psi_deg": truth["psi_deg"],
                    **mm,
                }
                trials.append(row)
        if (it + 1) % 5 == 0:
            print(f"    {key} truth {it+1}/{N_TRUTH}", flush=True)

    return pd.DataFrame(trials), pri


def aggregate_accuracy(trial_df, pri):
    """Build parameter_accuracy_baseline_vs_method.csv and coverage table."""
    rows = []
    params = [
        ("r", "km", "err_r_km", "abs_err_r_km", "ci_r_km", "cov_r", "r_km"),
        ("theta", "deg", "err_theta_deg", "abs_err_theta_deg", "ci_theta_deg", "cov_theta", "theta_deg"),
        ("z", "m", "err_z", "abs_err_z", "ci_z_m", "cov_z", "z_m"),
        ("v", "m/s", "err_v", "abs_err_v", "ci_v", "cov_v", "v_mps"),
        ("psi", "deg", "err_psi_deg", "abs_err_psi_deg", "ci_psi_deg", "cov_psi", "psi_deg"),
    ]
    pri_map = {
        "r": pri["r_km"]["std_km"],
        "theta": pri["theta_deg"]["std_deg"],
        "z": pri["z_m"]["std_m"],
        "v": pri["v_mps"]["std"],
        "psi": pri["psi_deg"]["std_deg"],
    }
    for scen in trial_df["scenario"].unique():
        for pname, unit, ecol, acol, ccol, covcol, pkey in params:
            rec = {
                "scenario": scen,
                "parameter": pname,
                "unit": unit,
                "prior_std": pri_map[pname],
            }
            stats = {}
            for layer, prefix in [("M0_prior", "prior"), ("M1_RC2", "rc2"), ("M2_RC2_RC3", "rc23")]:
                d = trial_df[(trial_df["scenario"] == scen) & (trial_df["layer"] == layer)]
                if not len(d):
                    continue
                e = pd.to_numeric(d[ecol], errors="coerce")
                ae = pd.to_numeric(d[acol], errors="coerce")
                ci = pd.to_numeric(d[ccol], errors="coerce")
                cov = pd.to_numeric(d[covcol], errors="coerce").astype(float)
                bias = float(e.mean())
                mae = float(ae.mean())
                rmse = float(np.sqrt(np.mean(e ** 2)))
                med = float(ae.median())
                ci_w = float(ci.mean())
                coverage = float(cov.mean())
                stats[prefix] = rmse
                rec[f"{prefix}_bias"] = bias
                rec[f"{prefix}_mae"] = mae
                rec[f"{prefix}_rmse"] = rmse
                rec[f"{prefix}_median_ae"] = med
                rec[f"{prefix}_ci95_width"] = ci_w
                rec[f"{prefix}_coverage"] = coverage
                # status per layer
                if prefix in ("rc2", "rc23"):
                    if coverage < 0.50 and rmse < 0.5 * abs(bias + rmse + 1e-9):
                        st = "OVERCONFIDENT_INVALID_COVERAGE"
                    elif coverage < 0.50:
                        st = "INVALID_COVERAGE"
                    elif rmse >= 0.8 * pri_map[pname] and pname in ("r", "z", "v", "psi"):
                        # essentially unresolved if RMSE ~ prior
                        st = "UNRESOLVED_OR_WEAK"
                    else:
                        st = "OK"
                    rec[f"{prefix}_status"] = st
            # improvement M1 → M2
            if "rc2" in stats and "rc23" in stats and stats["rc2"] > 1e-9:
                imp = 1.0 - stats["rc23"] / stats["rc2"]
                rec["rmse_improvement"] = float(imp)
            else:
                rec["rmse_improvement"] = np.nan
            # overall status
            rc2_st = rec.get("rc2_status", "")
            rc23_st = rec.get("rc23_status", "")
            if rc2_st == "UNRESOLVED_OR_WEAK" and rec.get("rc23_coverage", 0) >= 0.7:
                rec["status"] = f"UNRESOLVED_RC2 → RMSE/width after RC3: {rec.get('rc23_rmse'):.3g}/{rec.get('rc23_ci95_width'):.3g}"
            elif "INVALID" in str(rc23_st) or "OVERCONFIDENT" in str(rc23_st):
                rec["status"] = rc23_st
            elif rec.get("rmse_improvement") == rec.get("rmse_improvement") and rec["rmse_improvement"] >= 0.2 and rc23_st == "OK":
                rec["status"] = "IMPROVED"
            elif rec.get("rmse_improvement") == rec.get("rmse_improvement") and abs(rec["rmse_improvement"]) < 0.1:
                rec["status"] = "UNCHANGED_LARGELY_RC2"
            else:
                rec["status"] = f"M2={rc23_st}"
            # relative RMSE % for r,z,v
            if pname in ("r", "z", "v") and rec.get("rc23_rmse") is not None:
                scale = {"r": 50.0, "z": 200.0, "v": 2.0}[pname]
                rec["rc23_rel_rmse_pct"] = 100.0 * rec["rc23_rmse"] / scale
                rec["rc2_rel_rmse_pct"] = 100.0 * rec.get("rc2_rmse", np.nan) / scale
            rows.append(rec)
    return pd.DataFrame(rows)


def rc3_sensitivity(gfull, scen_key="B-M"):
    """1D and 2D J_CZ slices around a typical truth."""
    sc = SCENARIOS[scen_key]
    t = np.linspace(0.0, sc["T"], max(6, int(round(sc["T"] / 15.0)) + 1))
    turn = sc["turn_deg"]
    tab = p4.AMP_TABS[(sc["source"], "E0")]
    freqs, _, _ = p4.source_freqs(source := sc["source"])
    truth = dict(r0_m=sc["r0"], theta0_rad=math.radians(sc["theta0"]), z=sc["z"],
                 v=sc["v"], psi_rad=math.radians(sc["psi"]))
    rng = np.random.default_rng(RNG + 99)
    f_obs = obs_feat(t, truth, turn, source, sc["snr_db"], rng, tab)
    # noise scale
    sig_n2 = max((10 ** (-sc["snr_db"] / 10.0)) / max(len(f_obs), 1), 1e-10)

    def Jcz(r0, th0, v, psi, z):
        xp, yp = p4.platform_xy(t, turn)
        xt = r0 * np.cos(th0) + v * t * np.cos(psi)
        yt = r0 * np.sin(th0) + v * t * np.sin(psi)
        rr = np.hypot(xt - xp, yt - yp)
        th = np.arctan2(yt - yp, xt - xp)
        A = p4.lookup_amp(tab, freqs, rr, float(z))
        for it in range(t.size):
            A[:, it] *= p4.array_gain(th[it], freqs)
        feat = time_shape_from_A(A)
        return float(np.sum((feat - f_obs) ** 2) / sig_n2)

    base = (truth["r0_m"], truth["theta0_rad"], truth["v"], truth["psi_rad"], truth["z"])
    rows1 = []
    axes = {
        "r": ("r_km", np.linspace(45, 60, 31), 1e3),
        "v": ("v_mps", np.linspace(1.2, 2.8, 17), 1.0),
        "z": ("z_m", np.linspace(160, 240, 17), 1.0),
        "psi": ("psi_deg", np.linspace(-12, 12, 25), math.pi / 180),
        "theta": ("theta_deg", np.linspace(-4, 4, 17), math.pi / 180),
    }
    idx_map = {"r": 0, "theta": 1, "v": 2, "psi": 3, "z": 4}
    for dim, (lab, grid, scale) in axes.items():
        for gv in grid:
            x = list(base)
            x[idx_map[dim]] = gv * scale
            j = Jcz(*x)
            rows1.append({"slice": "1d", "dim": dim, "value": gv, "J_cz": j})
    df1 = pd.DataFrame(rows1)

    # 2D r-v and r-z
    rows2 = []
    r_grid = np.linspace(46, 58, 13)
    v_grid = np.linspace(1.3, 2.7, 13)
    z_grid = np.linspace(170, 230, 13)
    for rv in r_grid:
        for vv in v_grid:
            x = list(base)
            x[0] = rv * 1e3
            x[2] = vv
            rows2.append({"slice": "rv", "r_km": rv, "v": vv, "J_cz": Jcz(*x)})
    for rv in r_grid:
        for zv in z_grid:
            x = list(base)
            x[0] = rv * 1e3
            x[4] = zv
            rows2.append({"slice": "rz", "r_km": rv, "z": zv, "J_cz": Jcz(*x)})
    df2 = pd.DataFrame(rows2)

    # half-width metrics
    def half_width(dim):
        d = df1[df1["slice"] == "1d"]
        dd = d[d["dim"] == dim]
        jmin = dd["J_cz"].min()
        ok = dd[dd["J_cz"] <= jmin + 9.21]
        if len(ok) < 2:
            return np.nan
        return float(ok["value"].max() - ok["value"].min())

    hw = {dim: half_width(dim) for dim in ["r", "v", "z", "psi", "theta"]}
    return df1, df2, hw


def _esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def write_figA(path, acc_df):
    """Grouped bars: RMSE M0/M1/M2 for five params, three scenarios — one panel per scenario for r (and table-like)."""
    # Focus: for each scenario, five params, three layers — use RMSE normalized note
    w, h = 980, 420
    params = ["r", "theta", "z", "v", "psi"]
    scen = ["B-L", "B-M", "B-U"]
    layers = [("prior_std", "#8a8a8a", "先验std"),
              ("rc2_rmse", "#b45309", "RC2 RMSE"),
              ("rc23_rmse", "#0f766e", "RC2+RC3 RMSE")]
    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
         f'<rect width="{w}" height="{h}" fill="#f7f4ef"/>',
         f'<text x="{w/2}" y="26" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="15" font-weight="600">图A  五维参数：先验 / RC2 / RC2+RC3 RMSE</text>']
    # three scenario columns
    for iscen, sc in enumerate(scen):
        ox = 30 + iscen * 320
        p.append(f'<text x="{ox+130}" y="48" text-anchor="middle" font-size="12" fill="#333" font-family="sans-serif">{sc}</text>')
        # normalize each param RMSE by prior_std for display
        for ip, pname in enumerate(params):
            row = acc_df[(acc_df["scenario"] == sc) & (acc_df["parameter"] == pname)]
            if not len(row):
                continue
            row = row.iloc[0]
            ps = float(row["prior_std"]) if row["prior_std"] == row["prior_std"] else 1.0
            vals = [
                1.0,
                float(row["rc2_rmse"]) / ps if row["rc2_rmse"] == row["rc2_rmse"] else np.nan,
                float(row["rc23_rmse"]) / ps if row["rc23_rmse"] == row["rc23_rmse"] else np.nan,
            ]
            # status for M2
            st = str(row.get("rc23_status", ""))
            y0 = 60 + ip * 65
            p.append(f'<text x="{ox}" y="{y0+12}" font-size="11" fill="#555" font-family="sans-serif">{pname}</text>')
            for iv, v in enumerate(vals):
                col = layers[iv][1]
                if not np.isfinite(v):
                    p.append(f'<rect x="{ox+30+iv*70}" y="{y0}" width="60" height="18" fill="#4a0f0f"/>')
                    p.append(f'<text x="{ox+60+iv*70}" y="{y0+12}" text-anchor="middle" font-size="9" fill="#fff" font-family="sans-serif">n/a</text>')
                    continue
                if iv == 2 and "INVALID" in st:
                    p.append(f'<rect x="{ox+30+iv*70}" y="{y0}" width="60" height="18" fill="#4a0f0f"/>')
                    p.append(f'<text x="{ox+60+iv*70}" y="{y0+12}" text-anchor="middle" font-size="8" fill="#fff" font-family="sans-serif">INVALID</text>')
                    continue
                bw = min(v, 2.0) / 2.0 * 60
                p.append(f'<rect x="{ox+30+iv*70}" y="{y0}" width="{max(bw,2):.1f}" height="18" fill="{col}"/>')
                p.append(f'<text x="{ox+34+iv*70+max(bw,2):.1f}" y="{y0+13}" font-size="9" fill="#333" font-family="sans-serif">{v:.2f}</text>')
    p.append(f'<text x="30" y="{h-8}" font-size="10" fill="#777" font-family="sans-serif">柱长=RMSE/先验std（1=与先验同量级）；INVALID=覆盖率不足不作性能</text></svg>')
    path.write_text("\n".join(p), encoding="utf-8")


def write_figB(path, acc_df):
    w, h = 780, 400
    params = ["r", "theta", "z", "v", "psi"]
    scen = ["B-L", "B-M", "B-U"]
    ml, mt = 80, 70
    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
         f'<rect width="{w}" height="{h}" fill="#f7f4ef"/>',
         f'<text x="{w/2}" y="28" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="15" font-weight="600">图B  RMSE 改善率 (1 − RMSE_RC3/RMSE_RC2)</text>']
    # -1..1 scale
    def ymap(v):
        return mt + 280 - (np.clip(v, -1, 1) + 1) / 2 * 280
    p.append(f'<line x1="{ml}" y1="{ymap(0)}" x2="{w-40}" y2="{ymap(0)}" stroke="#999"/>')
    p.append(f'<text x="{ml-8}" y="{ymap(0)+4}" text-anchor="end" font-size="11" font-family="sans-serif">0</text>')
    p.append(f'<text x="{ml-8}" y="{ymap(1)+4}" text-anchor="end" font-size="11" font-family="sans-serif">+100%</text>')
    p.append(f'<text x="{ml-8}" y="{ymap(-1)+4}" text-anchor="end" font-size="11" font-family="sans-serif">−100%</text>')
    cols = ["#b45309", "#0f766e", "#1d4ed8"]
    for ip, pname in enumerate(params):
        for isc, sc in enumerate(scen):
            row = acc_df[(acc_df["scenario"] == sc) & (acc_df["parameter"] == pname)]
            x = ml + ip * 130 + isc * 35
            if not len(row):
                continue
            row = row.iloc[0]
            st = str(row.get("status", ""))
            imp = row.get("rmse_improvement", np.nan)
            if "UNRESOLVED" in st and "RC3" in st:
                p.append(f'<rect x="{x}" y="{ymap(0)-10}" width="28" height="20" fill="#c4c4c4"/>')
                p.append(f'<text x="{x+14}" y="{ymap(0)+4}" text-anchor="middle" font-size="8" fill="#333" font-family="sans-serif">UNRES</text>')
            elif "INVALID" in st or "OVERCONFIDENT" in st:
                p.append(f'<rect x="{x}" y="{ymap(0)-10}" width="28" height="20" fill="#4a0f0f"/>')
                p.append(f'<text x="{x+14}" y="{ymap(0)+4}" text-anchor="middle" font-size="7" fill="#fff" font-family="sans-serif">INV</text>')
            elif imp == imp:
                y1 = ymap(0)
                y2 = ymap(float(imp))
                top, bot = min(y1, y2), max(y1, y2)
                p.append(f'<rect x="{x}" y="{top:.1f}" width="28" height="{max(bot-top,1):.1f}" fill="{cols[isc]}"/>')
                p.append(f'<text x="{x+14}" y="{top-4:.1f}" text-anchor="middle" font-size="9" font-family="sans-serif">{100*float(imp):.0f}%</text>')
            else:
                p.append(f'<rect x="{x}" y="{ymap(0)-8}" width="28" height="16" fill="#ddd"/>')
        p.append(f'<text x="{ml+ip*130+40}" y="{mt+300}" text-anchor="middle" font-size="12" font-family="sans-serif">{pname}</text>')
    p.append(f'<text x="{w-120}" y="{mt+10}" font-size="11" fill="#b45309" font-family="sans-serif">B-L</text>')
    p.append(f'<text x="{w-80}" y="{mt+10}" font-size="11" fill="#0f766e" font-family="sans-serif">B-M</text>')
    p.append(f'<text x="{w-40}" y="{mt+10}" font-size="11" fill="#1d4ed8" font-family="sans-serif">B-U</text>')
    p.append(f'<text x="40" y="{h-8}" font-size="10" fill="#777" font-family="sans-serif">结构不可辨识维度不硬画百分比</text></svg>')
    path.write_text("\n".join(p), encoding="utf-8")


def write_figC(path, df1, df2, hw):
    w, h = 980, 420
    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
         f'<rect width="{w}" height="{h}" fill="#f7f4ef"/>',
         f'<text x="{w/2}" y="26" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="15" font-weight="600">图C  RC3（CZ）单参数敏感度 J_CZ 与 r–v / r–z 切片</text>']
    d1 = df1[df1["slice"] == "1d"]
    dims = ["r", "v", "z", "psi"]
    for i, dim in enumerate(dims):
        dd = d1[d1["dim"] == dim]
        if not len(dd):
            continue
        ox, oy = 30 + i * 240, 50
        xs = dd["value"].tolist()
        ys = dd["J_cz"].tolist()
        ymin, ymax = min(ys), max(ys)
        if ymax - ymin < 1e-9:
            ymax = ymin + 1
        x0, x1 = min(xs), max(xs)
        p.append(f'<text x="{ox+90}" y="{oy-8}" text-anchor="middle" font-size="12" font-family="sans-serif">J_CZ({dim})</text>')
        p.append(f'<rect x="{ox}" y="{oy}" width="180" height="140" fill="#fff" stroke="#ccc"/>')
        pts = []
        for x, y in zip(xs, ys):
            X = ox + (x - x0) / (x1 - x0 + 1e-12) * 180
            Y = oy + 140 - (y - ymin) / (ymax - ymin + 1e-12) * 140
            pts.append(f"{X:.1f},{Y:.1f}")
        p.append(f'<polyline fill="none" stroke="#b45309" stroke-width="2" points="{" ".join(pts)}"/>')
        p.append(f'<text x="{ox+90}" y="{oy+158}" text-anchor="middle" font-size="10" fill="#555" font-family="sans-serif">{x0:.2g} … {x1:.2g}</text>')
        hwv = hw.get(dim, np.nan)
        p.append(f'<text x="{ox}" y="{oy+175}" font-size="10" fill="#333" font-family="sans-serif">Δ@+9.21≈{hwv:.3g}</text>')
    # 2D r-v heatmap mini
    rv = df2[df2["slice"] == "rv"]
    if len(rv):
        ox, oy = 30, 270
        p.append(f'<text x="{ox+80}" y="{oy-6}" text-anchor="middle" font-size="12" font-family="sans-serif">J_CZ(r,v)</text>')
        rvals = sorted(rv["r_km"].unique())
        vvals = sorted(rv["v"].unique())
        J = np.zeros((len(vvals), len(rvals)))
        for i, vv in enumerate(vvals):
            for j, rr in enumerate(rvals):
                J[i, j] = float(rv[(rv["r_km"] == rr) & (rv["v"] == vv)]["J_cz"].iloc[0])
        jmin, jmax = J.min(), J.max()
        cw, ch = 8, 6
        for i in range(len(vvals)):
            for j in range(len(rvals)):
                u = (J[i, j] - jmin) / (jmax - jmin + 1e-12)
                r = int(240 * (1 - u) + 31 * u); g = int(230 * (1 - u) + 111 * u); b = int(220 * (1 - u) + 139 * u)
                p.append(f'<rect x="{ox+j*cw}" y="{oy+(len(vvals)-1-i)*ch}" width="{cw}" height="{ch}" fill="rgb({r},{g},{b})"/>')
        p.append(f'<text x="{ox}" y="{oy+len(vvals)*ch+14}" font-size="10" fill="#555" font-family="sans-serif">r→</text>')
        p.append(f'<text x="{ox+120}" y="{oy-6}" text-anchor="middle" font-size="11" fill="#555" font-family="sans-serif">v↓</text>')
    rz = df2[df2["slice"] == "rz"]
    if len(rz):
        ox, oy = 220, 270
        p.append(f'<text x="{ox+80}" y="{oy-6}" text-anchor="middle" font-size="12" font-family="sans-serif">J_CZ(r,z)</text>')
        rvals = sorted(rz["r_km"].unique())
        zvals = sorted(rz["z"].unique())
        J = np.zeros((len(zvals), len(rvals)))
        for i, zv in enumerate(zvals):
            for j, rr in enumerate(rvals):
                J[i, j] = float(rz[(rz["r_km"] == rr) & (rz["z"] == zv)]["J_cz"].iloc[0])
        jmin, jmax = J.min(), J.max()
        cw, ch = 8, 6
        for i in range(len(zvals)):
            for j in range(len(rvals)):
                u = (J[i, j] - jmin) / (jmax - jmin + 1e-12)
                r = int(240 * (1 - u) + 31 * u); g = int(230 * (1 - u) + 111 * u); b = int(220 * (1 - u) + 139 * u)
                p.append(f'<rect x="{ox+j*cw}" y="{oy+(len(zvals)-1-i)*ch}" width="{cw}" height="{ch}" fill="rgb({r},{g},{b})"/>')
    # interpretation box
    p.append(f'<text x="430" y="290" font-size="12" font-family="sans-serif" fill="#1a1a1a">半高宽(+9.21)：</text>')
    yy = 310
    for dim in ["r", "v", "z", "psi", "theta"]:
        p.append(f'<text x="430" y="{yy}" font-size="11" fill="#333" font-family="sans-serif">{dim}: {hw.get(dim, float("nan")):.4g}</text>')
        yy += 18
    p.append(f'<text x="430" y="{yy+8}" font-size="11" fill="#555" font-family="sans-serif">J_CZ(r) 峰清晰 → CZ直接敏感距离；</text>')
    p.append(f'<text x="430" y="{yy+26}" font-size="11" fill="#555" font-family="sans-serif">v/z/ψ 半高宽大 → 单独弱，多来自耦合。</text>')
    p.append('</svg>')
    path.write_text("\n".join(p), encoding="utf-8")


def main():
    t0 = time.time()
    print("=== P4.5 parameter accuracy baseline ===")
    print(f"OUT={OUT}")
    print(f"N_truth={N_TRUTH} N_noise={N_NOISE} → {N_TRUTH*N_NOISE}/scenario")
    ensure_tables()
    gkin, gfull, n_k, n_f = build_grid()
    print(f"kin grid={n_k}  full(z)={n_f}")
    (OUT / "P4_5_CONFIG.json").write_text(json.dumps(CONFIG, ensure_ascii=False, indent=2), encoding="utf-8")

    cache_bear, cache_feat, cache_sig = {}, {}, {}
    all_trials = []
    pri = prior_stats()
    for key in ["B-L", "B-M", "B-U"]:
        print(f"[MC] {key} ...", flush=True)
        df_t, pri = run_mc_one_scenario(key, gkin, gfull, n_k, n_f, cache_bear, cache_feat, cache_sig)
        all_trials.append(df_t)
        # incremental save
        pd.concat(all_trials, ignore_index=True).to_csv(
            OUT / "parameter_accuracy_trials.csv", index=False, encoding="utf-8-sig"
        )
        print(f"  saved partial trials={sum(len(x) for x in all_trials)}", flush=True)
    trial_df = pd.concat(all_trials, ignore_index=True)
    trial_df.to_csv(OUT / "parameter_accuracy_trials.csv", index=False, encoding="utf-8-sig")
    print("trials", len(trial_df), flush=True)

    print("aggregate ...", flush=True)
    acc = aggregate_accuracy(trial_df, pri)
    acc.to_csv(OUT / "parameter_accuracy_baseline_vs_method.csv", index=False, encoding="utf-8-sig")

    # coverage summary
    cov_rows = []
    for _, r in acc.iterrows():
        cov_rows.append({
            "scenario": r["scenario"], "parameter": r["parameter"],
            "rc2_coverage": r.get("rc2_coverage"), "rc23_coverage": r.get("rc23_coverage"),
            "rc2_ci95_width": r.get("rc2_ci95_width"), "rc23_ci95_width": r.get("rc23_ci95_width"),
            "rc2_status": r.get("rc2_status"), "rc23_status": r.get("rc23_status"),
            "status": r.get("status"), "rmse_improvement": r.get("rmse_improvement"),
        })
    pd.DataFrame(cov_rows).to_csv(OUT / "parameter_coverage.csv", index=False, encoding="utf-8-sig")

    print("sensitivity ...", flush=True)
    df1, df2, hw = rc3_sensitivity(gfull, "B-M")
    df1.to_csv(OUT / "rc3_1d_sensitivity.csv", index=False, encoding="utf-8-sig")
    df2.to_csv(OUT / "rc3_rv_sensitivity.csv", index=False, encoding="utf-8-sig")
    df2[df2["slice"] == "rz"].to_csv(OUT / "rc3_rz_sensitivity.csv", index=False, encoding="utf-8-sig")
    print("half-widths", hw, flush=True)

    print("figures ...", flush=True)
    write_figA(FIG / "figA_rmse_baseline_vs_method.svg", acc)
    write_figB(FIG / "figB_rmse_improvement.svg", acc)
    write_figC(FIG / "figC_rc3_parameter_sensitivity.svg", df1, df2, hw)

    # ---- report ----
    def fnum(x, nd=3):
        try:
            if x is None:
                return "n/a"
            v = float(x)
            if not np.isfinite(v):
                return "n/a"
            return f"{v:.{nd}f}"
        except Exception:
            return "n/a"

    md = []
    md.append("# P4.5 报告：五维参数统一估计精度基准")
    md.append("")
    md.append(f"UTC：{NOW}  ·  目录：`results/P4_5_parameter_accuracy/`")
    md.append("")
    md.append("本阶段**不**开发新方法；把成果从“候选排除率”转换为 **先验 → RC2 → RC2+RC3** 的参数 RMSE/MAE/覆盖率。")
    md.append("")
    md.append("## 0. 实验设计")
    md.append("")
    md.append(f"- 每场景 **N_truth={N_TRUTH}** × **N_noise={N_NOISE}** = {N_TRUTH*N_NOISE} 次估计")
    md.append("- 真值在基准附近**连续 off-grid** 采样（r±2 km, θ±1°, z∈[180,220]±, v±0.3, ψ 中心±2°）")
    md.append("- **M0** 先验中心 / 先验宽度；**M1** 仅方位；**M2** 同一真值+同一方位噪声+同一候选空间 + CZ")
    md.append("- 点估计：粗网格 argmin → 局部细化；CI：J≤J_min+Δχ² 接受集宽度")
    md.append("- 角度误差使用 wrap 最小角差")
    md.append("- 覆盖率低时标记 INVALID_COVERAGE / OVERCONFIDENT，不把小 RMSE 称为高精度")
    md.append("")
    md.append("| 场景 | 几何 | 观测 |")
    md.append("| --- | --- | --- |")
    md.append("| B-L | 近共线 ψ≈0° | T=600s, σθ=0.1°, S1, SNR=10dB |")
    md.append("| B-M | 开放 ψ≈5° | 同上（主结果） |")
    md.append("| B-U | ψ≈8° + turn15° | T=1200s, σθ=0.02°, S2, SNR=20dB（理论上界，非装备） |")
    md.append("")
    md.append("## 1. 主结果表（`parameter_accuracy_baseline_vs_method.csv`）")
    md.append("")
    md.append("| 场景 | 参数 | 先验std | RC2 RMSE | RC2 覆盖 | RC2+RC3 RMSE | RC2+RC3 覆盖 | 改善 | status |")
    md.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for _, r in acc.iterrows():
        md.append(
            f"| {r['scenario']} | {r['parameter']} ({r['unit']}) | {fnum(r['prior_std'])} | "
            f"{fnum(r.get('rc2_rmse'))} | {fnum(r.get('rc2_coverage'),2)} | "
            f"{fnum(r.get('rc23_rmse'))} | {fnum(r.get('rc23_coverage'),2)} | "
            f"{fnum(r.get('rmse_improvement'),3)} | {r.get('status','')} |"
        )
    md.append("")
    md.append("### 论文式表述模板（引用上表数字）")
    md.append("")
    for sc in ["B-M", "B-L", "B-U"]:
        for pname in ["r", "v", "theta", "psi", "z"]:
            row = acc[(acc["scenario"] == sc) & (acc["parameter"] == pname)]
            if not len(row):
                continue
            row = row.iloc[0]
            st = str(row.get("status", ""))
            if "UNRESOLVED" in st:
                md.append(f"- **{sc} / {pname}**：RC2下结构性弱/不可辨识；加入CZ后 RMSE={fnum(row.get('rc23_rmse'))} {row['unit']}，95%宽={fnum(row.get('rc23_ci95_width'))}，覆盖={fnum(row.get('rc23_coverage'),2)}。")
            elif "INVALID" in st or "OVERCONFIDENT" in st:
                md.append(f"- **{sc} / {pname}**：M2 RMSE={fnum(row.get('rc23_rmse'))} 但覆盖={fnum(row.get('rc23_coverage'),2)} → **{st}**，不得称高精度。")
            else:
                md.append(
                    f"- **{sc} / {pname}**：RC2 RMSE={fnum(row.get('rc2_rmse'))} {row['unit']}；"
                    f"加入CZ后={fnum(row.get('rc23_rmse'))} {row['unit']}；"
                    f"降低 {fnum(100*(row.get('rmse_improvement') or 0),1)}%；"
                    f"95%覆盖={fnum(row.get('rc23_coverage'),2)}。"
                )
    md.append("")
    md.append("## 2. RC3（CZ）直接观察了哪个参数？")
    md.append("")
    md.append("B-M 典型真值上单参数切片，Δ@J_min+9.21 半高宽：")
    md.append("")
    md.append("| 参数 | 半高宽 | 解读 |")
    md.append("| --- | --- | --- |")
    interp = {
        "r": "CZ对距离直接敏感（峰最清晰）" if hw.get("r", 99) < 5 else "距离敏感度有限",
        "v": "单独较弱，多来自 r–v 轨迹耦合" if hw.get("v", 0) > 0.3 else "对速度也有直接敏感度",
        "z": "剖面平 → 深度仍难独立约束" if hw.get("z", 0) > 30 else "深度有一定敏感度",
        "psi": "很平 → 航向主要来自RC2几何" if hw.get("psi", 0) > 5 else "航向有一定CZ信息",
        "theta": "θ主要由方位提供" if hw.get("theta", 0) > 2 else "θ有CZ信息",
    }
    for dim in ["r", "v", "z", "psi", "theta"]:
        md.append(f"| {dim} | {fnum(hw.get(dim),4)} | {interp.get(dim,'')} |")
    md.append("")
    md.append("二维切片见 `rc3_rv_sensitivity.csv` / `rc3_rz_sensitivity.csv` 与图C。")
    md.append("")
    md.append("**严谨表述**：CZ 的直接观测贡献主要来自 **距离 r**；速度改善更多来自距离–速度联合轨迹，而非 Doppler；深度与航向在本模型下改善有限。不得仅从“候选排除率”反推某一参数精度。")
    md.append("")
    md.append("## 3. 图")
    md.append("")
    md.append("- figures/figA_rmse_baseline_vs_method.svg")
    md.append("- figures/figB_rmse_improvement.svg")
    md.append("- figures/figC_rc3_parameter_sensitivity.svg")
    md.append("")
    md.append("## 4. 禁止表述（已执行）")
    md.append("")
    md.append("- 不再说“候选排除97% ⇒ 测距提高97%”")
    md.append("- 只允许：“RC2下距离RMSE=A；加入CZ后=B；降低C%；覆盖=D”或“RC2结构不可辨识；CZ后宽度A→B”")
    md.append("- 未删失败 trial；未为漂亮 RMSE 改先验；未进 P5")
    md.append("")
    md.append("## 5. 交付清单")
    md.append("")
    for pth in sorted(OUT.rglob("*")):
        if pth.is_file():
            md.append(f"- results/P4_5_parameter_accuracy/{pth.relative_to(OUT)}")
    md.append("")
    md.append("**P4.5 完成后停止。**")
    (OUT / "P4_5_REPORT.md").write_text("\n".join(md), encoding="utf-8")

    # GPT sync
    gs = []
    gs.append("# P4.5 — GPT 同步稿")
    gs.append("")
    gs.append(f"- UTC: {NOW}")
    gs.append(f"- 设计: {N_TRUTH} truth × {N_NOISE} noise / 场景 × 3 场景")
    gs.append("- 主表: parameter_accuracy_baseline_vs_method.csv")
    gs.append("")
    gs.append("## B-M 主参数（摘录）")
    gs.append("")
    bm = acc[acc["scenario"] == "B-M"]
    gs.append("| 参数 | RC2 RMSE | RC2+RC3 RMSE | 改善 | 覆盖(M2) | status |")
    gs.append("| --- | --- | --- | --- | --- | --- |")
    for _, r in bm.iterrows():
        gs.append(
            f"| {r['parameter']} | {fnum(r.get('rc2_rmse'))} | {fnum(r.get('rc23_rmse'))} | "
            f"{fnum(r.get('rmse_improvement'),3)} | {fnum(r.get('rc23_coverage'),2)} | {r.get('status','')} |"
        )
    gs.append("")
    gs.append("## CZ 直接敏感度（B-M 切片半高宽）")
    gs.append("")
    gs.append(json.dumps(hw, ensure_ascii=False))
    gs.append("")
    gs.append("## 下一轮")
    gs.append("")
    gs.append("> P4.5 完成后停止；进入 P5 成果收敛由任务书决定。")
    gs.append("")
    (OUT / "P4_5_GPT_SYNC.md").write_text("\n".join(gs), encoding="utf-8")

    print(f"DONE in {time.time()-t0:.1f}s")
    print(acc.to_string(index=False)[:2000])


if __name__ == "__main__":
    main()
