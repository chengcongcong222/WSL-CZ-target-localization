#!/usr/bin/env python3
"""R3-A1: WI/fringe method fidelity gate + HLA-MFP reference correction.

Outputs: results/R3_A1_fidelity/
No RC2 hard-pair MC, no R3-B/C, no P5.
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
OUT = ROOT / "results" / "R3_A1_fidelity"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

C_WATER = 1500.0

# Anchor paper-equivalent controlled scenario (PE-WI)
# Two-mode / ideal-waveguide WI with analytically known beta_true
PE = dict(
    name="PE-WI 理想波导两模态条纹（锚定文献等价受控场景）",
    literature_anchor=(
        "Chuprov waveguide-invariant line; practical range-frequency fringe ranging "
        "(e.g. waveguide-invariant passive ranging literature). "
        "Controlled reproduction uses standard two-mode interference: "
        "I(r,f)=|A1 e^{i k1 r}+A2 e^{i k2 r}|^2, "
        "beta_true = omega * d(k1-k2)/d(omega) / (k1-k2)."
    ),
    c0=C_WATER,
    # vertical wavenumbers — larger Δk → denser fringes for paper-condition fidelity
    gamma1=0.0005,
    gamma2=0.0500,
    A1=1.0,
    A2=0.90,
    r_true_km=25.0,
    f_lo=150.0,
    f_hi=375.0,
    n_f=51,
    r_lo_km=8.0,
    r_hi_km=45.0,
    n_r=151,
)

# E-STD
R_TRUE = 50e3
Z_TRUE = 200.0
Z_S = 200.0
Z_R = 200.0
R_ESTD = np.linspace(45e3, 60e3, 61)
Z_ESTD = np.linspace(150.0, 250.0, 41)
F_ESTD = np.linspace(150.0, 375.0, 32)

# HLA configs (theoretical reference, not equipment requirement)
HLA_CONFIGS = [
    dict(tag="single", n=1, d_m=0.0, note="单通道/当前R3-0参考"),
    dict(tag="hla_8x2m", n=8, d_m=2.0, note="基准水平阵 8×2m（孔径约14m）"),
    dict(tag="hla_diag_8x10m", n=8, d_m=10.0, note="大孔径诊断参考 8×10m（孔径约70m，非装备要求）"),
]

CONFIG = {
    "package": "R3_A1_fidelity",
    "created_utc": NOW,
    "anchor_paper": PE["literature_anchor"],
    "pe_scenario": {k: v for k, v in PE.items() if k != "literature_anchor"},
    "estd": {
        "r_km": [45, 60], "z_m": [150, 250], "r_true_km": 50.0, "z_true_m": 200.0,
        "band_hz": [150, 375], "propagation": "E-STD modal theory (P3-compatible), no Bellhop",
    },
    "hla_configs": HLA_CONFIGS,
    "metric_separation": [
        "fringe_score_width (alignment score peak width) — NOT localization accuracy",
        "range_error |r_hat-r_true|",
        "MC RMSE — only if A1-PASS (not this round unless pass)",
    ],
    "beta_judgment": ["CONSTANT", "LOCALLY_CONSTANT", "NON_CONSTANT"],
    "rc3a_decision": ["A1-PASS", "A1-CONDITIONAL", "A1-NOT-DIRECTLY-TRANSFERABLE"],
    "stop": "after R3-A1; no hard-pair MC unless A1-PASS; no R3-B/C; no P5",
    "r3_0_freeze_corrected": [
        "R3-0 MFP = frozen-observation matched-field reference ambiguity (not full HLA physical ceiling)",
        "R3-0 RC3-A = algorithm skeleton only, method fidelity NOT yet verified",
        "score width != localization RMSE",
    ],
}


# ---------------------------------------------------------------------------
# PE-WI controlled scenario: analytic beta + paper-like estimator
# ---------------------------------------------------------------------------
def pe_k(omega, gamma, c0=C_WATER):
    k2 = (omega / c0) ** 2 - gamma ** 2
    return np.sqrt(np.maximum(k2, 1e-12))


def pe_beta_true(omega, gamma1, gamma2, c0=C_WATER):
    """Fringe-slope invariant: beta = -omega * d(k1-k2)/domega / (k1-k2).

    Constant-phase fringes: (k1-k2) r = const ⇒ dr/df = -r (dΔk/df)/Δk
    ⇒ (f/r) dr/df = -f (dΔk/df)/Δk = -ω (dΔk/dω)/Δk  (literature beta for r(f)∝(f)^β local law)
    """
    d_omega = 1e-3
    k1 = pe_k(omega, gamma1, c0)
    k2 = pe_k(omega, gamma2, c0)
    k1p = pe_k(omega + d_omega, gamma1, c0)
    k2p = pe_k(omega + d_omega, gamma2, c0)
    dk = k1 - k2
    ddk = ((k1p - k2p) - dk) / d_omega
    beta_fringe = float(-omega * ddk / dk)
    return beta_fringe, float(dk), float(k1), float(k2)


def pe_intensity(freqs, r_grid, gamma1, gamma2, A1=1.0, A2=0.85, c0=C_WATER):
    """I(f,r) for two-mode interference (paper-condition field)."""
    I = np.zeros((len(freqs), len(r_grid)))
    for i, f in enumerate(freqs):
        w = 2 * np.pi * f
        k1 = pe_k(w, gamma1, c0)
        k2 = pe_k(w, gamma2, c0)
        p = A1 * np.exp(1j * k1 * r_grid) + A2 * np.exp(1j * k2 * r_grid)
        I[i, :] = np.abs(p) ** 2
    return I


def detrend_I(I, freqs):
    L = np.log(I + 1e-20)
    L0 = L - L.mean(axis=1, keepdims=True)
    L1 = L0.copy()
    for j in range(L0.shape[1]):
        c = np.polyfit(freqs, L0[:, j], 2)
        L1[:, j] = L0[:, j] - np.polyval(c, freqs)
    return L1


def fringe_beta_from_local_slope(L1, freqs, r_grid, r_center, r_halfwidth_km=3.0):
    """Estimate beta via local gradient slope in a range window: beta=(f/r)dr/df."""
    dL_df = np.gradient(L1, freqs, axis=0)
    dL_dr = np.gradient(L1, r_grid, axis=1)
    mask_r = np.abs(r_grid - r_center) <= r_halfwidth_km * 1e3
    if mask_r.sum() < 3:
        mask_r = np.ones_like(r_grid, dtype=bool)
    betas = []
    for i in range(len(freqs)):
        for j in np.where(mask_r)[0]:
            if not np.isfinite(dL_dr[i, j]) or abs(dL_dr[i, j]) < 1e-8:
                continue
            dr_df = -dL_df[i, j] / dL_dr[i, j]
            b = (freqs[i] / r_grid[j]) * dr_df
            if np.isfinite(b) and abs(b) < 5:
                betas.append(b)
    if not betas:
        return np.nan, np.nan
    betas = np.array(betas)
    return float(np.median(betas)), float(np.std(betas))


def fringe_beta_search(L1, freqs, r_grid, beta_grid):
    """Power-law alignment score (secondary check)."""
    f0 = float(np.mean(freqs))
    scores = []
    for beta in beta_grid:
        sc = 0.0
        n = 0
        for jref in range(0, len(r_grid), 3):
            r_ref = r_grid[jref]
            acc = np.zeros(len(freqs))
            for i, f in enumerate(freqs):
                r_f = r_ref * (f / f0) ** beta
                acc[i] = np.interp(r_f, r_grid, L1[i, :])
            sc += float(np.var(acc))
            n += 1
        scores.append((float(beta), sc / max(n, 1)))
    scores = np.array(scores)
    k = int(np.argmax(scores[:, 1]))
    return float(scores[k, 0]), float(scores[k, 1]), scores


def range_from_fringe_peaks(freqs, I_obs_1d, r_search, gamma1, gamma2, c0=C_WATER):
    """Paper-like ranging: match intensity fringe peaks in frequency to Δk(f)·r model.

    Observed: I(f) at true range (1D spectrum from two-mode field).
    Predicted fringe maxima for candidate r: where cos-phase peaks, I ∝ |A1|^2+|A2|^2+2A1A2 cos(Δk r)
    Score: correlation of detrended observed vs predicted cos(Δk(f) r).
    """
    L_obs = np.log(I_obs_1d + 1e-20)
    L_obs = L_obs - np.mean(L_obs)
    # optional light detrend
    c = np.polyfit(freqs, L_obs, 2)
    L_obs = L_obs - np.polyval(c, freqs)
    costs = []
    for r in r_search:
        w = 2 * np.pi * freqs
        k1 = pe_k(w, gamma1, c0)
        k2 = pe_k(w, gamma2, c0)
        phase = (k1 - k2) * r
        pred = np.cos(phase)
        pred = pred - np.mean(pred)
        cc = np.polyfit(freqs, pred, 2)
        pred = pred - np.polyval(cc, freqs)
        denom = (np.linalg.norm(pred) * np.linalg.norm(L_obs)) + 1e-30
        corr = float(np.dot(pred, L_obs) / denom)
        costs.append(1.0 - corr)
    costs = np.asarray(costs)
    k = int(np.argmin(costs))
    cmin, cmax = costs.min(), costs.max()
    # score width on cost (NOT RMSE)
    thr = cmin + 0.05 * (cmax - cmin + 1e-15)
    idx = np.where(costs <= thr)[0]
    score_width = float(r_search[idx.max()] - r_search[idx.min()]) / 1e3 if len(idx) >= 2 else 0.0
    return float(r_search[k]), costs, score_width


def range_from_beta(L1, freqs, r_grid, beta, r_search):
    """Secondary power-law alignment range search (score width ≠ RMSE)."""
    f0 = float(np.mean(freqs))
    costs = []
    for r_ref in r_search:
        acc = np.zeros(len(freqs))
        for i, f in enumerate(freqs):
            r_f = r_ref * (f / f0) ** float(beta)
            acc[i] = np.interp(r_f, r_grid, L1[i, :])
        costs.append(-float(np.var(acc)))
    costs = np.asarray(costs)
    k = int(np.argmin(costs))
    c = -costs
    cmin, cmax = c.min(), c.max()
    thr = cmin + 0.05 * (cmax - cmin + 1e-15)
    idx = np.where(c <= thr)[0]
    score_width = float(r_search[idx.max()] - r_search[idx.min()]) / 1e3 if len(idx) >= 2 else 0.0
    return float(r_search[k]), costs, score_width


def range_from_fringe_peaks(freqs, I_obs_1d, r_search, gamma1, gamma2, c0=C_WATER):
    """Paper-like ranging: match intensity fringe pattern in frequency to cos(Δk(f)·r)."""
    L_obs = np.log(I_obs_1d + 1e-20)
    c = np.polyfit(freqs, L_obs, 2)
    L_obs = L_obs - np.polyval(c, freqs)
    costs = []
    for r in r_search:
        w = 2 * np.pi * freqs
        k1 = pe_k(w, gamma1, c0)
        k2 = pe_k(w, gamma2, c0)
        phase = (k1 - k2) * r
        pred = np.cos(phase)
        cc = np.polyfit(freqs, pred, 2)
        pred = pred - np.polyval(cc, freqs)
        denom = (np.linalg.norm(pred) * np.linalg.norm(L_obs)) + 1e-30
        corr = float(np.dot(pred, L_obs) / denom)
        costs.append(1.0 - corr)
    costs = np.asarray(costs)
    k = int(np.argmin(costs))
    cmin, cmax = float(costs.min()), float(costs.max())
    thr = cmin + 0.05 * (cmax - cmin + 1e-15)
    idx = np.where(costs <= thr)[0]
    score_width = float(r_search[idx.max()] - r_search[idx.min()]) / 1e3 if len(idx) >= 2 else 0.0
    return float(r_search[k]), costs, score_width


def fringe_count_range(freqs, I_obs_1d, gamma1, gamma2, c0=C_WATER):
    """Analytic-style: number of spectral fringes × Δk band → r_hat."""
    L = np.log(I_obs_1d + 1e-20)
    c = np.polyfit(freqs, L, 2)
    L = L - np.polyval(c, freqs)
    # count zero crossings of detrended log-I as fringe count proxy
    s = np.sign(L)
    n_cross = int(np.sum(np.abs(np.diff(s)) > 0))
    # expected crossings ≈ 2 * |Δk(f_hi)-Δk(f_lo)| * r / (2π)  (up+down per fringe)
    w_hi = 2 * np.pi * freqs[-1]
    w_lo = 2 * np.pi * freqs[0]
    dk_hi = pe_k(w_hi, gamma1, c0) - pe_k(w_hi, gamma2, c0)
    dk_lo = pe_k(w_lo, gamma1, c0) - pe_k(w_lo, gamma2, c0)
    span = abs(dk_hi - dk_lo)
    if span < 1e-12:
        return np.nan, n_cross
    # n_cross ≈ 2 * span * r / (2π) = span * r / π  → r ≈ n_cross * π / span
    r_hat = n_cross * math.pi / span
    return float(r_hat), n_cross


def run_paper_reproduction():
    """A. PE-WI controlled: beta_true and r_true known."""
    freqs = np.linspace(PE["f_lo"], PE["f_hi"], PE["n_f"])
    r_grid = np.linspace(PE["r_lo_km"], PE["r_hi_km"], PE["n_r"]) * 1e3
    r_true = PE["r_true_km"] * 1e3
    I = pe_intensity(freqs, r_grid, PE["gamma1"], PE["gamma2"], PE["A1"], PE["A2"])
    L1 = detrend_I(I, freqs)

    betas = [pe_beta_true(2 * np.pi * f, PE["gamma1"], PE["gamma2"])[0] for f in freqs]
    beta_true_mean = float(np.mean(betas))
    beta_true_std = float(np.std(betas))

    beta_local, beta_local_std = fringe_beta_from_local_slope(L1, freqs, r_grid, r_true)
    beta_grid = np.linspace(0.2, 2.5, 47)
    beta_align, beta_align_score, score_tab = fringe_beta_search(L1, freqs, r_grid, beta_grid)

    j_true = int(np.argmin(np.abs(r_grid - r_true)))
    I_obs = I[:, j_true]
    r_hat, costs, score_w = range_from_fringe_peaks(
        freqs, I_obs, r_grid, PE["gamma1"], PE["gamma2"]
    )
    r_hat_cnt, n_cross = fringe_count_range(freqs, I_obs, PE["gamma1"], PE["gamma2"])
    r_hat_pl, costs_pl, score_w_pl = range_from_beta(
        L1, freqs, r_grid,
        beta_align if np.isfinite(beta_align) else beta_true_mean,
        r_grid,
    )

    err_peak = r_hat / 1e3 - PE["r_true_km"]
    err_cnt = (r_hat_cnt / 1e3 - PE["r_true_km"]) if np.isfinite(r_hat_cnt) else np.nan
    # pass if ANY primary paper-like estimator is within 2 km and beta recovered
    pass_beta = bool(
        abs(beta_true_mean) > 0.3 and (
            (np.isfinite(beta_local) and abs(beta_local - beta_true_mean) < max(0.15, 0.25 * abs(beta_true_mean)))
            or (np.isfinite(beta_align) and abs(beta_align - beta_true_mean) < max(0.15, 0.25 * abs(beta_true_mean)))
        )
    )
    pass_range = bool(
        abs(err_peak) < 2.0 or (np.isfinite(err_cnt) and abs(err_cnt) < 2.0)
        or abs(r_hat_pl / 1e3 - PE["r_true_km"]) < 2.0
    )

    rows = [{
        "scenario": "PE_WI_PAPER_CONDITION",
        "beta_true_mean": beta_true_mean,
        "beta_true_std": beta_true_std,
        "beta_local_slope": beta_local,
        "beta_local_slope_std": beta_local_std,
        "beta_align_search": beta_align,
        "r_true_km": PE["r_true_km"],
        "r_hat_peakmatch_km": r_hat / 1e3,
        "err_r_peakmatch_km": err_peak,
        "r_hat_fringe_count_km": r_hat_cnt / 1e3 if np.isfinite(r_hat_cnt) else np.nan,
        "err_r_fringe_count_km": err_cnt,
        "n_fringe_crossings": n_cross,
        "r_hat_powerlaw_km": r_hat_pl / 1e3,
        "err_r_powerlaw_km": r_hat_pl / 1e3 - PE["r_true_km"],
        "range_error_km": abs(err_peak),
        "fringe_score_width_km": score_w,
        "fringe_score_width_powerlaw_km": score_w_pl,
        "pass_beta": pass_beta,
        "pass_range": pass_range,
        "note": "beta_true uses fringe-slope convention; score_width is NOT localization RMSE",
    }]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "PAPER_REPRO.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(score_tab, columns=["beta", "score"]).to_csv(
        OUT / "PAPER_REPRO_beta_scores.csv", index=False, encoding="utf-8-sig"
    )
    print("PE-WI beta_true=", beta_true_mean,
          "beta_local=", beta_local,
          "beta_align=", beta_align,
          "r_peak=", r_hat / 1e3,
          "r_count=", r_hat_cnt / 1e3 if np.isfinite(r_hat_cnt) else None,
          "pass_beta", pass_beta, "pass_range", pass_range, flush=True)
    return df


# ---------------------------------------------------------------------------
# E-STD local beta_eff(r,f)
# ---------------------------------------------------------------------------
def estd_intensity(mode, freqs, r_grid, z_s=Z_S, z_r=Z_R):
    I = np.zeros((len(freqs), len(r_grid)))
    for i, f in enumerate(freqs):
        for j, r in enumerate(r_grid):
            I[i, j] = mode.amp(float(f), float(r), z_s=z_s, z_r=z_r) ** 2
    return I


def beta_eff_map(I, freqs, r_grid):
    """Local stripe slope via 2D gradient of detrended log-I.

    beta_eff = (f/r) * dr/df along intensity ridge / local phase of log residual.
    Practical estimator: for each (f_i, r_j), use local FFT or gradient:
      dI/df, dI/dr on high-pass residual L; stripe normal direction (dL/df, dL/dr);
      slope dr/df = -(dL/df)/(dL/dr); beta = (f/r)*dr/df
    """
    L = detrend_I(I, freqs)
    # gradients
    dL_df = np.gradient(L, freqs, axis=0)
    dL_dr = np.gradient(L, r_grid, axis=1)
    # avoid div0
    dL_dr_safe = np.where(np.abs(dL_dr) < 1e-8, np.nan, dL_dr)
    dr_df = -dL_df / dL_dr_safe
    F, R = np.meshgrid(freqs, r_grid, indexing="ij")
    with np.errstate(invalid="ignore", divide="ignore"):
        beta = (F / R) * dr_df
    # clip insane values
    beta = np.where(np.isfinite(beta) & (np.abs(beta) < 20), beta, np.nan)
    return beta, L, dr_df


def summarize_beta(beta, freqs, r_grid):
    b = beta[np.isfinite(beta)]
    if b.size == 0:
        return dict(status="NON_CONSTANT", note="no finite beta")
    med = float(np.median(b))
    q25, q75 = np.percentile(b, [25, 75])
    iqr = float(q75 - q25)
    p5, p95 = np.percentile(b, [5, 95])
    # local std over r at mid frequency
    i_f = len(freqs) // 2
    row = beta[i_f, :]
    row = row[np.isfinite(row)]
    local_std = float(np.std(row)) if row.size else np.nan
    # judgment
    # CONSTANT: IQR small and median stable across f slices
    med_by_f = []
    for i in range(len(freqs)):
        bi = beta[i, :]
        bi = bi[np.isfinite(bi)]
        if bi.size:
            med_by_f.append(float(np.median(bi)))
    med_by_f = np.array(med_by_f) if med_by_f else np.array([med])
    spread_f = float(np.std(med_by_f))
    if iqr < 0.25 and spread_f < 0.25 and np.all(np.abs(med_by_f - med) < 0.4):
        status = "CONSTANT"
    elif iqr < 0.8 and spread_f < 0.6:
        status = "LOCALLY_CONSTANT"
    else:
        status = "NON_CONSTANT"
    return dict(
        status=status,
        median=med,
        iqr=iqr,
        p5=float(p5),
        p95=float(p95),
        std_all=float(np.std(b)),
        median_by_f_std=spread_f,
        local_std_mid_f=local_std,
        frac_abs_lt_1p5=float(np.mean(np.abs(b) < 1.5)),
        frac_abs_gt_3=float(np.mean(np.abs(b) > 3.0)),
    )


# ---------------------------------------------------------------------------
# HLA-MFP reference
# ---------------------------------------------------------------------------
def hla_response(mode, freqs, r_target, z_s, z_r, n_el, d_m, theta=0.0):
    """Spatial snapshot y(f) ∈ C^M for HLA along y, target at range r on +x (theta=0)."""
    M = max(n_el, 1)
    ys = (np.arange(M) - (M - 1) / 2.0) * d_m
    Y = np.zeros((len(freqs), M), dtype=complex)
    for i, f in enumerate(freqs):
        for m in range(M):
            xt = r_target * math.cos(theta)
            yt = r_target * math.sin(theta)
            rm = math.hypot(xt - 0.0, yt - ys[m])
            # use amplitude with phase from modal sum via pressure if available else amp
            if hasattr(mode, "pressure"):
                Y[i, m] = mode.pressure(float(f), float(rm), z_s=z_s, z_r=z_r)
            else:
                # reconstruct complex phase from modal sum
                ms = mode.modes(float(f))
                acc = 0j
                rr = max(float(rm), 1.0)
                for kr, phi in ms:
                    a_s = float(np.interp(z_s, mode.z, phi))
                    a_r = float(np.interp(z_r, mode.z, phi))
                    att = np.exp(-2e-5 * (f / 200.0) * rr / 1000.0)
                    acc += (a_s * a_r / math.sqrt(kr * rr)) * np.exp(1j * kr * rr) * att
                Y[i, m] = acc
    return Y, ys


def hla_mfp_surface(mode, freqs, R, Z, Y_obs, z_s, n_el, d_m, theta=0.0):
    """MFP cost.

    M=1: broadband spectral match (spatial correlation is identically 1 — not usable).
    M>1: normalized spatial correlation per frequency (source S(f) eliminated), multi-freq mean.
    """
    M = max(n_el, 1)
    J = np.zeros((len(Z), len(R)))
    for iz, z in enumerate(Z):
        for ir, r in enumerate(R):
            H, _ = hla_response(mode, freqs, r, z_s, z, n_el, d_m, theta)
            if M == 1:
                # spectral amplitude match (unknown S(f))
                a_o = np.abs(Y_obs[:, 0])
                a_h = np.abs(H[:, 0])
                no = a_o / (np.linalg.norm(a_o) + 1e-30)
                nh = a_h / (np.linalg.norm(a_h) + 1e-30)
                J[iz, ir] = 1.0 - float(np.dot(no, nh))
            else:
                Bf = []
                for i in range(len(freqs)):
                    y = Y_obs[i]
                    h = H[i]
                    num = abs(np.vdot(h, y)) ** 2
                    den = (np.linalg.norm(h) ** 2) * (np.linalg.norm(y) ** 2) + 1e-30
                    Bf.append(num / den)
                J[iz, ir] = 1.0 - float(np.mean(Bf))
    return J


def mainlobe(J, R, Z, r_true, z_true):
    iz, ir = np.unravel_index(np.argmin(J), J.shape)
    r_hat, z_hat = float(R[ir]), float(Z[iz])
    jmin = float(J.min())

    def w1d(vec, grid, thr):
        idx = np.where(vec <= thr)[0]
        return float(grid[idx.max()] - grid[idx.min()]) / 1e3 if len(idx) >= 2 else 0.0  # km for r

    thr = jmin + 0.05
    wr = w1d(J[iz, :], R, thr)  # km
    wz = float(Z[np.where(J[:, ir] <= thr)[0]].max() - Z[np.where(J[:, ir] <= thr)[0]].min()) if np.any(J[:, ir] <= thr) else 0.0
    mask = J <= thr
    corr = np.nan
    if mask.sum() >= 5:
        idx = np.where(mask)
        rr = R[idx[1]] / 1e3
        zz = Z[idx[0]]
        if np.std(rr) > 1e-9 and np.std(zz) > 1e-9:
            corr = float(np.corrcoef(rr, zz)[0, 1])
    return dict(
        r_hat_km=r_hat / 1e3,
        z_hat_m=z_hat,
        err_r_km=r_hat / 1e3 - r_true / 1e3,
        err_z_m=z_hat - z_true,
        r_mainlobe_km=wr,
        z_mainlobe_m=wz,
        corr_rz=corr,
        truth_is_best=bool(iz == int(np.argmin(np.abs(Z - z_true))) and ir == int(np.argmin(np.abs(R - r_true)))),
        J_min=jmin,
    )


def run_hla_mfp():
    mode = p4.MODE_ENV["E0"]
    freqs = F_ESTD
    rows = []
    for cfg in HLA_CONFIGS:
        print(f"  HLA-MFP {cfg['tag']} ...", flush=True)
        Y_obs, _ = hla_response(mode, freqs, R_TRUE, Z_S, Z_TRUE, cfg["n"], cfg["d_m"], 0.0)
        J = hla_mfp_surface(mode, freqs, R_ESTD, Z_ESTD, Y_obs, Z_S, cfg["n"], cfg["d_m"], 0.0)
        m = mainlobe(J, R_ESTD, Z_ESTD, R_TRUE, Z_TRUE)
        rows.append({
            "config": cfg["tag"],
            "n_elements": cfg["n"],
            "spacing_m": cfg["d_m"],
            "aperture_m": cfg["d_m"] * max(cfg["n"] - 1, 0),
            "note": cfg["note"],
            **m,
        })
        # save S0-like broadband ambiguity for each config
        rec = []
        for iz, z in enumerate(Z_ESTD):
            for ir, r in enumerate(R_ESTD):
                rec.append({"r_km": r / 1e3, "z_m": z, "J": J[iz, ir]})
        pd.DataFrame(rec).to_csv(OUT / f"HLA_MFP_surface_{cfg['tag']}.csv", index=False, encoding="utf-8-sig")
        write_heatmap(FIG / f"HLA_MFP_{cfg['tag']}.svg", J, R_ESTD, Z_ESTD,
                      f"HLA-MFP J(r,z) — {cfg['tag']}", cfg["note"])
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "HLA_MFP_COMPARE.csv", index=False, encoding="utf-8-sig")
    print(df.to_string(index=False), flush=True)
    return df


def write_heatmap(path, J, R, Z, title, note):
    w, h = 760, 400
    ml, mr, mt, mb = 70, 90, 48, 56
    pw, ph = w - ml - mr, h - mt - mb
    S = 1.0 - J
    smin, smax = float(np.nanmin(S)), float(np.nanmax(S))
    nx, ny = len(R), len(Z)
    cw, ch = pw / nx, ph / ny
    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
         f'<rect width="{w}" height="{h}" fill="#f7f4ef"/>',
         f'<text x="{w/2}" y="26" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="14" font-weight="600">{title}</text>']
    for iz in range(ny):
        for ir in range(nx):
            u = (S[iz, ir] - smin) / (smax - smin + 1e-15)
            r = int(240 * (1 - u) + 31 * u)
            g = int(230 * (1 - u) + 111 * u)
            b = int(220 * (1 - u) + 139 * u)
            x = ml + ir * cw
            y = mt + (ny - 1 - iz) * ch
            p.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(cw,0.4):.1f}" height="{max(ch,0.4):.1f}" fill="rgb({r},{g},{b})"/>')
    ir = int(np.argmin(np.abs(R - 50e3)))
    iz = int(np.argmin(np.abs(Z - 200)))
    p.append(f'<circle cx="{ml+ir*cw+cw/2:.1f}" cy="{mt+(ny-1-iz)*ch+ch/2:.1f}" r="7" fill="none" stroke="#fff" stroke-width="2"/>')
    p.append(f'<text x="{ml+pw/2}" y="{h-16}" text-anchor="middle" font-size="11" font-family="sans-serif">r (km) 45→60</text>')
    p.append(f'<text x="{mt+ph/2}" y="16" font-size="10" fill="#777" font-family="sans-serif">{note}</text></svg>'.replace(
        f'y="16"', f'y="{h-36}"'))
    path.write_text("\n".join(p), encoding="utf-8")


def write_beta_figure(path, beta, freqs, r_grid, summary):
    w, h = 760, 400
    ml, mr, mt, mb = 70, 90, 48, 56
    pw, ph = w - ml - mr, h - mt - mb
    B = np.clip(beta, -5, 5)
    bmin, bmax = -3.0, 3.0
    nx, ny = len(r_grid), len(freqs)
    cw, ch = pw / nx, ph / ny
    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
         f'<rect width="{w}" height="{h}" fill="#f7f4ef"/>',
         f'<text x="{w/2}" y="26" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="14" font-weight="600">β_eff(r,f) — E-STD 模态场局部条纹斜率</text>']
    for i in range(ny):
        for j in range(nx):
            v = B[i, j]
            if not np.isfinite(beta[i, j]):
                fill = "#ddd"
            else:
                u = (np.clip(v, bmin, bmax) - bmin) / (bmax - bmin)
                r = int(240 * (1 - u) + 31 * u)
                g = int(230 * (1 - u) + 111 * u)
                b = int(220 * (1 - u) + 139 * u)
                fill = f"rgb({r},{g},{b})"
            x = ml + j * cw
            y = mt + (ny - 1 - i) * ch
            p.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(cw,0.4):.1f}" height="{max(ch,0.4):.1f}" fill="{fill}"/>')
    p.append(f'<text x="{ml+pw/2}" y="{h-28}" text-anchor="middle" font-size="12" font-family="sans-serif">r (km) 45→60</text>')
    p.append(f'<text x="16" y="{mt+ph/2}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90 16 {mt+ph/2})">f (Hz)</text>')
    p.append(f'<text x="{ml}" y="{h-8}" font-size="11" font-family="sans-serif">判定={summary.get("status")}  median={summary.get("median", float("nan")):.3g}  IQR={summary.get("iqr", float("nan")):.3g}</text></svg>')
    path.write_text("\n".join(p), encoding="utf-8")


def decide_rc3a(paper_df, beta_summary, hla_df):
    pr = paper_df.iloc[0]
    paper_ok = bool(pr["pass_beta"] and pr["pass_range"])
    status = beta_summary.get("status", "NON_CONSTANT")
    if not paper_ok:
        decision = "A1-IMPLEMENTATION_FIX_REQUIRED"
        why = (
            f"Paper-condition reproduction did not pass "
            f"(beta_local={pr.get('beta_local_slope')}, beta_true={pr.get('beta_true_mean')}, "
            f"err_r_peakmatch={pr.get('err_r_peakmatch_km')} km). "
            f"E-STD beta_eff status={status}. Do not transfer until PE-WI passes."
        )
        next_step = "fix WI estimator on PE-WI until beta_hat≈beta_true and r_hat≈r_true; re-run R3-A1 gate"
    elif status == "CONSTANT":
        decision = "A1-PASS"
        why = f"Paper repro passed; E-STD beta_eff CONSTANT (median={beta_summary.get('median')})."
        next_step = "proceed RC3-A to RC2 hard pairs + MC RMSE"
    elif status == "LOCALLY_CONSTANT":
        decision = "A1-CONDITIONAL"
        why = (
            f"Paper repro passed; E-STD beta_eff LOCALLY_CONSTANT "
            f"(median={beta_summary.get('median')}, IQR={beta_summary.get('iqr')}). "
            f"Only sub-regions may support constant-beta ranging."
        )
        next_step = "restrict RC3-A RMSE to locally-flat (r,f) windows; report applicability region as science result"
    else:
        decision = "A1-NOT-DIRECTLY-TRANSFERABLE"
        why = (
            f"Paper repro passed but E-STD beta_eff NON_CONSTANT "
            f"(median={beta_summary.get('median')}, IQR={beta_summary.get('iqr')}, "
            f"median_by_f_std={beta_summary.get('median_by_f_std')}). "
            f"Classical constant-beta fringe ranging does not directly transfer to deep first-CZ."
        )
        next_step = "STOP RC3-A parameter tuning; proceed to RC3-B (MMAC) next round"
    hla_note = ""
    if hla_df is not None and len(hla_df):
        hla_note = " | HLA-MFP: " + "; ".join(
            f"{r['config']}: r宽={r['r_mainlobe_km']:.2f}km z宽={r['z_mainlobe_m']:.1f}m"
            for _, r in hla_df.iterrows()
        )
    return {
        "rc3a_decision": decision,
        "why": why + hla_note,
        "next_step": next_step,
        "paper_condition_pass": paper_ok,
        "beta_eff_status": status,
        "beta_eff_summary": beta_summary,
        "paper_repro": pr.to_dict(),
        "hla_mfp": hla_df.to_dict(orient="records") if hla_df is not None else [],
        "metric_separation": {
            "fringe_score_width": "alignment/cost score peak width — NOT localization RMSE",
            "range_error": "single-trial |r_hat-r_true|",
            "mc_rmse": "not computed this round unless A1-PASS",
        },
        "created_utc": NOW,
        "stop": "R3-A1 only; no R3-B/C auto; no P5",
    }


def main():
    t0 = time.time()
    print("=== R3-A1 fidelity gate ===")
    print(f"OUT={OUT}")
    (OUT / "R3_A1_CONFIG.json").write_text(json.dumps(CONFIG, ensure_ascii=False, indent=2), encoding="utf-8")

    print("[1] PE-WI paper-condition reproduction ...", flush=True)
    paper_df = run_paper_reproduction()

    print("[2] E-STD beta_eff(r,f) ...", flush=True)
    mode = p4.MODE_ENV["E0"]
    I = estd_intensity(mode, F_ESTD, R_ESTD)
    beta, L1, dr_df = beta_eff_map(I, F_ESTD, R_ESTD)
    summary = summarize_beta(beta, F_ESTD, R_ESTD)
    # write maps
    rows = []
    for i, f in enumerate(F_ESTD):
        for j, r in enumerate(R_ESTD):
            rows.append({"f_hz": f, "r_km": r / 1e3, "beta_eff": beta[i, j], "L_detrend": L1[i, j]})
    pd.DataFrame(rows).to_csv(OUT / "beta_eff_map.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame([summary]).to_csv(OUT / "beta_eff_summary.csv", index=False, encoding="utf-8-sig")
    print("beta summary", summary, flush=True)
    write_beta_figure(FIG / "beta_eff_ESTD.svg", beta, F_ESTD, R_ESTD, summary)

    print("[3] HLA-MFP reference ...", flush=True)
    hla_df = run_hla_mfp()

    print("[4] reports + decision ...", flush=True)
    decision = decide_rc3a(paper_df, summary, hla_df)
    (OUT / "R3_A1_DECISION.json").write_text(
        json.dumps(decision, ensure_ascii=False, indent=2, default=str), encoding="utf-8"
    )

    # ANCHOR_PAPER_METHOD.md
    def fnum(x, nd=4):
        try:
            v = float(x)
            return "n/a" if not np.isfinite(v) else f"{v:.{nd}f}"
        except Exception:
            return "n/a"

    am = []
    am.append("# 锚定文献方法说明（R3-A1）")
    am.append("")
    am.append(f"UTC：{NOW}")
    am.append("")
    am.append("## 锚定选择")
    am.append("")
    am.append("项目库内无单篇 PDF 原文时，采用**可复算的文献等价受控场景 PE-WI**，锚定在：")
    am.append("")
    am.append(PE["literature_anchor"])
    am.append("")
    am.append("理由：任务书要求“明确 β 定义 + 条纹–距离关系 + 可复算场景”。两模态干涉是波导不变量文献中的标准受控模型，β_true 可解析写出，用于**方法保真**，不冒充深海实测。")
    am.append("")
    am.append("## 原公式（受控场景）")
    am.append("")
    am.append("- 观测量：\\(I(r,f)=|A_1 e^{ik_1 r}+A_2 e^{ik_2 r}|^2\\)")
    am.append("- 水平波数：\\(k_i=\\sqrt{(\\omega/c_0)^2-\\gamma_i^2}\\)")
    am.append("- 波导不变量：\\(\\beta=\\omega\\,\\dfrac{d(k_1-k_2)/d\\omega}{k_1-k_2}\\)")
    am.append("- 杼纹/距离：沿常相位 \\((k_1-k_2)r\\approx\\mathrm{const}\\) 的 \\(dr/df\\approx \\beta r/f\\)")
    am.append("- 场景参数：见 `R3_A1_CONFIG.json` → `pe_scenario`")
    am.append("")
    am.append("## 复现结果（PAPER_REPRO.csv）")
    am.append("")
    pr = paper_df.iloc[0]
    am.append(f"- β_true (band mean) = **{fnum(pr['beta_true_mean'])}** ± {fnum(pr['beta_true_std'])}")
    am.append(f"- β local-slope = **{fnum(pr.get('beta_local_slope'))}**；β align-search = **{fnum(pr.get('beta_align_search'))}**")
    am.append(f"- r_true = **{fnum(pr['r_true_km'],2)} km**")
    am.append(f"- r_hat peak-match = **{fnum(pr.get('r_hat_peakmatch_km'),3)} km**, 误差 = **{fnum(pr.get('err_r_peakmatch_km'),3)} km**")
    am.append(f"- r_hat fringe-count = **{fnum(pr.get('r_hat_fringe_count_km'),3)} km**, 误差 = **{fnum(pr.get('err_r_fringe_count_km'),3)} km**")
    am.append(f"- r_hat power-law = **{fnum(pr.get('r_hat_powerlaw_km'),3)} km**, 误差 = **{fnum(pr.get('err_r_powerlaw_km'),3)} km**")
    am.append(f"- 评分峰宽（≠定位精度）= **{fnum(pr.get('fringe_score_width_km'),3)} km**")
    am.append(f"- pass_beta={pr.get('pass_beta')}, pass_range={pr.get('pass_range')}")
    am.append("")
    am.append("## 指标三分离（强制）")
    am.append("")
    am.append("| 量 | 含义 | 能否当地位精度 |")
    am.append("| --- | --- | --- |")
    am.append("| fringe_score_width | 对齐评分峰宽 | **否** |")
    am.append("| range_error | 单次 \\(|\\hat r-r_{true}|\\) | 仅单次 |")
    am.append("| MC RMSE | 多真值×多噪声 | 是（A1-PASS 后才做） |")
    am.append("")
    (OUT / "ANCHOR_PAPER_METHOD.md").write_text("\n".join(am), encoding="utf-8")

    # Main report
    rp = []
    rp.append("# R3-A1 报告：条纹方法保真 + HLA-MFP 参考修正")
    rp.append("")
    rp.append(f"UTC：{NOW}  ·  输出：`results/R3_A1_fidelity/`")
    rp.append("")
    rp.append("## 0. R3-0 口径修正（冻结）")
    rp.append("")
    for line in CONFIG["r3_0_freeze_corrected"]:
        rp.append(f"- {line}")
    rp.append("")
    rp.append(f"## 1. RC3-A1 判定：`{decision['rc3a_decision']}`")
    rp.append("")
    rp.append(decision["why"])
    rp.append("")
    rp.append(f"**下一步**：{decision['next_step']}")
    rp.append("")
    rp.append("## 2. 论文条件复现（PE-WI）")
    rp.append("")
    rp.append("| β_true | β_local | β_align | r误差 peak/count/pl | 评分峰宽 | pass |")
    rp.append("| --- | --- | --- | --- | --- | --- |")
    rp.append(
        f"| {fnum(pr['beta_true_mean'])} | {fnum(pr.get('beta_local_slope'))} | {fnum(pr.get('beta_align_search'))} | "
        f"{fnum(pr.get('err_r_peakmatch_km'),3)} / {fnum(pr.get('err_r_fringe_count_km'),3)} / {fnum(pr.get('err_r_powerlaw_km'),3)} | "
        f"{fnum(pr.get('fringe_score_width_km'),3)} | β:{pr.get('pass_beta')} r:{pr.get('pass_range')} |"
    )
    rp.append("")
    rp.append("## 3. E-STD β_eff(r,f)")
    rp.append("")
    rp.append(f"- 判定：**{summary.get('status')}**")
    rp.append(f"- median={fnum(summary.get('median'))}, IQR={fnum(summary.get('iqr'))}, "
              f"P5–P95=[{fnum(summary.get('p5'))}, {fnum(summary.get('p95'))}]")
    rp.append(f"|β|<1.5 占比={fnum(summary.get('frac_abs_lt_1p5'),2)}，|β|>3 占比={fnum(summary.get('frac_abs_gt_3'),2)}")
    rp.append(f"|f| 方向 median 漂移 std={fnum(summary.get('median_by_f_std'))}")
    rp.append("")
    rp.append("**禁止**：仅扩大 β 搜索范围而不检验 β_eff(r,f) 是否常数。")
    rp.append("")
    if summary.get("status") == "NON_CONSTANT":
        rp.append("> **结论**：深海第一会聚区 E-STD 场中 β_eff 随 (r,f) 明显变化 → **经典常数波导不变量模型不能直接迁移**（迁移条件不成立，不是“方法无效”）。")
    elif summary.get("status") == "LOCALLY_CONSTANT":
        rp.append("> **结论**：β_eff 仅局部近似恒定 → RC3-A 仅能在局部 (r,f) 窗口做条件化验证。")
    else:
        rp.append("> **结论**：β_eff 近似常数，经典条纹模型在 E-STD 具备迁移条件（仍须通过论文复现门）。")
    rp.append("")
    rp.append("## 4. HLA-MFP 参考（空间孔径）")
    rp.append("")
    rp.append("对比：单通道 / 8×2m / 大孔径诊断（8×10m，非装备要求）。S0 宽带，无噪，环境匹配，未知源按频消去（归一化空间相关）。")
    rp.append("")
    if hla_df is not None and len(hla_df):
        rp.append("| 配置 | 孔径 m | r̂ km | ẑ m | r主瓣宽 km | z主瓣宽 m | corr(r,z) | 真值最优 |")
        rp.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
        for _, r in hla_df.iterrows():
            rp.append(
                f"| {r['config']} | {fnum(r['aperture_m'],1)} | {fnum(r['r_hat_km'],2)} | {fnum(r['z_hat_m'],1)} | "
                f"{fnum(r['r_mainlobe_km'],2)} | {fnum(r['z_mainlobe_m'],1)} | {fnum(r['corr_rz'])} | {r['truth_is_best']} |"
            )
        rp.append("")
        rp.append("解读：")
        rp.append("- R3-0 的“MFP”实为**冻结观测（单通道/多频）匹配场参考**，不是完整水平阵物理上限。")
        rp.append("- 若大孔径 HLA 的 r/z 主瓣仍宽 → 更接近“当前场景传播信息有限”；若明显变窄 → 空间孔径有潜在增益（理论诊断，非型号指标）。")
        rp.append("- **主瓣宽 ≠ 最终 RMSE 下限**（离网点估计误差可以远小于主瓣宽）。")
    rp.append("")
    rp.append("## 5. 图与数据")
    rp.append("")
    rp.append("- figures/beta_eff_ESTD.svg")
    rp.append("- figures/HLA_MFP_single.svg / HLA_MFP_hla_8x2m.svg / HLA_MFP_hla_diag_8x10m.svg")
    rp.append("- PAPER_REPRO.csv, beta_eff_map.csv, beta_eff_summary.csv, HLA_MFP_COMPARE.csv")
    rp.append("")
    rp.append("## 6. 停止条件")
    rp.append("")
    rp.append("- 未跑 RC2 困难候选 / Monte Carlo（除非 A1-PASS）")
    rp.append("- 未进入 RC3-B/C")
    rp.append("- 未进入 P5")
    rp.append(f"- **本轮判定：{decision['rc3a_decision']}**")
    rp.append("")
    (OUT / "R3_A1_REPORT.md").write_text("\n".join(rp), encoding="utf-8")

    gs = []
    gs.append("# R3-A1 — GPT 同步稿")
    gs.append("")
    gs.append(f"- UTC: {NOW}")
    gs.append(f"- **判定：{decision['rc3a_decision']}**")
    gs.append(f"- {decision['why']}")
    gs.append(f"- 下一步：{decision['next_step']}")
    gs.append("")
    gs.append("## PE-WI 复现")
    gs.append(paper_df.to_string(index=False))
    gs.append("")
    gs.append("## β_eff (E-STD)")
    gs.append(json.dumps(summary, ensure_ascii=False, default=str))
    gs.append("")
    gs.append("## HLA-MFP")
    if hla_df is not None:
        gs.append(hla_df.to_string(index=False))
    gs.append("")
    gs.append("## R3-0 口径")
    gs.append("- MFP=冻结观测参考，非完整HLA物理上限")
    gs.append("- RC3-A骨架≠论文复现完成")
    gs.append("- score width ≠ RMSE")
    gs.append("")
    gs.append("停止：无 hard-pair MC（除非PASS），无 R3-B/C，无 P5。")
    gs.append("")
    (OUT / "R3_A1_GPT_SYNC.md").write_text("\n".join(gs), encoding="utf-8")

    print(f"DONE in {time.time()-t0:.1f}s")
    print("DECISION", decision["rc3a_decision"])
    print(decision["why"])


if __name__ == "__main__":
    main()
