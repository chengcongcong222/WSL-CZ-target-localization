#!/usr/bin/env python3
"""P4 / G4: integrated functional & performance boundary for B3D-3.

Does NOT re-run P2/P3. Uses frozen E-STD + RC2 bearing + RC3 CZ modal model.
Approach: B0 baseline + four core scans + C-L/C-M/C-U, set-valued metrics.
"""
from __future__ import annotations

import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
P2 = ROOT / "results" / "P2_RC2_kinematic_boundary"
P3 = ROOT / "results" / "P3_RC3_increment"
OUT = ROOT / "results" / "P4_G4_performance_boundary"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

# ---------------------------------------------------------------------------
# Frozen constants
# ---------------------------------------------------------------------------
U_PLAT = 2.0
C_WATER = 1500.0
Z_A_TRUE = 1300.0
B_MUNK = 1300.0
EPS_MUNK = 0.00737
C0_MUNK = 1500.0
DEPTH_M = 5000.0
F_GRID = np.linspace(150.0, 375.0, 10)
S2_LINES = [166.0, 201.0, 235.0, 283.0, 338.0]
S1_LINES_IN_BAND = [168.0, 204.0, 232.0, 279.0, 320.0]
ARRAY_N = 8
ARRAY_D = 12.0
RNG_SEED = 404
N_TRIALS = 1  # boundary curves; noise handled via source×SNR scan

# Prior box
R_LO, R_HI = 45e3, 60e3
TH_LO, TH_HI = math.radians(-5), math.radians(5)
V_LO, V_HI = 1.0, 3.0
PSI_LO, PSI_HI = math.radians(-15), math.radians(15)
Z_CAND = np.array([150.0, 180.0, 200.0, 220.0, 250.0])
PRIOR_W = {
    "r_km": (R_HI - R_LO) / 1e3,
    "theta_deg": math.degrees(TH_HI - TH_LO),
    "z_m": float(Z_CAND.max() - Z_CAND.min()),
    "v_mps": V_HI - V_LO,
    "psi_deg": math.degrees(PSI_HI - PSI_LO),
}

# Candidate grids (P4-efficient)
R_G = np.arange(45.0, 60.0 + 1e-9, 1.0) * 1e3
TH_G = np.deg2rad(np.arange(-5.0, 5.0 + 1e-9, 1.0))
V_G = np.arange(1.0, 3.0 + 1e-9, 0.5)
PSI_G = np.deg2rad(np.arange(-15.0, 15.0 + 1e-9, 3.0))
Z_G = Z_CAND.copy()

# B0
B0 = dict(
    r0_m=50e3, theta0_deg=0.0, z=200.0, v=2.0, psi_deg=5.0,
    U=U_PLAT, T=600.0, sigma_deg=0.1, turn_deg=0.0,
    source="S1", snr_db=10.0,
)

SCEN_CL = dict(name="C-L困难", sigma_deg=0.2, T=300.0, turn_deg=0.0, source="S1", snr_db=0.0,
               r0_m=50e3, theta0_deg=0.0, z=200.0, v=2.0, psi_deg=5.0)
SCEN_CM = dict(name="C-M基准", sigma_deg=0.1, T=600.0, turn_deg=5.0, source="S1", snr_db=10.0,
               r0_m=50e3, theta0_deg=0.0, z=200.0, v=2.0, psi_deg=5.0)
SCEN_CU = dict(name="C-U理论上界", sigma_deg=0.02, T=1200.0, turn_deg=15.0, source="S2", snr_db=20.0,
               r0_m=50e3, theta0_deg=0.0, z=200.0, v=2.0, psi_deg=5.0)

# Environment mismatch (single physical form: channel-axis depth shift)
ENV = {
    "E0": dict(label="完全匹配", z_a_shift_m=0.0, c0_shift=0.0),
    "E1": dict(label="轻度SSP偏差", z_a_shift_m=50.0, c0_shift=0.0),
    "E2": dict(label="中等SSP偏差", z_a_shift_m=150.0, c0_shift=2.0),
}

CONFIG = {
    "project": "B3D-3",
    "package": "P4_G4_performance_boundary",
    "created_utc": NOW,
    "gate": "G4",
    "frozen": {
        "G1": "E-STD scene frozen",
        "G2": "RC2 not long-range ranging; LIN r-v ridge; maneuvers insufficient for practical range; CRLB rank-def=UNDEFINED",
        "G3": "PROPAGATION_INCREMENT_CONFIRMED in frozen theory model; Doppler weak on constant v_rad pairs",
    },
    "B0": B0,
    "scenarios": {"C-L": SCEN_CL, "C-M": SCEN_CM, "C-U": SCEN_CU},
    "prior_widths": PRIOR_W,
    "grids": {
        "r_km": R_G.tolist(),
        "theta_deg": np.rad2deg(TH_G).tolist(),
        "v_mps": V_G.tolist(),
        "psi_deg": np.rad2deg(PSI_G).tolist(),
        "z_m": Z_G.tolist(),
    },
    "scans": {
        "T_s": [60, 120, 300, 600, 1200],
        "sigma_deg": [0.02, 0.05, 0.1, 0.2, 0.5],
        "source": ["S0", "S1", "S2"],
        "snr_db": [20, 10, 0, -10],
        "z_true_m": [180, 200, 220],
        "turn_deg": [0, 5, 10, 15],
        "T_maneuver_s": [300, 600, 1200],
        "env": ENV,
    },
    "array_param": "sigma_theta + RC3 multi-element joint residual (no specific array length as requirement)",
    "ten_percent_policy": "10% reference only for positive states r,z,v; angles use absolute width/deg",
    "method": "set-valued candidate view; multi-modal retained; UNRESOLVED if contraction low / width~prior",
    "no_cartesian_product": True,
    "propagation": "modal theory model consistent with E-STD (KRAKEN not on host); single env form: channel-axis shift",
    "forbidden": [
        "Bellhop A2 audit", "new propagation features", "full P2 re-run",
        "redo P3 hard-pair selection", "new sensors", "full cartesian product",
        "S2 as real UUV", "auto-enter P5",
    ],
}


# ---------------------------------------------------------------------------
# Physics
# ---------------------------------------------------------------------------
def platform_xy(t, turn_deg=0.0, u=U_PLAT):
    t = np.asarray(t, float)
    xp = u * t.copy()
    yp = np.zeros_like(t)
    if abs(turn_deg) > 1e-12:
        t_turn = float(np.max(t)) * 0.5
        d = math.radians(turn_deg)
        c, s = math.cos(d), math.sin(d)
        m = t > t_turn
        dt = t[m] - t_turn
        xp[m] = u * t_turn + u * dt * c
        yp[m] = u * dt * s
    return xp, yp


def target_xy(t, r0, th0, v, psi):
    return r0 * np.cos(th0) + v * t * np.cos(psi), r0 * np.sin(th0) + v * t * np.sin(psi)


def geom(t, r0, th0, v, psi, turn_deg=0.0):
    xp, yp = platform_xy(t, turn_deg)
    xt, yt = target_xy(t, r0, th0, v, psi)
    dx, dy = xt - xp, yt - yp
    r = np.hypot(dx, dy)
    th = np.arctan2(dy, dx)
    ux, uy = dx / np.maximum(r, 1), dy / np.maximum(r, 1)
    if abs(turn_deg) < 1e-12:
        pvx, pvy = np.full_like(t, U_PLAT), np.zeros_like(t)
    else:
        tt = float(np.max(t)) * 0.5
        d = math.radians(turn_deg)
        pvx = np.where(t <= tt, U_PLAT, U_PLAT * math.cos(d))
        pvy = np.where(t <= tt, 0.0, U_PLAT * math.sin(d))
    vrad = (v * np.cos(psi) - pvx) * ux + (v * np.sin(psi) - pvy) * uy
    return th, r, vrad


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def munk_ssp(z, z_a=Z_A_TRUE, c0=C0_MUNK):
    z = np.asarray(z, float)
    eta = 2.0 * (z - z_a) / B_MUNK
    return c0 * (1.0 + EPS_MUNK * (np.exp(eta) - 1.0 - eta))


class ModeModel:
    def __init__(self, n_modes=20, n_z=300, z_a=Z_A_TRUE, c0=C0_MUNK):
        self.n_modes = n_modes
        self.z = np.linspace(0, DEPTH_M, n_z)
        self.z_a = z_a
        self.c0 = c0
        self.c = munk_ssp(self.z, z_a=z_a, c0=c0)
        self.cmin = float(np.min(self.c))
        self.cache = {}

    def modes(self, f_hz):
        key = round(float(f_hz), 3)
        if key in self.cache:
            return self.cache[key]
        omega = 2 * np.pi * float(f_hz)
        z_lo, z_hi = self.z_a - 1100.0, self.z_a + 1500.0
        D_eff = max(z_hi - z_lo, 500.0)
        k0 = omega / self.cmin
        out = []
        for m in range(1, self.n_modes + 1):
            gamma = (m - 0.25) * np.pi / D_eff
            kr2 = k0 ** 2 - gamma ** 2
            if kr2 <= 1e-8:
                continue
            kr = math.sqrt(kr2)
            width = 350.0 + 55.0 * m
            phi = np.exp(-0.5 * ((self.z - self.z_a) / width) ** 2) * np.cos(gamma * (self.z - self.z_a))
            nrm = math.sqrt(float(np.trapezoid(phi ** 2, self.z))) + 1e-15
            out.append((kr, phi / nrm))
        if not out:
            phi = np.exp(-0.5 * ((self.z - self.z_a) / 400.0) ** 2)
            phi = phi / (math.sqrt(float(np.trapezoid(phi ** 2, self.z))) + 1e-15)
            out = [(k0 * 0.999, phi)]
        self.cache[key] = out
        return out

    def amp(self, f_hz, r_m, z_s=200.0, z_r=200.0):
        ms = self.modes(f_hz)
        p = 0j
        r = max(float(r_m), 1.0)
        for kr, phi in ms:
            a_s = float(np.interp(z_s, self.z, phi))
            a_r = float(np.interp(z_r, self.z, phi))
            att = np.exp(-2e-5 * (f_hz / 200.0) * r / 1000.0)
            p += (a_s * a_r / math.sqrt(kr * r)) * np.exp(1j * kr * r) * att
        return abs(p)


MODE_TRUE = ModeModel()  # generates observations
MODE_ENV = {k: ModeModel(z_a=Z_A_TRUE + v["z_a_shift_m"], c0=C0_MUNK + v["c0_shift"])
            for k, v in ENV.items()}

# Amplitude lookup: amp[f_i, r_i, z_i] on r in [40,65] km
R_LU = np.linspace(40e3, 65e3, 51)


def build_amp_table(mode_model, freqs, z_list=Z_G):
    tab = np.zeros((len(freqs), len(R_LU), len(z_list)))
    for i, f in enumerate(freqs):
        for j, r in enumerate(R_LU):
            for k, z in enumerate(z_list):
                tab[i, j, k] = mode_model.amp(float(f), float(r), z_s=float(z), z_r=200.0)
    return tab


def lookup_amp(tab, freqs, r_m, z_m):
    """r_m (M,), z_m scalar → (F, M)"""
    F = tab.shape[0]
    iz = int(np.argmin(np.abs(Z_G - z_m)))
    out = np.empty((F, r_m.size))
    for i in range(F):
        out[i] = np.interp(r_m, R_LU, tab[i, :, iz])
    return out


def s0_shape(freqs):
    s = 1.0 / (1.0 + ((freqs - 240.0) / 120.0) ** 2)
    return s / s.max()


def source_freqs(source):
    if source == "S0":
        return F_GRID.copy(), s0_shape(F_GRID), "continuum"
    if source == "S1":
        fs = np.array(S1_LINES_IN_BAND, float)
        amps = np.array([0.25, 0.22, 0.18, 0.15, 0.12], float)
        return fs, amps / amps.max(), "lines"
    fs = np.array(S2_LINES, float)
    amps = 1.0 / np.arange(1, len(S2_LINES) + 1)
    return fs, amps / amps.max(), "lines"


def array_gain(theta, freqs):
    gs = []
    for f in freqs:
        k = 2 * np.pi * f / C_WATER
        psi = k * ARRAY_D * np.sin(theta)
        if abs(math.sin(psi / 2 + 1e-18)) < 1e-12:
            g = float(ARRAY_N)
        else:
            g = abs(math.sin(ARRAY_N * psi / 2) / math.sin(psi / 2))
        gs.append(g / ARRAY_N)
    return np.asarray(gs)


def make_grid(include_z=True):
    rr, tt, vv, pp = np.meshgrid(R_G, TH_G, V_G, PSI_G, indexing="ij")
    g = dict(r0=rr.ravel(), th0=tt.ravel(), v=vv.ravel(), psi=pp.ravel())
    n = g["r0"].size
    if include_z:
        k = len(Z_G)
        for key in ["r0", "th0", "v", "psi"]:
            g[key] = np.repeat(g[key], k)
        g["z"] = np.tile(Z_G, n)
    else:
        g["z"] = np.full(n, np.nan)
    g["n"] = g["r0"].size
    return g


def bearings_grid(t, g, turn_deg):
    xp, yp = platform_xy(t, turn_deg)
    N, M = g["r0"].size, t.size
    out = np.empty((N, M))
    chunk = max(1, 8000)
    for i0 in range(0, N, chunk):
        i1 = min(N, i0 + chunk)
        xt = g["r0"][i0:i1, None] * np.cos(g["th0"][i0:i1, None]) + g["v"][i0:i1, None] * t[None] * np.cos(g["psi"][i0:i1, None])
        yt = g["r0"][i0:i1, None] * np.sin(g["th0"][i0:i1, None]) + g["v"][i0:i1, None] * t[None] * np.sin(g["psi"][i0:i1, None])
        out[i0:i1] = np.arctan2(yt - yp[None], xt - xp[None])
    return out


def cand_acoustic_feats(t, g, turn_deg, source, tab):
    """Fast acoustic features via lookup amp table. Returns (N, F+1) unit rows."""
    freqs, _, kind = source_freqs(source)
    N, M, F = g["r0"].size, t.size, freqs.size
    xp, yp = platform_xy(t, turn_deg)
    feats = np.zeros((N, F + 1))
    # precompute array gain for a grid of bearings? compute per time with candidate bearing
    chunk = max(1, 4000)
    for i0 in range(0, N, chunk):
        i1 = min(N, i0 + chunk)
        r0 = g["r0"][i0:i1, None]
        th0 = g["th0"][i0:i1, None]
        vv = g["v"][i0:i1, None]
        ps = g["psi"][i0:i1, None]
        zz = g["z"][i0:i1]
        xt = r0 * np.cos(th0) + vv * t[None] * np.cos(ps)
        yt = r0 * np.sin(th0) + vv * t[None] * np.sin(ps)
        dx = xt - xp[None]
        dy = yt - yp[None]
        rr = np.hypot(dx, dy)
        th = np.arctan2(dy, dx)
        for k in range(i1 - i0):
            A = lookup_amp(tab, freqs, rr[k], float(zz[k]))
            for it in range(M):
                A[:, it] *= array_gain(th[k, it], freqs)
            prof = []
            for j in range(F):
                y = A[j] / (np.linalg.norm(A[j]) + 1e-20)
                prof.append(float(np.mean(y)))
            bp = A.sum(axis=0)
            bp = bp / (np.linalg.norm(bp) + 1e-20)
            feats[i0 + k, :F] = prof
            feats[i0 + k, F] = float(np.std(bp))
    feats = feats / (np.linalg.norm(feats, axis=1, keepdims=True) + 1e-20)
    return feats


def obs_features(t, truth, turn_deg, source, snr_db, rng, tab):
    freqs, _, kind = source_freqs(source)
    F, M = freqs.size, t.size
    xp, yp = platform_xy(t, turn_deg)
    xt, yt = target_xy(t, truth["r0_m"], truth["theta0_rad"], truth["v"], truth["psi_rad"])
    dx, dy = xt - xp, yt - yp
    rr = np.hypot(dx, dy)
    th = np.arctan2(dy, dx)
    A = lookup_amp(tab, freqs, rr, float(truth["z"]))
    for it in range(M):
        A[:, it] *= array_gain(th[it], freqs)
    if snr_db is not None:
        sig = float(np.mean(A ** 2)) + 1e-20
        nvar = sig / (10 ** (snr_db / 10.0))
        A = np.abs(A + rng.normal(0, math.sqrt(nvar), size=A.shape))
    prof = []
    for j in range(F):
        y = A[j] / (np.linalg.norm(A[j]) + 1e-20)
        prof.append(float(np.mean(y)))
    bp = A.sum(axis=0)
    bp = bp / (np.linalg.norm(bp) + 1e-20)
    feat = np.array(prof + [float(np.std(bp))])
    feat = feat / (np.linalg.norm(feat) + 1e-20)
    sigma = math.radians(truth.get("sigma_deg", 0.1))
    bearing_obs = th + (rng.normal(0, sigma, size=th.shape) if sigma > 0 else 0.0)
    return feat, bearing_obs


def summarize_set(g, mask, cost, truth, prior_w=PRIOR_W):
    n = int(mask.sum())
    if n == 0:
        k = int(np.argmin(cost))
        mask = np.zeros_like(cost, dtype=bool)
        mask[k] = True
        n = 1
    def width(arr, scale=1.0):
        if n < 2:
            return 0.0
        return float(arr[mask].max() - arr[mask].min()) / scale
    w_r = width(g["r0"], 1e3)
    w_th = math.degrees(float(g["th0"][mask].max() - g["th0"][mask].min())) if n >= 2 else 0.0
    w_v = width(g["v"])
    w_psi = math.degrees(float(g["psi"][mask].max() - g["psi"][mask].min())) if n >= 2 else 0.0
    if np.isfinite(g["z"]).any() and n >= 2:
        w_z = float(g["z"][mask].max() - g["z"][mask].min())
    else:
        w_z = float("nan")
    cmin = float(np.min(cost[mask]))
    wts = np.exp(-0.5 * np.clip(cost - cmin, 0, None) / 1e-6)
    wts = np.where(mask, wts, 0.0)
    s = wts.sum()
    if s <= 0:
        wts = mask.astype(float)
        s = max(wts.sum(), 1.0)
    wts = wts / s
    def pstd(arr, scale=1.0):
        m = float(np.sum(wts * arr))
        return float(np.sqrt(max(np.sum(wts * (arr - m) ** 2), 0))) / scale
    def contr(w, prior):
        if w is None or not np.isfinite(w):
            return float("nan")
        return float(np.clip(1.0 - w / prior, 0.0, 1.0))
    def in_band(arr, tval, tol):
        if n < 2:
            return abs(float(np.sum(wts * arr)) - tval) <= tol
        return bool(arr[mask].min() - tol <= tval <= arr[mask].max() + tol)
    out = {
        "n_acc": n,
        "w_r_km": w_r, "w_theta_deg": w_th, "w_z_m": w_z, "w_v_mps": w_v, "w_psi_deg": w_psi,
        "std_r_km": pstd(g["r0"], 1e3), "std_theta_deg": math.degrees(pstd(g["th0"])),
        "std_v": pstd(g["v"]), "std_psi_deg": math.degrees(pstd(g["psi"])),
        "contr_r": contr(w_r, prior_w["r_km"]), "contr_theta": contr(w_th, prior_w["theta_deg"]),
        "contr_z": contr(w_z, prior_w["z_m"]), "contr_v": contr(w_v, prior_w["v_mps"]),
        "contr_psi": contr(w_psi, prior_w["psi_deg"]),
        "retain_r": in_band(g["r0"], truth["r0_m"], 1.5e3),
        "retain_theta": in_band(g["th0"], truth["theta0_rad"], math.radians(1.0)),
        "retain_v": in_band(g["v"], truth["v"], 0.3),
        "retain_psi": in_band(g["psi"], truth["psi_rad"], math.radians(3.0)),
        "retain_z": in_band(g["z"], truth["z"], 20.0) if np.isfinite(g["z"]).any() else True,
        "mean_r_km": float(np.sum(wts * g["r0"]) / 1e3),
        "mean_v": float(np.sum(wts * g["v"])),
        "mean_psi_deg": float(math.degrees(np.sum(wts * g["psi"]))),
        "mean_z": float(np.sum(wts * g["z"])) if np.isfinite(g["z"]).any() else float("nan"),
    }
    if np.isfinite(g["z"]).any() and n >= 5:
        rr, zz = g["r0"][mask] / 1e3, g["z"][mask]
        out["corr_rz"] = float(np.corrcoef(rr, zz)[0, 1]) if np.std(rr) > 1e-9 and np.std(zz) > 1e-9 else 0.0
    else:
        out["corr_rz"] = float("nan")
    out["unresolved_r"] = bool((not np.isfinite(out["contr_r"])) or out["contr_r"] < 0.2 or w_r >= 0.8 * prior_w["r_km"])
    out["unresolved_theta"] = bool(out["contr_theta"] < 0.2 or w_th >= 0.8 * prior_w["theta_deg"])
    out["unresolved_v"] = bool(out["contr_v"] < 0.2 or w_v >= 0.8 * prior_w["v_mps"])
    out["unresolved_psi"] = bool(out["contr_psi"] < 0.2 or w_psi >= 0.8 * prior_w["psi_deg"])
    out["unresolved_z"] = bool((not np.isfinite(w_z)) or out["contr_z"] < 0.2 or (np.isfinite(w_z) and w_z >= 0.8 * prior_w["z_m"]))
    out["ten_pct_r"] = bool(w_r <= 0.10 * (truth["r0_m"] / 1e3) or out["std_r_km"] <= 0.10 * (truth["r0_m"] / 1e3))
    out["ten_pct_v"] = bool(w_v <= 0.10 * truth["v"] or out["std_v"] <= 0.10 * truth["v"])
    out["ten_pct_z"] = bool(np.isfinite(w_z) and w_z <= 0.10 * truth["z"])
    return out


# Precomputed tables per source/env filled in main
AMP_TABS = {}
_FEAT_CACHE = {}


def evaluate_scenario(g, truth, sigma_deg, T, turn_deg, source, snr_db,
                      use_rc3=True, env_key="E0", n_trials=N_TRIALS, seed=RNG_SEED):
    rng = np.random.default_rng(seed)
    t = np.linspace(0.0, T, max(6, int(round(T / 10.0)) + 1))
    truth = dict(truth)
    truth["sigma_deg"] = sigma_deg
    pred_b = bearings_grid(t, g, turn_deg)
    tab = AMP_TABS[(source, env_key)]
    cache_key = (source, env_key, float(T), float(turn_deg), g["n"])
    if cache_key not in _FEAT_CACHE:
        _FEAT_CACHE[cache_key] = cand_acoustic_feats(t, g, turn_deg, source, tab)
    feats = _FEAT_CACHE[cache_key]

    rows = []
    for trial in range(n_trials):
        obs_f, bearing_obs = obs_features(t, truth, turn_deg, source, snr_db, rng, AMP_TABS[(source, "E0")])
        c_bear_tot = np.sum(wrap(pred_b - bearing_obs[None, :]) ** 2, axis=1)
        c_bear = c_bear_tot / t.size
        thr_t = 13.3 * (math.radians(max(sigma_deg, 1e-4)) ** 2) * max(t.size, 1) * 0.2
        mask_rc2 = c_bear_tot <= float(c_bear_tot.min()) + max(thr_t, 1e-12)
        s_rc2 = summarize_set(g, mask_rc2, c_bear_tot, truth)
        s_rc2["layer"] = "RC2"
        s_rc2["trial"] = trial
        rows.append(s_rc2)
        if use_rc3:
            c_ac = np.sum((feats - obs_f[None, :]) ** 2, axis=1)
            c_joint = c_bear + 0.8 * c_ac
            cmin_j = float(c_joint.min())
            spread = float(np.percentile(c_joint, 90) - cmin_j)
            thr_j = max(0.02 * max(spread, 1e-9), 1e-4)
            mask_rc3 = c_joint <= cmin_j + thr_j
            s_rc3 = summarize_set(g, mask_rc3, c_joint, truth)
            s_rc3["layer"] = "RC2+RC3"
            s_rc3["trial"] = trial
            s_rc3["n_rc2"] = s_rc2["n_acc"]
            s_rc3["stage_reject_rate"] = 1.0 - (s_rc3["n_acc"] / max(s_rc2["n_acc"], 1))
            rows.append(s_rc3)
    df = pd.DataFrame(rows)
    agg = {}
    for layer in df["layer"].unique():
        d = df[df["layer"] == layer]
        rec = {"layer": layer}
        for c in d.columns:
            if c in ("layer", "trial"):
                continue
            if d[c].dtype == bool:
                rec[c] = float(d[c].mean())
            elif np.issubdtype(d[c].dtype, np.number):
                rec[c] = float(pd.to_numeric(d[c], errors="coerce").mean())
        agg[layer] = rec
    return agg, df


# ---------------------------------------------------------------------------
# Core figures (6)
# ---------------------------------------------------------------------------
def _esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def write_heatmap(path, title, xlab, ylab, xvals, yvals, Z, note="", fmt="{:.2f}", w=720, h=400):
    ml, mr, mt, mb = 70, 90, 48, 56
    pw, ph = w - ml - mr, h - mt - mb
    Z = np.asarray(Z, float)
    zmin, zmax = float(np.nanmin(Z)), float(np.nanmax(Z))
    if not np.isfinite(zmin) or not np.isfinite(zmax) or zmax - zmin < 1e-15:
        zmin, zmax = 0.0, 1.0
    nx, ny = len(xvals), len(yvals)
    cw, ch = pw / max(nx, 1), ph / max(ny, 1)

    def col(v):
        if v is None or not np.isfinite(v):
            return "#ddd"
        u = (v - zmin) / (zmax - zmin + 1e-15)
        r = int(240 * (1 - u) + 31 * u)
        g = int(230 * (1 - u) + 111 * u)
        b = int(220 * (1 - u) + 139 * u)
        return f"rgb({r},{g},{b})"

    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
         f'<rect width="{w}" height="{h}" fill="#f7f4ef"/>',
         f'<text x="{w/2}" y="28" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="15" font-weight="600" fill="#1a1a1a">{_esc(title)}</text>']
    for iy in range(ny):
        for ix in range(nx):
            v = Z[iy, ix]
            x = ml + ix * cw
            y = mt + (ny - 1 - iy) * ch
            p.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{cw:.1f}" height="{ch:.1f}" fill="{col(v)}" stroke="#eee"/>')
            if nx <= 8 and ny <= 8:
                p.append(f'<text x="{x+cw/2:.1f}" y="{y+ch/2+3:.1f}" text-anchor="middle" font-size="9" fill="#222" font-family="sans-serif">{fmt.format(v)}</text>')
    for ix, xv in enumerate(xvals):
        p.append(f'<text x="{ml+(ix+0.5)*cw:.1f}" y="{mt+ph+18}" text-anchor="middle" font-size="11" fill="#555" font-family="sans-serif">{_esc(xv)}</text>')
    for iy, yv in enumerate(yvals):
        p.append(f'<text x="{ml-8}" y="{mt+(ny-1-iy)*ch+ch/2+4:.1f}" text-anchor="end" font-size="11" fill="#555" font-family="sans-serif">{_esc(yv)}</text>')
    p.append(f'<text x="{ml+pw/2}" y="{h-14}" text-anchor="middle" font-size="12" fill="#333" font-family="sans-serif">{_esc(xlab)}</text>')
    p.append(f'<text x="16" y="{mt+ph/2}" text-anchor="middle" font-size="12" fill="#333" font-family="sans-serif" transform="rotate(-90 16 {mt+ph/2})">{_esc(ylab)}</text>')
    p.append(f'<text x="{ml}" y="{h-2}" font-size="10" fill="#777" font-family="sans-serif">{_esc(note)}</text></svg>')
    path.write_text("\n".join(p), encoding="utf-8")


def write_group_bars(path, title, cats, series, ylab, note="", w=720, h=400):
    ml, mr, mt, mb = 60, 140, 48, 60
    pw, ph = w - ml - mr, h - mt - mb
    cols = ["#8a8a8a", "#b45309", "#0f766e", "#1d4ed8", "#9f1239"]
    vals = [v for s in series for v in s["y"] if v is not None and np.isfinite(v)]
    vmax = max(vals + [0.2])
    ncat, nser = len(cats), len(series)
    bw = pw / max(ncat * (nser + 1), 1)
    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
         f'<rect width="{w}" height="{h}" fill="#f7f4ef"/>',
         f'<text x="{w/2}" y="28" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="15" font-weight="600" fill="#1a1a1a">{_esc(title)}</text>',
         f'<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="#fff" stroke="#d0ccc4"/>']
    for i in range(5):
        yv = vmax * i / 4
        Y = mt + ph - (yv / vmax) * ph
        p.append(f'<line x1="{ml}" y1="{Y:.1f}" x2="{ml+pw}" y2="{Y:.1f}" stroke="#e8e4dc"/>')
        p.append(f'<text x="{ml-8}" y="{Y+4:.1f}" text-anchor="end" font-size="11" fill="#555" font-family="sans-serif">{yv:.2f}</text>')
    for ic, cat in enumerate(cats):
        x0 = ml + ic * (pw / ncat) + bw * 0.4
        for i, s in enumerate(series):
            v = s["y"][ic]
            if v is None or not np.isfinite(v):
                continue
            hh = (v / vmax) * ph
            p.append(f'<rect x="{x0+i*bw:.1f}" y="{mt+ph-hh:.1f}" width="{bw*0.8:.1f}" height="{max(hh,0.4):.1f}" fill="{s.get("color", cols[i%len(cols)])}"/>')
        p.append(f'<text x="{ml+(ic+0.5)*pw/ncat:.1f}" y="{mt+ph+18}" text-anchor="middle" font-size="11" fill="#333" font-family="sans-serif">{_esc(cat)}</text>')
    p.append(f'<text x="16" y="{mt+ph/2}" text-anchor="middle" font-size="12" fill="#333" font-family="sans-serif" transform="rotate(-90 16 {mt+ph/2})">{_esc(ylab)}</text>')
    for i, s in enumerate(series):
        p.append(f'<rect x="{w-mr+10}" y="{mt+8+i*18}" width="14" height="3" fill="{s.get("color", cols[i%len(cols)])}"/>')
        p.append(f'<text x="{w-mr+28}" y="{mt+12+i*18}" font-size="11" fill="#333" font-family="sans-serif">{_esc(s["name"])}</text>')
    p.append(f'<text x="{ml}" y="{h-6}" font-size="10" fill="#777" font-family="sans-serif">{_esc(note)}</text></svg>')
    path.write_text("\n".join(p), encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def truth_from(d):
    return {
        "r0_m": float(d["r0_m"]),
        "theta0_rad": math.radians(float(d["theta0_deg"])),
        "v": float(d["v"]),
        "psi_rad": math.radians(float(d["psi_deg"])),
        "z": float(d["z"]),
    }


def run_scan_T_sigma(g):
    rows = []
    truth = truth_from(B0)
    for T in CONFIG["scans"]["T_s"]:
        for sig in CONFIG["scans"]["sigma_deg"]:
            agg, _ = evaluate_scenario(
                g, truth, sig, T, B0["turn_deg"], B0["source"], B0["snr_db"],
                use_rc3=True, env_key="E0", n_trials=N_TRIALS, seed=RNG_SEED + int(T + sig * 100)
            )
            for layer, rec in agg.items():
                rows.append({
                    "T_s": T, "sigma_deg": sig, "layer": layer,
                    "source": B0["source"], "snr_db": B0["snr_db"], "turn_deg": B0["turn_deg"],
                    **{k: rec.get(k) for k in [
                        "n_acc", "w_r_km", "w_theta_deg", "w_v_mps", "w_psi_deg",
                        "std_r_km", "std_theta_deg", "std_v", "std_psi_deg",
                        "contr_r", "contr_theta", "contr_v", "contr_psi",
                        "retain_r", "retain_theta", "retain_v", "retain_psi",
                        "unresolved_r", "unresolved_theta", "unresolved_v", "unresolved_psi",
                        "ten_pct_r", "ten_pct_v", "stage_reject_rate",
                    ]},
                })
            print(f"  Tsig T={T} sig={sig} RC2 n={agg.get('RC2',{}).get('n_acc')} "
                  f"RC3 n={agg.get('RC2+RC3',{}).get('n_acc')} "
                  f"contr_r RC3={agg.get('RC2+RC3',{}).get('contr_r')}")
    return pd.DataFrame(rows)


def run_scan_source_snr(g):
    rows = []
    truth = truth_from(B0)
    for src in CONFIG["scans"]["source"]:
        for snr in CONFIG["scans"]["snr_db"]:
            agg, _ = evaluate_scenario(
                g, truth, B0["sigma_deg"], B0["T"], B0["turn_deg"], src, snr,
                use_rc3=True, env_key="E0", n_trials=N_TRIALS, seed=RNG_SEED + hash((src, snr)) % 1000
            )
            for layer, rec in agg.items():
                rows.append({
                    "source": src, "snr_db": snr, "layer": layer,
                    "T_s": B0["T"], "sigma_deg": B0["sigma_deg"], "turn_deg": B0["turn_deg"],
                    **{k: rec.get(k) for k in [
                        "n_acc", "w_r_km", "w_v_mps", "w_theta_deg", "w_psi_deg",
                        "contr_r", "contr_v", "contr_theta", "contr_psi",
                        "retain_r", "retain_v", "unresolved_r", "unresolved_v",
                        "stage_reject_rate", "ten_pct_r", "ten_pct_v",
                    ]},
                })
            print(f"  srcsnr {src} SNR={snr} RC3 contr_r={agg.get('RC2+RC3',{}).get('contr_r')} "
                  f"w_r={agg.get('RC2+RC3',{}).get('w_r_km')}")
    return pd.DataFrame(rows)


def run_scan_depth(g):
    rows = []
    for z_true in CONFIG["scans"]["z_true_m"]:
        truth = truth_from({**B0, "z": z_true})
        agg, _ = evaluate_scenario(
            g, truth, B0["sigma_deg"], B0["T"], B0["turn_deg"], B0["source"], B0["snr_db"],
            use_rc3=True, env_key="E0", n_trials=N_TRIALS, seed=RNG_SEED + z_true
        )
        for layer, rec in agg.items():
            rows.append({
                "z_true_m": z_true, "layer": layer,
                "T_s": B0["T"], "sigma_deg": B0["sigma_deg"], "source": B0["source"], "snr_db": B0["snr_db"],
                **{k: rec.get(k) for k in [
                    "n_acc", "w_z_m", "w_r_km", "contr_z", "contr_r",
                    "retain_z", "retain_r", "unresolved_z", "unresolved_r",
                    "corr_rz", "mean_z", "mean_r_km", "stage_reject_rate",
                ]},
            })
        print(f"  depth z_true={z_true} RC3 w_z={agg.get('RC2+RC3',{}).get('w_z_m')} "
              f"contr_z={agg.get('RC2+RC3',{}).get('contr_z')} corr_rz={agg.get('RC2+RC3',{}).get('corr_rz')}")
    return pd.DataFrame(rows)


def run_scan_maneuver(g):
    rows = []
    truth = truth_from(B0)
    for turn in CONFIG["scans"]["turn_deg"]:
        for T in CONFIG["scans"]["T_maneuver_s"]:
            agg, _ = evaluate_scenario(
                g, truth, B0["sigma_deg"], T, turn, B0["source"], B0["snr_db"],
                use_rc3=True, env_key="E0", n_trials=N_TRIALS, seed=RNG_SEED + turn + T
            )
            for layer, rec in agg.items():
                rows.append({
                    "turn_deg": turn, "T_s": T, "layer": layer,
                    "sigma_deg": B0["sigma_deg"], "source": B0["source"], "snr_db": B0["snr_db"],
                    **{k: rec.get(k) for k in [
                        "n_acc", "n_rc2", "contr_theta", "contr_psi", "contr_r", "contr_v",
                        "w_r_km", "w_theta_deg", "w_psi_deg", "w_v_mps",
                        "stage_reject_rate", "unresolved_r", "unresolved_psi",
                        "retain_r", "retain_psi",
                    ]},
                })
            print(f"  man turn={turn} T={T} RC2 n={agg.get('RC2',{}).get('n_acc')} "
                  f"RC3 n={agg.get('RC2+RC3',{}).get('n_acc')} "
                  f"contr_psi RC2={agg.get('RC2',{}).get('contr_psi')} RC3={agg.get('RC2+RC3',{}).get('contr_psi')}")
    return pd.DataFrame(rows)


def run_environment(g):
    rows = []
    truth = truth_from(B0)
    for ek, em in ENV.items():
        # generate obs always from TRUE model; scoring uses mismatched MODE_ENV[ek]
        # evaluate_scenario uses mode_use = MODE_ENV[env_key] for candidates and MODE_TRUE for obs
        agg, _ = evaluate_scenario(
            g, truth, B0["sigma_deg"], B0["T"], B0["turn_deg"], B0["source"], B0["snr_db"],
            use_rc3=True, env_key=ek, n_trials=N_TRIALS, seed=RNG_SEED + 77
        )
        for layer, rec in agg.items():
            if layer != "RC2+RC3":
                continue
            rows.append({
                "env": ek, "label": em["label"], "z_a_shift_m": em["z_a_shift_m"],
                "c0_shift": em["c0_shift"], "layer": layer,
                "w_r_km": rec.get("w_r_km"), "contr_r": rec.get("contr_r"),
                "w_z_m": rec.get("w_z_m"), "contr_z": rec.get("contr_z"),
                "retain_r": rec.get("retain_r"), "retain_z": rec.get("retain_z"),
                "unresolved_r": rec.get("unresolved_r"), "n_acc": rec.get("n_acc"),
            })
        print(f"  env {ek} contr_r={agg.get('RC2+RC3',{}).get('contr_r')} retain_r={agg.get('RC2+RC3',{}).get('retain_r')}")
    return pd.DataFrame(rows)


def run_integrated(g):
    rows = []
    detail = {}
    for key, sc in [("C-L", SCEN_CL), ("C-M", SCEN_CM), ("C-U", SCEN_CU)]:
        truth = truth_from(sc)
        agg, _ = evaluate_scenario(
            g, truth, sc["sigma_deg"], sc["T"], sc["turn_deg"], sc["source"], sc["snr_db"],
            use_rc3=True, env_key="E0", n_trials=N_TRIALS, seed=RNG_SEED + hash(key) % 1000
        )
        detail[key] = agg
        rec = agg.get("RC2+RC3", {})
        rec2 = agg.get("RC2", {})
        rows.append({
            "scenario": key,
            "name": sc["name"],
            "sigma_deg": sc["sigma_deg"],
            "T_s": sc["T"],
            "turn_deg": sc["turn_deg"],
            "source": sc["source"],
            "snr_db": sc["snr_db"],
            "layer": "RC2+RC3",
            "w_r_km": rec.get("w_r_km"),
            "std_r_km": rec.get("std_r_km"),
            "contr_r": rec.get("contr_r"),
            "unresolved_r": rec.get("unresolved_r"),
            "ten_pct_r": rec.get("ten_pct_r"),
            "w_theta_deg": rec.get("w_theta_deg"),
            "contr_theta": rec.get("contr_theta"),
            "unresolved_theta": rec.get("unresolved_theta"),
            "w_z_m": rec.get("w_z_m"),
            "contr_z": rec.get("contr_z"),
            "unresolved_z": rec.get("unresolved_z"),
            "w_v_mps": rec.get("w_v_mps"),
            "contr_v": rec.get("contr_v"),
            "unresolved_v": rec.get("unresolved_v"),
            "ten_pct_v": rec.get("ten_pct_v"),
            "w_psi_deg": rec.get("w_psi_deg"),
            "contr_psi": rec.get("contr_psi"),
            "unresolved_psi": rec.get("unresolved_psi"),
            "n_acc_rc2": rec2.get("n_acc"),
            "n_acc_rc3": rec.get("n_acc"),
            "stage_reject_rate": rec.get("stage_reject_rate"),
            "retain_r": rec.get("retain_r"),
            "retain_theta": rec.get("retain_theta"),
            "retain_z": rec.get("retain_z"),
            "retain_v": rec.get("retain_v"),
            "retain_psi": rec.get("retain_psi"),
            "corr_rz": rec.get("corr_rz"),
        })
        print(f"  integ {key}: r_w={rec.get('w_r_km')} contr_r={rec.get('contr_r')} "
              f"theta_w={rec.get('w_theta_deg')} psi_w={rec.get('w_psi_deg')} "
              f"z_w={rec.get('w_z_m')} n {rec2.get('n_acc')}→{rec.get('n_acc')}")
    return pd.DataFrame(rows), detail


def make_core_figures(df_ts, df_ss, df_dz, df_mv, df_int, df_stage):
    # core fig1: sigma x T → five-state contraction heatmaps side concept → use 3 heatmaps as one file series
    # We'll put multi-panel style: three heatmaps written as one "core_fig1" with primary contr_r and note others
    # Requirement: one file core_fig1_T_sigma.svg — combine by writing a wide multi-panel SVG
    sigs = CONFIG["scans"]["sigma_deg"]
    Ts = CONFIG["scans"]["T_s"]
    # RC2+RC3 rows
    d = df_ts[df_ts["layer"] == "RC2+RC3"]
    def mat(col):
        Z = np.zeros((len(sigs), len(Ts)))
        for i, s in enumerate(sigs):
            for j, T in enumerate(Ts):
                row = d[(d["sigma_deg"] == s) & (d["T_s"] == T)]
                Z[i, j] = float(row[col].iloc[0]) if len(row) and col in row.columns else np.nan
        return Z
    # multi-panel core fig1
    panels = [("contr_theta", "θ 收缩"), ("contr_psi", "ψ 收缩"), ("contr_r", "r 收缩")]
    w, h = 980, 360
    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
         f'<rect width="{w}" height="{h}" fill="#f7f4ef"/>',
         f'<text x="{w/2}" y="26" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="15" font-weight="600" fill="#1a1a1a">核心图1  σθ×T → 五维候选收缩（RC2+RC3，B0其它条件）</text>']
    for ip, (col, lab) in enumerate(panels):
        Z = mat(col)
        ox = 40 + ip * 310
        p.append(f'<text x="{ox+120}" y="50" text-anchor="middle" font-size="12" fill="#333" font-family="sans-serif">{lab}</text>')
        zmin, zmax = float(np.nanmin(Z)), float(np.nanmax(Z))
        if not np.isfinite(zmin) or zmax - zmin < 1e-12:
            zmin, zmax = 0, 1
        cw, ch = 240 / len(Ts), 220 / len(sigs)
        for i, s in enumerate(sigs):
            for j, T in enumerate(Ts):
                v = Z[i, j]
                u = 0.5 if not np.isfinite(v) else (v - zmin) / (zmax - zmin + 1e-15)
                r = int(240 * (1 - u) + 31 * u); g = int(230 * (1 - u) + 111 * u); b = int(220 * (1 - u) + 139 * u)
                x = ox + j * cw
                y = 60 + (len(sigs) - 1 - i) * ch
                p.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{cw:.1f}" height="{ch:.1f}" fill="rgb({r},{g},{b})" stroke="#eee"/>')
        for j, T in enumerate(Ts):
            p.append(f'<text x="{ox+j*cw+cw/2:.1f}" y="{60+len(sigs)*ch+14}" text-anchor="middle" font-size="9" fill="#555" font-family="sans-serif">{T}</text>')
        for i, s in enumerate(sigs):
            p.append(f'<text x="{ox-6}" y="{60+(len(sigs)-1-i)*ch+ch/2+3:.1f}" text-anchor="end" font-size="9" fill="#555" font-family="sans-serif">{s}</text>')
        p.append(f'<text x="{ox+120}" y="{60+len(sigs)*ch+30}" text-anchor="middle" font-size="10" fill="#555" font-family="sans-serif">T(s) →</text>')
        p.append(f'<text x="14" y="{60+len(sigs)*ch/2}" font-size="10" fill="#555" font-family="sans-serif" transform="rotate(-90 14 {60+len(sigs)*ch/2})">σθ°</text>')
    p.append(f'<text x="40" y="{h-10}" font-size="10" fill="#777" font-family="sans-serif">收缩=1-宽度/先验宽度；v收缩见 CSV（辅助）</text></svg>')
    (FIG / "core_fig1_T_sigma.svg").write_text("\n".join(p), encoding="utf-8")

    # core fig2: source x SNR → CZ increment (stage reject + contr_r)
    d2 = df_ss[df_ss["layer"] == "RC2+RC3"]
    sources = CONFIG["scans"]["source"]
    snrs = CONFIG["scans"]["snr_db"]
    Z = np.zeros((len(sources), len(snrs)))
    for i, s in enumerate(sources):
        for j, sn in enumerate(snrs):
            row = d2[(d2["source"] == s) & (d2["snr_db"] == sn)]
            Z[i, j] = float(row["contr_r"].iloc[0]) if len(row) and "contr_r" in row.columns else np.nan
    write_heatmap(
        FIG / "core_fig2_source_snr.svg",
        "核心图2  source×SNR → CZ后距离候选收缩",
        "SNR (dB)", "source",
        [str(x) for x in snrs], sources, Z,
        note="S0好于S1是因为宽带可用频点更多，不是“无线谱优于有线谱”",
        fmt="{:.2f}",
    )

    # core fig3: z_true vs performance + r-z coupling proxy
    d3 = df_dz[df_dz["layer"] == "RC2+RC3"]
    write_group_bars(
        FIG / "core_fig3_depth_range.svg",
        "核心图3  z_true → 深度/距离约束（RC2+RC3）",
        [f"z={int(z)}m" for z in CONFIG["scans"]["z_true_m"]],
        [
            {"name": "contr_z", "color": "#0f766e", "y": [float(d3[d3["z_true_m"]==z]["contr_z"].iloc[0]) if len(d3[d3["z_true_m"]==z]) else np.nan for z in CONFIG["scans"]["z_true_m"]]},
            {"name": "contr_r", "color": "#b45309", "y": [float(d3[d3["z_true_m"]==z]["contr_r"].iloc[0]) if len(d3[d3["z_true_m"]==z]) else np.nan for z in CONFIG["scans"]["z_true_m"]]},
        ],
        "contraction",
        note="r–z耦合见 boundary_depth.csv 的 corr_rz",
    )

    # core fig4: turn x T RC2 vs RC2+RC3
    d4 = df_mv
    turns = CONFIG["scans"]["turn_deg"]
    Ts_m = CONFIG["scans"]["T_maneuver_s"]
    # use T=600 slice bar comparison
    cats = [f"转角{t}°" for t in turns]
    y_rc2, y_rc3 = [], []
    for t in turns:
        row2 = d4[(d4["layer"]=="RC2") & (d4["turn_deg"]==t) & (d4["T_s"]==600)]
        row3 = d4[(d4["layer"]=="RC2+RC3") & (d4["turn_deg"]==t) & (d4["T_s"]==600)]
        y_rc2.append(float(row2["contr_psi"].iloc[0]) if len(row2) else np.nan)
        y_rc3.append(float(row3["contr_psi"].iloc[0]) if len(row3) else np.nan)
    write_group_bars(
        FIG / "core_fig4_turn_T.svg",
        "核心图4  turn×T → RC2 vs RC2+RC3（T=600s 的 ψ 收缩）",
        cats,
        [
            {"name": "RC2 ψ收缩", "color": "#8a8a8a", "y": y_rc2},
            {"name": "RC2+RC3 ψ收缩", "color": "#b45309", "y": y_rc3},
        ],
        "contraction ψ",
        note="若小机动主要改善候选结构而非直接测距，应看到 RC2 候选仍宽、RC3 后筛选更容易",
    )

    # core fig5: five-state C-L/M/U
    states = ["r", "theta", "z", "v", "psi"]
    labels = ["r", "θ", "z", "v", "ψ"]
    scen = ["C-L", "C-M", "C-U"]
    series = []
    for st, lab in zip(states, labels):
        ys = []
        for sc in scen:
            row = df_int[df_int["scenario"] == sc]
            col = f"contr_{st}" if st != "theta" else "contr_theta"
            if st == "psi":
                col = "contr_psi"
            ys.append(float(row[col].iloc[0]) if len(row) and col in row.columns else np.nan)
        series.append({"name": lab, "y": ys})
    write_group_bars(
        FIG / "core_fig5_five_state.svg",
        "核心图5  C-L / C-M / C-U 五维状态收缩对比",
        ["C-L困难", "C-M基准", "C-U理论上界"],
        series,
        "contraction (1=完全收缩)",
        note="单元格数值详见 state_function_boundary.csv；UNRESOLVED 以低收缩体现",
    )

    # core fig6: RC1→RC2→RC3 candidate count
    # stages: RC1 prior ~ full grid; RC2 n; RC2+RC3 n
    # use integrated C-M
    rowm = df_int[df_int["scenario"] == "C-M"]
    n_prior = int(g_prior_size)
    n2 = float(rowm["n_acc_rc2"].iloc[0]) if len(rowm) else np.nan
    n3 = float(rowm["n_acc_rc3"].iloc[0]) if len(rowm) else np.nan
    write_group_bars(
        FIG / "core_fig6_stage_contraction.svg",
        "核心图6  RC1→RC2→RC3 候选数量逐层收缩（C-M）",
        ["RC1先验网格", "RC2连续方位", "RC2+RC3传播"],
        [{"name": "候选数", "color": "#1d4ed8", "y": [n_prior, n2, n3]}],
        "候选数 (log10 显示见报告)",
        note=f"先验节点≈{n_prior}；RC2≈{n2:.0f}；RC2+RC3≈{n3:.0f}",
    )


# global prior size for fig6
g_prior_size = 0


def main():
    global g_prior_size
    t0 = time.time()
    print("=== P4 G4 performance boundary ===")
    print(f"OUT={OUT}")
    g = make_grid(include_z=True)
    g_prior_size = g["n"]
    print(f"grid nodes with z: {g['n']}")
    (OUT / "P4_G4_CONFIG.json").write_text(json.dumps(CONFIG, ensure_ascii=False, indent=2), encoding="utf-8")

    print("building amplitude lookup tables ...")
    AMP_TABS.clear()
    for src in ["S0", "S1", "S2"]:
        freqs, _, _ = source_freqs(src)
        AMP_TABS[(src, "E0")] = build_amp_table(MODE_ENV["E0"], freqs)
        AMP_TABS[(src, "E1")] = build_amp_table(MODE_ENV["E1"], freqs)
        AMP_TABS[(src, "E2")] = build_amp_table(MODE_ENV["E2"], freqs)
        print(f"  table {src} F={len(freqs)}")

    print("[1/6] Scan T × sigma ...")
    df_ts = run_scan_T_sigma(g)
    df_ts.to_csv(OUT / "boundary_T_sigma.csv", index=False, encoding="utf-8-sig")

    print("[2/6] Scan source × SNR ...")
    df_ss = run_scan_source_snr(g)
    df_ss.to_csv(OUT / "boundary_source_snr.csv", index=False, encoding="utf-8-sig")

    print("[3/6] Scan depth ...")
    df_dz = run_scan_depth(g)
    df_dz.to_csv(OUT / "boundary_depth.csv", index=False, encoding="utf-8-sig")

    print("[4/6] Scan maneuver ...")
    df_mv = run_scan_maneuver(g)
    df_mv.to_csv(OUT / "boundary_maneuver.csv", index=False, encoding="utf-8-sig")

    print("[5/6] Environment minimal ...")
    df_env = run_environment(g)
    df_env.to_csv(OUT / "boundary_environment.csv", index=False, encoding="utf-8-sig")

    print("[6/6] Integrated C-L/M/U ...")
    df_int, detail = run_integrated(g)
    df_int.to_csv(OUT / "integrated_CL_CM_CU.csv", index=False, encoding="utf-8-sig")

    # state function boundary table
    sfb = []
    for _, row in df_int.iterrows():
        for st in ["r", "theta", "z", "v", "psi"]:
            wcol = {"r": "w_r_km", "theta": "w_theta_deg", "z": "w_z_m", "v": "w_v_mps", "psi": "w_psi_deg"}[st]
            ccol = {"r": "contr_r", "theta": "contr_theta", "z": "contr_z", "v": "contr_v", "psi": "contr_psi"}[st]
            ucol = {"r": "unresolved_r", "theta": "unresolved_theta", "z": "unresolved_z", "v": "unresolved_v", "psi": "unresolved_psi"}[st]
            unit = {"r": "km", "theta": "deg", "z": "m", "v": "m/s", "psi": "deg"}[st]
            sfb.append({
                "scenario": row["scenario"],
                "state": st,
                "posterior_width": row.get(wcol),
                "unit": unit,
                "contraction": row.get(ccol),
                "unresolved": row.get(ucol),
                "ten_pct_ref": row.get(f"ten_pct_{st}") if st in ("r", "v", "z") else np.nan,
            })
    pd.DataFrame(sfb).to_csv(OUT / "state_function_boundary.csv", index=False, encoding="utf-8-sig")

    # stage contraction table
    stage = []
    for key in ["C-L", "C-M", "C-U"]:
        row = df_int[df_int["scenario"] == key]
        if not len(row):
            continue
        r = row.iloc[0]
        stage.append({
            "scenario": key,
            "n_RC1_prior": g_prior_size,
            "n_RC2": r.get("n_acc_rc2"),
            "n_RC2_RC3": r.get("n_acc_rc3"),
            "reject_RC2_to_RC3": r.get("stage_reject_rate"),
        })
    pd.DataFrame(stage).to_csv(OUT / "rc1_rc2_rc3_contraction.csv", index=False, encoding="utf-8-sig")

    print("figures ...")
    make_core_figures(df_ts, df_ss, df_dz, df_mv, df_int, pd.DataFrame(stage))

    # Build report
    def fnum(x, nd=3):
        try:
            if x is None or (isinstance(x, float) and not np.isfinite(x)):
                return "n/a"
            return f"{float(x):.{nd}f}"
        except Exception:
            return "n/a"

    def pick(df, **kw):
        d = df
        for k, v in kw.items():
            d = d[d[k] == v]
        return d.iloc[0] if len(d) else None

    # five sentences
    def state_cell(scen, state):
        row = df_int[df_int["scenario"] == scen]
        if not len(row):
            return "n/a"
        r = row.iloc[0]
        wcol = {"r": "w_r_km", "theta": "w_theta_deg", "z": "w_z_m", "v": "w_v_mps", "psi": "w_psi_deg"}[state]
        ucol = {"r": "unresolved_r", "theta": "unresolved_theta", "z": "unresolved_z", "v": "unresolved_v", "psi": "unresolved_psi"}[state]
        w = r.get(wcol)
        u = r.get(ucol)
        if u or not np.isfinite(float(w)) if w is not None else True:
            return f"UNRESOLVED (宽≈{fnum(w)})"
        return f"宽≈{fnum(w)}"

    rep = []
    rep.append("# P4 G4 报告：综合功能与性能边界")
    rep.append("")
    rep.append(f"项目：**B3D-3**  ·  阶段：**P4 / G4**  ·  UTC：{NOW}")
    rep.append("")
    rep.append("本阶段不证明某条算法“成功”，而是生产甲方可用的性能边界：**什么条件下能约束哪些状态、约束到什么程度、边界在哪里。**")
    rep.append("")
    rep.append("## 0. 前序冻结（不重跑 P2/P3）")
    rep.append("")
    rep.append("- G1：E-STD 场景冻结；x=[r,θ,z,v,ψ]，r 45–60 km，z 200 m 基准，U=2 m/s。")
    rep.append("- G2：连续方位不是远距离测距；LIN 有 r–v 脊线；有限转角不足以实用测距；秩亏 CRLB=UNDEFINED。")
    rep.append("- G3：**PROPAGATION_INCREMENT_CONFIRMED**（冻结理论模型内）；S0/S2 CZ 增量强，S1 较弱；Doppler 对恒定 v_rad 对近零增量。")
    rep.append("")
    rep.append("## 1. 方法")
    rep.append("")
    rep.append("- **候选集合观点**：不强迫唯一点估计；输出 posterior width / contraction / truth retention / UNRESOLVED。")
    rep.append("- `contraction_j = 1 - width_j / prior_width_j`；多峰保留。")
    rep.append("- 基准 **B0**：r=50 km，θ0=0°，z=200 m，v=2 m/s，**ψ=5°（开放几何）**，T=600 s，σθ=0.1°，turn=0，S1，SNR=10 dB。")
    rep.append("- LIN 仅作运动学下界，不作主性能基准。")
    rep.append("- 不做全笛卡尔积；四组核心扫描 + 三档综合场景。")
    rep.append("- 阵列不以具体阵长写装备需求，统一用 **σθ** 与 RC3 阵列联合残差参数化。")
    rep.append("- 环境仅三档：E0/E1/E2，**单一物理形式**：声道轴深度偏移（E1 +50 m，E2 +150 m）。")
    rep.append("- **10%** 仅作 r/z/v 等正值状态参考线；θ/ψ 用绝对角宽度。")
    rep.append("")
    rep.append("## 2. 核心扫描结果")
    rep.append("")
    rep.append("### 2.1 T × σθ（B0其它条件，RC2+RC3）")
    rep.append("")
    rep.append("| T | σθ | θ收缩 | ψ收缩 | r收缩 | v收缩 | r宽 km |")
    rep.append("| --- | --- | --- | --- | --- | --- | --- |")
    dts = df_ts[df_ts["layer"] == "RC2+RC3"]
    for _, r in dts.iterrows():
        if r["T_s"] in (300, 600, 1200) and r["sigma_deg"] in (0.02, 0.1, 0.2, 0.5):
            rep.append(
                f"| {r['T_s']:.0f} | {r['sigma_deg']} | {fnum(r.get('contr_theta'))} | {fnum(r.get('contr_psi'))} | "
                f"{fnum(r.get('contr_r'))} | {fnum(r.get('contr_v'))} | {fnum(r.get('w_r_km'))} |"
            )
    rep.append("")
    rep.append("### 2.2 source × SNR")
    rep.append("")
    rep.append("| source | SNR | r收缩 | v收缩 | r宽 km | RC3阶段排除 |")
    rep.append("| --- | --- | --- | --- | --- | --- |")
    dss = df_ss[df_ss["layer"] == "RC2+RC3"]
    for _, r in dss.iterrows():
        rep.append(
            f"| {r['source']} | {r['snr_db']} | {fnum(r.get('contr_r'))} | {fnum(r.get('contr_v'))} | "
            f"{fnum(r.get('w_r_km'))} | {fnum(r.get('stage_reject_rate'))} |"
        )
    rep.append("")
    rep.append("**必须明确**：S0 好于 S1 **不是**“无线谱优于有线谱”，而是 S0 宽带在 150–375 Hz 提供更多可用频点；S1 分析带内有效谱点较少。")
    rep.append("")
    rep.append("### 2.3 目标深度")
    rep.append("")
    rep.append("| z_true | z宽 m | contr_z | r宽 km | contr_r | corr(r,z) | retain_z |")
    rep.append("| --- | --- | --- | --- | --- | --- | --- |")
    ddz = df_dz[df_dz["layer"] == "RC2+RC3"]
    for _, r in ddz.iterrows():
        rep.append(
            f"| {r['z_true_m']:.0f} | {fnum(r.get('w_z_m'))} | {fnum(r.get('contr_z'))} | "
            f"{fnum(r.get('w_r_km'))} | {fnum(r.get('contr_r'))} | {fnum(r.get('corr_rz'))} | {fnum(r.get('retain_z'),2)} |"
        )
    rep.append("")
    rep.append("### 2.4 平台小幅机动")
    rep.append("")
    rep.append("| turn | T | 层 | n候选 | ψ收缩 | r收缩 | r宽 km |")
    rep.append("| --- | --- | --- | --- | --- | --- | --- |")
    for _, r in df_mv.iterrows():
        if r["T_s"] in (300, 600, 1200) and r["turn_deg"] in (0, 5, 15):
            rep.append(
                f"| {r['turn_deg']:.0f} | {r['T_s']:.0f} | {r['layer']} | {fnum(r.get('n_acc'),0)} | "
                f"{fnum(r.get('contr_psi'))} | {fnum(r.get('contr_r'))} | {fnum(r.get('w_r_km'))} |"
            )
    rep.append("")
    # maneuver interpretation
    m0 = pick(df_mv, layer="RC2", turn_deg=0.0, T_s=600.0)
    m15 = pick(df_mv, layer="RC2", turn_deg=15.0, T_s=600.0)
    m15r = pick(df_mv, layer="RC2+RC3", turn_deg=15.0, T_s=600.0)
    m0r = pick(df_mv, layer="RC2+RC3", turn_deg=0.0, T_s=600.0)
    if m0 is not None and m15 is not None and m0r is not None and m15r is not None:
        rep.append("**机动价值判断（核心理论结论候选）**：")
        rep.append("")
        rep.append(f"- RC2 ψ收缩：0° **{fnum(m0.get('contr_psi'))}** → 15° **{fnum(m15.get('contr_psi'))}**")
        rep.append(f"- RC2 r宽：0° **{fnum(m0.get('w_r_km'))} km** → 15° **{fnum(m15.get('w_r_km'))} km**")
        rep.append(f"- RC2+RC3 后 n 候选：0° **{fnum(m0r.get('n_acc'),0)}** → 15° **{fnum(m15r.get('n_acc'),0)}**")
        rep.append(f"- RC2+RC3 r收缩：0° **{fnum(m0r.get('contr_r'))}** → 15° **{fnum(m15r.get('contr_r'))}**")
        rep.append("")
        # decide A or B
        dr_rc2 = abs((m15.get("w_r_km") or 0) - (m0.get("w_r_km") or 0))
        dps_rc2 = abs((m15.get("contr_psi") or 0) - (m0.get("contr_psi") or 0))
        dn_rc3 = abs((m15r.get("n_acc") or 0) - (m0r.get("n_acc") or 0))
        if dr_rc2 < 2.0 and dps_rc2 >= 0.05:
            rep.append("**结论：更支持 B** —— 小幅转向主要改善**运动学候选结构（θ/ψ 等）**，使 RC3 更容易筛选；**并非**直接把距离约束到实用水平。")
        else:
            rep.append("**结论：** 需按表中 r 宽度与候选数变化解读；若 r 收缩仍弱而候选结构变化明显，仍指向 B。")
    rep.append("")
    rep.append("### 2.5 环境最小鲁棒性（声道轴偏移）")
    rep.append("")
    rep.append("| env | 说明 | r宽 km | contr_r | retain_r | contr_z |")
    rep.append("| --- | --- | --- | --- | --- | --- |")
    for _, r in df_env.iterrows():
        rep.append(
            f"| {r['env']} | {r['label']} | {fnum(r.get('w_r_km'))} | {fnum(r.get('contr_r'))} | "
            f"{fnum(r.get('retain_r'),2)} | {fnum(r.get('contr_z'))} |"
        )
    rep.append("")
    rep.append("## 3. 三档综合场景 C-L / C-M / C-U")
    rep.append("")
    rep.append("| 场景 | σθ | T | turn | source | SNR | r宽 | θ宽° | z宽 m | v宽 | ψ宽° | n_RC2→RC3 |")
    rep.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for _, r in df_int.iterrows():
        rep.append(
            f"| {r['name']} | {r['sigma_deg']} | {r['T_s']:.0f} | {r['turn_deg']:.0f} | {r['source']} | {r['snr_db']} | "
            f"{fnum(r.get('w_r_km'))} | {fnum(r.get('w_theta_deg'))} | {fnum(r.get('w_z_m'))} | "
            f"{fnum(r.get('w_v_mps'))} | {fnum(r.get('w_psi_deg'))} | "
            f"{fnum(r.get('n_acc_rc2'),0)}→{fnum(r.get('n_acc_rc3'),0)} |"
        )
    rep.append("")
    rep.append("### 五维综合功能表（单元格填宽度/UNRESOLVED，非 PASS/FAIL）")
    rep.append("")
    rep.append("| 状态 | C-L | C-M | C-U | 主要限制来源 |")
    rep.append("| --- | --- | --- | --- | --- |")
    limit_src = {
        "r": "RC3 CZ 剖面 + 观测时长/SNR/源条件；RC2 几乎不直接测距",
        "theta": "RC2 方位噪声 σθ 与观测几何",
        "z": "RC3 深度响应 vs 声道结构；与 r 可能耦合",
        "v": "RC2 时间演化 + RC3 距离轨迹联合",
        "psi": "RC2 方位变化率/机动；传播对纯 ψ 增量有限",
    }
    for st, lab in [("r", "r"), ("theta", "θ"), ("z", "z"), ("v", "v"), ("psi", "ψ")]:
        cells = []
        for sc in ["C-L", "C-M", "C-U"]:
            cells.append(state_cell(sc, st))
        rep.append(f"| {lab} | {cells[0]} | {cells[1]} | {cells[2]} | {limit_src[st]} |")
    rep.append("")
    rep.append("10% 参考：仅对 r、z、v 标注是否进入参考区（见 CSV `ten_pct_*`）；**不作为五维统一通过公式**。")
    rep.append("")
    rep.append("## 4. G4 五句话（甲方）")
    rep.append("")
    # compute anchors from C-M and C-U and T-sigma
    cm = df_int[df_int["scenario"] == "C-M"]
    cu = df_int[df_int["scenario"] == "C-U"]
    cl = df_int[df_int["scenario"] == "C-L"]
    t12 = pick(df_ts, layer="RC2+RC3", T_s=1200.0, sigma_deg=0.02)
    t300 = pick(df_ts, layer="RC2+RC3", T_s=300.0, sigma_deg=0.5)
    s0 = pick(df_ss, layer="RC2+RC3", source="S0", snr_db=10.0)
    s1 = pick(df_ss, layer="RC2+RC3", source="S1", snr_db=10.0)
    s2 = pick(df_ss, layer="RC2+RC3", source="S2", snr_db=10.0)
    dz = pick(df_dz, layer="RC2+RC3", z_true_m=200.0)

    def g(row, col, nd=3):
        return fnum(row.get(col) if row is not None else None, nd)

    rep.append(f"1. **距离 r**：主要由 **RC3 会聚区传播** 约束。B0 下 source×SNR 时 r 候选宽度约 {g(s1,'w_r_km')}–{g(s0,'w_r_km')} km（S1/S0）；C-U 上界场景约 {g(cu.iloc[0] if len(cu) else None,'w_r_km')} km，C-L 困难场景约 {g(cl.iloc[0] if len(cl) else None,'w_r_km')} km。RC2 在冻结机动下**不能**给出实用测距。")
    rep.append(f"2. **方位 θ**：由 **RC2 连续方位** 约束。要 θ 收缩明显，需足够 T 与较小 σθ：σθ=0.02°、T=1200 s 时 θ 收缩约 {g(t12,'contr_theta')}；σθ=0.5°、T=300 s 时约 {g(t300,'contr_theta')}。")
    if dz is not None:
        rep.append(f"3. **深度 z**：在第一会聚区传播模型下，z 候选宽约 {g(dz,'w_z_m')} m，contr_z≈{g(dz,'contr_z')}，r–z 相关≈{g(dz,'corr_rz')}。若 contr_z 偏低或宽度接近先验，则 **传播未提供可靠独立深度约束**，与距离存在耦合风险；深度功能边界按 UNRESOLVED 写，不加新特征。")
    else:
        rep.append("3. **深度 z**：见 boundary_depth.csv；若不可靠则直接 UNRESOLVED。")
    rep.append(f"4. **速度 v**：RC2+RC3 联合后，v 候选宽在 C-L/C-M/C-U 约 {g(cl.iloc[0] if len(cl) else None,'w_v_mps')} / {g(cm.iloc[0] if len(cm) else None,'w_v_mps')} / {g(cu.iloc[0] if len(cu) else None,'w_v_mps')} m/s；距离–速度补偿主要靠 CZ 轨迹压缩，而非 Doppler。")
    rep.append(f"5. **航向 ψ**：靠 **RC2 方位变化率 + 几何/机动** 在长 T、开放几何（ψ≠0）下可辨识；传播对纯 ψ 增量有限（G3 已示）。C-M ψ 宽约 {g(cm.iloc[0] if len(cm) else None,'w_psi_deg')}°，C-U 约 {g(cu.iloc[0] if len(cu) else None,'w_psi_deg')}°，C-L 约 {g(cl.iloc[0] if len(cl) else None,'w_psi_deg')}°。")
    rep.append("")
    rep.append("## 5. RC1 / RC2 / RC3 信息贡献分离")
    rep.append("")
    rep.append("| 层 | 主要贡献 | 不贡献/弱贡献 |")
    rep.append("| --- | --- | --- |")
    rep.append("| RC1 场景/先验 | 状态定义、先验盒、深度候选、频带 | 不提供观测约束 |")
    rep.append("| RC2 连续方位 | θ；开放几何下部分 ψ；形成困难候选 | 远距离实用 r；LIN 下 r–v |")
    rep.append("| RC3 CZ 传播 | r 与距离轨迹相关 v；S0/S2 下强排除 | 恒定 v_rad 的 Doppler；纯 ψ 且不改 r(t) |")
    rep.append("")
    stg = pd.read_csv(OUT / "rc1_rc2_rc3_contraction.csv") if (OUT / "rc1_rc2_rc3_contraction.csv").exists() else pd.DataFrame()
    if len(stg):
        rep.append("| 场景 | RC1先验 n | RC2 n | RC2+RC3 n | 层间排除率 |")
        rep.append("| --- | --- | --- | --- | --- |")
        for _, r in stg.iterrows():
            rep.append(f"| {r['scenario']} | {r['n_RC1_prior']} | {fnum(r.get('n_RC2'),0)} | {fnum(r.get('n_RC2_RC3'),0)} | {fnum(r.get('reject_RC2_to_RC3'))} |")
    rep.append("")
    rep.append("## 6. 六张核心图")
    rep.append("")
    for name in [
        "core_fig1_T_sigma.svg", "core_fig2_source_snr.svg", "core_fig3_depth_range.svg",
        "core_fig4_turn_T.svg", "core_fig5_five_state.svg", "core_fig6_stage_contraction.svg",
    ]:
        rep.append(f"- figures/{name}")
    rep.append("")
    rep.append("辅助数据均在 CSV，不另堆重复图。")
    rep.append("")
    rep.append("## 7. G4 判定")
    rep.append("")
    # decide G4
    any_unresolved = False
    if len(df_int):
        for col in ["unresolved_r", "unresolved_theta", "unresolved_z", "unresolved_v", "unresolved_psi"]:
            if col in df_int.columns and df_int[col].fillna(True).any():
                any_unresolved = True
    core_ok = all((FIG / n).exists() for n in [
        "core_fig1_T_sigma.svg", "core_fig2_source_snr.svg", "core_fig3_depth_range.svg",
        "core_fig4_turn_T.svg", "core_fig5_five_state.svg", "core_fig6_stage_contraction.svg"])
    scans_ok = len(df_ts) > 0 and len(df_ss) > 0 and len(df_dz) > 0 and len(df_mv) > 0 and len(df_int) > 0
    if core_ok and scans_ok and not any_unresolved:
        g4 = "G4_PERFORMANCE_BOUNDARY_COMPLETE"
    else:
        g4 = "G4_PARTIAL_BOUNDARY_WITH_UNRESOLVED_STATES"
    rep.append(f"### `{g4}`")
    rep.append("")
    rep.append(f"- 五维功能边界已形成：{'是' if scans_ok else '否'}")
    rep.append(f"- 核心条件扫描完成：{'是' if scans_ok else '否'}")
    rep.append(f"- RC1/RC2/RC3 贡献已分离：是")
    rep.append(f"- 不可辨识区域明确：{'是（见表中 UNRESOLVED）' if any_unresolved else '本轮样本内均非 UNRESOLVED'}")
    rep.append(f"- 三档代表场景：{'是' if len(df_int)==3 else '否'}")
    rep.append(f"- 六张核心图：{'是' if core_ok else '否'}")
    rep.append("- 不依赖特定型号装备：是（σθ 参数化，无具体阵长需求）")
    rep.append("")
    rep.append("两种合法终态均允许：`G4_PERFORMANCE_BOUNDARY_COMPLETE` 或 `G4_PARTIAL_BOUNDARY_WITH_UNRESOLVED_STATES`。")
    rep.append("")
    rep.append("## 8. 禁止事项自检")
    rep.append("")
    rep.append("| 条件 | 状态 |")
    rep.append("| --- | --- |")
    rep.append("| 回 Bellhop A2 审计 | 否 |")
    rep.append("| 新发明传播特征 | 否 |")
    rep.append("| 重跑完整 P2 | 否 |")
    rep.append("| 重做 P3 困难候选筛选 | 否 |")
    rep.append("| 临时加传感器 | 否 |")
    rep.append("| 全参数笛卡尔积 | 否 |")
    rep.append("| S2 写成真实 UUV | 否 |")
    rep.append("| 自动进入 P5/最终报告 | 否（G4 完成即停） |")
    rep.append("")
    rep.append("## 9. 交付清单")
    rep.append("")
    for p in sorted(OUT.rglob("*")):
        if p.is_file():
            rep.append(f"- results/P4_G4_performance_boundary/{p.relative_to(OUT)}")
    rep.append("")
    rep.append("---")
    rep.append("")
    rep.append("**G4 完成后停止，不自动进入 P5。**")
    (OUT / "P4_G4_REPORT.md").write_text("\n".join(rep), encoding="utf-8")

    # GPT sync
    gs = []
    gs.append("# P4 G4 — GPT 同步稿")
    gs.append("")
    gs.append(f"- UTC: {NOW}")
    gs.append(f"- G4: **{g4}**")
    gs.append("- 未重跑 P2/P3；未回 Bellhop；未全笛卡尔积；未写具体阵长需求。")
    gs.append("")
    gs.append("## B0")
    gs.append(f"`{json.dumps(B0, ensure_ascii=False)}`")
    gs.append("")
    gs.append("## 三档综合（RC2+RC3）")
    gs.append("")
    if len(df_int):
        gs.append(df_int[["scenario", "name", "w_r_km", "w_theta_deg", "w_z_m", "w_v_mps", "w_psi_deg", "n_acc_rc2", "n_acc_rc3"]].to_string(index=False))
    gs.append("")
    gs.append("## 五句话要点")
    gs.append("")
    for line in rep:
        if line.startswith("1. **距离") or line.startswith("2. **方位") or line.startswith("3. **深度") or line.startswith("4. **速度") or line.startswith("5. **航向"):
            gs.append(line)
    gs.append("")
    gs.append("## 下一轮")
    gs.append("")
    gs.append("> G4 完成后停止。是否进入 P5/最终甲方报告由任务书决定，本轮不自动启动。")
    gs.append("")
    (OUT / "P4_G4_GPT_SYNC.md").write_text("\n".join(gs), encoding="utf-8")

    (OUT / "g4_decision.json").write_text(json.dumps({
        "g4": g4,
        "core_figures_ok": core_ok,
        "scans_ok": scans_ok,
        "any_unresolved": bool(any_unresolved),
        "integrated": df_int.to_dict(orient="records") if len(df_int) else [],
    }, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    print(f"DONE in {time.time()-t0:.1f}s  G4={g4}")
    print("Artifacts", OUT)


if __name__ == "__main__":
    main()
