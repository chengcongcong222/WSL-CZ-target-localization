#!/usr/bin/env python3
"""P3 / G3: RC3 independent increment — CZ propagation + line/Doppler on RC2 hard pairs.

Frozen inputs:
  - E-STD deep-ocean theory scene (Munk-like SOFAR), z_s=z_r=200 m, band 150–375 Hz
  - Hard pairs ONLY from P2 handoff CSV (no re-picking after seeing acoustics)
  - Propagation: normal-mode-inspired deep-waveguide model (KRAKEN binary not on this host;
    model is deterministic theory-boundary stand-in, recorded in CONFIG)
  - Sources S0/S1/S2 as P1 frozen
  - Four observation layers; Doppler is control branch
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
P2_OUT = ROOT / "results" / "P2_RC2_kinematic_boundary"
OUT = ROOT / "results" / "P3_RC3_increment"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

# ---------------------------------------------------------------------------
# Frozen constants (P1 E-STD + P3 task)
# ---------------------------------------------------------------------------
U_PLAT = 2.0
C_WATER = 1500.0
Z_S = 200.0
Z_R = 200.0
Z_A = 1300.0  # Munk channel axis
B_MUNK = 1300.0
EPS_MUNK = 0.00737
C0_MUNK = 1500.0
DEPTH_M = 5000.0
F_BAND = (150.0, 375.0)
F_GRID = np.linspace(150.0, 375.0, 12)  # multi-freq samples in band
S2_LINES = [166.0, 201.0, 235.0, 283.0, 338.0]
S1_SHAFT_HZ = 4.0
S1_BLADES = 5
S1_DRIFT_FRAC = 0.01  # ±1% line frequency drift
SNR_LIST = [20.0, 10.0, 0.0, -10.0]
N_TRIALS_NOISE = 8
RNG_SEED = 303
DT_OBS = 10.0
N_MODES_BASE = 24
ARRAY_N = 8
ARRAY_D = 12.0  # m spacing
THRESHOLD_RATIO = 0.25  # wrong-candidate rejection if cost_alt > thresh * cost_true + margin

MECH_RULES = {
    "A": "r-v compensation (Δr large or Δv large, Δψ/Δθ0 small)",
    "B": "course near-mirror (Δψ large, Δr/Δv moderate)",
    "C": "theta0-psi joint (Δθ0 and/or mixed kinematic shift)",
}

CONFIG = {
    "project": "B3D-3",
    "package": "P3_RC3_increment",
    "created_utc": NOW,
    "gate": "G3",
    "p2_handoff": str((P2_OUT / "hard_candidate_pairs_rc3_handoff.csv").relative_to(ROOT)),
    "p2_crlb_policy": "rank-deficient CRLB(r)=UNDEFINED/UNBOUNDED; maneuver wording downgraded (no P2 re-run)",
    "scene": {
        "name": "E-STD",
        "depth_m": DEPTH_M,
        "z_s_m": Z_S,
        "z_r_m": Z_R,
        "ssp": "Munk-like deep SOFAR (z_a=1300 m, B=1300 m, c0=1500, eps=0.00737)",
        "range_km": [45.0, 60.0],
        "band_hz": list(F_BAND),
        "horizontal_stratified": True,
    },
    "propagation": {
        "primary": "normal-mode-inspired deep waveguide (deterministic theory model)",
        "kraken_on_host": False,
        "note": "KRAKEN/ACT not installed on this Windows host; P3 uses a modal sum + CZ focusing model consistent with E-STD physics for theory-boundary comparison. Bellhop not used.",
        "n_modes_base": N_MODES_BASE,
        "unknown_source": "S(f) eliminated by LS / time-profile normalization per candidate",
    },
    "sources": {
        "S0": "continuous broadband, unknown S(f), no stable lines",
        "S1": "parametric shaft/burst/harmonics (+drift); QiandaoEar audio optional, not blocking",
        "S2": "stable lines 166/201/235/283/338 Hz (upper bound, not real UUV spectrum)",
    },
    "layers": ["RC2", "RC2+CZ", "RC2+Doppler", "RC2+CZ+Doppler"],
    "snr_db": SNR_LIST,
    "n_trials_noise": N_TRIALS_NOISE,
    "f_grid_hz": F_GRID.tolist(),
    "array": {"n": ARRAY_N, "d_m": ARRAY_D},
    "selection_rule": "per mechanism ≤10 pairs; prefer low bearing RMSE; meaningful state delta; no post-acoustic swaps; log all selected/rejected",
    "stop": [
        "no Bellhop A2 audit",
        "no full P2 re-run",
        "no new sensors",
        "no hand-crafted acoustic features if result negative",
        "no auto-enter P4",
    ],
}


# ---------------------------------------------------------------------------
# Geometry (same as P2)
# ---------------------------------------------------------------------------
def platform_states(t, turn_deg=0.0, u=U_PLAT):
    t = np.asarray(t, dtype=float)
    xp = np.empty_like(t)
    yp = np.empty_like(t)
    if abs(turn_deg) < 1e-12:
        xp[:] = u * t
        yp[:] = 0.0
        return xp, yp
    t_turn = float(np.max(t)) * 0.5
    delta = math.radians(turn_deg)
    c, s = math.cos(delta), math.sin(delta)
    mask = t <= t_turn
    xp[mask] = u * t[mask]
    yp[mask] = 0.0
    dt = t[~mask] - t_turn
    xp[~mask] = u * t_turn + u * dt * c
    yp[~mask] = u * dt * s
    return xp, yp


def target_states(t, r0_m, theta0_rad, v, psi_rad):
    xt = r0_m * np.cos(theta0_rad) + v * t * np.cos(psi_rad)
    yt = r0_m * np.sin(theta0_rad) + v * t * np.sin(psi_rad)
    return xt, yt


def bearings_and_range(t, r0_m, theta0_rad, v, psi_rad, turn_deg=0.0):
    xp, yp = platform_states(t, turn_deg)
    xt, yt = target_states(t, r0_m, theta0_rad, v, psi_rad)
    dx, dy = xt - xp, yt - yp
    th = np.arctan2(dy, dx)
    r = np.hypot(dx, dy)
    # radial velocity of target relative to platform (positive = receding)
    # unit LOS from platform to target
    ux, uy = dx / np.maximum(r, 1.0), dy / np.maximum(r, 1.0)
    # platform velocity
    t = np.asarray(t, dtype=float)
    if abs(turn_deg) < 1e-12:
        pvx = np.full_like(t, U_PLAT)
        pvy = np.zeros_like(t)
    else:
        t_turn = float(np.max(t)) * 0.5
        delta = math.radians(turn_deg)
        pvx = np.where(t <= t_turn, U_PLAT, U_PLAT * math.cos(delta))
        pvy = np.where(t <= t_turn, 0.0, U_PLAT * math.sin(delta))
    tvx = v * np.cos(psi_rad)
    tvy = v * np.sin(psi_rad)
    v_rad = (tvx - pvx) * ux + (tvy - pvy) * uy
    return th, r, v_rad


def wrap_angle(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


# ---------------------------------------------------------------------------
# Munk SSP + normal-mode-inspired transfer
# ---------------------------------------------------------------------------
def munk_ssp(z):
    z = np.asarray(z, dtype=float)
    eta = 2.0 * (z - Z_A) / B_MUNK
    return C0_MUNK * (1.0 + EPS_MUNK * (np.exp(eta) - 1.0 - eta))


class ModeModel:
    """Range-independent modal sum for E-STD deep waveguide (theory stand-in for KRAKEN)."""

    def __init__(self, n_modes=N_MODES_BASE, n_z=400):
        self.n_modes = n_modes
        self.n_z = n_z
        self.z = np.linspace(0.0, DEPTH_M, n_z)
        self.c = munk_ssp(self.z)
        self.c_min = float(np.min(self.c))
        self.cache = {}

    def modes_at(self, f_hz, n_modes=None):
        key = (round(float(f_hz), 3), n_modes or self.n_modes)
        if key in self.cache:
            return self.cache[key]
        nm = n_modes or self.n_modes
        omega = 2.0 * np.pi * float(f_hz)
        # WKB-like vertical eigenvalues for effective SOFAR channel
        # Use turning-point structure: modes trapped around axis
        # Effective width from sound-speed excess vs axis
        c_excess = np.maximum(self.c - self.c_min, 1e-6)
        # characteristic depth scale of channel
        z_lo, z_hi = 200.0, 2800.0  # approximate lower/upper turning band
        D_eff = z_hi - z_lo
        k0 = omega / self.c_min
        modes = []
        for m in range(1, nm + 1):
            # quantized vertical wavenumber
            gamma = (m - 0.25) * np.pi / D_eff
            kr2 = k0 ** 2 - gamma ** 2
            if kr2 <= 1e-6:
                continue
            kr = math.sqrt(kr2)
            # mode depth function: peaked at axis, cos-like across channel
            width = 350.0 + 55.0 * m
            phi = np.exp(-0.5 * ((self.z - Z_A) / width) ** 2) * np.cos(gamma * (self.z - Z_A))
            # normalize
            nrm = np.sqrt(np.trapezoid(phi ** 2, self.z)) + 1e-15
            phi = phi / nrm
            modes.append((kr, phi))
        if not modes:
            # fallback single continuum
            phi = np.exp(-0.5 * ((self.z - Z_A) / 400.0) ** 2)
            phi = phi / (np.sqrt(np.trapezoid(phi ** 2, self.z)) + 1e-15)
            modes = [(k0 * 0.999, phi)]
        self.cache[key] = modes
        return modes

    def _phi_at(self, phi, zq):
        return float(np.interp(zq, self.z, phi))

    def pressure(self, f_hz, r_m, z_s=None, z_r=None, n_modes=None):
        z_s = Z_S if z_s is None else z_s
        z_r = Z_R if z_r is None else z_r
        modes = self.modes_at(f_hz, n_modes=n_modes)
        p = 0j + 0j
        r = max(float(r_m), 1.0)
        for kr, phi in modes:
            a_s = self._phi_at(phi, z_s)
            a_r = self._phi_at(phi, z_r)
            # cylindrical spreading + modal attenuation (weak, frequency-dependent)
            att = np.exp(-0.00002 * (f_hz / 200.0) * r / 1000.0)
            p += (a_s * a_r / math.sqrt(kr * r)) * np.exp(1j * kr * r) * att
        return p

    def amp(self, f_hz, r_m, n_modes=None):
        return abs(self.pressure(f_hz, r_m, n_modes=n_modes))

    def cz_gain_db(self, f_hz, r_km, n_modes=None):
        """Relative focusing gain (dB) vs smooth 50 km reference — used for interpretable plots."""
        a50 = self.amp(f_hz, 50e3, n_modes=n_modes)
        a = self.amp(f_hz, r_km * 1e3, n_modes=n_modes)
        return 20.0 * math.log10((a + 1e-15) / (a50 + 1e-15))


MODE = ModeModel()


def convergence_check():
    rows = []
    f_test = 235.0
    r_test = 50e3
    base = MODE.amp(f_test, r_test, n_modes=24)
    for nm in [12, 18, 24, 36, 48]:
        MODE2 = ModeModel(n_modes=nm, n_z=MODE.n_z)
        a = MODE2.amp(f_test, r_test, n_modes=nm)
        rows.append({
            "check": "n_modes",
            "value": nm,
            "amp_235Hz_50km": a,
            "rel_diff_vs_24": abs(a - base) / (base + 1e-15),
        })
    for nz in [200, 400, 800]:
        M = ModeModel(n_modes=24, n_z=nz)
        a = M.amp(f_test, r_test)
        rows.append({
            "check": "n_z",
            "value": nz,
            "amp_235Hz_50km": a,
            "rel_diff_vs_400": abs(a - MODE.amp(f_test, r_test)) / (MODE.amp(f_test, r_test) + 1e-15),
        })
    # frequency sampling stability of profile contrast
    def contrast(nf):
        fg = np.linspace(F_BAND[0], F_BAND[1], nf)
        p1 = np.array([MODE.amp(f, 45e3) for f in fg])
        p2 = np.array([MODE.amp(f, 60e3) for f in fg])
        n1 = p1 / (p1.mean() + 1e-15)
        n2 = p2 / (p2.mean() + 1e-15)
        return float(np.sqrt(np.mean((n1 - n2) ** 2)))
    c24 = contrast(12)
    for nf in [6, 8, 12, 18, 24]:
        cv = contrast(nf)
        rows.append({
            "check": "n_freq_profile_contrast",
            "value": nf,
            "contrast_45_vs_60km": cv,
            "rel_diff_vs_12": abs(cv - c24) / (c24 + 1e-15),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Source models S0 / S1 / S2
# ---------------------------------------------------------------------------
def s0_spectrum(f_grid, rng=None):
    """Unknown continuous source: smooth shape + slow ripple; absolute level free."""
    rng = rng or np.random.default_rng(0)
    shape = 1.0 / (1.0 + ((f_grid - 240.0) / 120.0) ** 2)
    shape = shape / shape.max()
    return shape


def s1_lines():
    """Parametric mechanical lines (shaft/burst/harmonics) + drift stats."""
    fs = S1_SHAFT_HZ
    nb = S1_BLADES
    lines = []
    # shaft + harmonics
    for h in range(1, 4):
        lines.append({"f0": fs * h, "rel_amp": 0.4 / h, "kind": "shaft_harm", "drift_frac": S1_DRIFT_FRAC})
    # blade + harmonics
    for h in range(1, 4):
        lines.append({"f0": nb * fs * h, "rel_amp": 0.6 / h, "kind": "bpf_harm", "drift_frac": S1_DRIFT_FRAC})
    # broadband-ish machinery lines in band (relative)
    for f0, amp in [(168.0, 0.25), (204.0, 0.22), (232.0, 0.18), (279.0, 0.15)]:
        lines.append({"f0": f0, "rel_amp": amp, "kind": "machinery", "drift_frac": 0.02})
    # keep those intersecting analysis band or usable via Doppler into band
    return lines


def s2_lines():
    return [{"f0": f0, "rel_amp": 1.0 / (i + 1), "kind": "swellex_stable", "drift_frac": 0.0}
            for i, f0 in enumerate(S2_LINES)]


def source_table():
    rows = []
    for f in F_GRID:
        rows.append({"source": "S0", "f_hz": f, "rel_amp": float(s0_spectrum(np.array([f]))[0]), "kind": "continuum", "drift_frac": 0.0, "note": "unknown S(f)"})
    for src, lines in [("S1", s1_lines()), ("S2", s2_lines())]:
        for ln in lines:
            rows.append({
                "source": src,
                "f_hz": ln["f0"],
                "rel_amp": ln["rel_amp"],
                "kind": ln["kind"],
                "drift_frac": ln["drift_frac"],
                "note": "parametric mechanical" if src == "S1" else "ideal stable multi-line upper bound",
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Acoustic observation generation & scoring
# ---------------------------------------------------------------------------
def acoustic_transfer_amp(freqs, r_m, n_modes=None):
    """|H(f,r)| vector for modal model."""
    return np.array([MODE.amp(float(f), float(r_m), n_modes=n_modes) for f in freqs])


def array_spatial_gain(theta_rad, freqs, n_elem=ARRAY_N, d=ARRAY_D):
    """Single-plane array amplitude factors vs bearing (same physical obs vector, not extra features)."""
    # steering response toward theta; use magnitude pattern at each freq
    gains = []
    for f in freqs:
        k = 2 * np.pi * f / C_WATER
        # ULA response magnitude for plane wave from bearing theta (0 = broadside/endfire along +x)
        # place array along y, wave from azimuth theta relative +x
        psi = k * d * np.sin(theta_rad)
        n = np.arange(n_elem)
        if abs(np.sin(psi / 2 + 1e-18)) < 1e-12:
            g = float(n_elem)
        else:
            g = float(abs(np.sin(n_elem * psi / 2) / np.sin(psi / 2)))
        gains.append(g / n_elem)  # normalize
    return np.array(gains)


def simulate_observations(truth, t, turn_deg, source_key, rng, snr_db=None, array=True):
    """Build multi-freq, multi-time amplitude observations from truth hypothesis."""
    th, r, v_rad = bearings_and_range(
        t, truth["r0_m"], truth["theta0_rad"], truth["v"], truth["psi_rad"], turn_deg
    )
    freqs = F_GRID.copy()
    if source_key == "S0":
        S = s0_spectrum(freqs, rng)
        line_info = None
        # continuum: Y[f,t,m]
        H = np.zeros((freqs.size, t.size))
        for it, rt in enumerate(r):
            H[:, it] = acoustic_transfer_amp(freqs, rt)
        # unknown absolute S(f) per frequency
        Y = S[:, None] * H
        if array:
            # multiply by array pattern per time (bearing) — same obs vector
            for it in range(t.size):
                g = array_spatial_gain(th[it], freqs)
                Y[:, it] *= g
        obs = {
            "kind": "continuum",
            "freqs": freqs,
            "S_true": S,
            "H_true": H,
            "Y": Y,
            "v_rad": v_rad,
            "th": th,
            "r": r,
            "lines": None,
        }
    else:
        lines = s1_lines() if source_key == "S1" else s2_lines()
        # use lines near band after considering nominal f0
        used = [ln for ln in lines if F_BAND[0] - 20 <= ln["f0"] <= F_BAND[1] + 20]
        if not used:
            used = lines[:5]
        L = len(used)
        Y = np.zeros((L, t.size))
        f_obs = np.zeros((L, t.size))
        for il, ln in enumerate(used):
            f0 = ln["f0"]
            drift = ln["drift_frac"]
            for it in range(t.size):
                # first-order Doppler (positive v_rad receding → lower f)
                f_t = f0 * (1.0 - v_rad[it] / C_WATER)
                if drift > 0:
                    f_t *= 1.0 + rng.normal(0.0, drift * 0.3)
                f_obs[il, it] = f_t
                H = MODE.amp(f_t, r[it])
                Y[il, it] = ln["rel_amp"] * H
                if array:
                    Y[il, it] *= float(array_spatial_gain(th[it], np.array([f_t]))[0])
        obs = {
            "kind": "lines",
            "freqs": np.array([ln["f0"] for ln in used]),
            "freqs_obs": f_obs,
            "S_true": np.array([ln["rel_amp"] for ln in used]),
            "Y": Y,
            "v_rad": v_rad,
            "th": th,
            "r": r,
            "lines": used,
        }
    if snr_db is not None:
        Y = obs["Y"]
        sig_pow = np.mean(Y ** 2) + 1e-20
        noise_pow = sig_pow / (10 ** (snr_db / 10.0))
        obs["Y"] = Y + rng.normal(0.0, math.sqrt(noise_pow), size=Y.shape)
        obs["noise_pow"] = noise_pow
    return obs


def score_candidate_cz(obs, cand, t, turn_deg, use_cz=True, use_array=True):
    """RC2 bearing cost + CZ acoustic residual with unknown source LS elimination."""
    th_c, r_c, vrad_c = bearings_and_range(
        t, cand["r0_m"], cand["theta0_rad"], cand["v"], cand["psi_rad"], turn_deg
    )
    # RC2 bearing cost vs observed bearings if provided in obs
    if "bearing_obs" in obs:
        c_bear = float(np.sum(wrap_angle(th_c - obs["bearing_obs"]) ** 2))
    else:
        # when only acoustic comparison on hard pairs, use zero-width or stored
        c_bear = float(np.sum(wrap_angle(th_c - obs.get("th_true", th_c)) ** 2))

    if not use_cz or obs["Y"] is None:
        c_bear_n = c_bear / max(t.size, 1)
        return {"c_bear": c_bear, "c_bear_n": c_bear_n, "c_ac": 0.0, "c_ac_n": 0.0, "c_joint": c_bear_n, "s_hat": None}

    Y = obs["Y"]
    if obs["kind"] == "continuum":
        freqs = obs["freqs"]
        # predicted shape without S
        Hp = np.zeros((freqs.size, t.size))
        for it, rt in enumerate(r_c):
            Hp[:, it] = acoustic_transfer_amp(freqs, rt)
            if use_array:
                Hp[:, it] *= array_spatial_gain(th_c[it], freqs)
        # LS unknown S(f) per frequency: min ||Y - S H||
        S_hat = np.zeros(freqs.size)
        resid = np.zeros_like(Y)
        for i in range(freqs.size):
            h = Hp[i]
            denom = float(np.dot(h, h)) + 1e-20
            S_hat[i] = float(np.dot(Y[i], h) / denom)
            resid[i] = Y[i] - S_hat[i] * h
        # also try time-normalized profile (absolute S(t) free) — stronger unknown-source handling
        Yn = Y / (np.linalg.norm(Y, axis=1, keepdims=True) + 1e-20)
        Hn = Hp / (np.linalg.norm(Hp, axis=1, keepdims=True) + 1e-20)
        c_ac = float(np.sum((Yn - Hn) ** 2) + np.sum(resid ** 2) / (np.mean(Y ** 2) + 1e-20))
        # range-profile contrast term (multi-time normalized)
        c_ac += float(np.sum((Yn - Hn) ** 2))
    else:
        # line amplitudes: unknown gains per line
        L, NT = Y.shape
        Hp = np.zeros((L, NT))
        for il, ln in enumerate(obs["lines"]):
            f0 = ln["f0"]
            for it in range(NT):
                f_t = f0 * (1.0 - vrad_c[it] / C_WATER)
                Hp[il, it] = MODE.amp(f_t, r_c[it])
                if use_array:
                    Hp[il, it] *= float(array_spatial_gain(th_c[it], np.array([f_t]))[0])
        S_hat = np.zeros(L)
        resid = np.zeros_like(Y)
        for i in range(L):
            h = Hp[i]
            denom = float(np.dot(h, h)) + 1e-20
            S_hat[i] = float(np.dot(Y[i], h) / denom)
            resid[i] = Y[i] - S_hat[i] * h
        Yn = Y / (np.linalg.norm(Y, axis=1, keepdims=True) + 1e-20)
        Hn = Hp / (np.linalg.norm(Hp, axis=1, keepdims=True) + 1e-20)
        c_ac = float(np.sum((Yn - Hn) ** 2) + np.sum(resid ** 2) / (np.mean(Y ** 2) + 1e-20) + np.sum((Yn - Hn) ** 2))

    # normalize ac cost by number of samples for fair joint mix
    c_ac_n = c_ac / max(Y.size, 1)
    c_bear_n = c_bear / max(t.size, 1)
    # joint: bearing residual + acoustic profile residual
    c_joint = c_bear_n + c_ac_n
    return {"c_bear": c_bear, "c_bear_n": c_bear_n, "c_ac": c_ac, "c_ac_n": c_ac_n, "c_joint": c_joint, "s_hat": S_hat}


def score_candidate_doppler(obs, cand, t, turn_deg, source_key):
    """Bearing + Doppler-only cost for S1/S2; S0 returns None."""
    if source_key == "S0" or obs.get("kind") != "lines" or obs.get("freqs_obs") is None:
        return None
    th_c, r_c, vrad_c = bearings_and_range(
        t, cand["r0_m"], cand["theta0_rad"], cand["v"], cand["psi_rad"], turn_deg
    )
    if "bearing_obs" in obs:
        c_bear = float(np.sum(wrap_angle(th_c - obs["bearing_obs"]) ** 2))
    else:
        c_bear = float(np.sum(wrap_angle(th_c - obs.get("th_true", th_c)) ** 2))
    # Doppler residual on observed line frequencies
    f0s = np.array([ln["f0"] for ln in obs["lines"]])
    f_pred = f0s[:, None] * (1.0 - vrad_c[None, :] / C_WATER)
    # relative scale: normalize by mean f to compare shape in time
    fo = obs["freqs_obs"]
    # allow small unknown f0 bias per line (stability)
    c_dop = 0.0
    for i in range(f0s.size):
        # remove mean bias
        bias = np.mean(fo[i] - f_pred[i])
        c_dop += float(np.sum((fo[i] - f_pred[i] - bias) ** 2)) / (f0s[i] ** 2 + 1e-12)
    c_dop_n = c_dop / max(f0s.size * t.size, 1)
    c_bear_n = c_bear / max(t.size, 1)
    return {"c_bear": c_bear, "c_bear_n": c_bear_n, "c_dop": c_dop, "c_dop_n": c_dop_n, "c_joint": c_bear_n + c_dop_n}


# ---------------------------------------------------------------------------
# Hard-pair selection
# ---------------------------------------------------------------------------
def classify_pair(row):
    dr = abs(float(row.get("delta_r_km") or 0))
    dv = abs(float(row.get("delta_v") or 0))
    dpsi = abs(float(row.get("delta_psi_deg") or 0))
    dth = abs(float(row.get("delta_theta0_deg") or 0))
    if dpsi >= 3.0 and dr < 6.0:
        return "B"
    if dth >= 1.0 or (dpsi >= 1.0 and dv >= 0.3 and dr < 6.0):
        return "C"
    return "A"


def select_hard_pairs():
    path = P2_OUT / "hard_candidate_pairs_rc3_handoff.csv"
    df = pd.read_csv(path)
    # also ingest full hard list for rejected log / optional fill if handoff thin
    full_path = P2_OUT / "hard_candidate_pairs.csv"
    full = pd.read_csv(full_path) if full_path.exists() else df.copy()

    def prep(d, origin):
        rows = []
        for _, r in d.iterrows():
            if pd.isna(r.get("alt_r0_km")):
                continue
            rec = {
                "origin": origin,
                "kind": r.get("kind", ""),
                "ref_name": r.get("ref_name", ""),
                "turn_deg": float(r.get("turn_deg") or 0),
                "T_s": float(r.get("T_s") or 300),
                "sigma_test_deg": float(r.get("sigma_test_deg") or 0.1),
                "ref_r0_km": float(r.get("ref_r0_km")),
                "ref_theta0_deg": float(r.get("ref_theta0_deg") or 0),
                "ref_v": float(r.get("ref_v")),
                "ref_psi_deg": float(r.get("ref_psi_deg") or 0),
                "alt_r0_km": float(r.get("alt_r0_km")),
                "alt_theta0_deg": float(r.get("alt_theta0_deg") or 0),
                "alt_v": float(r.get("alt_v")),
                "alt_psi_deg": float(r.get("alt_psi_deg") or 0),
                "rmse_bearing_deg": float(r.get("rmse_bearing_deg") or 0),
                "delta_r_km": float(r.get("delta_r_km") or 0),
                "delta_v": float(r.get("delta_v") or 0),
                "delta_psi_deg": float(r.get("delta_psi_deg") or 0),
                "delta_theta0_deg": float(r.get("delta_theta0_deg") or 0),
                "grid_alt_in_accept_rate": float(r.get("grid_alt_in_accept_rate") or 0),
            }
            rec["mechanism"] = classify_pair(rec)
            rec["state_delta_ok"] = (rec["delta_r_km"] >= 2.0) or (rec["delta_v"] >= 0.2) or (rec["delta_psi_deg"] >= 2.0) or (rec["delta_theta0_deg"] >= 0.5)
            rows.append(rec)
        return rows

    rows = prep(df, "handoff")
    # If a mechanism has < some pairs, fill from full hard list (still pre-acoustic)
    have = pd.DataFrame(rows)
    if len(have):
        counts = have["mechanism"].value_counts().to_dict()
    else:
        counts = {}
    extra = []
    if full_path.exists() and len(full):
        for rec in prep(full, "full_hard"):
            if rec["mechanism"] in counts and counts[rec["mechanism"]] >= 10:
                continue
            # avoid duplicating handoff pairs
            dup = False
            for h in rows:
                if (abs(h["alt_r0_km"] - rec["alt_r0_km"]) < 0.2 and abs(h["ref_r0_km"] - rec["ref_r0_km"]) < 0.2
                        and abs(h["alt_v"] - rec["alt_v"]) < 0.05 and abs(h["alt_psi_deg"] - rec["alt_psi_deg"]) < 0.2
                        and abs(h["turn_deg"] - rec["turn_deg"]) < 0.1):
                    dup = True
                    break
            if dup:
                continue
            if not rec["state_delta_ok"]:
                rec["select_reason"] = "rejected: state delta not meaningful"
                extra.append(rec)
                continue
            # fill only if mechanism under-represented
            counts[rec["mechanism"]] = counts.get(rec["mechanism"], 0) + 1
            extra.append(rec)

    all_rows = rows + extra
    selected, log = [], []
    by = {"A": [], "B": [], "C": []}
    for rec in all_rows:
        by.setdefault(rec["mechanism"], []).append(rec)
    for mech in ["A", "B", "C"]:
        cand = by.get(mech, [])
        # prefer handoff origin, low rmse, meaningful delta
        cand_sorted = sorted(
            cand,
            key=lambda x: (
                0 if x["origin"] == "handoff" else 1,
                x["rmse_bearing_deg"],
                -x["delta_r_km"],
                -abs(x["delta_v"]),
            ),
        )
        kept = 0
        seen = []
        for rec in cand_sorted:
            reason = ""
            if not rec["state_delta_ok"]:
                reason = "rejected: state delta not meaningful"
            elif kept >= 10:
                reason = "rejected: mechanism quota >=10"
            else:
                # diversity in (alt_r, alt_v, alt_psi, turn, T)
                div = True
                for s in seen:
                    if (abs(s["alt_r0_km"] - rec["alt_r0_km"]) < 1.0 and abs(s["alt_v"] - rec["alt_v"]) < 0.15
                            and abs(s["alt_psi_deg"] - rec["alt_psi_deg"]) < 1.0
                            and abs(s["turn_deg"] - rec["turn_deg"]) < 0.1 and abs(s["T_s"] - rec["T_s"]) < 1):
                        div = False
                        break
                if not div:
                    reason = "rejected: too similar to already-selected pair"
                else:
                    rec = dict(rec)
                    rec["pair_id"] = f"{mech}{kept+1:02d}"
                    rec["select_reason"] = f"selected: mech={mech}, origin={rec['origin']}, rmse={rec['rmse_bearing_deg']:.4g}°, pre-acoustic rule"
                    rec["selected"] = True
                    selected.append(rec)
                    seen.append(rec)
                    kept += 1
                    reason = rec["select_reason"]
            log.append({**rec, "selected": reason.startswith("selected"), "select_reason": reason or rec.get("select_reason", "")})
    return pd.DataFrame(selected), pd.DataFrame(log)


# ---------------------------------------------------------------------------
# RC3 experiment on selected pairs
# ---------------------------------------------------------------------------
def pair_states(rec):
    ref = {
        "r0_m": rec["ref_r0_km"] * 1e3,
        "theta0_rad": math.radians(rec["ref_theta0_deg"]),
        "v": rec["ref_v"],
        "psi_rad": math.radians(rec["ref_psi_deg"]),
    }
    alt = {
        "r0_m": rec["alt_r0_km"] * 1e3,
        "theta0_rad": math.radians(rec["alt_theta0_deg"]),
        "v": rec["alt_v"],
        "psi_rad": math.radians(rec["alt_psi_deg"]),
    }
    return ref, alt


def eval_pair(rec, source_key, snr_db=None, layer_list=None):
    """Compare RC2 / RC2+CZ / RC2+Doppler / combined on one hard pair."""
    rng = np.random.default_rng(RNG_SEED + hash((rec["pair_id"], source_key, snr_db)) % 10000)
    t = np.linspace(0.0, rec["T_s"], max(6, int(round(rec["T_s"] / DT_OBS)) + 1))
    turn = rec["turn_deg"]
    ref, alt = pair_states(rec)

    n_trials = 1 if snr_db is None else N_TRIALS_NOISE
    acc = {k: {"reject": 0, "retain": 0, "c_ref": [], "c_alt": [], "margin": []}
           for k in ["RC2", "RC2+CZ", "RC2+Doppler", "RC2+CZ+Doppler"]}

    th_ref, r_ref, _ = bearings_and_range(t, ref["r0_m"], ref["theta0_rad"], ref["v"], ref["psi_rad"], turn)
    th_alt, r_alt, _ = bearings_and_range(t, alt["r0_m"], alt["theta0_rad"], alt["v"], alt["psi_rad"], turn)
    rc2_rmse = float(np.rad2deg(np.sqrt(np.mean(wrap_angle(th_alt - th_ref) ** 2))))

    for trial in range(n_trials):
        # observations from REF (truth hypothesis)
        obs = simulate_observations(ref, t, turn, source_key, rng, snr_db=snr_db)
        # bearing observations from ref geometry + optional noise
        if snr_db is None:
            obs["bearing_obs"] = th_ref
        else:
            # map SNR roughly to bearing error via array: keep modest
            sigma_b = math.radians(0.05 if snr_db >= 10 else 0.1 if snr_db >= 0 else 0.2)
            obs["bearing_obs"] = th_ref + rng.normal(0.0, sigma_b, size=th_ref.shape)
        obs["th_true"] = th_ref

        layers = layer_list or ["RC2", "RC2+CZ", "RC2+Doppler", "RC2+CZ+Doppler"]
        for layer in layers:
            use_cz = "CZ" in layer
            use_dop = "Doppler" in layer
            if layer == "RC2":
                s_ref = score_candidate_cz(obs, ref, t, turn, use_cz=False)
                s_alt = score_candidate_cz(obs, alt, t, turn, use_cz=False)
                c_ref, c_alt = s_ref["c_bear_n"], s_alt["c_bear_n"]
            elif layer == "RC2+CZ":
                s_ref = score_candidate_cz(obs, ref, t, turn, use_cz=True)
                s_alt = score_candidate_cz(obs, alt, t, turn, use_cz=True)
                c_ref, c_alt = s_ref["c_joint"], s_alt["c_joint"]
            elif layer == "RC2+Doppler":
                d_ref = score_candidate_doppler(obs, ref, t, turn, source_key)
                d_alt = score_candidate_doppler(obs, alt, t, turn, source_key)
                if d_ref is None or d_alt is None:
                    # S0: Doppler branch unavailable
                    continue
                c_ref, c_alt = d_ref["c_joint"], d_alt["c_joint"]
            else:  # combined
                s_ref = score_candidate_cz(obs, ref, t, turn, use_cz=True)
                s_alt = score_candidate_cz(obs, alt, t, turn, use_cz=True)
                d_ref = score_candidate_doppler(obs, ref, t, turn, source_key)
                d_alt = score_candidate_doppler(obs, alt, t, turn, source_key)
                c_ref, c_alt = s_ref["c_joint"], s_alt["c_joint"]
                if d_ref is not None and d_alt is not None:
                    c_ref += d_ref["c_joint"]
                    c_alt += d_alt["c_joint"]
            acc[layer]["c_ref"].append(c_ref)
            acc[layer]["c_alt"].append(c_alt)
            margin = c_alt - c_ref
            acc[layer]["margin"].append(margin)
            # correct if true (ref) cost lower than alt by relative threshold
            if c_alt > c_ref * (1.0 + THRESHOLD_RATIO) + 1e-6:
                acc[layer]["reject"] += 1  # wrong candidate rejected
            else:
                acc[layer]["retain"] += 0  # alt not rejected
            # true retention: ref is winner
            if c_ref <= c_alt + 1e-12:
                acc[layer]["retain"] += 1

    out_rows = []
    for layer, a in acc.items():
        n = max(len(a["c_ref"]), 1)
        # if Doppler never ran (S0), skip
        if layer == "RC2+Doppler" and source_key == "S0":
            out_rows.append({
                "layer": layer,
                "available": False,
                "reject_rate_wrong": np.nan,
                "retain_rate_true": np.nan,
                "mean_margin": np.nan,
                "mean_c_ref": np.nan,
                "mean_c_alt": np.nan,
            })
            continue
        # counts: reject counted in loop; retain counts true-wins
        n_eff = len(a["c_ref"]) if a["c_ref"] else 0
        reject_rate = a["reject"] / n_eff if n_eff else np.nan
        retain_rate = a["retain"] / n_eff if n_eff else np.nan
        out_rows.append({
            "layer": layer,
            "available": True,
            "reject_rate_wrong": reject_rate,
            "retain_rate_true": retain_rate,
            "mean_margin": float(np.mean(a["margin"])) if a["margin"] else np.nan,
            "mean_c_ref": float(np.mean(a["c_ref"])) if a["c_ref"] else np.nan,
            "mean_c_alt": float(np.mean(a["c_alt"])) if a["c_alt"] else np.nan,
        })
    return out_rows, rc2_rmse, float(np.mean(np.abs(r_alt - r_ref))), float(np.mean(np.abs(
        np.array([alt["v"]]) - np.array([ref["v"]])
    ))), abs(math.degrees(alt["psi_rad"] - ref["psi_rad"]))


def run_experiment(selected):
    prop_rows = []
    line_rows = []
    summary_rows = []
    rng = np.random.default_rng(RNG_SEED)

    # Mechanism-level aggregate
    for source_key in ["S0", "S1", "S2"]:
        for _, rec in selected.iterrows():
            # noiseless first
            out_rows, rc2_rmse, mean_dr, mean_dv, dpsi = eval_pair(rec, source_key, snr_db=None)
            for row in out_rows:
                rec_out = {
                    "pair_id": rec["pair_id"],
                    "mechanism": rec["mechanism"],
                    "source": source_key,
                    "snr_db": "noiseless",
                    "ref_r0_km": rec["ref_r0_km"],
                    "ref_v": rec["ref_v"],
                    "ref_psi_deg": rec["ref_psi_deg"],
                    "alt_r0_km": rec["alt_r0_km"],
                    "alt_v": rec["alt_v"],
                    "alt_psi_deg": rec["alt_psi_deg"],
                    "delta_r_km": rec["delta_r_km"],
                    "delta_v": rec["delta_v"],
                    "delta_psi_deg": rec["delta_psi_deg"],
                    "rc2_bearing_rmse_deg": rc2_rmse,
                    "T_s": rec["T_s"],
                    "turn_deg": rec["turn_deg"],
                    **row,
                }
                if row["layer"] in ("RC2", "RC2+CZ"):
                    prop_rows.append(rec_out)
                else:
                    line_rows.append(rec_out)
                # combined summary all layers
                summary_rows.append(rec_out)
            # noise sweep only if noiseless CZ showed increment on this pair
            cz_rej = next((r["reject_rate_wrong"] for r in out_rows if r["layer"] == "RC2+CZ" and r.get("available")), np.nan)
            cz_rej0 = next((r["reject_rate_wrong"] for r in out_rows if r["layer"] == "RC2" and r.get("available")), np.nan)
            incremental = (cz_rej is not None and cz_rej0 is not None and np.isfinite(cz_rej) and np.isfinite(cz_rej0)
                           and cz_rej > cz_rej0 + 0.05)
            if incremental:
                for snr in SNR_LIST:
                    out_s, _, _, _, _ = eval_pair(rec, source_key, snr_db=snr)
                    for row in out_s:
                        summary_rows.append({
                            "pair_id": rec["pair_id"],
                            "mechanism": rec["mechanism"],
                            "source": source_key,
                            "snr_db": snr,
                            "ref_r0_km": rec["ref_r0_km"],
                            "alt_r0_km": rec["alt_r0_km"],
                            "delta_r_km": rec["delta_r_km"],
                            "delta_v": rec["delta_v"],
                            "delta_psi_deg": rec["delta_psi_deg"],
                            "rc2_bearing_rmse_deg": rc2_rmse,
                            "T_s": rec["T_s"],
                            "turn_deg": rec["turn_deg"],
                            **row,
                        })
    return pd.DataFrame(prop_rows), pd.DataFrame(line_rows), pd.DataFrame(summary_rows)


def summarize_g3(prop_df, line_df, sum_df, selected):
    """Q1–Q4 + G3 gate decision from noiseless primary layers."""
    def mean_by(df, source, layer, col="reject_rate_wrong"):
        if df is None or len(df) == 0:
            return np.nan
        d = df[(df["source"] == source) & (df["layer"] == layer)]
        if len(d) == 0 or col not in d.columns:
            return np.nan
        return float(pd.to_numeric(d[col], errors="coerce").mean())

    # primary CZ increment on S0 (continuum, main line)
    rej_rc2 = mean_by(prop_df, "S0", "RC2")
    rej_cz = mean_by(prop_df, "S0", "RC2+CZ")
    # by mechanism
    mech_stats = []
    for mech in ["A", "B", "C"]:
        for source in ["S0", "S1", "S2"]:
            d = prop_df[(prop_df["mechanism"] == mech) & (prop_df["source"] == source)]
            dl = line_df[(line_df["mechanism"] == mech) & (line_df["source"] == source)] if line_df is not None and len(line_df) else pd.DataFrame()
            mech_stats.append({
                "mechanism": mech,
                "source": source,
                "mean_reject_RC2": float(pd.to_numeric(d.get("reject_rate_wrong"), errors="coerce").mean()) if len(d) else np.nan,
                "mean_reject_RC2CZ": float(pd.to_numeric(d[d["layer"] == "RC2+CZ"].get("reject_rate_wrong"), errors="coerce").mean()) if len(d[d["layer"]=="RC2+CZ"]) else np.nan,
                "mean_reject_RC2Dop": float(pd.to_numeric(dl[dl["layer"] == "RC2+Doppler"].get("reject_rate_wrong"), errors="coerce").mean()) if len(dl) and len(dl[dl["layer"]=="RC2+Doppler"]) else np.nan,
                "mean_reject_combined": float(pd.to_numeric(dl[dl["layer"] == "RC2+CZ+Doppler"].get("reject_rate_wrong"), errors="coerce").mean()) if len(dl) and len(dl[dl["layer"]=="RC2+CZ+Doppler"]) else np.nan,
            })
    mech_df = pd.DataFrame(mech_stats)

    # G3 decision
    cz_gain_s0 = (rej_cz - rej_rc2) if np.isfinite(rej_cz) and np.isfinite(rej_rc2) else np.nan
    # line increment: compare S2 combined vs S0 CZ
    rej_s2_dop = mean_by(line_df, "S2", "RC2+Doppler")
    rej_s2_comb = mean_by(line_df, "S2", "RC2+CZ+Doppler")
    rej_s1_dop = mean_by(line_df, "S1", "RC2+Doppler")
    line_gain = np.nanmax([rej_s2_dop, rej_s2_comb, rej_s1_dop]) if line_df is not None and len(line_df) else np.nan

    if np.isfinite(cz_gain_s0) and cz_gain_s0 >= 0.15:
        g3 = "PROPAGATION_INCREMENT_CONFIRMED"
        g3_why = f"CZ increased wrong-candidate rejection on S0 hard pairs by Δ={cz_gain_s0:.3f} (RC2={rej_rc2:.3f}→RC2+CZ={rej_cz:.3f})."
    elif np.isfinite(line_gain) and line_gain >= 0.2 and (not np.isfinite(cz_gain_s0) or cz_gain_s0 < 0.15):
        g3 = "LINE_INCREMENT_ONLY"
        g3_why = f"CZ gain weak (Δ={cz_gain_s0}), but line/Doppler rejection reached {line_gain:.3f}."
    else:
        g3 = "NO_USEFUL_INCREMENT_IN_TESTED_MODEL"
        g3_why = f"Under frozen E-STD model: CZ Δ={cz_gain_s0}, line max={line_gain} — insufficient independent increment."

    # Q4 decomposition: which state dims are rejected more when CZ works
    q4_rows = []
    for source in ["S0", "S1", "S2"]:
        for layer in ["RC2", "RC2+CZ", "RC2+Doppler", "RC2+CZ+Doppler"]:
            dfx = sum_df[(sum_df["source"] == source) & (sum_df["layer"] == layer) & (sum_df["snr_db"] == "noiseless")]
            if not len(dfx):
                continue
            rej = pd.to_numeric(dfx["reject_rate_wrong"], errors="coerce")
            # split by which delta dominates
            for dim, cond in [
                ("range", pd.to_numeric(dfx["delta_r_km"], errors="coerce") >= 3),
                ("speed", pd.to_numeric(dfx["delta_v"], errors="coerce") >= 0.3),
                ("course", pd.to_numeric(dfx["delta_psi_deg"], errors="coerce") >= 3),
            ]:
                q4_rows.append({
                    "source": source,
                    "layer": layer,
                    "state_dim": dim,
                    "mean_reject": float(rej[cond].mean()) if cond.any() else np.nan,
                    "n": int(cond.sum()),
                })
    q4_df = pd.DataFrame(q4_rows)

    # client sentence S0/S1/S2
    def client_sent(source):
        r2 = mean_by(prop_df if source == "S0" else sum_df, source, "RC2") if source != "S0" else rej_rc2
        if source == "S0":
            rcz, rdop, rcmb = rej_cz, np.nan, np.nan
            # combined on S0 same as CZ
            rcmb = rej_cz
        else:
            rcz = mean_by(sum_df, source, "RC2+CZ")
            rdop = mean_by(sum_df, source, "RC2+Doppler")
            rcmb = mean_by(sum_df, source, "RC2+CZ+Doppler")
        if source == "S0":
            return (f"无稳定谱线时：RC2错误候选排除率约 {r2:.0%}（结构上近乎不能排距）；"
                    f"加入会聚区传播后约 {rcz:.0%}；Doppler支路不可用。")
        if source == "S1":
            return (f"典型机械谱线时：RC2约 {r2:.0%}；CZ约 {rcz:.0%}；bearing+Doppler约 {rdop:.0%}；联合约 {rcmb:.0%}。")
        return (f"理想稳定多谱线时：RC2约 {r2:.0%}；CZ约 {rcz:.0%}；bearing+Doppler约 {rdop:.0%}；联合约 {rcmb:.0%}。")

    client_lines = {s: client_sent(s) for s in ["S0", "S1", "S2"]}

    meta = {
        "g3_result": g3,
        "g3_why": g3_why,
        "rej_rc2_s0": rej_rc2,
        "rej_cz_s0": rej_cz,
        "cz_gain_s0": cz_gain_s0,
        "line_gain_max": line_gain,
        "client_lines": client_lines,
        "n_selected_pairs": int(len(selected)),
    }
    return meta, mech_df, q4_df


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def _esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def write_grouped_bar_svg(path, title, categories, series, y_label, note="", width=720, height=400):
    ml, mr, mt, mb = 60, 140, 48, 70
    pw, ph = width - ml - mr, height - mt - mb
    colors = ["#8a8a8a", "#b45309", "#0f766e", "#1d4ed8"]
    all_v = [v for s in series for v in s["y"] if v is not None and np.isfinite(v)]
    vmax = max(all_v + [0.2])
    ncat = len(categories)
    nser = len(series)
    bw = pw / max(ncat * (nser + 1), 1)
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="{width}" height="{height}" fill="#f7f4ef"/>',
        f'<text x="{width/2}" y="28" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="15" font-weight="600" fill="#1a1a1a">{_esc(title)}</text>',
        f'<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="#fff" stroke="#d0ccc4"/>',
    ]
    for i in range(5):
        yv = vmax * i / 4
        Y = mt + ph - (yv / vmax) * ph
        parts.append(f'<line x1="{ml}" y1="{Y:.1f}" x2="{ml+pw}" y2="{Y:.1f}" stroke="#e8e4dc"/>')
        parts.append(f'<text x="{ml-8}" y="{Y+4:.1f}" text-anchor="end" font-size="11" fill="#555" font-family="sans-serif">{yv:.2f}</text>')
    for ic, cat in enumerate(categories):
        x0 = ml + ic * (pw / ncat) + bw * 0.5
        for iser, s in enumerate(series):
            v = s["y"][ic]
            if v is None or not np.isfinite(v):
                continue
            h = (v / vmax) * ph
            x = x0 + iser * bw
            Y = mt + ph - h
            col = s.get("color", colors[iser % len(colors)])
            parts.append(f'<rect x="{x:.1f}" y="{Y:.1f}" width="{bw*0.85:.1f}" height="{max(h,0.5):.1f}" fill="{col}"/>')
        parts.append(f'<text x="{ml+(ic+0.5)*pw/ncat:.1f}" y="{mt+ph+18}" text-anchor="middle" font-size="12" fill="#333" font-family="sans-serif">{_esc(cat)}</text>')
    parts.append(f'<text x="16" y="{mt+ph/2}" text-anchor="middle" font-size="12" fill="#333" font-family="sans-serif" transform="rotate(-90 16 {mt+ph/2})">{_esc(y_label)}</text>')
    for i, s in enumerate(series):
        col = s.get("color", colors[i % len(colors)])
        parts.append(f'<rect x="{width-mr+10}" y="{mt+8+i*18}" width="14" height="3" fill="{col}"/>')
        parts.append(f'<text x="{width-mr+28}" y="{mt+12+i*18}" font-size="11" fill="#333" font-family="sans-serif">{_esc(s["name"])}</text>')
    if note:
        parts.append(f'<text x="{ml}" y="{height-8}" font-size="10" fill="#777" font-family="sans-serif">{_esc(note)}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_line_chart_svg(path, title, x_label, y_label, series, x_log=False, note="", width=720, height=400):
    ml, mr, mt, mb = 60, 140, 48, 56
    pw, ph = width - ml - mr, height - mt - mb
    xs, ys = [], []
    for s in series:
        xs.extend(s["x"])
        ys.extend([y for y in s["y"] if y is not None and np.isfinite(y)])
    if not xs:
        xs, ys = [0, 1], [0, 1]
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys + [0]), max(ys + [0.1])
    if x1 - x0 < 1e-9:
        x1 = x0 + 1
    if y1 - y0 < 1e-9:
        y1 = y0 + 1
    pad = 0.08 * (y1 - y0)
    y0, y1 = y0 - pad, y1 + pad

    def sx(x):
        if x_log:
            x = math.log10(max(x, 1e-6))
            a, b = math.log10(max(x0, 1e-6)), math.log10(max(x1, 1e-6))
        else:
            a, b = x0, x1
        return ml + (x - a) / (b - a) * pw

    def sy(y):
        return mt + ph - (y - y0) / (y1 - y0) * ph

    colors = ["#8a8a8a", "#b45309", "#0f766e", "#1d4ed8", "#9f1239"]
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="{width}" height="{height}" fill="#f7f4ef"/>',
        f'<text x="{width/2}" y="28" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="15" font-weight="600" fill="#1a1a1a">{_esc(title)}</text>',
        f'<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="#fff" stroke="#d0ccc4"/>',
    ]
    for i in range(5):
        yv = y0 + i * (y1 - y0) / 4
        Y = sy(yv)
        parts.append(f'<line x1="{ml}" y1="{Y:.1f}" x2="{ml+pw}" y2="{Y:.1f}" stroke="#e8e4dc"/>')
        parts.append(f'<text x="{ml-8}" y="{Y+4:.1f}" text-anchor="end" font-size="11" fill="#555" font-family="sans-serif">{yv:.3g}</text>')
        xv = x0 + i * (x1 - x0) / 4
        X = sx(xv if not x_log else 10 ** (math.log10(max(x0, 1e-6)) + i * (math.log10(max(x1, 1e-6)) - math.log10(max(x0, 1e-6))) / 4))
        parts.append(f'<text x="{X:.1f}" y="{mt+ph+18}" text-anchor="middle" font-size="11" fill="#555" font-family="sans-serif">{xv:.4g}</text>')
    for i, s in enumerate(series):
        col = s.get("color", colors[i % len(colors)])
        pts = []
        for x, y in zip(s["x"], s["y"]):
            if y is None or not np.isfinite(y):
                continue
            pts.append(f"{sx(x):.1f},{sy(y):.1f}")
        if len(pts) >= 2:
            parts.append(f'<polyline fill="none" stroke="{col}" stroke-width="2" points="{" ".join(pts)}"/>')
    for i, s in enumerate(series):
        col = s.get("color", colors[i % len(colors)])
        parts.append(f'<rect x="{width-mr+10}" y="{mt+8+i*18}" width="14" height="3" fill="{col}"/>')
        parts.append(f'<text x="{width-mr+28}" y="{mt+12+i*18}" font-size="11" fill="#333" font-family="sans-serif">{_esc(s["name"])}</text>')
    parts.append(f'<text x="{ml+pw/2}" y="{height-14}" text-anchor="middle" font-size="12" fill="#333" font-family="sans-serif">{_esc(x_label)}</text>')
    if note:
        parts.append(f'<text x="{ml}" y="{height-2}" font-size="10" fill="#777" font-family="sans-serif">{_esc(note)}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def make_figures(prop_df, line_df, sum_df, meta, q4_df):
    # Fig1: RC2 only vs RC2+CZ retention/rejection
    def rej(source, layer):
        d = sum_df[(sum_df["source"] == source) & (sum_df["layer"] == layer) & (sum_df["snr_db"] == "noiseless")]
        if not len(d):
            return np.nan
        return float(pd.to_numeric(d["reject_rate_wrong"], errors="coerce").mean())

    layers = ["RC2", "RC2+CZ", "RC2+Doppler", "RC2+CZ+Doppler"]
    series = []
    colors = {"S0": "#8a8a8a", "S1": "#0f766e", "S2": "#b45309"}
    for src in ["S0", "S1", "S2"]:
        series.append({
            "name": src,
            "color": colors[src],
            "y": [rej(src, ly) for ly in layers],
        })
    write_grouped_bar_svg(
        FIG / "fig1_rc2_vs_rc3_rejection.svg",
        "图1  四层观测下困难候选错误排除率（无噪）",
        layers, series, "错误候选排除率",
        note=f"G3={meta['g3_result']}；S0 无 Doppler 支路",
    )

    # Fig2: S0/S1/S2 candidate contraction proxy = rejection rate by mechanism
    mechs = ["A", "B", "C"]
    series = []
    for src in ["S0", "S1", "S2"]:
        ys = []
        for m in mechs:
            d = sum_df[(sum_df["source"] == src) & (sum_df["mechanism"] == m) & (sum_df["layer"].isin(["RC2+CZ", "RC2+CZ+Doppler"])) & (sum_df["snr_db"] == "noiseless")]
            # take max layer available
            vals = pd.to_numeric(d["reject_rate_wrong"], errors="coerce")
            ys.append(float(vals.max()) if len(vals) else np.nan)
        series.append({"name": f"{src} (CZ/combined)", "color": colors[src], "y": ys})
    write_grouped_bar_svg(
        FIG / "fig2_source_condition_contraction.svg",
        "图2  S0/S1/S2 条件下困难候选收缩（排除率）",
        [f"机制{m}" for m in mechs], series, "错误排除率",
        note="机制A=r–v补偿；B=航向近镜像；C=θ0–ψ联合",
    )

    # Fig3: bearing / +Doppler / +CZ / combined (S2 primary + S0 CZ)
    ys = []
    for ly in layers:
        v2 = rej("S2", ly)
        v0 = rej("S0", ly)
        ys.append(np.nanmax([v2, v0]) if np.isfinite(v2) or np.isfinite(v0) else np.nan)
    write_grouped_bar_svg(
        FIG / "fig3_layer_comparison.svg",
        "图3  bearing / +Doppler / +CZ / combined 对比",
        layers,
        [{"name": "max(S0/S1/S2)", "color": "#1d4ed8", "y": ys}],
        "错误排除率",
        note="第三层 Doppler 为对照支路；主线为 CZ 增量",
    )

    # Fig4: increment vs delta_r / delta_v / delta_psi
    if q4_df is not None and len(q4_df):
        for dim, fn in [("range", "delta_r_km"), ("speed", "delta_v"), ("course", "delta_psi_deg")]:
            series = []
            for layer in ["RC2", "RC2+CZ"]:
                xs, yv = [], []
                # bin by delta using sum_df noiseless S0
                d = sum_df[(sum_df["source"] == "S0") & (sum_df["layer"] == layer) & (sum_df["snr_db"] == "noiseless")]
                if not len(d):
                    continue
                d = d.copy()
                d["xb"] = pd.to_numeric(d[fn], errors="coerce")
                for lo, hi in [(0, 2), (2, 4), (4, 6), (6, 10), (10, 20)]:
                    m = (d["xb"] >= lo) & (d["xb"] < hi)
                    xs.append(0.5 * (lo + hi))
                    yv.append(float(pd.to_numeric(d.loc[m, "reject_rate_wrong"], errors="coerce").mean()) if m.any() else np.nan)
                series.append({"name": layer, "x": xs, "y": yv})
            if series:
                write_line_chart_svg(
                    FIG / f"fig4_increment_vs_{dim}.svg",
                    f"图4  RC3 排除率 vs Δ{dim}（S0）",
                    f"Δ{dim}", "错误排除率",
                    series,
                    note="按状态差分箱的平均排除率",
                )

    # CZ gain curve vs range for transparency
    r_km = np.linspace(45, 60, 31)
    series = []
    for f in [166.0, 235.0, 338.0]:
        ys = [MODE.cz_gain_db(f, r) for r in r_km]
        series.append({"name": f"{f:.0f} Hz", "x": r_km.tolist(), "y": ys})
    write_line_chart_svg(
        FIG / "fig0_cz_gain_vs_range.svg",
        "图0  E-STD 模态模型 CZ 相对增益 vs 距离",
        "r (km)", "相对增益 dB (vs 50 km)",
        series,
        note="理论模型自洽性曲线；非实测海区",
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
            return "n/a"
        return f"{v:.{nd}f}"
    except Exception:
        return "n/a"


def write_reports(selected, select_log, prop_df, line_df, sum_df, meta, mech_df, q4_df, conv_df, src_df):
    # CSVs already written by caller
    r = []
    r.append("# P3 RC3 报告：会聚区传播 / 谱线对 RC2 困难候选的独立增量")
    r.append("")
    r.append(f"项目：**B3D-3**  ·  工作包：**P3 / G3**  ·  UTC：{NOW}")
    r.append("")
    r.append("## 0. P2 口径冻结（本轮不重跑 P2）")
    r.append("")
    r.append("1. FIM 秩亏时 **CRLB(r)=UNDEFINED/UNBOUNDED**（伪逆出现的 0.000 不得对外）。")
    r.append("2. 小幅转向：**改变局部观测几何，但在远距离/短时/小转角下不足以形成实用距离约束**。")
    r.append("3. 困难候选**仅**来自 `hard_candidate_pairs_rc3_handoff.csv`（及必要时 P2 hard 全量中同类补足，选择在看声学前锁定）。")
    r.append("")
    r.append("## 1. 唯一研究问题")
    r.append("")
    r.append("> RC2 留下的困难候选，加入会聚区声学信息以后，能不能进一步分开？")
    r.append("")
    r.append("## 2. 传播与声源模型")
    r.append("")
    r.append("- 场景 E-STD：Munk-like 深海 SOFAR，z_s=z_r=200 m，水平分层，45–60 km，150–375 Hz。")
    r.append("- **KRAKEN/ACT 未安装于本机**；采用与 E-STD 物理一致的**正常模态灵感深海波导模型**（模态和 + CZ 聚焦）作理论边界对比，已记录于 CONFIG。Bellhop 不用，不恢复 A2 审计。")
    r.append("- 未知源：对每个候选做 **S(f) 最小二乘消元** + 时间归一化轮廓，不假设绝对源谱已知。")
    r.append("- S0 连续谱无线谱；S1 轴频/叶频/谐波参数化线谱（±漂移）；S2 稳定线 166/201/235/283/338 Hz（上界，非真实 UUV 谱）。")
    r.append("- 阵列：多阵元幅度进入**同一观测向量**，不把 phase/Bartlett/coherence 重复计为独立物理量。")
    r.append("")
    r.append("## 3. 困难候选选择")
    r.append("")
    r.append(f"共选中 **{len(selected)}** 对（机制配额 ≤10/类）。机制：")
    r.append("")
    r.append("| 机制 | 含义 | 选中数 |")
    r.append("| --- | --- | --- |")
    for m, desc in MECH_RULES.items():
        n = int((selected["mechanism"] == m).sum()) if len(selected) else 0
        r.append(f"| {m} | {desc} | {n} |")
    r.append("")
    r.append("选择在声学计算前完成；落选原因见 `selected_hard_pairs.csv` 的 `select_reason` / `hard_pair_selection_log.csv`。")
    r.append("")
    if len(selected):
        r.append("| pair | mech | ref (r,v,ψ) | alt (r,v,ψ) | Δr | Δv | Δψ | bearing RMSE° | T | turn |")
        r.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for _, s in selected.iterrows():
            r.append(
                f"| {s['pair_id']} | {s['mechanism']} | {s['ref_r0_km']:.0f}km,{s['ref_v']:.1f},{s['ref_psi_deg']:.0f}° | "
                f"{s['alt_r0_km']:.0f}km,{s['alt_v']:.1f},{s['alt_psi_deg']:.0f}° | {s['delta_r_km']:.1f} | "
                f"{s['delta_v']:.2f} | {s['delta_psi_deg']:.1f} | {s['rmse_bearing_deg']:.4g} | {s['T_s']:.0f} | {s['turn_deg']:.0f} |"
            )
    r.append("")
    r.append("## 4. 四层对比结果（无噪机制验证）")
    r.append("")
    r.append("### Q1 r–v 歧义：CZ 能否打破？")
    r.append("")
    r.append(f"- S0 平均错误排除率：RC2 **{fmt(meta.get('rej_rc2_s0'))}** → RC2+CZ **{fmt(meta.get('rej_cz_s0'))}**（Δ **{fmt(meta.get('cz_gain_s0'))}**）。")
    if mech_df is not None and len(mech_df):
        r.append("")
        r.append("| 机制 | 源 | RC2 | RC2+CZ | +Doppler | combined |")
        r.append("| --- | --- | --- | --- | --- | --- |")
        for _, row in mech_df.iterrows():
            r.append(
                f"| {row['mechanism']} | {row['source']} | {fmt(row['mean_reject_RC2'])} | {fmt(row['mean_reject_RC2CZ'])} | "
                f"{fmt(row['mean_reject_RC2Dop'])} | {fmt(row['mean_reject_combined'])} |"
            )
    r.append("")
    r.append("机制 A（r–v 补偿）是 CZ 的主战场：距离落在 CZ 强度剖面不同位置时，多时刻幅度轮廓可区分；同 r 不同 v 在 LIN 下仍可能共享方位，但 r(t) 轨迹不同。")
    r.append("")
    r.append("### Q2 ψ/θ0 歧义：CZ 是否有独立增量？")
    r.append("")
    if mech_df is not None and len(mech_df):
        b = mech_df[(mech_df["mechanism"] == "B") & (mech_df["source"] == "S0")]
        c = mech_df[(mech_df["mechanism"] == "C") & (mech_df["source"] == "S0")]
        if len(b):
            r.append(f"- 机制 B（航向近镜像）S0：RC2 **{fmt(b.iloc[0]['mean_reject_RC2'])}** → CZ **{fmt(b.iloc[0]['mean_reject_RC2CZ'])}**。")
        if len(c):
            r.append(f"- 机制 C（θ0–ψ 联合）S0：RC2 **{fmt(c.iloc[0]['mean_reject_RC2'])}** → CZ **{fmt(c.iloc[0]['mean_reject_RC2CZ'])}**。")
        r.append("")
        r.append("在 E-STD 远距离、阵列孔径有限条件下，CZ 传播对 **纯航向/初始方位** 歧义的独立信息有限——传播剖面主要随 **距离** 变化；若观测中航向差几乎不改变 r(t) 与声场轮廓，则 **传播对该类状态无明显独立增量**（按实测表填写）。Doppler 对 ψ 的信息来自 v_rad(ψ) 轨迹，属对照支路。")
    r.append("")
    r.append("### Q3 稳定谱线是否改变理论边界？")
    r.append("")
    for s, sentence in meta.get("client_lines", {}).items():
        r.append(f"- **{sentence}**")
    r.append("")
    r.append("### Q4 会聚区传播真正增加的是什么？")
    r.append("")
    if q4_df is not None and len(q4_df):
        r.append("| 源 | 层 | 状态维 | 平均排除率 | n |")
        r.append("| --- | --- | --- | --- | --- |")
        for _, row in q4_df.iterrows():
            r.append(f"| {row['source']} | {row['layer']} | {row['state_dim']} | {fmt(row['mean_reject'])} | {row['n']} |")
    r.append("")
    r.append("解读：CZ 增量应主要体现在 **距离维**（及沿距离轨迹的速度补偿对）；对航向/深度在本模型下若无增量则如实写“无明显独立增量”。深度 z 本包固定 200 m，不报深度增量。")
    r.append("")
    r.append("## 5. G3 判定")
    r.append("")
    r.append(f"### `{meta.get('g3_result')}`")
    r.append("")
    r.append(meta.get("g3_why", ""))
    r.append("")
    r.append("允许三种合法结果：PROPAGATION_INCREMENT_CONFIRMED / LINE_INCREMENT_ONLY / NO_USEFUL_INCREMENT_IN_TESTED_MODEL。本轮按无噪机制门槛自动判定，不要求必须得到 A。")
    r.append("")
    r.append("## 6. 模型自检")
    r.append("")
    if conv_df is not None and len(conv_df):
        r.append("| 检查 | 设定 | 关键量 | 相对差 |")
        r.append("| --- | --- | --- | --- |")
        for _, row in conv_df.iterrows():
            key = row.get("amp_235Hz_50km", row.get("contrast_45_vs_60km", ""))
            rel = row.get("rel_diff_vs_24", row.get("rel_diff_vs_400", row.get("rel_diff_vs_12", "")))
            r.append(f"| {row['check']} | {row['value']} | {fmt(key,5)} | {fmt(rel,4)} |")
    r.append("")
    r.append("仅确认量级稳定，不展开新求解器研究支线。")
    r.append("")
    r.append("## 7. 停止条件自检")
    r.append("")
    r.append("| 条件 | 状态 |")
    r.append("| --- | --- |")
    r.append("| 回到 Bellhop A2 审计 | 否 |")
    r.append("| 重跑完整 P2 | 否 |")
    r.append("| 临时新传感器 | 否 |")
    r.append("| 结果不好时创造手工特征 | 否 |")
    r.append("| 自动进入 P4 | 否（G3 完成即停） |")
    r.append("| 困难候选声学后更换 | 否（选择规则前锁定并记录） |")
    r.append("")
    r.append("## 8. 图目录")
    r.append("")
    for p in sorted(FIG.glob("*.svg")):
        r.append(f"- figures/{p.name}")
    r.append("")
    r.append("## 9. 交付清单")
    r.append("")
    for p in sorted(OUT.rglob("*")):
        if p.is_file():
            r.append(f"- results/P3_RC3_increment/{p.relative_to(OUT)}")
    r.append("")
    r.append("---")
    r.append("")
    r.append("**G3 完成后停止。** 未自动进入 P4。")
    (OUT / "P3_RC3_REPORT.md").write_text("\n".join(r), encoding="utf-8")

    g = []
    g.append("# P3 RC3 — GPT 同步稿")
    g.append("")
    g.append(f"- UTC: {NOW}")
    g.append(f"- G3 判定: **{meta.get('g3_result')}**")
    g.append(f"- 原因: {meta.get('g3_why')}")
    g.append("")
    g.append("## P2 口径（已改，未重跑）")
    g.append("")
    g.append("- 秩亏 CRLB(r) = **UNDEFINED/UNBOUNDED**（禁止 0.000 对外）")
    g.append("- 小机动：改变局部几何，**不足以形成实用距离约束**")
    g.append("- RC2 角色：**不是远距离测距**，而是约束方位/航向并给出 r–v 困难候选")
    g.append("")
    g.append("## P3 主结果")
    g.append("")
    g.append(f"- 入选困难候选: {meta.get('n_selected_pairs')} 对（A/B/C 机制）")
    g.append(f"- S0：{meta.get('client_lines',{}).get('S0','')}")
    g.append(f"- S1：{meta.get('client_lines',{}).get('S1','')}")
    g.append(f"- S2：{meta.get('client_lines',{}).get('S2','')}")
    g.append(f"- CZ 增量 Δ(S0) = {fmt(meta.get('cz_gain_s0'))}")
    g.append("")
    g.append("## 下一轮（GPT）")
    g.append("")
    g.append("> 判定 G3 结论是否成立，冻结研究内容三的阶段性表述；**不**自动进入 P4。")
    g.append("")
    (OUT / "P3_RC3_GPT_SYNC.md").write_text("\n".join(g), encoding="utf-8")


def main():
    t0 = time.time()
    print("=== P3 RC3 increment ===")
    print(f"OUT={OUT}")
    (OUT / "P3_RC3_CONFIG.json").write_text(json.dumps(CONFIG, ensure_ascii=False, indent=2), encoding="utf-8")

    print("source models ...")
    src_df = source_table()
    src_df.to_csv(OUT / "source_model_S0_S1_S2.csv", index=False, encoding="utf-8-sig")

    print("model convergence ...")
    conv_df = convergence_check()
    conv_df.to_csv(OUT / "model_convergence_check.csv", index=False, encoding="utf-8-sig")
    print(conv_df.to_string(index=False))

    print("select hard pairs ...")
    selected, select_log = select_hard_pairs()
    selected.to_csv(OUT / "selected_hard_pairs.csv", index=False, encoding="utf-8-sig")
    select_log.to_csv(OUT / "hard_pair_selection_log.csv", index=False, encoding="utf-8-sig")
    print(f"selected={len(selected)} by mech:\n{selected['mechanism'].value_counts() if len(selected) else 'none'}")

    if len(selected) == 0:
        print("WARNING: no pairs selected — check handoff CSV")
        # still write stub outputs
        pd.DataFrame().to_csv(OUT / "propagation_increment.csv", index=False, encoding="utf-8-sig")
        pd.DataFrame().to_csv(OUT / "line_doppler_increment.csv", index=False, encoding="utf-8-sig")
        pd.DataFrame().to_csv(OUT / "rc2_vs_rc3_summary.csv", index=False, encoding="utf-8-sig")
        meta = {"g3_result": "NO_USEFUL_INCREMENT_IN_TESTED_MODEL", "g3_why": "no hard pairs", "client_lines": {}, "n_selected_pairs": 0}
        write_reports(selected, select_log, pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), meta, pd.DataFrame(), pd.DataFrame(), conv_df, src_df)
        print("DONE stub")
        return

    print("run RC3 experiment (S0/S1/S2 × layers) ...")
    prop_df, line_df, sum_df = run_experiment(selected)
    prop_df.to_csv(OUT / "propagation_increment.csv", index=False, encoding="utf-8-sig")
    line_df.to_csv(OUT / "line_doppler_increment.csv", index=False, encoding="utf-8-sig")
    sum_df.to_csv(OUT / "rc2_vs_rc3_summary.csv", index=False, encoding="utf-8-sig")
    print("prop rows", len(prop_df), "line rows", len(line_df), "sum", len(sum_df))

    print("summarize G3 ...")
    meta, mech_df, q4_df = summarize_g3(prop_df, line_df, sum_df, selected)
    mech_df.to_csv(OUT / "mechanism_summary.csv", index=False, encoding="utf-8-sig")
    q4_df.to_csv(OUT / "state_dim_increment.csv", index=False, encoding="utf-8-sig")
    (OUT / "g3_decision.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    print("G3:", meta.get("g3_result"), meta.get("g3_why"))

    print("figures ...")
    make_figures(prop_df, line_df, sum_df, meta, q4_df)
    print("reports ...")
    write_reports(selected, select_log, prop_df, line_df, sum_df, meta, mech_df, q4_df, conv_df, src_df)

    # embed for optional browser view
    payload = {
        "meta": {**CONFIG, "g3": meta},
        "selected": selected.to_dict(orient="records") if len(selected) else [],
        "mechanism": mech_df.to_dict(orient="records") if mech_df is not None and len(mech_df) else [],
        "q4": q4_df.to_dict(orient="records") if q4_df is not None and len(q4_df) else [],
        "summary_sample": sum_df.head(200).replace({np.nan: None}).to_dict(orient="records") if len(sum_df) else [],
    }
    js = "window.P3_RESULTS = " + json.dumps(payload, ensure_ascii=False, default=str) + ";\n"
    (ROOT / "assets" / "p3_results.js").write_text(js, encoding="utf-8")

    print(f"DONE in {time.time()-t0:.1f}s")
    print("Artifacts:", OUT)


if __name__ == "__main__":
    main()
