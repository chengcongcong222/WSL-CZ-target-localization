#!/usr/bin/env python3
"""P4-G4 revision round (limited scope).

Fixes only:
  1) Nested candidates: C_RC2+RC3 = C_RC2 ∩ C_RC3  (n must not increase)
  2) SNR-consistent standardized residual J = sum |y-yhat|^2 / sigma_n^2, fixed threshold
  3) retain > width: INVALID_NOT_COVERED when truth not retained
Does NOT re-run depth scan, P2, P3, Bellhop.
Outputs: P4_G4_REVISION.md + 3 CSVs + core figs 2/5/6
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
OUT = ROOT / "results" / "P4_G4_performance_boundary"
FIG = OUT / "figures"
NOW = datetime.now(timezone.utc).isoformat()

# Fixed acceptance: chi-square-like on standardized residual, same for all SNR
# 4 effective state dims (r,θ,v,ψ) × mild factor on acoustic part
DF_BEAR = 4
DF_AC = 4
THR_BEAR = 13.28  # ~chi2 0.99, 4 dof
THR_AC = 13.28
W_AC = 1.0


def ensure_tables():
    if not p4.AMP_TABS:
        for src in ["S0", "S1", "S2"]:
            freqs, _, _ = p4.source_freqs(src)
            p4.AMP_TABS[(src, "E0")] = p4.build_amp_table(p4.MODE_ENV["E0"], freqs)
            p4.AMP_TABS[(src, "E1")] = p4.build_amp_table(p4.MODE_ENV["E1"], freqs)
            p4.AMP_TABS[(src, "E2")] = p4.build_amp_table(p4.MODE_ENV["E2"], freqs)


def time_shape_feats_from_A(A):
    """S(f)-invariant multi-freq temporal shape features.
    For each frequency: remove mean (unknown source level), L2-normalize time series.
    """
    F, M = A.shape
    parts = []
    for j in range(F):
        y = A[j] - np.mean(A[j])
        parts.append(y / (np.linalg.norm(y) + 1e-20))
    feat = np.concatenate(parts)
    return feat / (np.linalg.norm(feat) + 1e-20)


def cand_time_shape_feats(t, g, turn_deg, source, tab):
    freqs, _, _ = p4.source_freqs(source)
    N, M, F = g["r0"].size, t.size, freqs.size
    xp, yp = p4.platform_xy(t, turn_deg)
    feats = np.zeros((N, F * M))
    chunk = max(1, 2500)
    for i0 in range(0, N, chunk):
        i1 = min(N, i0 + chunk)
        r0 = g["r0"][i0:i1, None]
        th0 = g["th0"][i0:i1, None]
        vv = g["v"][i0:i1, None]
        ps = g["psi"][i0:i1, None]
        zz = g["z"][i0:i1]
        xt = r0 * np.cos(th0) + vv * t[None] * np.cos(ps)
        yt = r0 * np.sin(th0) + vv * t[None] * np.sin(ps)
        rr = np.hypot(xt - xp[None], yt - yp[None])
        th = np.arctan2(yt - yp[None], xt - xp[None])
        for k in range(i1 - i0):
            A = p4.lookup_amp(tab, freqs, rr[k], float(zz[k]))
            for it in range(M):
                A[:, it] *= p4.array_gain(th[k, it], freqs)
            feats[i0 + k] = time_shape_feats_from_A(A)
    return feats


def obs_time_shape_feats(t, truth, turn_deg, source, snr_db, rng, tab):
    freqs, _, _ = p4.source_freqs(source)
    F, M = freqs.size, t.size
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
    feat = time_shape_feats_from_A(A)
    sigma = math.radians(truth.get("sigma_deg", 0.1))
    bearing_obs = th + (rng.normal(0, sigma, size=th.shape) if sigma > 0 else 0.0)
    return feat, bearing_obs


def evaluate_nested(g, truth, sigma_deg, T, turn_deg, source, snr_db,
                    env_key="E0", seed=p4.RNG_SEED):
    """RC2 then RC3∩RC2 with standardized residuals; n_RC23 <= n_RC2."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0.0, T, max(6, int(round(T / 10.0)) + 1))
    truth = dict(truth)
    truth["sigma_deg"] = sigma_deg
    pred_b = p4.bearings_grid(t, g, turn_deg)
    cache_key = ("tshape", source, env_key, float(T), float(turn_deg), g["n"])
    if cache_key not in p4._FEAT_CACHE:
        p4._FEAT_CACHE[cache_key] = cand_time_shape_feats(
            t, g, turn_deg, source, p4.AMP_TABS[(source, env_key)]
        )
    feats = p4._FEAT_CACHE[cache_key]

    rows = []
    obs_f, bearing_obs = obs_time_shape_feats(
        t, truth, turn_deg, source, snr_db, rng, p4.AMP_TABS[(source, "E0")]
    )
    sig2 = (math.radians(max(sigma_deg, 1e-6))) ** 2
    dth = p4.wrap(pred_b - bearing_obs[None, :])
    J_bear = np.sum(dth ** 2, axis=1) / sig2
    Jb_min = float(J_bear.min())
    mask_rc2 = J_bear <= Jb_min + THR_BEAR
    s_rc2 = p4.summarize_set(g, mask_rc2, J_bear, truth)
    s_rc2["layer"] = "RC2"
    s_rc2["n_acc"] = int(mask_rc2.sum())
    rows.append(s_rc2)

    # acoustic standardized residual (same confidence rule all SNR)
    # residual energy of unit shape vectors ~ chi2-like; scale by feature noise
    d_ac = feats - obs_f[None, :]
    e_ac = np.sum(d_ac ** 2, axis=1)
    if snr_db is None:
        sig_n2 = 1e-8
    else:
        # temporal-shape noise scales with SNR; fixed functional form
        sig_n2 = (10 ** (-snr_db / 10.0)) / max(feats.shape[1], 1)
        sig_n2 = max(sig_n2, 1e-10)
    J_ac = e_ac / sig_n2
    # nested filter only inside RC2, fixed Delta chi2
    if mask_rc2.any():
        Jmin = float(J_ac[mask_rc2].min())
        mask_rc3 = mask_rc2 & (J_ac <= Jmin + THR_AC)
        if not mask_rc3.any():
            mask_rc3 = np.zeros_like(mask_rc2)
            mask_rc3[int(np.argmin(np.where(mask_rc2, J_ac, np.inf)))] = True
    else:
        mask_rc3 = mask_rc2.copy()

    s_rc3 = p4.summarize_set(g, mask_rc3, J_ac, truth)
    s_rc3["layer"] = "RC2+RC3"
    s_rc3["n_acc"] = int(mask_rc3.sum())
    s_rc3["n_rc2"] = int(mask_rc2.sum())
    s_rc3["stage_reject_rate"] = 1.0 - (s_rc3["n_acc"] / max(s_rc2["n_acc"], 1))
    s_rc3["nested_ok"] = bool(s_rc3["n_acc"] <= s_rc2["n_acc"] and np.all(~mask_rc3 | mask_rc2))
    rows.append(s_rc3)

    return {
        "RC2": s_rc2,
        "RC2+RC3": s_rc3,
        "nested_ok": s_rc3["nested_ok"],
        "n_RC2": int(mask_rc2.sum()),
        "n_RC2_RC3": int(mask_rc3.sum()),
        "J_bear_min": Jb_min,
        "sig_n2": sig_n2,
        "ac_sep": float(np.std(J_ac[mask_rc2])) if mask_rc2.any() else np.nan,
    }


def display_metrics(rec, truth):
    """retain-first display metrics."""
    out = {
        "layer": rec.get("layer"),
        "n_acc": rec.get("n_acc"),
        "retain_r": rec.get("retain_r"),
        "retain_theta": rec.get("retain_theta"),
        "retain_v": rec.get("retain_v"),
        "retain_psi": rec.get("retain_psi"),
        "retain_z": rec.get("retain_z"),
        "nested_ok": rec.get("nested_ok"),
        "n_rc2": rec.get("n_rc2"),
        "stage_reject_rate": rec.get("stage_reject_rate"),
    }
    mapping = [
        ("r", "w_r_km", "contr_r", "retain_r", "km"),
        ("theta", "w_theta_deg", "contr_theta", "retain_theta", "deg"),
        ("z", "w_z_m", "contr_z", "retain_z", "m"),
        ("v", "w_v_mps", "contr_v", "retain_v", "m/s"),
        ("psi", "w_psi_deg", "contr_psi", "retain_psi", "deg"),
    ]
    for state, wcol, ccol, rcol, unit in mapping:
        w = rec.get(wcol)
        c = rec.get(ccol)
        ret = rec.get(rcol)
        if ret is False or ret == 0.0:
            out[f"{state}_status"] = "INVALID_NOT_COVERED"
            out[f"{state}_width"] = np.nan
            out[f"{state}_contraction"] = np.nan
            out[f"{state}_unit"] = unit
        elif w is None or (isinstance(w, float) and not np.isfinite(w)):
            out[f"{state}_status"] = "UNRESOLVED"
            out[f"{state}_width"] = np.nan
            out[f"{state}_contraction"] = np.nan
            out[f"{state}_unit"] = unit
        elif c is not None and np.isfinite(c) and c < 0.2:
            out[f"{state}_status"] = "UNRESOLVED"
            out[f"{state}_width"] = w
            out[f"{state}_contraction"] = c
            out[f"{state}_unit"] = unit
        else:
            out[f"{state}_status"] = "OK"
            out[f"{state}_width"] = w
            out[f"{state}_contraction"] = c
            out[f"{state}_unit"] = unit
    return out


def run_source_snr(g):
    rows = []
    truth = p4.truth_from(p4.B0)
    nested_flags = []
    for src in ["S0", "S1", "S2"]:
        for snr in [20, 10, 0, -10]:
            res = evaluate_nested(
                g, truth, p4.B0["sigma_deg"], p4.B0["T"], p4.B0["turn_deg"], src, snr,
                env_key="E0", seed=p4.RNG_SEED + abs(hash((src, snr))) % 997
            )
            nested_flags.append(res["nested_ok"])
            for layer in ["RC2", "RC2+RC3"]:
                rec = dict(res[layer])
                disp = display_metrics(rec, truth)
                rows.append({
                    "source": src,
                    "snr_db": snr,
                    "layer": layer,
                    "T_s": p4.B0["T"],
                    "sigma_deg": p4.B0["sigma_deg"],
                    "turn_deg": p4.B0["turn_deg"],
                    "n_RC2": res["n_RC2"],
                    "n_RC2_RC3": res["n_RC2_RC3"],
                    "n_acc": rec.get("n_acc"),
                    "nested_ok": res["nested_ok"],
                    "stage_reject_rate": rec.get("stage_reject_rate"),
                    "J_bear_min": res["J_bear_min"],
                    "sig_n2": res["sig_n2"],
                    **{k: disp[k] for k in disp if k.startswith(("r_", "theta_", "v_", "psi_", "z_"))},
                    # raw retain
                    "retain_r": rec.get("retain_r"),
                    "retain_v": rec.get("retain_v"),
                    "retain_theta": rec.get("retain_theta"),
                    "retain_psi": rec.get("retain_psi"),
                    "raw_w_r_km": rec.get("w_r_km"),
                    "raw_contr_r": rec.get("contr_r"),
                    "raw_w_v_mps": rec.get("w_v_mps"),
                    "raw_contr_v": rec.get("contr_v"),
                })
            print(f"  snr {src} {snr}dB  n {res['n_RC2']}→{res['n_RC2_RC3']}  "
                  f"nested={res['nested_ok']}  "
                  f"r: {display_metrics(res['RC2+RC3'], truth)['r_status']} "
                  f"w={display_metrics(res['RC2+RC3'], truth)['r_width']}")
    return pd.DataFrame(rows), all(nested_flags)


def run_integrated(g):
    rows = []
    nested_all = []
    detail = {}
    for key, sc in [("C-L", p4.SCEN_CL), ("C-M", p4.SCEN_CM), ("C-U", p4.SCEN_CU)]:
        truth = p4.truth_from(sc)
        res = evaluate_nested(
            g, truth, sc["sigma_deg"], sc["T"], sc["turn_deg"], sc["source"], sc["snr_db"],
            env_key="E0", seed=p4.RNG_SEED + abs(hash(key)) % 997
        )
        nested_all.append(res["nested_ok"])
        detail[key] = res
        for layer in ["RC2", "RC2+RC3"]:
            rec = dict(res[layer])
            disp = display_metrics(rec, truth)
            rows.append({
                "scenario": key,
                "name": sc["name"],
                "sigma_deg": sc["sigma_deg"],
                "T_s": sc["T"],
                "turn_deg": sc["turn_deg"],
                "source": sc["source"],
                "snr_db": sc["snr_db"],
                "layer": layer,
                "n_RC2": res["n_RC2"],
                "n_RC2_RC3": res["n_RC2_RC3"],
                "n_acc": rec.get("n_acc"),
                "nested_ok": res["nested_ok"],
                "stage_reject_rate": rec.get("stage_reject_rate") if layer == "RC2+RC3" else np.nan,
                "r_status": disp["r_status"],
                "r_width": disp["r_width"],
                "r_contraction": disp["r_contraction"],
                "theta_status": disp["theta_status"],
                "theta_width": disp["theta_width"],
                "theta_contraction": disp["theta_contraction"],
                "z_status": disp["z_status"],
                "z_width": disp["z_width"],
                "z_contraction": disp["z_contraction"],
                "v_status": disp["v_status"],
                "v_width": disp["v_width"],
                "v_contraction": disp["v_contraction"],
                "psi_status": disp["psi_status"],
                "psi_width": disp["psi_width"],
                "psi_contraction": disp["psi_contraction"],
                "retain_r": rec.get("retain_r"),
                "retain_theta": rec.get("retain_theta"),
                "retain_z": rec.get("retain_z"),
                "retain_v": rec.get("retain_v"),
                "retain_psi": rec.get("retain_psi"),
                "raw_w_r_km": rec.get("w_r_km"),
                "raw_contr_r": rec.get("contr_r"),
            })
        print(f"  integ {key}: n {res['n_RC2']}→{res['n_RC2_RC3']} nested={res['nested_ok']}")
        for st in ["r", "theta", "z", "v", "psi"]:
            print(f"    {st}: {display_metrics(res['RC2+RC3'], truth)[st+'_status']} "
                  f"w={display_metrics(res['RC2+RC3'], truth)[st+'_width']}")
    return pd.DataFrame(rows), all(nested_all), detail


def write_stage_csv(detail, g_n):
    rows = []
    for key, res in detail.items():
        rows.append({
            "scenario": key,
            "n_RC1_prior": g_n,
            "n_RC2": res["n_RC2"],
            "n_RC2_RC3": res["n_RC2_RC3"],
            "reject_RC2_to_RC3": 1.0 - res["n_RC2_RC3"] / max(res["n_RC2"], 1),
            "nested_ok": res["nested_ok"],
        })
    return pd.DataFrame(rows)


# --- core figures 2, 5, 6 with retain-first display ---

def _esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def write_heatmap_status(path, title, xlab, ylab, xvals, yvals, Z, status, note="", w=760, h=400):
    """Heatmap where INVALID cells are hatched/dark labeled."""
    ml, mr, mt, mb = 70, 100, 48, 56
    pw, ph = w - ml - mr, h - mt - mb
    Z = np.asarray(Z, float)
    finite = Z[np.isfinite(Z)]
    zmin = float(finite.min()) if finite.size else 0.0
    zmax = float(finite.max()) if finite.size else 1.0
    if zmax - zmin < 1e-12:
        zmax = zmin + 1
    nx, ny = len(xvals), len(yvals)
    cw, ch = pw / max(nx, 1), ph / max(ny, 1)
    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
         f'<rect width="{w}" height="{h}" fill="#f7f4ef"/>',
         f'<text x="{w/2}" y="28" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="15" font-weight="600" fill="#1a1a1a">{_esc(title)}</text>']
    for iy in range(ny):
        for ix in range(nx):
            st = status[iy][ix]
            v = Z[iy, ix]
            x = ml + ix * cw
            y = mt + (ny - 1 - iy) * ch
            if st == "INVALID_NOT_COVERED":
                fill = "#4a0f0f"
                txt = "INVALID"
            elif st == "UNRESOLVED" or not np.isfinite(v):
                fill = "#c4c4c4"
                txt = "UNRES"
            else:
                u = (v - zmin) / (zmax - zmin + 1e-15)
                r = int(240 * (1 - u) + 31 * u)
                g = int(230 * (1 - u) + 111 * u)
                b = int(220 * (1 - u) + 139 * u)
                fill = f"rgb({r},{g},{b})"
                txt = f"{v:.2f}"
            p.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{cw:.1f}" height="{ch:.1f}" fill="{fill}" stroke="#eee"/>')
            p.append(f'<text x="{x+cw/2:.1f}" y="{y+ch/2+3:.1f}" text-anchor="middle" font-size="9" fill="#fff" font-family="sans-serif">{txt}</text>')
    for ix, xv in enumerate(xvals):
        p.append(f'<text x="{ml+(ix+0.5)*cw:.1f}" y="{mt+ph+18}" text-anchor="middle" font-size="11" fill="#555" font-family="sans-serif">{_esc(xv)}</text>')
    for iy, yv in enumerate(yvals):
        p.append(f'<text x="{ml-8}" y="{mt+(ny-1-iy)*ch+ch/2+4:.1f}" text-anchor="end" font-size="11" fill="#555" font-family="sans-serif">{_esc(yv)}</text>')
    p.append(f'<text x="{ml+pw/2}" y="{h-14}" text-anchor="middle" font-size="12" fill="#333" font-family="sans-serif">{_esc(xlab)}</text>')
    p.append(f'<text x="16" y="{mt+ph/2}" text-anchor="middle" font-size="12" fill="#333" font-family="sans-serif" transform="rotate(-90 16 {mt+ph/2})">{_esc(ylab)}</text>')
    p.append(f'<text x="{ml}" y="{h-2}" font-size="10" fill="#777" font-family="sans-serif">{_esc(note)}  深红=INVALID_NOT_COVERED；灰=UNRESOLVED</text></svg>')
    path.write_text("\n".join(p), encoding="utf-8")


def write_group_status_bars(path, title, cats, series, ylab, note="", w=780, h=420):
    """Bars where invalid is a special marker not a high value."""
    ml, mr, mt, mb = 60, 150, 48, 70
    pw, ph = w - ml - mr, h - mt - mb
    cols = ["#8a8a8a", "#b45309", "#0f766e", "#1d4ed8", "#9f1239"]
    # series: [{name, y:[values or nan], status:[OK|INVALID|UNRES]}]
    all_v = []
    for s in series:
        for v, st in zip(s["y"], s.get("status", ["OK"] * len(s["y"]))):
            if st == "OK" and v is not None and np.isfinite(v):
                all_v.append(v)
    vmax = max(all_v + [0.2])
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
            st = s.get("status", ["OK"] * len(s["y"]))[ic]
            col = s.get("color", cols[i % len(cols)])
            x = x0 + i * bw
            if st == "INVALID_NOT_COVERED":
                p.append(f'<rect x="{x:.1f}" y="{mt+ph*0.55:.1f}" width="{bw*0.8:.1f}" height="{ph*0.45:.1f}" fill="#4a0f0f"/>')
                p.append(f'<text x="{x+bw*0.4:.1f}" y="{mt+ph*0.72:.1f}" text-anchor="middle" font-size="8" fill="#fff" font-family="sans-serif">INV</text>')
            elif st == "UNRESOLVED" or v is None or not np.isfinite(v):
                p.append(f'<rect x="{x:.1f}" y="{mt+ph*0.55:.1f}" width="{bw*0.8:.1f}" height="{ph*0.45:.1f}" fill="#c4c4c4"/>')
                p.append(f'<text x="{x+bw*0.4:.1f}" y="{mt+ph*0.72:.1f}" text-anchor="middle" font-size="8" fill="#333" font-family="sans-serif">UNRES</text>')
            else:
                hh = (v / vmax) * ph
                p.append(f'<rect x="{x:.1f}" y="{mt+ph-hh:.1f}" width="{bw*0.8:.1f}" height="{max(hh,0.4):.1f}" fill="{col}"/>')
        p.append(f'<text x="{ml+(ic+0.5)*pw/ncat:.1f}" y="{mt+ph+18}" text-anchor="middle" font-size="12" fill="#333" font-family="sans-serif">{_esc(cat)}</text>')
    p.append(f'<text x="16" y="{mt+ph/2}" text-anchor="middle" font-size="12" fill="#333" font-family="sans-serif" transform="rotate(-90 16 {mt+ph/2})">{_esc(ylab)}</text>')
    for i, s in enumerate(series):
        p.append(f'<rect x="{w-mr+10}" y="{mt+8+i*18}" width="14" height="3" fill="{s.get("color", cols[i%len(cols)])}"/>')
        p.append(f'<text x="{w-mr+28}" y="{mt+12+i*18}" font-size="11" fill="#333" font-family="sans-serif">{_esc(s["name"])}</text>')
    p.append(f'<text x="{ml}" y="{h-6}" font-size="10" fill="#777" font-family="sans-serif">{_esc(note)}</text></svg>')
    path.write_text("\n".join(p), encoding="utf-8")


def write_stage_bars(path, stage_df, w=720, h=400):
    cats = stage_df["scenario"].tolist()
    y = stage_df["n_RC2_RC3"].tolist()  # nested endpoint
    y2 = stage_df["n_RC2"].tolist()
    y1 = stage_df["n_RC1_prior"].tolist()
    # log-like display via status bars using actual counts
    ml, mr, mt, mb = 70, 140, 48, 60
    pw, ph = w - ml - mr, h - mt - mb
    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
         f'<rect width="{w}" height="{h}" fill="#f7f4ef"/>',
         f'<text x="{w/2}" y="28" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="15" font-weight="600" fill="#1a1a1a">核心图6  RC1→RC2→RC3 候选数（强制嵌套 n3≤n2）</text>',
         f'<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="#fff" stroke="#d0ccc4"/>']
    vmax = max(y1 + y2 + y + [1])
    ncat = len(cats)
    bw = pw / (ncat * 4)
    colors = ["#8a8a8a", "#b45309", "#0f766e"]
    series = [("RC1先验", y1), ("RC2", y2), ("RC2∩RC3", y)]
    for ic, cat in enumerate(cats):
        for i, (nm, arr) in enumerate(series):
            v = arr[ic]
            hh = (v / vmax) * ph
            x = ml + ic * (pw / ncat) + i * bw + bw * 0.2
            p.append(f'<rect x="{x:.1f}" y="{mt+ph-hh:.1f}" width="{bw*0.75:.1f}" height="{max(hh,0.5):.1f}" fill="{colors[i]}"/>')
            p.append(f'<text x="{x+bw*0.37:.1f}" y="{mt+ph-hh-4:.1f}" text-anchor="middle" font-size="9" fill="#333" font-family="sans-serif">{int(v)}</text>')
        p.append(f'<text x="{ml+(ic+0.5)*pw/ncat:.1f}" y="{mt+ph+18}" text-anchor="middle" font-size="12" fill="#333" font-family="sans-serif">{_esc(cat)}</text>')
    p.append(f'<text x="16" y="{mt+ph/2}" text-anchor="middle" font-size="12" fill="#333" font-family="sans-serif" transform="rotate(-90 16 {mt+ph/2})">候选数</text>')
    for i, (nm, _) in enumerate(series):
        p.append(f'<rect x="{w-mr+10}" y="{mt+8+i*18}" width="14" height="3" fill="{colors[i]}"/>')
        p.append(f'<text x="{w-mr+28}" y="{mt+12+i*18}" font-size="11" fill="#333" font-family="sans-serif">{nm}</text>')
    note = "；".join([f"{r['scenario']}: {int(r['n_RC1_prior'])}→{int(r['n_RC2'])}→{int(r['n_RC2_RC3'])}" for _, r in stage_df.iterrows()])
    p.append(f'<text x="{ml}" y="{h-6}" font-size="10" fill="#777" font-family="sans-serif">{_esc(note)}</text></svg>')
    path.write_text("\n".join(p), encoding="utf-8")


def main():
    t0 = time.time()
    print("=== P4-G4 REVISION (limited) ===")
    print(f"OUT={OUT}")
    ensure_tables()
    g = p4.make_grid(include_z=True)
    print(f"grid={g['n']}")

    print("[1] source × SNR (nested + standardized J) ...")
    df_ss, nest_ss_flag = run_source_snr(g)
    nest_ss = [bool(nest_ss_flag)]
    df_ss.to_csv(OUT / "boundary_source_snr_revised.csv", index=False, encoding="utf-8-sig")

    print("[2] C-L / C-M / C-U ...")
    df_int, nest_int_flag, detail = run_integrated(g)
    nest_int = [bool(nest_int_flag)]
    df_int.to_csv(OUT / "integrated_CL_CM_CU_revised.csv", index=False, encoding="utf-8-sig")

    print("[3] stage contraction ...")
    df_stage = write_stage_csv(detail, g["n"])
    df_stage.to_csv(OUT / "rc1_rc2_rc3_contraction_revised.csv", index=False, encoding="utf-8-sig")

    print("[4] core figures 2/5/6 ...")
    # fig2: source x SNR → r contraction with status
    sources = ["S0", "S1", "S2"]
    snrs = [20, 10, 0, -10]
    Z = np.zeros((len(sources), len(snrs)))
    ST = []
    d3 = df_ss[df_ss["layer"] == "RC2+RC3"]
    for i, s in enumerate(sources):
        row_st = []
        for j, sn in enumerate(snrs):
            r = d3[(d3["source"] == s) & (d3["snr_db"] == sn)]
            if not len(r):
                Z[i, j] = np.nan
                row_st.append("UNRESOLVED")
            else:
                r = r.iloc[0]
                Z[i, j] = r["r_contraction"] if r["r_contraction"] == r["r_contraction"] else np.nan
                row_st.append(r["r_status"])
        ST.append(row_st)
    write_heatmap_status(
        FIG / "core_fig2_source_snr.svg",
        "核心图2（修正） source×SNR → r 候选收缩（retain优先）",
        "SNR (dB)", "source",
        [str(x) for x in snrs], sources, Z, ST,
        note="嵌套筛选 C_RC2∩C_RC3；标准化残差 J=Σ|y-ŷ|²/σ² 固定阈值",
    )

    # fig5: five-state C-L/M/U with status
    states = [("r", "r"), ("theta", "θ"), ("z", "z"), ("v", "v"), ("psi", "ψ")]
    scen = ["C-L", "C-M", "C-U"]
    d5 = df_int[df_int["layer"] == "RC2+RC3"]
    series = []
    for key, lab in states:
        ys, sts = [], []
        for sc in scen:
            r = d5[d5["scenario"] == sc]
            if not len(r):
                ys.append(np.nan)
                sts.append("UNRESOLVED")
            else:
                r = r.iloc[0]
                ys.append(r[f"{key}_contraction"] if r[f"{key}_contraction"] == r[f"{key}_contraction"] else np.nan)
                sts.append(r[f"{key}_status"])
        series.append({"name": lab, "y": ys, "status": sts})
    write_group_status_bars(
        FIG / "core_fig5_five_state.svg",
        "核心图5（修正） C-L/C-M/C-U 五维状态（INVALID≠高收缩）",
        ["C-L困难", "C-M基准", "C-U理论上界"],
        series,
        "contraction (OK only)",
        note="INVALID_NOT_COVERED=真值未覆盖，不画成收缩=1",
    )

    # fig6 stage
    write_stage_bars(FIG / "core_fig6_stage_contraction.svg", df_stage)

    # G4 revision decision
    all_nested = all(nest_ss) and all(nest_int)
    # SNR monotonicity soft check: for each source, contr_r at -10 should not systematically exceed 20 by a lot among OK cells
    snr_issue = []
    for s in sources:
        vals = {}
        for sn in snrs:
            r = d3[(d3["source"] == s) & (d3["snr_db"] == sn)]
            if len(r) and r.iloc[0]["r_status"] == "OK":
                vals[sn] = float(r.iloc[0]["r_contraction"])
        if 20 in vals and -10 in vals:
            if vals[-10] > vals[20] + 0.3:
                snr_issue.append((s, vals))
    # validity: if any OK cells have nested and non-empty RC2/RC3 and some OK performance metrics
    ok_r = int(((df_int["layer"] == "RC2+RC3") & (df_int["r_status"] == "OK")).sum())
    invalid_n = int(((df_int["layer"] == "RC2+RC3") & (df_int["r_status"] == "INVALID_NOT_COVERED")).sum())
    unres_n = int(((df_int["layer"] == "RC2+RC3") & (df_int["r_status"] == "UNRESOLVED")).sum())

    if all_nested and not snr_issue:
        g4 = "G4_PARTIAL_BOUNDARY_VALIDATED"
        why = (
            f"Nested filtering holds for all revised cells (n_RC2+RC3<=n_RC2). "
            f"Standardized residual J uses fixed threshold across SNR. "
            f"Display is retain-first. Integrated RC2+RC3 statuses: OK r={ok_r}, "
            f"INVALID={invalid_n}, UNRESOLVED={unres_n}."
        )
    elif all_nested and snr_issue:
        g4 = "G4_PARTIAL_BOUNDARY_VALIDATED"
        why = (
            f"Nested filtering OK. Residual scale fixed. Some source×SNR cells remain "
            f"non-monotone in r contraction among OK cells (factors={snr_issue}); "
            f"treated as model/threshold sensitivity notes, not equipment claims. "
            f"Statuses retained-first."
        )
    else:
        g4 = "G4_BOUNDARY_LOGIC_NOT_RESOLVED"
        why = "Nesting condition failed in one or more cells; cannot claim layered contraction."

    # Frozen conclusions (not overturned)
    frozen = [
        "RC2: long-range pure bearing cannot provide practical ranging; collinear r-v ridge is intrinsic.",
        "RC3: under frozen theory propagation + pre-locked hard pairs, CZ provides independent increment on range/range-trajectory candidates.",
        "Depth: z=200/220 m essentially unconstrained in current model; depth remains UNRESOLVED.",
        "Maneuver wording (no re-run): 0–15° small turns change kinematic candidate structure, but no stable monotonic range/course performance gain in this far-range limited-T study.",
        "Environment (no re-run): mismatch can change candidate structure; three points E0/E1/E2 insufficient for a monotone robustness boundary; E1 spikes are not evidence that mild mismatch improves performance.",
    ]

    decision = {
        "g4_revision": g4,
        "why": why,
        "nested_all_ok": all_nested,
        "snr_issues": snr_issue,
        "integrated_status_counts": {"OK_r": ok_r, "INVALID_r": invalid_n, "UNRESOLVED_r": unres_n},
        "frozen_unaffected_conclusions": frozen,
        "outputs": [
            "P4_G4_REVISION.md",
            "boundary_source_snr_revised.csv",
            "integrated_CL_CM_CU_revised.csv",
            "rc1_rc2_rc3_contraction_revised.csv",
            "figures/core_fig2_source_snr.svg",
            "figures/core_fig5_five_state.svg",
            "figures/core_fig6_stage_contraction.svg",
        ],
        "created_utc": NOW,
        "not_rerun": ["depth scan", "P2", "P3", "Bellhop", "maneuver scan", "T×sigma scan"],
        "next": "stop; do not auto-enter P5 in this revision round",
    }
    (OUT / "g4_revision_decision.json").write_text(
        json.dumps(decision, ensure_ascii=False, indent=2, default=str), encoding="utf-8"
    )

    # Write revision markdown
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
    md.append("# P4-G4 修正回合报告")
    md.append("")
    md.append(f"UTC：{NOW}")
    md.append("")
    md.append("**范围**：仅修正 RC2→RC3 候选关系、SNR 统计尺度与展示口径；**不**重跑深度扫描、P2/P3/Bellhop；**不**启动 P5。")
    md.append("")
    md.append(f"## 修正后 G4 状态：`{g4}`")
    md.append("")
    md.append(why)
    md.append("")
    md.append("## 1. 结构修正：强制候选嵌套")
    md.append("")
    md.append("$$\\mathcal{C}_{RC2+RC3}=\\mathcal{C}_{RC2}\\cap\\mathcal{C}_{RC3}$$")
    md.append("")
    md.append("实现：先在标准化方位残差上得到 $\\mathcal{C}_{RC2}$，再**仅在该集合内**用标准化声学残差继续筛选。")
    md.append("")
    md.append(f"- source×SNR 全部格点 nested_ok：**{all(nest_ss)}**")
    md.append(f"- C-L/M/U nested_ok：**{all(nest_int)}**")
    md.append("")
    md.append("### 层间候选数（修正后）")
    md.append("")
    md.append("| 场景 | RC1先验 | RC2 | RC2∩RC3 | 排除率 | nested |")
    md.append("| --- | --- | --- | --- | --- | --- |")
    for _, r in df_stage.iterrows():
        md.append(
            f"| {r['scenario']} | {int(r['n_RC1_prior'])} | {int(r['n_RC2'])} | {int(r['n_RC2_RC3'])} | "
            f"{fnum(r['reject_RC2_to_RC3'])} | {bool(r['nested_ok'])} |"
        )
    md.append("")
    md.append("旧版 C-U 中 445→5122 的反向增加已消除：修正后 RC3 只能减少或保持候选数。")
    md.append("")
    md.append("## 2. SNR 统一标准化残差")
    md.append("")
    md.append("固定同一形式，**不**按 SNR 重调阈值：")
    md.append("")
    md.append("$$J=\\sum\\frac{|y-\\hat{y}|^2}{\\sigma_n^2}$$")
    md.append("")
    md.append(f"- 方位：$\\sigma_n=\\sigma_\\theta$，阈值 $\\Delta\\chi^2\\approx{THR_BEAR}$（固定）")
    md.append(r"- 声学特征：$\sigma_n^2 \propto 10^{-\mathrm{SNR}/10}/F$，阈值固定；只改变测量噪声尺度，不改接受规则")
    md.append("")
    md.append("### source×SNR（RC2+RC3，retain 优先）")
    md.append("")
    md.append("| source | SNR | n RC2→RC3 | r状态 | r宽度 | r收缩 | retain_r |")
    md.append("| --- | --- | --- | --- | --- | --- | --- |")
    for _, r in d3.iterrows():
        md.append(
            f"| {r['source']} | {r['snr_db']} | {int(r['n_RC2'])}→{int(r['n_RC2_RC3'])} | "
            f"{r['r_status']} | {fnum(r['r_width'])} | {fnum(r['r_contraction'])} | {fnum(r['retain_r'],2)} |"
        )
    md.append("")
    if snr_issue:
        md.append("仍有个别 source 在 OK 格点上出现低 SNR 收缩不低于高 SNR 的情况，已记录为模型/阈值敏感性，**不作装备结论**。")
    else:
        md.append("在 OK 格点上，SNR 反常现象已通过统一 $J$ 尺度得到抑制或标注为 UNRESOLVED/INVALID。")
    md.append("")
    md.append("## 3. retain 优先于 width")
    md.append("")
    md.append("若 `truth retention=0`：")
    md.append("")
    md.append("- 状态记 **`INVALID_NOT_COVERED`**")
    md.append("- **不报告**“宽 0”")
    md.append("- **不**把 contraction 画成 1")
    md.append("- 核心图 5 使用深红 INVALID 块，而不是高性能柱")
    md.append("")
    md.append("### 三档综合（RC2+RC3，修正展示）")
    md.append("")
    md.append("| 场景 | r | θ | z | v | ψ | n RC2→RC3 |")
    md.append("| --- | --- | --- | --- | --- | --- | --- |")
    d5i = df_int[df_int["layer"] == "RC2+RC3"]
    for _, r in d5i.iterrows():
        def cell(st):
            stt = r[f"{st}_status"]
            if stt == "INVALID_NOT_COVERED":
                return "**INVALID_NOT_COVERED**"
            if stt == "UNRESOLVED":
                return f"UNRESOLVED (w={fnum(r[f'{st}_width'])})"
            return f"w={fnum(r[f'{st}_width'])}, c={fnum(r[f'{st}_contraction'])}"
        md.append(
            f"| {r['scenario']} | {cell('r')} | {cell('theta')} | {cell('z')} | {cell('v')} | {cell('psi')} | "
            f"{int(r['n_RC2'])}→{int(r['n_RC2_RC3'])} |"
        )
    md.append("")
    md.append("## 4. 仅改口径、不重跑的结论")
    md.append("")
    md.append("### 小机动（冻结表述）")
    md.append("")
    md.append("> 0–15° 小转向改变运动学候选结构，但在本研究的远距离、有限观测时间条件下，**未观察到稳定单调的距离或航向性能提升**；其作用更多体现为候选几何变化，而非直接高精度测距。")
    md.append("")
    md.append("### 环境 E0/E1/E2（非核心性能）")
    md.append("")
    md.append("> 环境扰动可显著改变候选结构；当前三点**不足以**形成单调鲁棒性边界。E1 出现的 r/z 零宽尖峰**不得**解释为“轻度失配反而提升性能”。")
    md.append("")
    md.append("## 5. 不受本次修正影响的项目结论（冻结）")
    md.append("")
    for line in frozen:
        md.append(f"- {line}")
    md.append("")
    md.append("## 6. 本回合输出")
    md.append("")
    for item in decision["outputs"]:
        md.append(f"- `results/P4_G4_performance_boundary/{item}`")
    md.append("")
    md.append("## 7. 停止条件")
    md.append("")
    md.append("- 未启动 P5")
    md.append("- 未重跑深度 / P2 / P3 / Bellhop / T×σθ / 机动扫描")
    md.append(f"- 允许终态：`G4_PARTIAL_BOUNDARY_VALIDATED` 或 `G4_BOUNDARY_LOGIC_NOT_RESOLVED` → 本轮为 **{g4}**")
    md.append("")
    (OUT / "P4_G4_REVISION.md").write_text("\n".join(md), encoding="utf-8")

    print(f"DONE in {time.time()-t0:.1f}s")
    print("G4_REVISION =", g4)
    print("nested all", all_nested)


if __name__ == "__main__":
    main()
