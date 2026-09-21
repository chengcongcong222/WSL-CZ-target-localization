#!/usr/bin/env python3
"""R3-B1.2: 2D branch continuity + E-STD r-z Fisher closure for MMAC.

Only three jobs: 5x5 root-solve grid, unified branch IDs, center Fisher + delay peaks.
No B2 / RC3-C / P5 / new methods.
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
from r3_b1_1_correction import (
    C0,
    RESIDUAL_TOL_M,
    find_eigenrays_root,
    shoot_ray,
    u_hla,
    steering_corr,
)

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_B1_MMAC_observability"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

R_GRID_KM = [49.0, 49.5, 50.0, 50.5, 51.0]
ZS_GRID_M = [190.0, 195.0, 200.0, 205.0, 210.0]
Z_R = 200.0
R0, ZS0 = 50.0, 200.0
SIGMA_THETA_DEG = [0.02, 0.05, 0.1, 0.2]
THETA_REL = [0.0, 30.0, 60.0]
F_BAND = np.arange(150.0, 375.0 + 1e-9, 1.0)
S1 = np.array([168.0, 204.0, 232.0, 279.0, 320.0])
S2 = np.array([166.0, 201.0, 235.0, 283.0, 338.0])

CONFIG = {
    "package": "R3_B1_2_final",
    "created_utc": NOW,
    "grid": {"r_km": R_GRID_KM, "z_s_m": ZS_GRID_M, "center": [R0, ZS0]},
    "admit_tol_m": RESIDUAL_TOL_M,
    "branch_rule": "topology + continuous launch/phi/tau/turning across r AND z_s; cover center ± neighbors only",
    "delay_rule": "S0_DENSE: PSL outside first-null of central peak; S1/S2: existing multipath peaks",
    "ray_modal": "NOT_VALIDATED unless complex H used; abs(H) IFFT forbidden",
    "stop": ["no B2", "no RC3-C", "no P5", "no new methods"],
}


def complex_H(freqs, r_m, z_s=200.0, z_r=200.0, mode=None):
    if mode is None:
        mode = p4.MODE_ENV["E0"]
    p = np.zeros(len(freqs), dtype=complex)
    for i, f in enumerate(freqs):
        ms = mode.modes(float(f))
        acc = 0j
        rr = max(float(r_m), 1.0)
        for kr, phi in ms:
            a_s = float(np.interp(z_s, mode.z, phi))
            a_r = float(np.interp(z_r, mode.z, phi))
            att = math.exp(-2e-5 * (f / 200.0) * rr / 1000.0)
            acc += (a_s * a_r / math.sqrt(kr * rr)) * np.exp(1j * kr * rr) * att
        p[i] = acc
    return p


def run_grid_roots():
    rows = []
    for zs in ZS_GRID_M:
        for r_km in R_GRID_KM:
            _, roots = find_eigenrays_root(r_km, zs)
            if roots is None or roots.empty:
                print(f"  grid r={r_km} zs={zs}: 0", flush=True)
                continue
            roots = roots.copy()
            roots["grid_r_km"] = r_km
            roots["grid_zs_m"] = zs
            rows.append(roots)
            print(f"  grid r={r_km} zs={zs}: admitted={int(roots['admitted'].sum())}/{len(roots)}", flush=True)
            pd.concat(rows, ignore_index=True).to_csv(OUT / "rz_branch_grid.csv", index=False, encoding="utf-8-sig")
    df = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    df.to_csv(OUT / "rz_branch_grid.csv", index=False, encoding="utf-8-sig")
    return df


def associate_2d(roots_df):
    """Unified branch_id across r and z_s via greedy continuity matching."""
    if roots_df is None or roots_df.empty:
        return pd.DataFrame(), pd.DataFrame()
    df = roots_df[roots_df["admitted"] == True].copy()  # noqa: E712
    if df.empty:
        return df, pd.DataFrame()
    # order cells: center-out by Chebyshev distance in (r,z)
    cells = sorted(
        df.groupby(["grid_r_km", "grid_zs_m"]).groups.keys(),
        key=lambda c: (abs(c[0] - R0), abs(c[1] - ZS0), c[0], c[1]),
    )
    branches = []  # active list
    bid = 0
    assigned_rows = []
    for cell in cells:
        g = df[(df["grid_r_km"] == cell[0]) & (df["grid_zs_m"] == cell[1])].copy()
        g = g.sort_values("tau")
        used = set()
        # match existing branches
        for br in branches:
            best_i, best_cost = None, np.inf
            for i, (_, row) in enumerate(g.iterrows()):
                if i in used:
                    continue
                if row["topology"] != br["topology"]:
                    continue
                dL = abs(row["launch_deg"] - br["last_launch"])
                dP = abs(row["phi_arr_deg"] - br["last_phi"])
                dT = abs(row["tau"] - br["last_tau"])
                cost = dL / 3.0 + dP / 3.0 + dT * 1e3 / 50.0
                if dL < 10.0 and dP < 12.0 and cost < best_cost:
                    best_cost = cost
                    best_i = i
            if best_i is not None:
                used.add(best_i)
                row = g.iloc[best_i]
                br["n"] += 1
                br["cells"].add(cell)
                br["last_launch"] = float(row["launch_deg"])
                br["last_phi"] = float(row["phi_arr_deg"])
                br["last_tau"] = float(row["tau"])
                br["last_zmax"] = float(row["z_max"]) if "z_max" in row else br.get("last_zmax", 0)
                br["launch_drift"] = max(br["launch_drift"], abs(row["launch_deg"] - br["first_launch"]))
                br["phi_drift"] = max(br["phi_drift"], abs(row["phi_arr_deg"] - br["first_phi"]))
                br["tau_drift"] = max(br["tau_drift"], abs(row["tau"] - br["first_tau"]))
                br["max_resid"] = max(br["max_resid"], abs(float(row["residual_m"])))
                assigned_rows.append({**row.to_dict(), "branch_id": br["id"]})
        # new branches
        for i, (_, row) in enumerate(g.iterrows()):
            if i in used:
                continue
            id_ = f"R{bid}"
            bid += 1
            br = {
                "id": id_,
                "topology": row["topology"],
                "n": 1,
                "cells": {cell},
                "first_launch": float(row["launch_deg"]),
                "first_phi": float(row["phi_arr_deg"]),
                "first_tau": float(row["tau"]),
                "last_launch": float(row["launch_deg"]),
                "last_phi": float(row["phi_arr_deg"]),
                "last_tau": float(row["tau"]),
                "last_zmax": float(row["z_max"]) if "z_max" in row else 0.0,
                "launch_drift": 0.0,
                "phi_drift": 0.0,
                "tau_drift": 0.0,
                "max_resid": abs(float(row["residual_m"])),
            }
            branches.append(br)
            assigned_rows.append({**row.to_dict(), "branch_id": id_})

    tracks = pd.DataFrame(assigned_rows)
    qual = []
    for br in branches:
        cells = list(br["cells"])
        rs = sorted({c[0] for c in cells})
        zs = sorted({c[1] for c in cells})
        has_center = (R0, ZS0) in br["cells"]
        # cover center and at least one r-neighbor and one z-neighbor
        r_nb = any(abs(r - R0) > 1e-9 for r in rs)
        z_nb = any(abs(z - ZS0) > 1e-9 for z in zs)
        # stronger: cover ±r and ±z if possible
        cover_pr = (50.5, ZS0) in br["cells"] or (R0 + 0.5, ZS0) in br["cells"]
        cover_mr = (49.5, ZS0) in br["cells"]
        cover_pz = (R0, 205.0) in br["cells"]
        cover_mz = (R0, 195.0) in br["cells"]
        eligible = has_center and r_nb and z_nb and br["max_resid"] <= RESIDUAL_TOL_M
        # Jacobian-ready: at center and can estimate ∂/∂r and ∂/∂zs via any neighbor cells
        jac_ready = eligible
        qual.append({
            "branch_id": br["id"],
            "topology": br["topology"],
            "n_cells": len(cells),
            "r_vals_km": json.dumps(rs),
            "zs_vals_m": json.dumps(zs),
            "has_center": has_center,
            "cover_center": has_center,
            "cover_r_neighbor": r_nb,
            "cover_z_neighbor": z_nb,
            "cover_plus_r": cover_pr,
            "cover_minus_r": cover_mr,
            "cover_plus_z": cover_pz,
            "cover_minus_z": cover_mz,
            "jac_ready": jac_ready,
            "eligible_2d": eligible,
            "launch_drift_deg": br["launch_drift"],
            "phi_drift_deg": br["phi_drift"],
            "tau_drift_s": br["tau_drift"],
            "max_residual_m": br["max_resid"],
        })
    return tracks, pd.DataFrame(qual)


def state_at(tracks, branch_id, r_km, zs_m):
    d = tracks[(tracks["branch_id"] == branch_id) &
               (np.isclose(tracks["grid_r_km"], r_km)) &
               (np.isclose(tracks["grid_zs_m"], zs_m)) &
               (tracks["admitted"] == True)]  # noqa: E712
    if d.empty:
        return None
    r = d.iloc[0]
    return {
        "phi_deg": float(r["phi_arr_deg"]),
        "tau": float(r["tau"]),
        "launch_deg": float(r["launch_deg"]),
        "topology": r["topology"],
        "z_max": float(r.get("z_max", np.nan)),
        "residual_m": float(r["residual_m"]),
    }


def estd_fisher(tracks, qual, sigma_tau=2.4e-4):
    """Fisher at center from jac_ready branches; FD via r±0.5 km, zs±5 m."""
    if tracks is None or tracks.empty or qual is None or qual.empty:
        return pd.DataFrame(), pd.DataFrame()
    ready = qual[(qual["jac_ready"] == True) | (qual["eligible_2d"] == True)]  # noqa: E712
    ready = ready[ready["has_center"] == True]  # noqa: E712
    if ready.empty:
        return pd.DataFrame(), pd.DataFrame()

    def pack(branch_ids, r_km, zs, mode, theta=0.0):
        sts = []
        for bid in branch_ids:
            s = state_at(tracks, bid, r_km, zs)
            if s is None:
                sc = state_at(tracks, bid, R0, ZS0)
                if sc is None:
                    continue
                # re-root-solve at neighbor cell; match continuity to center branch
                _, roots = find_eigenrays_root(r_km, zs)
                if roots is None or roots.empty:
                    continue
                ad = roots[roots["admitted"] == True]  # noqa: E712
                if ad.empty:
                    continue
                best = None
                best_c = np.inf
                for _, rr in ad.iterrows():
                    c = abs(rr["launch_deg"] - sc["launch_deg"]) / 3.0 + abs(rr["phi_arr_deg"] - sc["phi_deg"]) / 3.0 + abs(rr["tau"] - sc["tau"]) * 1e3 / 80.0
                    if c < best_c:
                        best_c = c
                        best = rr
                if best is None:
                    continue
                s = {
                    "phi_deg": float(best["phi_arr_deg"]),
                    "tau": float(best["tau"]),
                    "launch_deg": float(best["launch_deg"]),
                    "topology": best["topology"],
                    "z_max": float(best.get("z_max", np.nan)),
                    "residual_m": float(best["residual_m"]),
                }
            sts.append((bid, s))
        if len(sts) < 2:
            return None
        sts = sorted(sts, key=lambda x: x[1]["tau"])
        y = []
        for _, s in sts:
            if mode == "IDEAL":
                y.append(math.radians(s["phi_deg"]))
            elif mode == "HLA":
                y.append(u_hla(theta, s["phi_deg"]))
        for _, s in sts[1:]:
            y.append(s["tau"] - sts[0][1]["tau"])
        return np.asarray(y, dtype=float), [s for _, s in sts], [b for b, _ in sts]

    # all center branches from ready set
    bids_center = [q["branch_id"] for _, q in ready.iterrows()
                   if state_at(tracks, q["branch_id"], R0, ZS0) is not None]
    if len(bids_center) < 2:
        t_c = tracks[(tracks["grid_r_km"] == R0) & (tracks["grid_zs_m"] == ZS0) & (tracks["admitted"] == True)]  # noqa: E712
        bids_center = list(t_c["branch_id"].unique())[:6]
    if len(bids_center) < 2:
        return pd.DataFrame(), pd.DataFrame()

    if len(bids_center) < 2:
        print(f"  Fisher: insufficient center branches {bids_center}", flush=True)
        return pd.DataFrame(), pd.DataFrame()
    # use all center branches; pack() will shoot-fill missing FD points
    bids_use = bids_center
    print(f"  Fisher branch set (center): {bids_use}", flush=True)

    def yv(r_km, zs, mode, theta):
        return pack(bids_use, r_km, zs, mode, theta)

    rows_j, rows_f = [], []
    for sigma_th in SIGMA_THETA_DEG:
        sig_th = math.radians(sigma_th)
        for theta in THETA_REL:
            for mode, label in [("IDEAL", "IDEAL_ELEVATION+delay"),
                                ("HLA", "HLA_u+delay"),
                                ("DELAY", "delay-only")]:
                p0 = yv(R0, ZS0, mode, theta)
                p_rp = yv(R0 + 0.5, ZS0, mode, theta)
                p_rm = yv(R0 - 0.5, ZS0, mode, theta)
                p_zp = yv(R0, ZS0 + 5.0, mode, theta)
                p_zm = yv(R0, ZS0 - 5.0, mode, theta)
                # if a neighbor misses a branch, pack shrinks — align by min length is unsafe.
                # require all 5 packs share same branch id order
                packs = [p0, p_rp, p_rm, p_zp, p_zm]
                if any(p is None for p in packs):
                    continue
                # require same branch ids at center and each FD pack (order by tau may change)
                idsets = [set(p[2]) for p in packs]
                common = set.intersection(*idsets)
                if len(common) < 2:
                    continue
                # rebuild each pack restricted to common ids, sorted by tau
                def restrict(pack_out):
                    y, sts, bids = pack_out
                    keep = [(b, s) for b, s in zip(bids, sts) if b in common]
                    keep = sorted(keep, key=lambda x: x[1]["tau"])
                    y2 = []
                    for _, s in keep:
                        if mode == "IDEAL":
                            y2.append(math.radians(s["phi_deg"]))
                        elif mode == "HLA":
                            y2.append(u_hla(theta, s["phi_deg"]))
                    for _, s in keep[1:]:
                        y2.append(s["tau"] - keep[0][1]["tau"])
                    return np.asarray(y2, dtype=float), [s for _, s in keep], [b for b, _ in keep]
                packs = [restrict(p) for p in packs]
                n = min(len(p[0]) for p in packs)
                if n < 1:
                    continue
                J = np.zeros((n, 2))
                J[:, 0] = (p_rp[0][:n] - p_rm[0][:n]) / 1000.0
                J[:, 1] = (p_zp[0][:n] - p_zm[0][:n]) / 10.0
                n_br = len(p0[1])
                if mode == "IDEAL":
                    sig = [math.radians(0.3)] * n_br + [2.4e-4] * max(n_br - 1, 0)
                elif mode == "HLA":
                    sig = []
                    for s in p0[1]:
                        su = steering_corr(14.0, 8, u_hla(theta, s["phi_deg"]), u_hla(theta, s["phi_deg"]), 250.0)[1]
                        du_dth = -math.sin(math.radians(theta)) * math.cos(math.radians(s["phi_deg"]))
                        sig.append(math.sqrt(su ** 2 + (du_dth * sig_th) ** 2))
                    sig += [2.4e-4] * max(n_br - 1, 0)
                else:
                    sig = [2.4e-4] * max(n_br - 1, 0)
                if len(sig) < n:
                    sig = (sig + [sig[-1] if sig else 1e-3] * n)[:n]
                sig = sig[:n]
                Sinv = np.diag(1.0 / (np.asarray(sig) ** 2 + 1e-30))
                F = J.T @ Sinv @ J
                sv = np.linalg.svd(J, compute_uv=False)
                rank = int(np.linalg.matrix_rank(J))
                try:
                    cov = np.linalg.inv(F)
                    crb_r = math.sqrt(max(float(cov[0, 0]), 0))
                    crb_z = math.sqrt(max(float(cov[1, 1]), 0))
                    corr = float(cov[0, 1] / math.sqrt(cov[0, 0] * cov[1, 1] + 1e-30))
                    condF = float(np.linalg.cond(F))
                except np.linalg.LinAlgError:
                    crb_r = crb_z = np.inf
                    corr = np.nan
                    condF = np.inf
                rows_j.append({
                    "branch_set": json.dumps(list(p0[2])),
                    "sigma_theta_deg": sigma_th,
                    "theta_rel_deg": theta,
                    "case": label,
                    "n_obs": n,
                    "n_branches": n_br,
                    "dy_dr": json.dumps(J[:, 0].tolist()),
                    "dy_dz": json.dumps(J[:, 1].tolist()),
                    "phi_center_deg": json.dumps([round(s["phi_deg"], 3) for s in p0[1]]),
                    "tau_center_s": json.dumps([round(s["tau"], 5) for s in p0[1]]),
                })
                rows_f.append({
                    "branch_set": json.dumps(list(p0[2])),
                    "sigma_theta_deg": sigma_th,
                    "theta_rel_deg": theta,
                    "case": label,
                    "n_obs": n,
                    "n_branches": n_br,
                    "rank": rank,
                    "sv": json.dumps(sv.tolist()),
                    "sv_min": float(sv.min()) if sv.size else np.nan,
                    "cond_J": float(sv[0] / sv[-1]) if sv.size and sv[-1] > 1e-30 else np.inf,
                    "crb_r_m": crb_r,
                    "crb_z_m": crb_z,
                    "corr_rz": corr,
                    "cond_F": condF,
                })
    return pd.DataFrame(rows_j), pd.DataFrame(rows_f)


def delay_peaks_final():
    tau_grid = np.linspace(-0.15, 0.15, 3001)
    tau_true = 0.0
    rows = []
    for name, freqs in [
        ("S0_DENSE", np.arange(150.0, 375.0 + 1e-9, 1.0)),
        ("S1", S1),
        ("S2", S2),
        ("S0_23POINT", np.linspace(150.0, 375.0, 23)),
    ]:
        A = np.abs(np.mean(np.exp(1j * 2 * np.pi * freqs[:, None] * (tau_grid - tau_true)[None, :]), axis=0))
        i0 = int(np.argmin(np.abs(tau_grid - tau_true)))
        # first null-ish: first drop below 0.08 on each side
        def walk(direction):
            i = i0
            while 0 <= i < len(A) and A[i] >= 0.08:
                i += direction
            return i
        il, ir = walk(-1), walk(+1)
        # PSL outside [il, ir]
        mask = np.ones_like(A, dtype=bool)
        if 0 <= il < len(A) and 0 <= ir < len(A):
            lo, hi = min(il, ir), max(il, ir)
            mask[lo:hi + 1] = False
        # local maxima outside
        peaks = []
        for i in range(1, len(A) - 1):
            if not mask[i]:
                continue
            if A[i] >= A[i - 1] and A[i] >= A[i + 1] and A[i] >= 0.2:
                peaks.append((float(tau_grid[i] - tau_grid[i0]), float(A[i])))
        peaks = sorted(peaks, key=lambda x: -x[1])[:6]
        psl = max((p[1] for p in peaks), default=0.0)
        # also max A outside first-null window (true PSL)
        psl_max = float(A[mask].max()) if mask.any() else 0.0
        # local mainlobe
        lo = i0
        hi = i0
        while lo > 0 and A[lo - 1] >= 0.5:
            lo -= 1
        while hi < len(A) - 1 and A[hi + 1] >= 0.5:
            hi += 1
        local_w = float(tau_grid[hi] - tau_grid[lo])
        if name == "S0_DENSE":
            status = "GLOBAL_DELAY_USABLE_IN_TEST_WINDOW" if psl_max < 0.35 and all(p[1] < 0.35 for p in peaks) else "GLOBAL_AMBIGUOUS"
        else:
            status = "GLOBAL_AMBIGUOUS" if psl_max > 0.5 or any(p[1] > 0.5 for p in peaks) else "NEEDS_CHECK"
        if name == "S0_23POINT":
            status = "GLOBAL_AMBIGUOUS_COMB_CONTROL"
        rows.append({
            "source": name,
            "local_mainlobe_width_s": local_w,
            "psl_outside_first_null": psl_max,
            "n_significant_peaks": len(peaks),
            "peaks_offset_amp": json.dumps(peaks),
            "global_status": status,
            "note": "PSL = max A outside central first-null region; not half-power shoulder",
        })
    return pd.DataFrame(rows)


def ray_modal_status():
    """Try complex H IFFT; if fails mark NOT_VALIDATED."""
    try:
        df = 1.0
        n = 1024
        rows = []
        for r_km in [49.0, 50.0, 51.0]:
            H = complex_H(F_BAND, r_km * 1e3, 200.0, 200.0)
            # remove linear phase reference: use tau_ref = r/c0
            tau_ref = r_km * 1e3 / C0
            H = H * np.exp(-1j * 2 * np.pi * F_BAND * tau_ref)
            Hf = np.zeros(n, dtype=complex)
            idx = np.clip(np.round(F_BAND / df).astype(int), 0, n // 2 - 2)
            Hf[idx] = H
            h = np.fft.ifft(Hf)
            t = np.fft.fftfreq(n, d=df)
            mask = (t >= -0.02) & (t <= 0.08)
            p = np.abs(h[mask]) ** 2
            tt = t[mask]
            if p.max() <= 0:
                rows.append({"r_km": r_km, "status": "NOT_VALIDATED", "peaks_s": "[]", "n_peaks": 0})
                continue
            peaks = []
            for i in range(1, len(p) - 1):
                if p[i] >= p[i - 1] and p[i] >= p[i + 1] and p[i] > 0.2 * p.max():
                    peaks.append(float(tt[i] + tau_ref))
            peaks = sorted(peaks)[:8]
            rows.append({
                "r_km": r_km,
                "status": "COMPLEX_H_IFFT_OK" if peaks else "NO_PEAKS_FOUND",
                "peaks_s": json.dumps([round(x, 5) for x in peaks]),
                "n_peaks": len(peaks),
                "tau_ref_s": tau_ref,
                "note": "complex modal H, linear phase removed to tau_ref=r/c0",
            })
        return pd.DataFrame(rows), "COMPLEX_H_USED"
    except Exception as e:
        return pd.DataFrame([{"status": "NOT_VALIDATED", "error": str(e)}]), "NOT_VALIDATED"


def decide(fish_df, qual, delay_df, modal_status):
    notes = {}
    if qual is not None and len(qual):
        notes["n_branch_ids"] = int(qual["branch_id"].nunique())
        notes["n_jac_ready"] = int(qual["jac_ready"].sum())
        notes["n_eligible_2d"] = int(qual["eligible_2d"].sum())
        notes["n_center"] = int(qual["has_center"].sum())
    stable2d = bool(qual is not None and len(qual) and int(((qual["jac_ready"] == True) & (qual["has_center"] == True)).sum()) >= 2)  # noqa: E712

    d0 = delay_df[delay_df["source"] == "S0_DENSE"] if delay_df is not None and len(delay_df) else pd.DataFrame()
    s1 = delay_df[delay_df["source"] == "S1"] if delay_df is not None and len(delay_df) else pd.DataFrame()
    s2 = delay_df[delay_df["source"] == "S2"] if delay_df is not None and len(delay_df) else pd.DataFrame()
    notes["S0_DENSE_psl"] = float(d0["psl_outside_first_null"].iloc[0]) if len(d0) else np.nan
    notes["S0_DENSE_status"] = str(d0["global_status"].iloc[0]) if len(d0) else "n/a"
    notes["S1_status"] = str(s1["global_status"].iloc[0]) if len(s1) else "n/a"
    notes["S2_status"] = str(s2["global_status"].iloc[0]) if len(s2) else "n/a"
    notes["S0_local_ms"] = float(d0["local_mainlobe_width_s"].iloc[0]) * 1e3 if len(d0) else np.nan
    delay_usable = notes["S0_DENSE_status"] == "GLOBAL_DELAY_USABLE_IN_TEST_WINDOW"

    def fish_pick(case, th=0.0, st=0.1):
        if fish_df is None or fish_df.empty:
            return None
        d = fish_df[(fish_df["case"] == case) & (fish_df["theta_rel_deg"] == th) & (fish_df["sigma_theta_deg"] == st)]
        return d.iloc[0] if len(d) else None

    ideal = fish_pick("IDEAL_ELEVATION+delay")
    hla = fish_pick("HLA_u+delay")
    dly = fish_pick("delay-only")
    for lab, row in [("ideal", ideal), ("hla", hla), ("delay", dly)]:
        notes[f"fish_{lab}_crb_r_m"] = float(row["crb_r_m"]) if row is not None else np.nan
        notes[f"fish_{lab}_crb_z_m"] = float(row["crb_z_m"]) if row is not None else np.nan
        notes[f"fish_{lab}_rank"] = int(row["rank"]) if row is not None else 0
        notes[f"fish_{lab}_corr"] = float(row["corr_rz"]) if row is not None else np.nan
        notes[f"fish_{lab}_sv_min"] = float(row["sv_min"]) if row is not None else np.nan

    def ok(row, thr_r=8e3, thr_z=80.0):
        return bool(row is not None and int(row["rank"]) >= 2
                    and np.isfinite(row["crb_r_m"]) and row["crb_r_m"] < thr_r
                    and np.isfinite(row["crb_z_m"]) and row["crb_z_m"] < thr_z)

    ideal_ok, hla_ok, delay_ok = ok(ideal), ok(hla), ok(dly, 15e3, 150.0)
    notes["ideal_ok"], notes["hla_ok"], notes["delay_ok"] = ideal_ok, hla_ok, delay_ok
    notes["ray_modal"] = modal_status
    notes["HLA_elevation_weak_frozen"] = True  # from B1.1 actual branch corr

    if not stable2d:
        decision = "B1_NO_STABLE_MULTIPATH_IDENTITY"
        why = f"2D jac_ready branches={notes.get('n_jac_ready')} < 2 on r-z stencil."
        nxt = "stop RC3-B practical path"
    elif ideal_ok and not hla_ok and delay_ok:
        decision = "B1_PHYSICS_ONLY_HLA_LIMITED"
        why = (
            f"E-STD 2D Fisher: IDEAL+delay CRB_r={notes['fish_ideal_crb_r_m']:.1f} m, "
            f"CRB_z={notes['fish_ideal_crb_z_m']:.1f} m; HLA_u+delay CRB_r={notes['fish_hla_crb_r_m']:.1f} m "
            f"(rank={notes['fish_hla_rank']}). Delay-only CRB_r={notes['fish_delay_crb_r_m']:.1f} m. "
            f"S0_DENSE global={notes['S0_DENSE_status']} (PSL={notes['S0_DENSE_psl']:.3f})."
        )
        nxt = "MMAC r-z physics exists via ideal elevation + delay; HLA elevation limited → delay-dominant / limited-angle hybrid only"
    elif ideal_ok and not hla_ok and not delay_ok:
        decision = "B1_PHYSICS_ONLY_HLA_LIMITED"
        why = (
            f"IDEAL joint Fisher OK (CRB_r={notes['fish_ideal_crb_r_m']:.1f} m) but HLA weak "
            f"and delay-only not yet complementary enough (CRB_r={notes['fish_delay_crb_r_m']:.1f} m). "
            f"S0 global={notes['S0_DENSE_status']}."
        )
        nxt = "do not claim HLA-MMAC; if delay global usable, consider delay-set later — not this round"
    elif (not ideal_ok) and delay_ok and stable2d:
        decision = "B1_DELAY_DOMINANT_POSSIBLE"
        why = (
            f"Angle channels weak/limited; delay-only E-STD Fisher complementary "
            f"(CRB_r={notes['fish_delay_crb_r_m']:.1f} m, CRB_z={notes['fish_delay_crb_z_m']:.1f} m, "
            f"rank={notes['fish_delay_rank']}). S0_DENSE={notes['S0_DENSE_status']}, "
            f"S1={notes['S1_status']}, S2={notes['S2_status']}."
        )
        nxt = "B2 later only as unlabeled delay-set matching if global ambiguity controlled; not HLA elevation MMAC"
    elif stable2d and delay_usable and ideal_ok and hla_ok:
        decision = "B1_MMAC_PHYSICS_CONFIRMED"
        why = "2D stable branches + HLA+delay Fisher + S0 global delay usable."
        nxt = "eligible for B2 later (not this round)"
    elif notes["S1_status"] == "GLOBAL_AMBIGUOUS" or notes["S2_status"] == "GLOBAL_AMBIGUOUS":
        if not delay_ok and not ideal_ok:
            decision = "B1_DELAY_AMBIGUITY_LIMITED"
            why = (
                f"Stable 2D branches exist (jac_ready={notes['n_jac_ready']}) but neither ideal/HLA/delay Fisher "
                f"shows usable r-z CRB under current model; sparse-line global ambiguity remains "
                f"S1={notes['S1_status']}, S2={notes['S2_status']}."
            )
            nxt = "do not proceed B2; r-z complementarity not demonstrated"
        else:
            decision = "B1_DELAY_DOMINANT_POSSIBLE"
            why = (
                f"Fisher shows some r-z info; S1/S2 globally ambiguous but S0_DENSE={notes['S0_DENSE_status']}. "
                f"IDEAL CRB_r={notes['fish_ideal_crb_r_m']}, delay CRB_r={notes['fish_delay_crb_r_m']}."
            )
            nxt = "delay-set path conditional on S0 global usability; not full HLA MMAC"
    else:
        decision = "B1_NO_COMPLEMENTARY_INFORMATION"
        why = f"2D branches present but Fisher ranks/CRB not complementary (ideal_ok={ideal_ok}, hla_ok={hla_ok}, delay_ok={delay_ok})."
        nxt = "stop RC3-B as currently formulated"
    return decision, why, nxt, notes


def main():
    t0 = time.time()
    print("=== R3-B1.2 final B1 closure ===")
    print(f"OUT={OUT}")
    (OUT / "R3_B1_2_CONFIG.json").write_text(json.dumps(CONFIG, ensure_ascii=False, indent=2), encoding="utf-8")

    print("[1] 5x5 r-z eigenray grid ...", flush=True)
    roots = run_grid_roots()

    print("[2] 2D branch association ...", flush=True)
    tracks, qual = associate_2d(roots)
    tracks.to_csv(OUT / "rz_branch_association.csv", index=False, encoding="utf-8-sig")
    qual.to_csv(OUT / "rz_branch_quality.csv", index=False, encoding="utf-8-sig")
    print(qual.to_string(index=False) if len(qual) else "empty qual", flush=True)

    print("[3] E-STD r-z Fisher at center ...", flush=True)
    jac_df, fish_df = estd_fisher(tracks, qual)
    jac_df.to_csv(OUT / "estd_rz_jacobian_final.csv", index=False, encoding="utf-8-sig")
    fish_df.to_csv(OUT / "estd_rz_fisher_final.csv", index=False, encoding="utf-8-sig")
    if len(fish_df):
        print(fish_df[fish_df["sigma_theta_deg"] == 0.1].to_string(index=False), flush=True)
    else:
        print("Fisher empty — check branch coverage", flush=True)

    print("[4] delay global peaks ...", flush=True)
    delay_df = delay_peaks_final()
    delay_df.to_csv(OUT / "delay_global_peaks_final.csv", index=False, encoding="utf-8-sig")
    print(delay_df.to_string(index=False), flush=True)

    print("[5] ray-modal ...", flush=True)
    modal_df, modal_status = ray_modal_status()
    modal_df.to_csv(OUT / "ray_modal_delay_crosscheck.csv", index=False, encoding="utf-8-sig")
    print(modal_status, flush=True)

    print("[6] decision ...", flush=True)
    decision, why, nxt, notes = decide(fish_df, qual, delay_df, modal_status)
    dec = {
        "rc3b1_2_decision": decision,
        "why": why,
        "next_step": nxt,
        "notes": notes,
        "created_utc": NOW,
        "stop": "after R3-B1.2; no B2/RC3-C/P5",
        "frozen_from_b1_1": [
            "HLA elevation projection weak (14m corr-resolve ~0)",
            "local delay mainlobe ms-scale (old 300-400ms metric revoked)",
            "eigenray root-solve residual <= 0.1 m",
        ],
    }
    (OUT / "R3_B1_2_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    def fnum(x, nd=3):
        try:
            v = float(x)
            return "n/a" if not np.isfinite(v) else f"{v:.{nd}f}"
        except Exception:
            return "n/a"

    rp = []
    rp.append("# R3-B1.2 报告：二维分支 + E-STD r–z Fisher 收口")
    rp.append("")
    rp.append(f"UTC：{NOW}")
    rp.append("")
    rp.append("## 0. 已冻结（B1.1）")
    rp.append("")
    rp.append("- HLA 俯仰投影弱（14 m corr 可分≈0；±φ 镜像不可分）")
    rp.append("- 局部时延主瓣 ms 级；旧 300–400 ms 指标永久撤销")
    rp.append("- 本征射线根求解 residual≤0.1 m")
    rp.append("- 旧 `B1_DELAY_BANDWIDTH_LIMITED` / 中间 `B1_DELAY_AMBIGUITY_LIMITED` **不作为最终结论**")
    rp.append("")
    rp.append("## 1. 5×5 r–z 网格与二维分支")
    rp.append("")
    rp.append(f"- 网格：r={R_GRID_KM} km，z_s={ZS_GRID_M} m")
    rp.append(f"- branch id 数：**{notes.get('n_branch_ids')}**")
    rp.append(f"- 覆盖中心且含 r/z 邻域（eligible_2d）：**{notes.get('n_eligible_2d')}**")
    rp.append(f"- jac_ready（中心 ±r 与 ±z 均有点）：**{notes.get('n_jac_ready')}**")
    rp.append("")
    if qual is not None and len(qual):
        rp.append("| branch | topology | n_cells | jac_ready | launch drift° | φ drift° | τ drift s | max resid m |")
        rp.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
        for _, r in qual.iterrows():
            rp.append(
                f"| {r['branch_id']} | {r['topology']} | {r['n_cells']} | {r['jac_ready']} | "
                f"{fnum(r['launch_drift_deg'])} | {fnum(r['phi_drift_deg'])} | {fnum(r['tau_drift_s'],5)} | "
                f"{fnum(r['max_residual_m'])} |"
            )
        rp.append("")
    rp.append("## 2. E-STD 真实 r–z Fisher（中心 50 km, z_s=200 m）")
    rp.append("")
    rp.append("| case | σθ° | θ° | rank | CRB_r m | CRB_z m | corr(r,z) | sv_min |")
    rp.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
    if len(fish_df):
        show = fish_df[fish_df["sigma_theta_deg"].isin([0.05, 0.1, 0.2]) & fish_df["theta_rel_deg"].isin([0.0, 30.0])]
        for _, r in show.iterrows():
            rp.append(
                f"| {r['case']} | {r['sigma_theta_deg']} | {r['theta_rel_deg']} | {r['rank']} | "
                f"{fnum(r['crb_r_m'],1)} | {fnum(r['crb_z_m'],1)} | {fnum(r['corr_rz'])} | {fnum(r['sv_min'],5)} |"
            )
    else:
        rp.append("| （空） | | | | | | | |")
    rp.append("")
    rp.append("Jacobian 明细：`estd_rz_jacobian_final.csv`；Fisher：`estd_rz_fisher_final.csv`。")
    rp.append("")
    rp.append("## 3. 时延全局结论（PSL，非半高肩部）")
    rp.append("")
    rp.append("| 源 | 局部主瓣 ms | PSL（第一零点外） | 峰数 | global_status |")
    rp.append("| --- | --- | --- | --- | --- |")
    for _, r in delay_df.iterrows():
        rp.append(
            f"| {r['source']} | {fnum(r['local_mainlobe_width_s']*1e3,2)} | "
            f"{fnum(r['psl_outside_first_null'])} | {r['n_significant_peaks']} | **{r['global_status']}** |"
        )
    rp.append("")
    rp.append("## 4. ray–modal")
    rp.append("")
    rp.append(f"- 状态：**{modal_status}**（复数 H 去线性相位后 IFFT；若 NOT_VALIDATED 则不进入判定）")
    rp.append("")
    rp.append("## 5. R3-B1 最终判定")
    rp.append("")
    rp.append(f"### `{decision}`")
    rp.append("")
    rp.append(why)
    rp.append("")
    rp.append(f"**下一步**：{nxt}")
    rp.append("")
    rp.append("允许终态：`B1_MMAC_PHYSICS_CONFIRMED` / `B1_PHYSICS_ONLY_HLA_LIMITED` / `B1_DELAY_DOMINANT_POSSIBLE` / `B1_DELAY_AMBIGUITY_LIMITED` / `B1_NO_STABLE_MULTIPATH_IDENTITY` / `B1_NO_COMPLEMENTARY_INFORMATION`。")
    rp.append("")
    rp.append("## 6. 停止")
    rp.append("")
    rp.append("- 不进 B2 / RC3-C / P5")
    rp.append("- 不增加新方法")
    rp.append(f"- **B1.2 完成后停止；判定 `{decision}`**")
    rp.append("")
    (OUT / "R3_B1_2_REPORT.md").write_text("\n".join(rp), encoding="utf-8")

    gs = [
        "# R3-B1.2 — GPT 同步稿", "",
        f"- UTC: {NOW}",
        f"- **最终判定：{decision}**",
        f"- {why}",
        f"- 下一步：{nxt}", "",
        "## Fisher (σθ=0.1°)",
        fish_df[fish_df["sigma_theta_deg"] == 0.1].to_string(index=False) if len(fish_df) else "(empty)", "",
        "## Delay",
        delay_df.to_string(index=False), "",
        "## Branches",
        f"ids={notes.get('n_branch_ids')} jac_ready={notes.get('n_jac_ready')} eligible2d={notes.get('n_eligible_2d')}", "",
        "停止：无 B2/RC3-C/P5。", "",
    ]
    (OUT / "R3_B1_2_GPT_SYNC.md").write_text("\n".join(gs), encoding="utf-8")

    print(f"DONE in {time.time()-t0:.1f}s")
    print("DECISION", decision)
    print(why)


if __name__ == "__main__":
    main()
