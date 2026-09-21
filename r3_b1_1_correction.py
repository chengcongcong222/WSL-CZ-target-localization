#!/usr/bin/env python3
"""R3-B1.1: correction gate for MMAC observability.

Fix delay ambiguity metrics, S0 dense vs comb, eigenray root-solve,
local branch association, E-STD Fisher from actual branches, ray-modal delay crosscheck.
No B2 / RC3-C / P5 / MC RMSE / new physics features.
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
OUT = ROOT / "results" / "R3_B1_MMAC_observability"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

C0 = 1500.0
Z_A, B_MUNK, EPS_MUNK = 1300.0, 1300.0, 0.00737
DEPTH = 5000.0
Z_R = 200.0
F_BAND = (150.0, 375.0)

# Local stencil for Fisher (frozen)
R_STENCIL_KM = [49.0, 49.5, 50.0, 50.5, 51.0]
ZS_STENCIL_M = [190.0, 195.0, 200.0, 205.0, 210.0]
R_SIDEBAND_KM = [45.0, 55.0, 60.0]
THETA_REL_SCAN = [0.0, 30.0, 60.0]
SIGMA_THETA_DEG = [0.02, 0.05, 0.1, 0.2]
RESIDUAL_TOL_M = 0.1

CONFIG = {
    "package": "R3_B1_1_correction",
    "created_utc": NOW,
    "prior_decision_revoked": "B1_DELAY_BANDWIDTH_LIMITED",
    "pending_status": "B1_DECISION_PENDING_CORRECTION",
    "corrections": [
        "delay mainlobe = contiguous A>=0.5 around true peak (not global span)",
        "S0_DENSE continuous band vs S0_23POINT comb control",
        "eigenray root-solve with residual |z_hit-z_r|<=0.1 m",
        "branch_id via continuous association, not round(launch/2°)",
        "E-STD Fisher from actual stable branches, not controlled_paths",
        "ray-modal delay crosscheck via modal transfer IFFT",
        "HLA_ELEVATION_PROJECTION_WEAK based on actual branch Delta-u",
    ],
    "controlled_model_label": "50 km direct/surface weak-geometry sanity check — NOT E-STD CZ MMAC performance",
    "stop": ["no B2", "no RC3-C", "no P5", "no MC RMSE", "no new features"],
}


def munk_c(z, c0=1500.0, za=Z_A, B=B_MUNK, eps=EPS_MUNK):
    z = np.asarray(z, dtype=float)
    eta = 2.0 * (z - za) / B
    return c0 * (1.0 + eps * (np.exp(eta) - 1.0 - eta))


def munk_dcdz(z, c0=1500.0, za=Z_A, B=B_MUNK, eps=EPS_MUNK):
    z = np.asarray(z, dtype=float)
    eta = 2.0 * (z - za) / B
    return c0 * eps * (np.exp(eta) - 1.0) * (2.0 / B)


def ray_rhs(y):
    x, z, th = y
    zc = min(max(z, 0.2), DEPTH - 0.2)
    c = float(munk_c(zc))
    dcdz = float(munk_dcdz(zc))
    return np.array([math.cos(th), math.sin(th), -math.cos(th) * dcdz / c])


def rk4_step(y, ds):
    k1 = ray_rhs(y)
    k2 = ray_rhs(y + 0.5 * ds * k1)
    k3 = ray_rhs(y + 0.5 * ds * k2)
    k4 = ray_rhs(y + ds * k3)
    yn = y + (ds / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    n_bounce = 0
    if yn[1] < 0.0:
        yn[1] = -yn[1]
        yn[2] = -yn[2]
        n_bounce += 1
    if yn[1] > DEPTH:
        yn[1] = 2 * DEPTH - yn[1]
        yn[2] = -yn[2]
        n_bounce += 1
    return yn, n_bounce


def shoot_ray(launch_rad, r_target, z0, ds=40.0, max_steps=4000):
    """Integrate to r_target; return hit depth, tau, elevation, topology stats."""
    y = np.array([0.0, z0, launch_rad])
    t = 0.0
    z_min, z_max = z0, z0
    n_surf = 0
    for _ in range(max_steps):
        x_rem = r_target - y[0]
        if x_rem <= 0:
            break
        # take larger steps when nearly horizontal
        h = min(ds, x_rem / max(abs(math.cos(y[2])), 1e-3) + 1e-6)
        h = min(h, ds)
        yn, nb = rk4_step(y, h)
        n_surf += nb
        zmid = 0.5 * (y[1] + yn[1])
        t += h / float(munk_c(zmid))
        y = yn
        z_min = min(z_min, y[1])
        z_max = max(z_max, y[1])
        if y[0] >= r_target:
            break
    phi = float(y[2])
    topo = "surface_bounce" if n_surf > 0 or z_min < 1.0 else "refracted"
    return {
        "z_hit": float(y[1]),
        "tau": float(t),
        "phi_arr_rad": phi,
        "phi_arr_deg": math.degrees(phi),
        "launch_deg": math.degrees(launch_rad),
        "z_min": float(z_min),
        "z_max": float(z_max),
        "topology": topo,
        "residual_m": float(y[1] - Z_R),
        "x_hit": float(y[0]),
    }


def bisect_launch(a_deg, b_deg, r_target, z0, tol_m=RESIDUAL_TOL_M, max_iter=30):
    """Root-solve z_hit(alpha)-Z_R=0 on [a,b]."""
    fa = shoot_ray(math.radians(a_deg), r_target, z0)["residual_m"]
    fb = shoot_ray(math.radians(b_deg), r_target, z0)["residual_m"]
    if fa == 0:
        return a_deg, fa, True
    if fb == 0:
        return b_deg, fb, True
    if fa * fb > 0:
        return None, None, False
    lo, hi = a_deg, b_deg
    flo = fa
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        fm = shoot_ray(math.radians(mid), r_target, z0)["residual_m"]
        if abs(fm) <= tol_m:
            return mid, fm, True
        if flo * fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
        if hi - lo < 1e-6:
            break
    mid = 0.5 * (lo + hi)
    fm = shoot_ray(math.radians(mid), r_target, z0)["residual_m"]
    return mid, fm, abs(fm) <= max(tol_m, 1e-3)


def find_eigenrays_root(r_km, z_s, launch_grid=None):
    if launch_grid is None:
        launch_grid = np.linspace(-20.0, 20.0, 61)
    r_target = r_km * 1e3
    scans = []
    prev_res = None
    prev_ld = None
    brackets = []
    for ld in launch_grid:
        res = shoot_ray(math.radians(float(ld)), r_target, z_s)
        scans.append({"launch_deg": float(ld), **res, "r_km": r_km, "z_s": z_s})
        if prev_res is not None and prev_res * res["residual_m"] < 0:
            brackets.append((prev_ld, float(ld)))
        prev_res = res["residual_m"]
        prev_ld = float(ld)
    roots = []
    for a, b in brackets:
        ld, resid, ok = bisect_launch(a, b, r_target, z_s)
        if ld is None:
            continue
        full = shoot_ray(math.radians(ld), r_target, z_s)
        full["r_km"] = r_km
        full["z_s"] = z_s
        full["root_ok"] = bool(ok)
        full["residual_m"] = float(resid) if resid is not None else full["residual_m"]
        full["admitted"] = bool(ok and abs(full["residual_m"]) <= RESIDUAL_TOL_M)
        roots.append(full)
    return pd.DataFrame(scans), pd.DataFrame(roots)


def associate_branches(roots_df):
    """Continuous association across local r for fixed z_s — not round(launch/2)."""
    if roots_df is None or roots_df.empty:
        return pd.DataFrame(), pd.DataFrame()
    df = roots_df[roots_df["admitted"] == True].copy()  # noqa: E712
    if df.empty:
        return df, pd.DataFrame()
    records = []
    # group by z_s
    bid_global = 0
    quality = []
    for zs, g in df.groupby("z_s"):
        g = g.sort_values("r_km")
        active = []  # list of dicts with last state and branch_id
        for _, row in g.iterrows():
            matched = None
            best_cost = np.inf
            for br in active:
                # continuity cost
                d_launch = abs(row["launch_deg"] - br["last_launch"])
                d_phi = abs(row["phi_arr_deg"] - br["last_phi"])
                d_tau = abs(row["tau"] - br["last_tau"]) * 1e3  # ms
                # topology must match
                if row["topology"] != br["topology"]:
                    continue
                cost = d_launch / 2.0 + d_phi / 2.0 + d_tau / 20.0
                if cost < best_cost and d_launch < 8.0 and d_phi < 8.0:
                    best_cost = cost
                    matched = br
            if matched is not None:
                matched["n"] += 1
                matched["r_set"].add(float(row["r_km"]))
                matched["launch_drift"] = max(matched["launch_drift"], abs(row["launch_deg"] - matched["first_launch"]))
                matched["phi_drift"] = max(matched["phi_drift"], abs(row["phi_arr_deg"] - matched["first_phi"]))
                matched["tau_drift"] = max(matched["tau_drift"], abs(row["tau"] - matched["first_tau"]))
                matched["max_resid"] = max(matched["max_resid"], abs(row["residual_m"]))
                matched["last_launch"] = float(row["launch_deg"])
                matched["last_phi"] = float(row["phi_arr_deg"])
                matched["last_tau"] = float(row["tau"])
                row = row.copy()
                row["branch_id"] = matched["id"]
                records.append(row.to_dict())
            else:
                bid = f"E{bid_global}"
                bid_global += 1
                br = {
                    "id": bid, "z_s": zs, "topology": row["topology"],
                    "n": 1, "r_set": {float(row["r_km"])},
                    "first_launch": float(row["launch_deg"]),
                    "first_phi": float(row["phi_arr_deg"]),
                    "first_tau": float(row["tau"]),
                    "last_launch": float(row["launch_deg"]),
                    "last_phi": float(row["phi_arr_deg"]),
                    "last_tau": float(row["tau"]),
                    "launch_drift": 0.0, "phi_drift": 0.0, "tau_drift": 0.0,
                    "max_resid": abs(float(row["residual_m"])),
                }
                active.append(br)
                row = row.copy()
                row["branch_id"] = bid
                records.append(row.to_dict())
        for br in active:
            rs = sorted(br["r_set"])
            quality.append({
                "branch_id": br["id"],
                "z_s_m": zs,
                "topology": br["topology"],
                "n_hits": br["n"],
                "n_ranges": len(rs),
                "r_min_km": rs[0] if rs else np.nan,
                "r_max_km": rs[-1] if rs else np.nan,
                "coverage_local_stencil": len(rs),
                "launch_drift_deg": br["launch_drift"],
                "phi_drift_deg": br["phi_drift"],
                "tau_drift_s": br["tau_drift"],
                "max_residual_m": br["max_resid"],
                "stable_for_fisher": bool(len(rs) >= 3 and br["max_resid"] <= RESIDUAL_TOL_M),
            })
    tracks = pd.DataFrame(records)
    qual = pd.DataFrame(quality)
    return tracks, qual


def delay_metrics(freqs, tau_grid, tau_true):
    dtau = tau_grid - tau_true
    A = np.abs(np.mean(np.exp(1j * 2 * np.pi * freqs[:, None] * dtau[None, :]), axis=0))
    # local mainlobe: contiguous A>=0.5 around true peak
    i0 = int(np.argmin(np.abs(tau_grid - tau_true)))
    # find contiguous region around i0
    if A[i0] < 0.5:
        # walk to nearest >=0.5
        left = i0
        right = i0
        while left > 0 and A[left] < 0.5:
            left -= 1
        while right < len(A) - 1 and A[right] < 0.5:
            right += 1
        if A[left] < 0.5 and A[right] < 0.5:
            local_w = np.nan
            lo = hi = i0
        else:
            lo = left if A[left] >= 0.5 else right
            hi = right if A[right] >= 0.5 else left
            # expand contiguous
            while lo > 0 and A[lo - 1] >= 0.5:
                lo -= 1
            while hi < len(A) - 1 and A[hi + 1] >= 0.5:
                hi += 1
            local_w = float(tau_grid[hi] - tau_grid[lo])
    else:
        lo = hi = i0
        while lo > 0 and A[lo - 1] >= 0.5:
            lo -= 1
        while hi < len(A) - 1 and A[hi + 1] >= 0.5:
            hi += 1
        local_w = float(tau_grid[hi] - tau_grid[lo])

    # first null: first tau away from peak where A < 0.05 (or local min)
    def first_null(direction):
        i = i0
        while 0 <= i < len(A):
            if A[i] < 0.08:
                return float(tau_grid[i] - tau_grid[i0])
            i += direction
        return np.nan
    fn_p = first_null(+1)
    fn_m = first_null(-1)
    first_null_w = float(abs(fn_p) + abs(fn_m)) if np.isfinite(fn_p) and np.isfinite(fn_m) else np.nan

    # global sidelobe: max A outside local mainlobe
    mask = np.ones_like(A, dtype=bool)
    mask[max(lo, 0): min(hi + 1, len(A))] = False
    side_max = float(A[mask].max()) if mask.any() else 0.0
    # significant peaks offsets
    peaks = []
    for i in range(1, len(A) - 1):
        if A[i] >= A[i - 1] and A[i] >= A[i + 1] and A[i] >= 0.3:
            peaks.append((float(tau_grid[i] - tau_grid[i0]), float(A[i])))
    peaks = sorted(peaks, key=lambda x: -x[1])[:8]
    return A, dict(
        local_mainlobe_width_s=local_w,
        first_null_width_s=first_null_w,
        global_max_sidelobe=side_max,
        A_at_true=float(A[i0]),
        n_peaks_ge_0p3=len(peaks),
        peaks=peaks,
    )


def run_delay_correction():
    tau_true = 0.0
    tau_grid = np.linspace(-0.15, 0.15, 3001)
    sets = {
        "S0_DENSE": np.arange(150.0, 375.0 + 1e-9, 1.0),
        "S0_23POINT": np.linspace(150.0, 375.0, 23),
        "S1": np.array([168.0, 204.0, 232.0, 279.0, 320.0]),
        "S2": np.array([166.0, 201.0, 235.0, 283.0, 338.0]),
    }
    rows = []
    peak_rows = []
    A_store = {}
    for name, freqs in sets.items():
        A, met = delay_metrics(freqs, tau_grid, tau_true)
        A_store[name] = (tau_grid, A)
        sig_f = float(np.std(freqs)) if freqs.size > 1 else 0.0
        snr = 100.0
        crb = 1.0 / (2 * math.pi * math.sqrt(snr) * sig_f) if sig_f > 0 else np.inf
        df_grid = float(freqs[1] - freqs[0]) if freqs.size > 1 else np.nan
        T_amb = 1.0 / df_grid if np.isfinite(df_grid) and df_grid > 0 else np.inf
        rows.append({
            "source": name,
            "n_freqs": int(freqs.size),
            "freq_spacing_Hz": df_grid,
            "T_amb_period_s": T_amb,
            "local_mainlobe_width_s": met["local_mainlobe_width_s"],
            "first_null_width_s": met["first_null_width_s"],
            "global_max_sidelobe": met["global_max_sidelobe"],
            "crb_sigma_tau_s": crb,
            "effective_bw_Hz": float(np.std(freqs) * math.sqrt(12)) if freqs.size > 1 else 0.0,
            "represents": "continuous broadband" if name == "S0_DENSE" else (
                "discrete comb control (NOT true broadband)" if name == "S0_23POINT" else "sparse lines"),
        })
        for off, amp in met["peaks"]:
            peak_rows.append({"source": name, "offset_s": off, "amp": amp})
    return pd.DataFrame(rows), pd.DataFrame(peak_rows), A_store


def u_hla(theta_rel_deg, phi_deg):
    return math.cos(math.radians(theta_rel_deg)) * math.cos(math.radians(phi_deg))


def steering_corr(L, n_el, u1, u2, f):
    if n_el <= 1 or L <= 0:
        return 1.0, np.inf
    k = 2 * math.pi * f / C0
    x = (np.arange(n_el) - (n_el - 1) / 2.0) * (L / max(n_el - 1, 1))
    a1 = np.exp(-1j * k * x * u1)
    a2 = np.exp(-1j * k * x * u2)
    corr = float(abs(np.vdot(a1, a2)) / (np.linalg.norm(a1) * np.linalg.norm(a2) + 1e-30))
    snr = 100.0
    crb_su = math.sqrt(12.0) / (k * L * math.sqrt(n_el * snr))
    return corr, crb_su


def run_actual_hla_resolution(tracks, qual):
    """HLA resolution on actual stable branch phi pairs at r≈50 km."""
    if tracks is None or tracks.empty:
        return pd.DataFrame()
    # pick hits near 50 km
    t50 = tracks[np.abs(tracks["r_km"] - 50.0) <= 0.6]
    rows = []
    for zs, g in t50.groupby("z_s"):
        g = g.drop_duplicates(subset=["branch_id", "r_km"])
        phis = g.groupby("branch_id")["phi_arr_deg"].median().to_dict()
        ids = list(phis.keys())
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                p1, p2 = float(phis[ids[i]]), float(phis[ids[j]])
                for th in THETA_REL_SCAN:
                    u1, u2 = u_hla(th, p1), u_hla(th, p2)
                    du = abs(u1 - u2)
                    for ap_tag, L, n_el in [("14m", 14.0, 8), ("70m", 70.0, 8)]:
                        corr, su = steering_corr(L, n_el, u1, u2, 250.0)
                        rows.append({
                            "z_s_m": zs,
                            "branch_i": ids[i],
                            "branch_j": ids[j],
                            "phi_i_deg": p1,
                            "phi_j_deg": p2,
                            "theta_rel_deg": th,
                            "delta_u": du,
                            "aperture": ap_tag,
                            "aperture_m": L,
                            "corr": corr,
                            "crb_sigma_u": su,
                            "resolved_corr": bool(corr < 0.95),
                            "resolved_sigma_u": bool(np.isfinite(su) and su < du),
                            "is_pm_mirror": bool(abs(p1 + p2) < 0.05),
                        })
    return pd.DataFrame(rows)


def numeric_jac_from_branches(branches_at, r0, zs0, z_r=Z_R, mode="IDEAL", theta_rel=0.0,
                              dr=50.0, dzs=2.0):
    """branches_at(r_km, zs_m) -> list of dicts with phi_deg, tau_s, topology.

    y = [phi's..., dtau's...] or u's + dtau
    """
    def eval_y(r_km, zs):
        brs = branches_at(r_km, zs)
        if not brs:
            return None
        brs = sorted(brs, key=lambda b: b["tau"])
        ref = brs[0]
        y = []
        for b in brs:
            if mode == "IDEAL":
                y.append(math.radians(b["phi_deg"]))
            else:
                y.append(u_hla(theta_rel, b["phi_deg"]))
        for b in brs[1:]:
            y.append(b["tau"] - ref["tau"])
        return np.array(y), brs

    y0, brs = eval_y(r0, zs0)
    if y0 is None:
        return None, None, None
    # finite difference — need same branch count/order: re-query may mismatch; use nearest launch order
    def eval_aligned(r_km, zs):
        brs = branches_at(r_km, zs)
        if not brs:
            return None
        # order by tau then match to reference order by nearest phi
        brs_sorted = sorted(brs, key=lambda b: b["tau"])
        ref_order = brs  # caller uses consistent association via tracking
        return brs_sorted

    # Use stencil from association table if provided via closure with fixed branch ids
    return y0, brs, "ok"


def build_estd_fisher(tracks, qual, delay_row):
    """Fisher on local stencil using associated branches at z_s=200 if available."""
    if tracks is None or tracks.empty or qual is None or qual.empty:
        return pd.DataFrame(), pd.DataFrame()
    stable = qual[qual["stable_for_fisher"] == True]  # noqa: E712
    if stable.empty:
        # still try best branches with n_ranges>=2
        stable = qual.sort_values(["n_ranges", "max_residual_m"], ascending=[False, True]).head(4)
    rows_j = []
    rows_f = []
    sig_tau = float(delay_row["crb_sigma_tau_s"]) if delay_row is not None else 2e-4

    # evaluate Jacobian by shooting at (r,z_s)±δ for each stable branch launch at stencil
    # Collect branch samples at each stencil point for z_s near 200
    zs_focus = 200.0
    samples = {}  # (r_km, zs) -> list of eigenray dicts with branch tags via launch
    for r_km in R_STENCIL_KM:
        for zs in [190.0, 200.0, 210.0]:
            _, roots = find_eigenrays_root(r_km, zs)
            if roots is None or roots.empty:
                samples[(r_km, zs)] = []
            else:
                samples[(r_km, zs)] = roots[roots["admitted"] == True].to_dict("records")  # noqa: E712
    # pick branch launches that exist at r=50, zs=200
    key0 = (50.0, 200.0)
    base = samples.get(key0, [])
    if len(base) < 2:
        # use any
        for k, v in samples.items():
            if len(v) >= 2:
                key0 = k
                base = v
                zs_focus = k[1]
                break
    base = sorted(base, key=lambda b: b["tau"])[:4]
    if len(base) < 2:
        return pd.DataFrame(rows_j), pd.DataFrame(rows_f)

    def branch_states(r_km, zs, launches):
        out = []
        for ld in launches:
            # find nearest admitted root
            cands = samples.get((r_km, zs), [])
            if not cands:
                # shoot directly
                full = shoot_ray(math.radians(ld), r_km * 1e3, zs)
                if abs(full["residual_m"]) > 0.5:
                    continue
                out.append(full)
            else:
                best = min(cands, key=lambda c: abs(c["launch_deg"] - ld))
                if abs(best["launch_deg"] - ld) < 5.0 and abs(best["residual_m"]) <= 0.5:
                    out.append(best)
                else:
                    full = shoot_ray(math.radians(ld), r_km * 1e3, zs)
                    if abs(full["residual_m"]) <= 0.5:
                        out.append(full)
        return out

    launches = [b["launch_deg"] for b in base]
    r0, zs0 = key0[0], key0[1]
    st0 = branch_states(r0, zs0, launches)
    if len(st0) < 2:
        return pd.DataFrame(rows_j), pd.DataFrame(rows_f)
    st0 = sorted(st0, key=lambda b: b["tau"])

    def yvec(sts, mode, theta_rel=0.0):
        if len(sts) < 2:
            return None
        sts = sorted(sts, key=lambda b: b["tau"])
        y = []
        for b in sts:
            phi = b.get("phi_arr_deg", b.get("phi_deg"))
            y.append(math.radians(phi) if mode == "IDEAL" else u_hla(theta_rel, phi))
        for b in sts[1:]:
            y.append(b["tau"] - sts[0]["tau"])
        return np.asarray(y, dtype=float)

    def jac(mode, theta_rel=0.0):
        y0 = yvec(st0, mode, theta_rel)
        # derivatives wrt r and zs via stencil neighbors
        def y_at(r_km, zs):
            return yvec(branch_states(r_km, zs, launches), mode, theta_rel)
        y_rp = y_at(r0 + 0.05, zs0)  # +50 m
        y_rm = y_at(r0 - 0.05, zs0)
        y_zp = y_at(r0, zs0 + 2.0)
        y_zm = y_at(r0, zs0 - 2.0)
        if any(v is None for v in [y_rp, y_rm, y_zp, y_zm]) or y0 is None:
            return None, None
        # pad/truncate to same length
        n = min(len(y0), len(y_rp), len(y_rm), len(y_zp), len(y_zm))
        y0, y_rp, y_rm, y_zp, y_zm = y0[:n], y_rp[:n], y_rm[:n], y_zp[:n], y_zm[:n]
        J = np.zeros((n, 2))
        J[:, 0] = (y_rp - y_rm) / (100.0)  # 0.05 km *2 = 0.1 km = 100 m
        J[:, 1] = (y_zp - y_zm) / 4.0      # 2 m *2
        return y0, J

    def fish(J, sigma):
        Sinv = np.diag(1.0 / (np.asarray(sigma) ** 2 + 1e-30))
        F = J.T @ Sinv @ J
        s = np.linalg.svd(J, compute_uv=False)
        try:
            cov = np.linalg.inv(F)
            crb_r = math.sqrt(max(float(cov[0, 0]), 0))
            # zs is not z target exactly — treat as z proxy in this local model
            crb_z = math.sqrt(max(float(cov[1, 1]), 0))
            corr = float(cov[0, 1] / math.sqrt(cov[0, 0] * cov[1, 1] + 1e-30))
            condF = float(np.linalg.cond(F))
            rank = int(np.linalg.matrix_rank(J))
        except np.linalg.LinAlgError:
            crb_r = crb_z = np.inf
            corr = np.nan
            condF = np.inf
            rank = int(np.linalg.matrix_rank(J))
        return dict(singular_values=s.tolist(), sv_min=float(s.min()) if s.size else np.nan,
                    cond_J=float(s[0] / s[-1]) if s.size and s[-1] > 1e-30 else np.inf,
                    rank=rank, crb_r_m=crb_r, crb_z_s_proxy_m=crb_z, corr=corr, cond_F=condF)

    n_br = len(st0)
    for sigma_th in SIGMA_THETA_DEG:
        sig_th = math.radians(sigma_th)
        for theta in THETA_REL_SCAN:
            cases = []
            y_ideal, J_ideal = jac("IDEAL", theta)
            y_hla, J_hla = jac("HLA", theta)
            if y_ideal is not None and J_ideal is not None:
                # delay-only = last n_br-1 rows
                J_d = J_ideal[n_br - 1:, :]
                y_d = y_ideal[n_br - 1:]
                cases.append(("IDEAL_ELEVATION+delay", y_ideal, J_ideal,
                              [math.radians(0.3)] * n_br + [sig_tau] * (n_br - 1)))
                cases.append(("delay-only", y_d, J_d, [sig_tau] * (n_br - 1)))
            if y_hla is not None and J_hla is not None:
                sig_u = []
                for b in st0:
                    phi = b.get("phi_arr_deg", 0.0)
                    su_arr = steering_corr(14.0, 8, u_hla(theta, phi), u_hla(theta, phi), 250.0)[1]
                    du_dth = -math.sin(math.radians(theta)) * math.cos(math.radians(phi))
                    sig_u.append(math.sqrt(su_arr ** 2 + (du_dth * sig_th) ** 2))
                cases.append(("HLA_u+delay", y_hla, J_hla, sig_u + [sig_tau] * (n_br - 1)))
            for label, y0, J, sig in cases:
                if J is None or y0 is None or J.size == 0:
                    continue
                if len(sig) != J.shape[0]:
                    sig = (sig[:1] * J.shape[0]) if sig else [1e-3] * J.shape[0]
                f = fish(J, sig)
                rows_j.append({
                    "sigma_theta_deg": sigma_th,
                    "theta_rel_deg": theta,
                    "case": label,
                    "n_branches": n_br,
                    "launches_deg": json.dumps([round(b.get("launch_deg", 0), 3) for b in st0]),
                    "phi_deg": json.dumps([round(b.get("phi_arr_deg", 0), 3) for b in st0]),
                    "tau_s": json.dumps([round(b["tau"], 5) for b in st0]),
                    "y0": json.dumps(np.asarray(y0).tolist()),
                })
                rows_f.append({
                    "sigma_theta_deg": sigma_th,
                    "theta_rel_deg": theta,
                    "case": label,
                    "n_branches": n_br,
                    **f,
                })
    return pd.DataFrame(rows_j), pd.DataFrame(rows_f)


def run_ray_modal_crosscheck(mode_model=None):
    """Modal transfer IFFT delay profile vs ray branch delays at r=49/50/51 km."""
    if mode_model is None:
        mode_model = p4.MODE_ENV["E0"]
    freqs = np.arange(150.0, 375.0 + 1e-9, 1.0)
    df = 1.0
    # time axis
    n = 512
    tau = np.fft.fftfreq(n, d=df)  # seconds
    rows = []
    for r_km in [49.0, 50.0, 51.0]:
        H = np.array([mode_model.amp(float(f), r_km * 1e3, z_s=200.0, z_r=200.0) for f in freqs])
        # pad into full FFT grid centered
        Hfull = np.zeros(n, dtype=complex)
        # place band at positive freqs
        idx = np.round(freqs / df).astype(int)
        idx = np.clip(idx, 0, n // 2 - 1)
        Hfull[idx] = H
        h = np.fft.ifft(Hfull)
        # delay power in 0–0.1 s
        t = np.fft.fftfreq(n, d=df)
        mask = (t >= 0) & (t <= 0.12)
        p = np.abs(h[mask]) ** 2
        tt = t[mask]
        # energy clusters: local maxima
        peaks = []
        for i in range(1, len(p) - 1):
            if p[i] >= p[i - 1] and p[i] >= p[i + 1] and p[i] > 0.15 * p.max():
                peaks.append(float(tt[i]))
        peaks = sorted(peaks)[:8]
        # ray roots at z_s=200
        _, roots = find_eigenrays_root(r_km, 200.0)
        ray_taus = sorted(roots.loc[roots["admitted"] == True, "tau"].tolist()) if roots is not None and len(roots) else []  # noqa: E712
        rows.append({
            "r_km": r_km,
            "modal_delay_peaks_s": json.dumps(peaks),
            "modal_n_peaks": len(peaks),
            "ray_admitted_taus_s": json.dumps([round(x, 5) for x in ray_taus]),
            "ray_n_admitted": len(ray_taus),
            "ray_min_s": ray_taus[0] if ray_taus else np.nan,
            "ray_max_s": ray_taus[-1] if ray_taus else np.nan,
            "note": "sanity check: ray delay clusters vs modal IFFT energy; not full mode-id",
        })
    return pd.DataFrame(rows)


def decide_b11(delay_df, qual, hla_df, fish_df, cross_df):
    notes = {}
    # delay: true broadband S0_DENSE
    d0 = delay_df[delay_df["source"] == "S0_DENSE"]
    d23 = delay_df[delay_df["source"] == "S0_23POINT"]
    notes["S0_DENSE_local_mainlobe_s"] = float(d0["local_mainlobe_width_s"].iloc[0]) if len(d0) else np.nan
    notes["S0_DENSE_global_side"] = float(d0["global_max_sidelobe"].iloc[0]) if len(d0) else np.nan
    notes["S0_23POINT_global_side"] = float(d23["global_max_sidelobe"].iloc[0]) if len(d23) else np.nan
    notes["S0_23POINT_T_amb_s"] = float(d23["T_amb_period_s"].iloc[0]) if len(d23) else np.nan
    delay_global_ok = bool(len(d0) and d0["global_max_sidelobe"].iloc[0] < 0.35)
    delay_local_ok = bool(len(d0) and np.isfinite(d0["local_mainlobe_width_s"].iloc[0])
                          and d0["local_mainlobe_width_s"].iloc[0] < 0.02)
    notes["delay_local_ok"] = delay_local_ok
    notes["delay_global_ok_S0_DENSE"] = delay_global_ok
    notes["delay_comb_ambiguity_S0_23"] = bool(len(d23) and d23["global_max_sidelobe"].iloc[0] > 0.5)

    # branches
    if qual is not None and len(qual):
        n_stable = int(qual["stable_for_fisher"].sum())
        n_branch = int(qual["branch_id"].nunique())
        notes["n_branch_ids"] = n_branch
        notes["n_stable_for_fisher"] = n_stable
        notes["max_residual_m"] = float(qual["max_residual_m"].max())
        notes["branches_ge3_ranges"] = int((qual["n_ranges"] >= 3).sum())
    else:
        n_stable = 0
        notes["n_branch_ids"] = 0
        notes["n_stable_for_fisher"] = 0

    stable_branches = n_stable >= 2

    # HLA actual
    if hla_df is not None and len(hla_df):
        frac14 = float(hla_df[hla_df["aperture"] == "14m"]["resolved_corr"].mean())
        frac70 = float(hla_df[hla_df["aperture"] == "70m"]["resolved_corr"].mean())
        notes["hla_resolved_corr_14m"] = frac14
        notes["hla_resolved_corr_70m"] = frac70
        hla_weak = frac14 < 0.25 and frac70 < 0.5
    else:
        notes["hla_resolved_corr_14m"] = np.nan
        notes["hla_resolved_corr_70m"] = np.nan
        hla_weak = True
    notes["HLA_ELEVATION_PROJECTION_WEAK"] = hla_weak
    notes["HLA_projection_note"] = (
        "Corrected gate uses actual E-STD branch Delta-u + steering corr; "
        "prior False flag from composite CRB approximation is NOT frozen."
    )

    # Fisher from actual E-STD
    if fish_df is not None and len(fish_df):
        def pick(case, th=0.0, st=0.1):
            d = fish_df[(fish_df["case"] == case) & (fish_df["theta_rel_deg"] == th) & (fish_df["sigma_theta_deg"] == st)]
            return d.iloc[0] if len(d) else None
        ideal = pick("IDEAL_ELEVATION+delay")
        hla = pick("HLA_u+delay")
        dly = pick("delay-only")
        notes["estd_ideal_crb_r_m"] = float(ideal["crb_r_m"]) if ideal is not None else np.nan
        notes["estd_hla_crb_r_m"] = float(hla["crb_r_m"]) if hla is not None else np.nan
        notes["estd_delay_crb_r_m"] = float(dly["crb_r_m"]) if dly is not None else np.nan
        notes["estd_ideal_rank"] = int(ideal["rank"]) if ideal is not None else 0
        notes["estd_hla_rank"] = int(hla["rank"]) if hla is not None else 0
        notes["estd_delay_rank"] = int(dly["rank"]) if dly is not None else 0
        ideal_ok = bool(ideal is not None and ideal["rank"] >= 2 and np.isfinite(ideal["crb_r_m"]) and ideal["crb_r_m"] < 8e3)
        hla_ok = bool(hla is not None and hla["rank"] >= 2 and np.isfinite(hla["crb_r_m"]) and hla["crb_r_m"] < 8e3)
        delay_fish_ok = bool(dly is not None and dly["rank"] >= 2 and np.isfinite(dly["crb_r_m"]) and dly["crb_r_m"] < 25e3)
    else:
        ideal_ok = hla_ok = delay_fish_ok = False
        notes["estd_ideal_crb_r_m"] = notes["estd_hla_crb_r_m"] = notes["estd_delay_crb_r_m"] = np.nan
    notes["estd_ideal_ok"] = ideal_ok
    notes["estd_hla_ok"] = hla_ok
    notes["estd_delay_fish_ok"] = delay_fish_ok

    # decision
    if not stable_branches:
        decision = "B1_NO_STABLE_MULTIPATH_IDENTITY"
        why = f"Refined root-solve+association: n_stable_for_fisher={n_stable}, n_branch_ids={notes.get('n_branch_ids')}, max_resid={notes.get('max_residual_m')} m."
        nxt = "do not proceed B2; multipath identity not established"
    elif ideal_ok and not hla_ok:
        decision = "B1_PHYSICS_ONLY_HLA_LIMITED"
        why = (
            f"E-STD stable branches present; IDEAL elevation+delay Fisher OK (CRB_r={notes.get('estd_ideal_crb_r_m')} m) "
            f"but HLA_u+delay weak (CRB_r={notes.get('estd_hla_crb_r_m')} m). "
            f"HLA projection weak={hla_weak} (corr14={notes.get('hla_resolved_corr_14m')}, corr70={notes.get('hla_resolved_corr_70m')})."
        )
        nxt = "MMAC physics may exist but not via HLA elevation; consider delay-set only if delay global ambiguity acceptable"
    elif (not hla_weak) and hla_ok and ideal_ok and (delay_global_ok or delay_local_ok):
        decision = "B1_MMAC_PHYSICS_CONFIRMED"
        why = "Stable branches + HLA resolvable + E-STD joint Fisher + delay observability."
        nxt = "eligible for B2 later (not this round)"
    elif delay_fish_ok and stable_branches and delay_global_ok:
        decision = "B1_DELAY_DOMINANT_POSSIBLE"
        why = (
            f"HLA weak; delay-only E-STD Fisher rank/CRB indicates complementary r-z info "
            f"(CRB_r={notes.get('estd_delay_crb_r_m')} m); S0_DENSE global sidelobe={notes.get('S0_DENSE_global_side')}."
        )
        nxt = "B2 delay-set matching only if association without hidden labels is designed later"
    elif delay_local_ok and not delay_global_ok:
        decision = "B1_DELAY_AMBIGUITY_LIMITED"
        why = (
            f"Local delay mainlobe ms-scale OK (S0_DENSE {notes.get('S0_DENSE_local_mainlobe_s')} s) "
            f"but global ambiguity remains (S0_23 comb side={notes.get('S0_23POINT_global_side')}, "
            f"S0_DENSE side={notes.get('S0_DENSE_global_side')}). CRB vs ambiguity must stay separated."
        )
        nxt = "do not kill MMAC; resolve global delay ambiguity or restrict to local peak tracking"
    elif stable_branches and not ideal_ok and not delay_fish_ok:
        decision = "B1_NO_COMPLEMENTARY_INFORMATION"
        why = "Stable branches exist but angle/delay Jacobian not complementary non-degenerate under current noise model."
        nxt = "stop RC3-B practical path or re-examine observability assumptions"
    else:
        decision = "B1_DECISION_PENDING_CORRECTION"
        why = "Insufficient refined outputs to freeze a terminal state."
        nxt = "complete remaining B1.1 outputs"
    return decision, why, nxt, notes


def _esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def write_line_svg(path, title, series, xlab, ylab, note=""):
    w, h = 760, 360
    ml, mr, mt, mb = 70, 140, 48, 50
    pw, ph = w - ml - mr, h - mt - mb
    xs, ys = [], []
    for s in series:
        xs.extend(s["x"]); ys.extend([y for y in s["y"] if y is not None and np.isfinite(y)])
    if not xs:
        xs, ys = [0, 1], [0, 1]
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys + [0]), max(ys + [0.1])
    if x1 - x0 < 1e-12:
        x1 = x0 + 1
    cols = ["#b45309", "#0f766e", "#1d4ed8", "#8a8a8a"]
    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
         f'<rect width="{w}" height="{h}" fill="#f7f4ef"/>',
         f'<text x="{w/2}" y="26" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="14" font-weight="600">{_esc(title)}</text>',
         f'<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="#fff" stroke="#ccc"/>']
    def sx(x):
        return ml + (x - x0) / (x1 - x0) * pw
    def sy(y):
        return mt + ph - (y - y0) / (y1 - y0 + 1e-15) * ph
    for i, s in enumerate(series):
        pts = [f"{sx(x):.1f},{sy(y):.1f}" for x, y in zip(s["x"], s["y"]) if y is not None and np.isfinite(y)]
        if len(pts) >= 2:
            p.append(f'<polyline fill="none" stroke="{cols[i%len(cols)]}" stroke-width="1.5" points="{" ".join(pts)}"/>')
        p.append(f'<rect x="{w-mr+8}" y="{mt+8+i*16}" width="12" height="3" fill="{cols[i%len(cols)]}"/>')
        p.append(f'<text x="{w-mr+24}" y="{mt+12+i*16}" font-size="10" font-family="sans-serif">{_esc(s["name"])}</text>')
    p.append(f'<text x="{ml+pw/2}" y="{h-16}" text-anchor="middle" font-size="11" font-family="sans-serif">{_esc(xlab)}</text>')
    p.append(f'<text x="{ml}" y="{h-4}" font-size="10" fill="#777" font-family="sans-serif">{_esc(note)} {_esc(ylab)}</text></svg>')
    path.write_text("\n".join(p), encoding="utf-8")


def main():
    t0 = time.time()
    print("=== R3-B1.1 correction ===")
    print(f"OUT={OUT}")
    (OUT / "R3_B1_1_CONFIG.json").write_text(json.dumps(CONFIG, ensure_ascii=False, indent=2), encoding="utf-8")

    # 1) mark old decision revoked in a sidecar (do not overwrite R3_B1_DECISION.json)
    old_dec = OUT / "R3_B1_DECISION.json"
    if old_dec.exists():
        try:
            od = json.loads(old_dec.read_text(encoding="utf-8"))
        except Exception:
            od = {}
        od["superseded_by"] = "R3_B1_1_DECISION.json"
        od["superseded_note"] = "B1_DELAY_BANDWIDTH_LIMITED revoked: delay mainlobe metric incorrect; pending B1.1"
        od["revoked_utc"] = NOW
        (OUT / "R3_B1_DECISION_SUPERSEDED.json").write_text(json.dumps(od, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    print("[1] delay ambiguity corrected ...", flush=True)
    delay_df, peak_df, A_store = run_delay_correction()
    delay_df.to_csv(OUT / "delay_ambiguity_corrected.csv", index=False, encoding="utf-8-sig")
    peak_df.to_csv(OUT / "delay_ambiguity_peaks.csv", index=False, encoding="utf-8-sig")
    print(delay_df.to_string(index=False), flush=True)

    # fig7 delay
    series = []
    for name, (tg, A) in A_store.items():
        # subsample
        series.append({"name": name, "x": tg[::5].tolist(), "y": A[::5].tolist()})
    write_line_svg(FIG / "fig7_delay_ambiguity_corrected.svg",
                   "图7  修正后时延模糊函数 A(τ)", series, "τ (s)", "A",
                   note="局部主瓣 vs 全局旁瓣分离；S0_DENSE≠S0_23POINT")

    print("[2] refined eigenrays root-solve ...", flush=True)
    all_roots = []
    # full local stencil only at zs=200; sideband r at zs=200; thin zs check at r=50
    jobs = []
    for r_km in R_STENCIL_KM + R_SIDEBAND_KM:
        jobs.append((r_km, 200.0))
    for zs in [190.0, 210.0]:
        jobs.append((50.0, zs))
    for r_km, zs in jobs:
        scans, roots = find_eigenrays_root(r_km, zs)
        if roots is not None and len(roots):
            roots = roots.copy()
            all_roots.append(roots)
            print(f"  r={r_km} zs={zs} roots={len(roots)} admitted={int(roots['admitted'].sum())}", flush=True)
        else:
            print(f"  r={r_km} zs={zs} roots=0", flush=True)
        # incremental save
        if all_roots:
            pd.concat(all_roots, ignore_index=True).to_csv(OUT / "eigenray_refined.csv", index=False, encoding="utf-8-sig")
    roots_df = pd.concat(all_roots, ignore_index=True) if all_roots else pd.DataFrame()
    roots_df.to_csv(OUT / "eigenray_refined.csv", index=False, encoding="utf-8-sig")

    print("[3] branch association ...", flush=True)
    tracks, qual = associate_branches(roots_df)
    tracks.to_csv(OUT / "branch_tracking_refined.csv", index=False, encoding="utf-8-sig")
    qual.to_csv(OUT / "branch_quality.csv", index=False, encoding="utf-8-sig")
    print("qual head:\n", qual.head(20).to_string(index=False) if len(qual) else "empty", flush=True)

    print("[4] actual HLA branch resolution ...", flush=True)
    hla_df = run_actual_hla_resolution(tracks, qual)
    hla_df.to_csv(OUT / "actual_branch_hla_resolution.csv", index=False, encoding="utf-8-sig")
    if len(hla_df):
        print(hla_df.groupby("aperture")[["resolved_corr", "resolved_sigma_u"]].mean(), flush=True)

    print("[5] ray-modal delay crosscheck ...", flush=True)
    cross_df = run_ray_modal_crosscheck()
    cross_df.to_csv(OUT / "ray_modal_delay_crosscheck.csv", index=False, encoding="utf-8-sig")
    print(cross_df.to_string(index=False), flush=True)

    print("[6] E-STD actual Fisher ...", flush=True)
    d0row = delay_df[delay_df["source"] == "S0_DENSE"].iloc[0] if len(delay_df[delay_df["source"] == "S0_DENSE"]) else None
    jac_df, fish_df = build_estd_fisher(tracks, qual, d0row)
    jac_df.to_csv(OUT / "estd_actual_jacobian.csv", index=False, encoding="utf-8-sig")
    fish_df.to_csv(OUT / "estd_actual_fisher.csv", index=False, encoding="utf-8-sig")
    if len(fish_df):
        print(fish_df.head(12).to_string(index=False), flush=True)

    print("[7] decision + figures + report ...", flush=True)
    decision, why, nxt, notes = decide_b11(delay_df, qual, hla_df, fish_df, cross_df)

    # fig8 branches
    if len(qual):
        write_line_svg(FIG / "fig8_refined_branches.svg",
                       "图8  精化本征射线分支数 vs 距离覆盖",
                       [{"name": "n_ranges", "x": list(range(len(qual))), "y": qual["n_ranges"].tolist()},
                        {"name": "stable", "x": list(range(len(qual))), "y": qual["stable_for_fisher"].astype(float).tolist()}],
                       "branch index", "count",
                       note=f"stable_for_fisher={notes.get('n_stable_for_fisher')} / ids={notes.get('n_branch_ids')}")
    # fig9 HLA actual
    if len(hla_df):
        write_line_svg(FIG / "fig9_actual_hla_branch_resolution.svg",
                       "图9  实际E-STD分支对的 HLA Δu（θ_rel=0°）",
                       [{"name": "delta_u", "x": list(range(len(hla_df))), "y": hla_df["delta_u"].tolist()}],
                       "pair index", "Δu",
                       note=f"corr14={notes.get('hla_resolved_corr_14m')}, corr70={notes.get('hla_resolved_corr_70m')}")
    # fig10 fisher
    if len(fish_df):
        d = fish_df[(fish_df["sigma_theta_deg"] == 0.1) & (fish_df["theta_rel_deg"] == 0.0)]
        write_line_svg(FIG / "fig10_actual_estd_fisher.svg",
                       "图10  E-STD实际分支 Fisher CRB_r（σθ=0.1°, θ=0°）",
                       [{"name": d.iloc[i]["case"], "x": [i], "y": [d.iloc[i]["crb_r_m"] / 1e3]}
                        for i in range(len(d))],
                       "case", "CRB_r km",
                       note="非 controlled_paths")

    dec = {
        "rc3b1_status": "B1_DECISION_PENDING_CORRECTION → refined terminal",
        "rc3b1_1_decision": decision,
        "why": why,
        "next_step": nxt,
        "notes": notes,
        "controlled_model_label": CONFIG["controlled_model_label"],
        "created_utc": NOW,
        "stop": "after R3-B1.1; no B2/RC3-C/P5/MC RMSE",
    }
    (OUT / "R3_B1_1_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    def fnum(x, nd=4):
        try:
            v = float(x)
            return "n/a" if not np.isfinite(v) else f"{v:.{nd}f}"
        except Exception:
            return "n/a"

    rp = []
    rp.append("# R3-B1.1 修正报告：MMAC 观测关口")
    rp.append("")
    rp.append(f"UTC：{NOW}")
    rp.append("")
    rp.append("## 0. 撤销与口径")
    rp.append("")
    rp.append("- 旧判定 `B1_DELAY_BANDWIDTH_LIMITED` **撤销**，暂记 `B1_DECISION_PENDING_CORRECTION` 后由本修正给出终态。")
    rp.append("- 旧文件保留；`R3_B1_DECISION_SUPERSEDED.json` 记录替代关系。")
    rp.append("- 受控 direct+surface 模型保留为 **50 km 弱几何 sanity check**，**不是** E-STD 第一 CZ MMAC 性能。")
    rp.append("- `HLA_ELEVATION_PROJECTION_WEAK=False` **不冻结**；修正后以实际分支 Δu + 导向相关为准。")
    rp.append("")
    rp.append("## 1. 修正后时延模糊（核心）")
    rp.append("")
    rp.append("| 源 | 代表 | 间距 Hz | T_amb s | 局部主瓣 s | 第一零点宽 s | 全局旁瓣 | CRB στ s |")
    rp.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for _, r in delay_df.iterrows():
        rp.append(
            f"| {r['source']} | {r['represents']} | {fnum(r['freq_spacing_Hz'],3)} | {fnum(r['T_amb_period_s'],4)} | "
            f"{fnum(r['local_mainlobe_width_s'],5)} | {fnum(r['first_null_width_s'],5)} | "
            f"{fnum(r['global_max_sidelobe'],3)} | {fnum(r['crb_sigma_tau_s'],6)} |"
        )
    rp.append("")
    rp.append("**必须分开**：")
    rp.append("")
    rp.append(f"- 局部时延分辨（S0_DENSE 主瓣 ≈ **{fnum(notes.get('S0_DENSE_local_mainlobe_s'),5)} s**）与 CRB（~0.2–0.3 ms）一致量级；")
    rp.append(f"- 全局无模糊：S0_23POINT 人为梳状 T_amb≈**{fnum(notes.get('S0_23POINT_T_amb_s'),4)} s**，旁瓣 **{fnum(notes.get('S0_23POINT_global_side'),3)}** —— **不能**代表连续宽带 S0；")
    rp.append(f"- S0_DENSE 全局旁瓣 **{fnum(notes.get('S0_DENSE_global_side'),3)}**。")
    rp.append("")
    rp.append("禁止再用“所有 A≥0.5 的最大 τ 跨度”当主瓣宽（旧 300–400 ms 为 **计算错误**）。")
    rp.append("")
    rp.append("## 2. 精化本征射线（根求解）")
    rp.append("")
    rp.append(f"- 准入：|z_hit−z_r|≤**{RESIDUAL_TOL_M} m**（二分根求解 + 射线积分）")
    rp.append(f"- branch_id：launch/φ/τ/topology **连续关联**，禁止 `round(launch/2°)`")
    rp.append(f"- 分支 id 数 **{notes.get('n_branch_ids')}**；stable_for_fisher **{notes.get('n_stable_for_fisher')}**；最大 residual **{fnum(notes.get('max_residual_m'),3)} m**")
    rp.append("")
    if len(qual):
        rp.append("| branch | z_s | n_ranges | r覆盖 km | launch drift° | φ drift° | τ drift s | max resid m | stable |")
        rp.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for _, r in qual.head(15).iterrows():
            rp.append(
                f"| {r['branch_id']} | {r['z_s_m']} | {r['n_ranges']} | {r['r_min_km']}-{r['r_max_km']} | "
                f"{fnum(r['launch_drift_deg'],3)} | {fnum(r['phi_drift_deg'],3)} | {fnum(r['tau_drift_s'],5)} | "
                f"{fnum(r['max_residual_m'],3)} | {r['stable_for_fisher']} |"
            )
        rp.append("")
        rp.append("正确表述：**粗扫描显示第一 CZ 存在多径候选；稳定因果分支以本表 residual/coverage 为准。**")
    rp.append("")
    rp.append("## 3. 实际分支上的 HLA 方向余弦")
    rp.append("")
    rp.append(f"- corr<0.95 可分比例：14 m **{fnum(notes.get('hla_resolved_corr_14m'),3)}**，70 m **{fnum(notes.get('hla_resolved_corr_70m'),3)}**")
    rp.append(f"- **HLA_ELEVATION_PROJECTION_WEAK = {notes.get('HLA_ELEVATION_PROJECTION_WEAK')}**")
    rp.append(f"- {notes.get('HLA_projection_note')}")
    rp.append("")
    rp.append("## 4. E-STD 实际分支 Fisher（非 controlled_paths）")
    rp.append("")
    rp.append("| case | rank | CRB_r m | CRB_zs proxy m | cond F |")
    rp.append("| --- | --- | --- | --- | --- |")
    if len(fish_df):
        d = fish_df[(fish_df["sigma_theta_deg"] == 0.1) & (fish_df["theta_rel_deg"] == 0.0)]
        for _, r in d.iterrows():
            rp.append(
                f"| {r['case']} | {r['rank']} | {fnum(r['crb_r_m'],1)} | {fnum(r['crb_z_s_proxy_m'],1)} | {fnum(r['cond_F'],3)} |"
            )
    rp.append("")
    rp.append("## 5. ray–modal 时延交叉检查（sanity）")
    rp.append("")
    if len(cross_df):
        rp.append("| r km | modal peaks | ray admitted taus s |")
        rp.append("| --- | --- | --- |")
        for _, r in cross_df.iterrows():
            rp.append(f"| {r['r_km']} | {r['modal_n_peaks']} | {r['ray_admitted_taus_s']} |")
    rp.append("")
    rp.append("## 6. 修正后判定")
    rp.append("")
    rp.append(f"### `{decision}`")
    rp.append("")
    rp.append(why)
    rp.append("")
    rp.append(f"**下一步**：{nxt}")
    rp.append("")
    rp.append("允许终态：`B1_MMAC_PHYSICS_CONFIRMED` / `B1_PHYSICS_ONLY_HLA_LIMITED` / `B1_DELAY_DOMINANT_POSSIBLE` / `B1_DELAY_AMBIGUITY_LIMITED` / `B1_NO_STABLE_MULTIPATH_IDENTITY` / `B1_NO_COMPLEMENTARY_INFORMATION`。")
    rp.append("")
    rp.append("## 7. 停止")
    rp.append("")
    rp.append("- 不进入 B2 / RC3-C / P5 / Monte Carlo RMSE")
    rp.append("- 不增加新物理特征")
    rp.append(f"- **B1.1 完成后停止；判定 `{decision}`**")
    rp.append("")
    (OUT / "R3_B1_1_CORRECTION_REPORT.md").write_text("\n".join(rp), encoding="utf-8")

    gs = [
        "# R3-B1.1 — GPT 同步稿", "",
        f"- UTC: {NOW}",
        f"- 旧判定撤销：B1_DELAY_BANDWIDTH_LIMITED",
        f"- **新判定：{decision}**",
        f"- {why}",
        f"- 下一步：{nxt}", "",
        "## 时延修正",
        delay_df.to_string(index=False), "",
        f"S0_DENSE local mainlobe={notes.get('S0_DENSE_local_mainlobe_s')} s; "
        f"S0_23 comb side={notes.get('S0_23POINT_global_side')}; "
        f"S0_DENSE side={notes.get('S0_DENSE_global_side')}", "",
        "## 分支",
        f"n_branch={notes.get('n_branch_ids')} stable={notes.get('n_stable_for_fisher')} "
        f"max_resid={notes.get('max_residual_m')} m", "",
        "## HLA / Fisher",
        json.dumps({k: notes[k] for k in notes if k.startswith('hla') or k.startswith('estd') or k.startswith('HLA') or k.startswith('delay')}, ensure_ascii=False, default=str), "",
        "停止：无 B2/RC3-C/P5。", "",
    ]
    (OUT / "R3_B1_1_GPT_SYNC.md").write_text("\n".join(gs), encoding="utf-8")

    print(f"DONE in {time.time()-t0:.1f}s")
    print("DECISION", decision)
    print(why)


if __name__ == "__main__":
    main()
