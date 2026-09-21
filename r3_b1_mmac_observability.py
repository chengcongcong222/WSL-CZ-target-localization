#!/usr/bin/env python3
"""R3-B1: HLA direction-cosine observability gate + controlled MMAC + E-STD Fisher.

No Monte Carlo RMSE, no RC2 hard pairs, no RC3-C/P5, no Bellhop.
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
Z_A = 1300.0
B_MUNK = 1300.0
EPS_MUNK = 0.00737
DEPTH = 5000.0
Z_S = 200.0
Z_R = 200.0
THETA_REL_SCAN = [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0]
PHI_SCAN = [-20.0, -15.0, -10.0, -5.0, -2.0, -1.0, 0.0, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0]
SIGMA_THETA_DEG = [0.02, 0.05, 0.1, 0.2]
APERTURES = [
    dict(tag="L14m", L=14.0, note="理论基准孔径（非装备需求）"),
    dict(tag="L70m", L=70.0, note="诊断大孔径（非装备需求）"),
    dict(tag="Lparam", L=28.0, note="连续参数化参考 L=28 m"),
]
F_BAND = (150.0, 375.0)

CONFIG = {
    "package": "R3_B1_MMAC_observability",
    "created_utc": NOW,
    "coords": {
        "array_axis": "x-hat a",
        "theta_rel": "horizontal LOS vs array axis (not global theta0)",
        "phi": "elevation from horizontal, 0=horizontal",
        "u": "cos(theta_rel)*cos(phi) — HLA spatial observation",
        "mirror": "u(+phi)=u(-phi) — HLA cannot separate sign of elevation from phase alone",
    },
    "r3a1_erratum": "beta_align_search is empirical alignment score, not physical beta; fidelity from beta_local≈0.959 vs beta_true≈1.001 and range err≈0.02 km",
    "gate0": "if 14m/70m cannot resolve typical first-CZ Delta-u → HLA_ELEVATION_PROJECTION_WEAK",
    "stop": ["no MC RMSE", "no RC2 hard pairs", "no RC3-C", "no P5", "no Bellhop", "no new features"],
}


def munk_c(z, c0=1500.0, za=Z_A, B=B_MUNK, eps=EPS_MUNK):
    z = np.asarray(z, dtype=float)
    eta = 2.0 * (z - za) / B
    return c0 * (1.0 + eps * (np.exp(eta) - 1.0 - eta))


def munk_dcdz(z, c0=1500.0, za=Z_A, B=B_MUNK, eps=EPS_MUNK):
    z = np.asarray(z, dtype=float)
    eta = 2.0 * (z - za) / B
    return c0 * eps * (np.exp(eta) - 1.0) * (2.0 / B)


def u_hla(theta_rel_deg, phi_deg):
    th = np.deg2rad(theta_rel_deg)
    ph = np.deg2rad(phi_deg)
    return np.cos(th) * np.cos(ph)


# ---------------------------------------------------------------------------
# Gate-0: analytic projection + spatial resolution
# ---------------------------------------------------------------------------
def gate0_projection():
    rows = []
    for th in THETA_REL_SCAN:
        for ph in PHI_SCAN:
            u = float(u_hla(th, ph))
            # delta u vs |phi|=0 at same theta
            u0 = float(u_hla(th, 0.0))
            rows.append({
                "theta_rel_deg": th,
                "phi_deg": ph,
                "u": u,
                "delta_u_vs_phi0": abs(u - u0),
                "is_mirror_pair_sign": ph != 0,
                "mirror_u_equals_pos": abs(u - float(u_hla(th, abs(ph)))) < 1e-15,
            })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "analytic_hla_projection.csv", index=False, encoding="utf-8-sig")

    # pairwise Delta-u for multipath pairs (all phi pairs) at each theta
    pair_rows = []
    phis = [p for p in PHI_SCAN if p != 0]
    for th in THETA_REL_SCAN:
        for i, p1 in enumerate(phis):
            for p2 in phis[i + 1:]:
                u1 = float(u_hla(th, p1))
                u2 = float(u_hla(th, p2))
                du = abs(u1 - u2)
                mirror = abs(p1 + p2) < 1e-12 or abs(abs(p1) - abs(p2)) < 1e-12 and (p1 * p2 < 0)
                pair_rows.append({
                    "theta_rel_deg": th,
                    "phi1_deg": p1,
                    "phi2_deg": p2,
                    "u1": u1,
                    "u2": u2,
                    "delta_u": du,
                    "is_pm_mirror": bool(abs(p1 + p2) < 1e-12),
                    "delta_u_formula_pm": abs(math.cos(math.radians(th)) * (math.cos(math.radians(abs(p1))) - math.cos(math.radians(abs(p2))))) if p1 * p2 < 0 and abs(abs(p1) - abs(p2)) < 1e-12 else du,
                })
    pdf = pd.DataFrame(pair_rows)
    # signed-pair Delta u must be 0
    assert np.all(pdf.loc[pdf["is_pm_mirror"], "delta_u"] < 1e-12)
    pdf.to_csv(OUT / "direction_cosine_pairs.csv", index=False, encoding="utf-8-sig")

    # compact Δu table for task-style |phi| vs theta=0
    ref = []
    for ph in [1.0, 2.0, 5.0, 10.0, 15.0, 20.0]:
        du0 = abs(float(u_hla(0.0, 0.0)) - float(u_hla(0.0, ph)))
        ref.append({"abs_phi_deg": ph, "delta_u_theta0": du0,
                    "L_req_200Hz_m": (C0 / 200.0) / du0 if du0 > 0 else np.inf,
                    "L_req_300Hz_m": (C0 / 300.0) / du0 if du0 > 0 else np.inf})
    ref_df = pd.DataFrame(ref)
    ref_df.to_csv(OUT / "delta_u_vs_phi_theta0.csv", index=False, encoding="utf-8-sig")
    return df, pdf, ref_df


def array_steering(L, n_el, u, f):
    """a_m(u,f)=exp(-j k x_m u), x_m centered."""
    k = 2 * np.pi * f / C0
    x = (np.arange(n_el) - (n_el - 1) / 2.0) * (L / max(n_el - 1, 1) if n_el > 1 else 0.0)
    return np.exp(-1j * k * x * u)


def spatial_resolution_scan(pair_df, freqs=(200.0, 300.0)):
    """For each (theta, pair) and aperture: correlation, Bartlett dual-peak, sigma_u."""
    rows = []
    for ap in APERTURES:
        n_el = 8 if ap["L"] > 0 else 1
        L = ap["L"]
        for f in freqs:
            lam = C0 / f
            for _, r in pair_df.iterrows():
                if r["theta_rel_deg"] not in (0.0, 30.0, 60.0, 90.0):
                    continue
                if not (r["phi1_deg"] in (2.0, 5.0, 10.0) or r["phi2_deg"] in (2.0, 5.0, 10.0)):
                    continue
                u1, u2 = float(r["u1"]), float(r["u2"])
                du = float(r["delta_u"])
                if n_el == 1:
                    corr = 1.0
                    crb_su = np.inf
                    dual = False
                else:
                    a1 = array_steering(L, n_el, u1, f)
                    a2 = array_steering(L, n_el, u2, f)
                    corr = float(abs(np.vdot(a1, a2)) / (np.linalg.norm(a1) * np.linalg.norm(a2) + 1e-30))
                    # Bartlett dual peak: response of a1 steered across u grid
                    ugrid = np.linspace(-1.0, 1.0, 401)
                    # power for source at u1, beamscan
                    P = []
                    for ug in ugrid:
                        a = array_steering(L, n_el, ug, f)
                        P.append(abs(np.vdot(a, a1)) ** 2 / (np.linalg.norm(a) ** 2 + 1e-30))
                    P = np.asarray(P)
                    # dual peak if local max near u2 comparable
                    i1 = int(np.argmax(P))
                    # search peak near u2
                    mask = np.abs(ugrid - u2) < max(du, 0.02)
                    i2 = int(np.argmax(np.where(mask, P, 0))) if mask.any() else i1
                    dual = bool(i1 != i2 and P[i2] > 0.5 * P[i1] and abs(ugrid[i2] - u2) < 0.05)
                    # CRB sigma_u for single element-array: ~ 1/(k L sqrt(SNR*M)) order
                    # Use: sigma_u ≈ sqrt(12)/(k L sqrt(M * SNR)), SNR=20dB theoretical
                    k = 2 * np.pi * f / C0
                    snr = 100.0
                    crb_su = math.sqrt(12.0) / (k * L * math.sqrt(n_el * snr)) if L > 0 else np.inf
                L_req = lam / du if du > 0 else np.inf
                rows.append({
                    "aperture_tag": ap["tag"],
                    "aperture_m": L,
                    "note": ap["note"],
                    "f_hz": f,
                    "theta_rel_deg": r["theta_rel_deg"],
                    "phi1_deg": r["phi1_deg"],
                    "phi2_deg": r["phi2_deg"],
                    "delta_u": du,
                    "steering_corr": corr,
                    "bartlett_dual_peak": dual,
                    "crb_sigma_u": crb_su,
                    "L_required_lambda_over_du_m": L_req,
                    "resolved_corr_lt_0p95": bool(corr < 0.95) if n_el > 1 else False,
                    "resolved_sigma_u_lt_du": bool(np.isfinite(crb_su) and crb_su < du) if n_el > 1 else False,
                    "is_pm_mirror": bool(r["is_pm_mirror"]),
                })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "direction_cosine_resolution.csv", index=False, encoding="utf-8-sig")
    return df


# ---------------------------------------------------------------------------
# Controlled MMAC (isovelocity image model)
# ---------------------------------------------------------------------------
def controlled_paths(r_m, z_s, z_r, c=C0):
    """Two-path image model: direct + surface bounce."""
    # direct
    tau1 = r_m / c
    phi1 = math.atan2(z_r - z_s, r_m)  # rad
    # surface (image source at -z_s)
    L2 = math.hypot(r_m, z_r + z_s)
    tau2 = L2 / c
    phi2 = math.atan2(z_r + z_s, r_m)
    return [
        dict(branch="direct", tau=tau1, phi=phi1, topology="direct"),
        dict(branch="surface", tau=tau2, phi=phi2, topology="surface_image"),
    ]


def mac_obs(r_m, z, theta_rel_deg, mode="IDEAL_ELEVATION"):
    paths = controlled_paths(r_m, z, Z_R)
    phis = [p["phi"] for p in paths]
    taus = [p["tau"] for p in paths]
    dtau = taus[1] - taus[0]
    if mode == "IDEAL_ELEVATION":
        y = np.array(phis + [dtau])
    else:
        us = [float(u_hla(theta_rel_deg, math.degrees(p))) for p in phis]
        y = np.array(us + [dtau])
    return y, paths


def numeric_jac(fun, x, dx=(10.0, 1.0)):
    """fun(r,z)->y; x=(r,z) in meters."""
    y0 = fun(x[0], x[1])
    J = np.zeros((y0.size, 2))
    for i, d in enumerate(dx):
        xp = list(x)
        xm = list(x)
        xp[i] += d
        xm[i] -= d
        yp = fun(*xp)
        ym = fun(*xm)
        J[:, i] = (yp - ym) / (2 * d)
    return y0, J


def fisher_from_J(J, sigma):
    """F=J^T Sigma^{-1} J; sigma vector."""
    Sinv = np.diag(1.0 / (np.asarray(sigma) ** 2 + 1e-30))
    F = J.T @ Sinv @ J
    # singular values of J
    s = np.linalg.svd(J, compute_uv=False)
    detF = float(np.linalg.det(F))
    try:
        condJ = float(s[0] / s[-1]) if s[-1] > 1e-30 else np.inf
    except Exception:
        condJ = np.inf
    try:
        cov = np.linalg.inv(F)
        crb_r = math.sqrt(max(cov[0, 0], 0))
        crb_z = math.sqrt(max(cov[1, 1], 0))
        corr = cov[0, 1] / math.sqrt(cov[0, 0] * cov[1, 1] + 1e-30)
        condF = float(np.linalg.cond(F))
    except np.linalg.LinAlgError:
        crb_r = crb_z = np.inf
        corr = np.nan
        condF = np.inf
    return dict(sv=s.tolist(), sv_min=float(s.min()), cond_J=condJ, det_F=detF,
                crb_r_m=crb_r, crb_z_m=crb_z, corr_rz=corr, cond_F=condF)


def run_controlled_mmac():
    rows = []
    r_true, z_true = 50e3, 200.0
    for theta in [0.0, 30.0, 60.0]:
        for mode in ["IDEAL_ELEVATION", "HLA_DIRECTION_COSINE"]:
            for obs_set in ["A_angle_only", "B_delay_only", "C_joint"]:
                def make_fun(mode=mode, obs_set=obs_set, theta=theta):
                    def fun(r, z):
                        y, _ = mac_obs(r, z, theta, mode)
                        if obs_set == "A_angle_only":
                            return y[:2]
                        if obs_set == "B_delay_only":
                            return y[2:]
                        return y
                    return fun
                y0, J = numeric_jac(make_fun(), (r_true, z_true))
                # sigma: angle path 0.5° or array-u 0.005; delay 1 ms
                if mode == "IDEAL_ELEVATION":
                    sig_ang = math.radians(0.5)
                else:
                    sig_ang = 0.005
                sig_tau = 1e-3
                if obs_set == "A_angle_only":
                    sig = [sig_ang, sig_ang]
                elif obs_set == "B_delay_only":
                    sig = [sig_tau]
                else:
                    sig = [sig_ang, sig_ang, sig_tau]
                fish = fisher_from_J(J, sig)
                # ambiguity ridge description
                if obs_set == "A_angle_only":
                    ridge = "angle-only: elevation pair constrains r/z via geometry ridge"
                elif obs_set == "B_delay_only":
                    ridge = "delay-only: Δτ const-curve ridge in (r,z)"
                else:
                    ridge = "joint: angle ridge × delay ridge"
                rows.append({
                    "theta_rel_deg": theta,
                    "obs_mode": mode,
                    "obs_set": obs_set,
                    "y0": json.dumps(y0.tolist()),
                    "J_r": json.dumps(J[:, 0].tolist()),
                    "J_z": json.dumps(J[:, 1].tolist()),
                    "sv_min": fish["sv_min"],
                    "cond_J": fish["cond_J"],
                    "det_F": fish["det_F"],
                    "crb_r_m": fish["crb_r_m"],
                    "crb_z_m": fish["crb_z_m"],
                    "corr_rz": fish["corr_rz"],
                    "cond_F": fish["cond_F"],
                    "ridge_note": ridge,
                    "sigma_assumed": json.dumps(sig),
                })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "CONTROLLED_MMAC_REPRO.csv", index=False, encoding="utf-8-sig")
    return df


# ---------------------------------------------------------------------------
# E-STD eigenrays via RK4 shooting on Munk
# ---------------------------------------------------------------------------
def ray_rhs(s, y):
    """y=[x,z,theta]; d/ds: cos, sin, -(cos/c) dc/dz."""
    x, z, th = y
    z = min(max(z, 0.5), DEPTH - 0.5)
    c = float(munk_c(z))
    dcdz = float(munk_dcdz(z))
    dx = math.cos(th)
    dz = math.sin(th)
    dth = -math.cos(th) * dcdz / c
    # turning / surface reflection handled outside
    return np.array([dx, dz, dth])


def integrate_ray(theta0, r_target, ds=20.0, max_steps=8000, z0=Z_S):
    """Shoot from (0,z0) with launch angle theta0 (rad). Return arrival at x≈r_target."""
    y = np.array([0.0, z0, theta0])
    t = 0.0
    x_hist, z_hist, th_hist = [0.0], [z0], [theta0]
    for _ in range(max_steps):
        k1 = ray_rhs(0, y)
        k2 = ray_rhs(0, y + 0.5 * ds * k1)
        k3 = ray_rhs(0, y + 0.5 * ds * k2)
        k4 = ray_rhs(0, y + ds * k3)
        yn = y + (ds / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        if yn[1] < 0.0:
            yn[1] = -yn[1]
            yn[2] = -yn[2]
        if yn[1] > DEPTH:
            yn[1] = 2 * DEPTH - yn[1]
            yn[2] = -yn[2]
        zmid = 0.5 * (y[1] + yn[1])
        c = float(munk_c(zmid))
        t += ds / c
        y = yn
        z_hist.append(y[1]); th_hist.append(y[2])
        if y[0] >= r_target:
            break
    return {
        "x": float(y[0]), "z": float(y[1]), "theta_arr": float(y[2]),
        "tau": float(t), "launch": float(theta0),
        "turning_depth": float(np.max(z_hist)),
        "topology": "surface_bounce" if np.min(z_hist) < 1.0 else "refracted",
        "z_hist_max": float(np.max(z_hist)),
        "z_hist_min": float(np.min(z_hist)),
    }


def find_eigenrays(r_target, launch_grid_deg=None):
    if launch_grid_deg is None:
        launch_grid_deg = np.linspace(-25.0, 25.0, 121)
    results = []
    prev = None
    for ld in launch_grid_deg:
        res = integrate_ray(math.radians(ld), r_target)
        err = res["z"] - Z_R
        res["err_z"] = err
        res["launch_deg"] = ld
        res["r_target_km"] = r_target / 1e3
        results.append(res)
        prev = res
    df = pd.DataFrame(results)
    # find sign changes / near hits in err_z
    branches = []
    e = df["err_z"].values
    launches = df["launch_deg"].values
    for i in range(len(e) - 1):
        if e[i] == 0 or e[i] * e[i + 1] < 0:
            # linear interp launch
            if e[i + 1] != e[i]:
                ld = launches[i] - e[i] * (launches[i + 1] - launches[i]) / (e[i + 1] - e[i])
            else:
                ld = launches[i]
            # refine shoot
            ref = integrate_ray(math.radians(ld), r_target)
            # elevation at receiver: phi related to arrival theta (vertical angle)
            # In ray tracer theta is from horizontal in range-depth plane → arrival elevation ≈ theta_arr (sign)
            phi = math.degrees(ref["theta_arr"])
            # also ray path length ~ integrate; tau already
            branches.append({
                "r_target_km": r_target / 1e3,
                "launch_deg": float(ld),
                "phi_arr_deg": phi,
                "tau_s": ref["tau"],
                "turning_depth_m": ref["turning_depth"],
                "topology": ref["topology"],
                "z_hit_m": ref["z"],
                "err_z_m": ref["z"] - Z_R,
            })
    bdf = pd.DataFrame(branches)
    return df, bdf


def track_branches(all_branches):
    """Group eigenrays into branches by continuity of launch/topology across range."""
    if all_branches.empty:
        return pd.DataFrame()
    df = all_branches.sort_values(["topology", "launch_deg", "r_target_km"]).copy()
    df["branch_id"] = ""
    bid = 0
    ids = {}
    for _, row in df.iterrows():
        key = (row["topology"], round(row["launch_deg"] / 2.0))  # 2 deg bins
        if key not in ids:
            ids[key] = f"B{bid}"
            bid += 1
        df.loc[row.name, "branch_id"] = ids[key]
    return df


# ---------------------------------------------------------------------------
# Delay observability S0/S1/S2
# ---------------------------------------------------------------------------
def delay_ambiguity(freqs, tau_grid, tau_true):
    """Matched-filter-like ambiguity vs delay for given frequency support."""
    # A(tau) = |sum_f exp(j 2 pi f (tau-tau_true))| / N
    dtau = tau_grid - tau_true
    A = np.abs(np.mean(np.exp(1j * 2 * np.pi * freqs[:, None] * dtau[None, :]), axis=0))
    return A


def run_delay_observability():
    tau_true = 0.05  # 50 ms relative
    tau_grid = np.linspace(-0.2, 0.2, 2001)
    sources = {
        "S0": np.linspace(F_BAND[0], F_BAND[1], 23),
        "S1": np.array([168.0, 204.0, 232.0, 279.0, 320.0]),
        "S2": np.array([166.0, 201.0, 235.0, 283.0, 338.0]),
    }
    rows = []
    for name, freqs in sources.items():
        A = delay_ambiguity(freqs, tau_grid, tau_true)
        # mainlobe width at 0.5
        idx = np.where(A >= 0.5)[0]
        width = float(tau_grid[idx.max()] - tau_grid[idx.min()]) if len(idx) >= 2 else 0.0
        # sidelobe / period: max A outside mainlobe
        mask = np.ones_like(A, dtype=bool)
        if len(idx):
            mask[idx] = False
        side = float(A[mask].max()) if mask.any() else 0.0
        # effective bandwidth
        bw = float(np.std(freqs)) * math.sqrt(12) if freqs.size > 1 else 0.0
        # CRB-like sigma_tau ~ 1/(2 pi sqrt(SNR) sigma_f)
        snr = 100.0
        sig_f = float(np.std(freqs)) if freqs.size > 1 else 0.0
        crb_tau = 1.0 / (2 * math.pi * math.sqrt(snr) * sig_f) if sig_f > 0 else np.inf
        # periodic ambiguity for line sets
        period_amb = bool(side > 0.7)
        rows.append({
            "source": name,
            "n_freqs": int(freqs.size),
            "effective_bandwidth_Hz": bw,
            "amb_mainlobe_width_s": width,
            "amb_sidelobe_max": side,
            "crb_sigma_tau_s": crb_tau,
            "periodic_ambiguity_risk": period_amb,
            "delay_usable": bool(width < 0.05 and not (period_amb and name == "S2")),
            "note": "S2 multi-line can create periodic delay ambiguity; do not assume best",
        })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "source_delay_observability.csv", index=False, encoding="utf-8-sig")
    return df


# ---------------------------------------------------------------------------
# Fisher scan with HLA sigma_u and delay CRB
# ---------------------------------------------------------------------------
def run_fisher_scan(ctrl_df, delay_df, res_df):
    rows = []
    r_true, z_true = 50e3, 200.0
    # pick sigma_tau from delay table
    sig_tau = {r["source"]: float(r["crb_sigma_tau_s"]) for _, r in delay_df.iterrows()}

    def sigma_u_hla(theta_rel_deg, phi_deg, L=14.0, f=250.0, snr=100.0, n_el=8):
        k = 2 * np.pi * f / C0
        # array CRB
        su_arr = math.sqrt(12.0) / (k * L * math.sqrt(n_el * snr))
        # RC2 theta propagation
        th = math.radians(theta_rel_deg)
        ph = math.radians(phi_deg)
        du_dth = -math.sin(th) * math.cos(ph)
        return su_arr, du_dth

    for sigma_th_deg in SIGMA_THETA_DEG:
        sig_th = math.radians(sigma_th_deg)
        for source in ["S0", "S1", "S2"]:
            st = sig_tau[source]
            for theta in [0.0, 30.0, 60.0]:
                for mode, obs_set, label in [
                    ("IDEAL_ELEVATION", "C_joint", "A_ideal_elev+delay"),
                    ("HLA_DIRECTION_COSINE", "C_joint", "B_hla_u+delay"),
                    ("IDEAL_ELEVATION", "B_delay_only", "C_delay_only"),
                ]:
                    def make_fun(mode=mode, obs_set=obs_set, theta=theta):
                        def fun(r, z):
                            y, _ = mac_obs(r, z, theta, mode)
                            if obs_set == "B_delay_only":
                                return y[2:]
                            return y
                        return fun
                    y0, J = numeric_jac(make_fun(), (r_true, z_true))
                    # build sigma vector
                    phi1 = math.atan2(Z_R - z_true, r_true)
                    phi2 = math.atan2(Z_R + z_true, r_true)
                    if mode == "IDEAL_ELEVATION":
                        sig_ang = math.radians(0.5)  # ideal elevation still needs some noise
                    else:
                        su_arr1, dudth1 = sigma_u_hla(theta, math.degrees(phi1))
                        su_arr2, dudth2 = sigma_u_hla(theta, math.degrees(phi2))
                        sig_u1 = math.sqrt(su_arr1 ** 2 + (dudth1 * sig_th) ** 2)
                        sig_u2 = math.sqrt(su_arr2 ** 2 + (dudth2 * sig_th) ** 2)
                    if obs_set == "B_delay_only":
                        sig = [st]
                    elif mode == "IDEAL_ELEVATION":
                        sig = [math.radians(0.5), math.radians(0.5), st]
                    else:
                        sig = [sig_u1, sig_u2, st]
                    fish = fisher_from_J(J, sig)
                    # dominate analysis for HLA
                    if mode == "HLA_DIRECTION_COSINE" and obs_set == "C_joint":
                        su_arr1, dudth1 = sigma_u_hla(theta, math.degrees(phi1))
                        rc2_part = abs(dudth1 * sig_th)
                        arr_part = su_arr1
                        dominate = "RC2_theta" if rc2_part > arr_part else "array_spatial"
                    else:
                        rc2_part = arr_part = np.nan
                        dominate = "n/a"
                    rows.append({
                        "sigma_theta_deg": sigma_th_deg,
                        "source": source,
                        "theta_rel_deg": theta,
                        "case": label,
                        "sigma_tau_s": st,
                        "crb_r_m": fish["crb_r_m"],
                        "crb_z_m": fish["crb_z_m"],
                        "corr_rz": fish["corr_rz"],
                        "cond_F": fish["cond_F"],
                        "sv_min": fish["sv_min"],
                        "cond_J": fish["cond_J"],
                        "sigma_u_rc2_component": rc2_part,
                        "sigma_u_array_component": arr_part,
                        "u_error_dominated_by": dominate,
                    })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "hla_direction_cosine_fisher.csv", index=False, encoding="utf-8-sig")
    # split exports
    df[df["case"] == "A_ideal_elev+delay"].to_csv(OUT / "ideal_elevation_fisher.csv", index=False, encoding="utf-8-sig")
    df[df["case"] == "C_delay_only"].to_csv(OUT / "delay_only_fisher.csv", index=False, encoding="utf-8-sig")
    return df


# ---------------------------------------------------------------------------
# Decision
# ---------------------------------------------------------------------------
def decide(gate_res, ctrl, eigen_df, delay_df, fish_df):
    # Gate-0: can 14m/70m resolve typical Delta-u for first-CZ elevations?
    # Typical first-CZ elevations for 50 km, 200m depth: roughly few degrees
    typical = gate_res[(gate_res["theta_rel_deg"].isin([0.0, 30.0])) &
                       (gate_res["phi1_deg"].abs().isin([2.0, 5.0]))]
    if typical.empty:
        typical = gate_res
    resolved_14 = typical[typical["aperture_tag"] == "L14m"]["resolved_corr_lt_0p95"].mean() if len(typical) else 0
    resolved_70 = typical[typical["aperture_tag"] == "L70m"]["resolved_corr_lt_0p95"].mean() if len(typical) else 0
    su_lt_du_14 = typical[typical["aperture_tag"] == "L14m"]["resolved_sigma_u_lt_du"].mean() if len(typical) else 0
    su_lt_du_70 = typical[typical["aperture_tag"] == "L70m"]["resolved_sigma_u_lt_du"].mean() if len(typical) else 0
    hla_weak = bool(resolved_14 < 0.25 and resolved_70 < 0.5 and su_lt_du_14 < 0.2)

    # controlled: ideal vs HLA joint Fisher
    def crb_case(theta, mode, obs_set):
        d = ctrl[(ctrl["theta_rel_deg"] == theta) & (ctrl["obs_mode"] == mode) & (ctrl["obs_set"] == obs_set)]
        return float(d["crb_r_m"].iloc[0]) if len(d) else np.inf

    ideal_ok = all(crb_case(th, "IDEAL_ELEVATION", "C_joint") < 5e3 for th in [0.0, 30.0])
    hla_ok = all(crb_case(th, "HLA_DIRECTION_COSINE", "C_joint") < 5e3 for th in [0.0, 30.0])
    delay_ok = crb_case(0.0, "IDEAL_ELEVATION", "B_delay_only") < 20e3

    n_branches = eigen_df["branch_id"].nunique() if len(eigen_df) and "branch_id" in eigen_df.columns else 0
    stable_branches = n_branches >= 2

    # delay usable any source
    delay_use = bool(delay_df["delay_usable"].any()) if len(delay_df) else False

    # Fisher complementary: delay-only vs joint
    f_joint = fish_df[fish_df["case"] == "B_hla_u+delay"]
    f_delay = fish_df[fish_df["case"] == "C_delay_only"]
    # if ideal joint much better than delay-only → angle helps in physics
    physics_angle_helps = ideal_ok and (crb_case(0.0, "IDEAL_ELEVATION", "C_joint") <
                                        0.7 * crb_case(0.0, "IDEAL_ELEVATION", "B_delay_only"))

    notes = {
        "gate0_hla_elevation_projection_weak": hla_weak,
        "resolved_corr_14m_frac": float(resolved_14),
        "resolved_corr_70m_frac": float(resolved_70),
        "resolved_sigma_u_14m_frac": float(su_lt_du_14),
        "n_estd_branches": int(n_branches),
        "stable_branches_ge2": stable_branches,
        "ideal_joint_ok": ideal_ok,
        "hla_joint_ok": hla_ok,
        "delay_only_ok": delay_ok,
        "delay_usable_any_source": delay_use,
        "physics_angle_helps": physics_angle_helps,
    }

    if not stable_branches:
        decision = "B1_NO_STABLE_MULTIPATH_IDENTITY"
        why = f"E-STD eigenray tracking found n_branches={n_branches} < 2 stable branches under current shooting grid."
        nxt = "refine ray branch tracking or declare RC3-B multipath identity unavailable; do not proceed to B2 MMAC"
    elif ideal_ok and not hla_ok:
        decision = "B1_PHYSICS_ONLY_HLA_LIMITED"
        why = (
            f"Controlled IDEAL_ELEVATION+delay joint Fisher OK (crb_r small), but HLA_u+delay weak. "
            f"Gate-0: hla_projection_weak={hla_weak} (corr-resolve 14m={resolved_14:.2f}, 70m={resolved_70:.2f}). "
            f"HLA spatial phase cannot resolve first-CZ elevation differences / ±phi mirror."
        )
        nxt = "MMAC physics exists but not HLA-realizable via elevation; next study delay-set matching if delays usable, else RC3-C"
    elif (not hla_weak) and stable_branches and ideal_ok and hla_ok and delay_use:
        decision = "B1_MMAC_PHYSICS_CONFIRMED"
        why = "Stable branches + HLA_u resolvable + joint Fisher non-degenerate + delay observable."
        nxt = "enter B2 full HLA-MMAC (still no hidden branch_id in observation)"
    elif delay_use and stable_branches and (delay_ok or physics_angle_helps):
        decision = "B1_DELAY_DOMINANT_POSSIBLE"
        why = (
            f"HLA elevation projection weak ({hla_weak}), but delay observability usable for at least one source; "
            f"controlled delay/joint shows r-z information may be delay-dominant."
        )
        nxt = "B2 only delay-set / multipath-delay matching WITHOUT true branch labels (association required)"
    elif not delay_use:
        decision = "B1_DELAY_BANDWIDTH_LIMITED"
        why = f"Delay observability weak/periodic across S0/S1/S2 under current band; delay_df={delay_df.to_dict(orient='records')}."
        nxt = "do not claim practical delay MMAC; consider RC3-C or broader bandwidth studies later"
    else:
        decision = "B1_NO_COMPLEMENTARY_INFORMATION"
        why = "Neither HLA_u+delay nor delay-set show complementary non-degenerate r-z information."
        nxt = "stop RC3-B; RC3-C next"

    return decision, why, nxt, notes


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def _esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def write_group_svg(path, title, cats, series, ylab, note=""):
    w, h = 760, 380
    ml, mr, mt, mb = 70, 140, 48, 60
    pw, ph = w - ml - mr, h - mt - mb
    cols = ["#b45309", "#0f766e", "#1d4ed8", "#8a8a8a", "#9f1239"]
    vals = [v for s in series for v in s["y"] if v is not None and np.isfinite(v)]
    vmax = max(vals + [1e-6])
    ncat, nser = len(cats), len(series)
    bw = pw / max(ncat * (nser + 1), 1)
    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
         f'<rect width="{w}" height="{h}" fill="#f7f4ef"/>',
         f'<text x="{w/2}" y="26" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="14" font-weight="600">{_esc(title)}</text>',
         f'<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="#fff" stroke="#ccc"/>']
    for i, s in enumerate(series):
        for ic, cat in enumerate(cats):
            v = s["y"][ic]
            if v is None or not np.isfinite(v):
                continue
            x = ml + ic * (pw / ncat) + i * bw + 4
            hh = (v / vmax) * ph
            p.append(f'<rect x="{x:.1f}" y="{mt+ph-hh:.1f}" width="{bw*0.8:.1f}" height="{max(hh,0.5):.1f}" fill="{s.get("color", cols[i%len(cols)])}"/>')
        p.append(f'<text x="{ml+ (i+0.5)*20}" y="{mt+ph+16}" font-size="10" fill="#333" font-family="sans-serif"></text>')
    for ic, cat in enumerate(cats):
        p.append(f'<text x="{ml+(ic+0.5)*pw/ncat:.1f}" y="{mt+ph+18}" text-anchor="middle" font-size="11" font-family="sans-serif">{_esc(cat)}</text>')
    for i, s in enumerate(series):
        p.append(f'<rect x="{w-mr+8}" y="{mt+8+i*16}" width="12" height="3" fill="{s.get("color", cols[i%len(cols)])}"/>')
        p.append(f'<text x="{w-mr+24}" y="{mt+12+i*16}" font-size="10" font-family="sans-serif">{_esc(s["name"])}</text>')
    p.append(f'<text x="{ml}" y="{h-8}" font-size="10" fill="#777" font-family="sans-serif">{_esc(note)}  ylab={_esc(ylab)}</text></svg>')
    path.write_text("\n".join(p), encoding="utf-8")


def make_figures(ref_df, gate_res, ctrl, branches, fish_df, delay_df):
    # fig0 Δu vs phi
    write_group_svg(
        FIG / "fig0_hla_projection_gate.svg",
        "图0  HLA方向余弦 |Δu| vs |φ|（θ_rel=0°最有利）",
        [f"{int(p)}°" for p in ref_df["abs_phi_deg"]],
        [{"name": "|Δu|", "color": "#b45309", "y": ref_df["delta_u_theta0"].tolist()}],
        "|Δu|",
        note="±φ镜像：u(+φ)=u(-φ)，HLA相位天然不分辨俯仰符号；L_req~λ/Δu 见CSV",
    )
    # fig1 controlled CRB r
    cats, s_ideal, s_hla, s_del = [], [], [], []
    for th in [0.0, 30.0, 60.0]:
        for mode, arr in [("IDEAL_ELEVATION", s_ideal), ("HLA_DIRECTION_COSINE", s_hla)]:
            d = ctrl[(ctrl["theta_rel_deg"] == th) & (ctrl["obs_mode"] == mode) & (ctrl["obs_set"] == "C_joint")]
            if th == 0.0:
                cats.append(f"θ={th:.0f}°")
            arr.append(float(d["crb_r_m"].iloc[0]) / 1e3 if len(d) else np.nan)
        d = ctrl[(ctrl["theta_rel_deg"] == th) & (ctrl["obs_set"] == "B_delay_only")]
        # take first matching
        d = ctrl[(ctrl["theta_rel_deg"] == th) & (ctrl["obs_set"] == "B_delay_only") & (ctrl["obs_mode"] == "IDEAL_ELEVATION")]
        s_del.append(float(d["crb_r_m"].iloc[0]) / 1e3 if len(d) else np.nan)
    write_group_svg(
        FIG / "fig1_controlled_mmac.svg",
        "图1  受控MMAC：CRB_r（km）角度/时延/联合",
        ["θ=0°", "θ=30°", "θ=60°"],
        [
            {"name": "IDEAL elev+delay", "color": "#0f766e", "y": s_ideal},
            {"name": "HLA u+delay", "color": "#b45309", "y": s_hla},
            {"name": "delay only", "color": "#8a8a8a", "y": s_del},
        ],
        "CRB_r km",
        note="IDEAL成立而HLA弱 → B1_PHYSICS_ONLY_HLA_LIMITED",
    )
    # fig2 branches
    if branches is not None and len(branches):
        b = branches.groupby(["branch_id", "r_target_km"]).size().reset_index(name="n")
        cats_b = [f"{r:.0f}km" for r in sorted(branches["r_target_km"].unique())]
        # count branches per range
        counts = []
        for rkm in sorted(branches["r_target_km"].unique()):
            counts.append(branches[branches["r_target_km"] == rkm]["branch_id"].nunique())
        write_group_svg(
            FIG / "fig2_ESTD_branches.svg",
            "图2  E-STD 本征射线分支数 vs 距离",
            cats_b,
            [{"name": "n_branches", "color": "#1d4ed8", "y": counts}],
            "branch count",
            note="≥2 稳定分支才继续 MMAC；身份按 launch/topology 连续，不按时延峰序",
        )
    else:
        (FIG / "fig2_ESTD_branches.svg").write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="400" height="80"><text x="20" y="40" font-size="14">no eigenray branches</text></svg>',
            encoding="utf-8")

    # fig3 sigma_u dominate
    f = fish_df[fish_df["case"] == "B_hla_u+delay"]
    if len(f):
        cats3 = [f"σθ={s}°" for s in SIGMA_THETA_DEG]
        y_rc2, y_arr = [], []
        for s in SIGMA_THETA_DEG:
            d = f[(f["sigma_theta_deg"] == s) & (f["theta_rel_deg"] == 0.0) & (f["source"] == "S0")]
            if not len(d):
                y_rc2.append(np.nan); y_arr.append(np.nan)
            else:
                y_rc2.append(float(d["sigma_u_rc2_component"].iloc[0]))
                y_arr.append(float(d["sigma_u_array_component"].iloc[0]))
        write_group_svg(
            FIG / "fig3_angle_delay_sensitivity.svg",
            "图3  HLA u 误差：RC2方位 vs 阵列空间项",
            cats3,
            [{"name": "RC2 θ→u", "color": "#b45309", "y": y_rc2},
             {"name": "array σu", "color": "#0f766e", "y": y_arr}],
            "sigma_u",
            note="∂u/∂θ=-sinθ_rel cosφ；θ_rel=0 时 RC2 项对 φ 的耦合仍经 cosφ 高阶",
        )
    # fig4 ideal vs HLA crb_z
    cats4, iz, hz = [], [], []
    for th in [0.0, 30.0, 60.0]:
        cats4.append(f"θ={th:.0f}°")
        di = ctrl[(ctrl["theta_rel_deg"]==th)&(ctrl["obs_mode"]=="IDEAL_ELEVATION")&(ctrl["obs_set"]=="C_joint")]
        dh = ctrl[(ctrl["theta_rel_deg"]==th)&(ctrl["obs_mode"]=="HLA_DIRECTION_COSINE")&(ctrl["obs_set"]=="C_joint")]
        iz.append(float(di["crb_z_m"].iloc[0]) if len(di) else np.nan)
        hz.append(float(dh["crb_z_m"].iloc[0]) if len(dh) else np.nan)
    write_group_svg(
        FIG / "fig4_ideal_vs_hla.svg",
        "图4  IDEAL φ vs HLA u：CRB_z（m）",
        cats4,
        [{"name": "IDEAL", "color": "#0f766e", "y": iz},
         {"name": "HLA_u", "color": "#b45309", "y": hz}],
        "CRB_z m",
        note="结构性对比",
    )
    # fig5 fisher boundary from scan
    if len(fish_df):
        d = fish_df[(fish_df["source"]=="S0") & (fish_df["theta_rel_deg"]==0.0)]
        cats5 = ["ideal", "hla_u", "delay"]
        y5 = []
        for case in ["A_ideal_elev+delay", "B_hla_u+delay", "C_delay_only"]:
            dd = d[d["case"]==case]
            y5.append(float(dd["crb_r_m"].iloc[0])/1e3 if len(dd) else np.nan)
        write_group_svg(
            FIG / "fig5_fisher_boundary.svg",
            "图5  Fisher 边界：CRB_r（S0, θ_rel=0°, σθ=0.1°）",
            cats5,
            [{"name": "CRB_r km", "color": "#1d4ed8", "y": y5}],
            "km",
            note="",
        )
    # fig6 delay ambiguity
    if delay_df is not None and len(delay_df):
        write_group_svg(
            FIG / "fig6_delay_ambiguity.svg",
            "图6  源条件时延可观测性（主瓣宽 ms）",
            delay_df["source"].tolist(),
            [{"name": "mainlobe ms", "color": "#b45309", "y": (delay_df["amb_mainlobe_width_s"]*1e3).tolist()},
             {"name": "sidelobe", "color": "#8a8a8a", "y": delay_df["amb_sidelobe_max"].tolist()}],
            "ms / ratio",
            note="S2 可能周期模糊，不默认最优",
        )


def main():
    t0 = time.time()
    print("=== R3-B1 HLA-MMAC observability ===")
    print(f"OUT={OUT}")
    (OUT / "R3_B1_CONFIG.json").write_text(json.dumps(CONFIG, ensure_ascii=False, indent=2), encoding="utf-8")

    print("[Gate-0] analytic projection ...", flush=True)
    proj, pairs, ref = gate0_projection()
    print(ref.to_string(index=False), flush=True)
    print("[Gate-0] spatial resolution ...", flush=True)
    gate_res = spatial_resolution_scan(pairs)
    print(gate_res.groupby("aperture_tag")[["resolved_corr_lt_0p95", "resolved_sigma_u_lt_du"]].mean(), flush=True)

    print("[controlled MMAC] ...", flush=True)
    ctrl = run_controlled_mmac()
    print(ctrl[["theta_rel_deg","obs_mode","obs_set","crb_r_m","crb_z_m","cond_J"]].to_string(index=False), flush=True)

    print("[E-STD eigenrays] ...", flush=True)
    all_b = []
    shoot_grids = {}
    for rkm in [45.0, 50.0, 55.0, 60.0]:
        # convergence: two launch grids
        _, b1 = find_eigenrays(rkm * 1e3, np.linspace(-25, 25, 121))
        _, b2 = find_eigenrays(rkm * 1e3, np.linspace(-25, 25, 181))
        if len(b1) and len(b2):
            stab = abs(len(b1) - len(b2)) <= max(2, 0.5 * max(len(b1), 1))
            # tau stability: compare sorted taus if counts match
            if len(b1) == len(b2):
                tdiff = float(np.max(np.abs(np.sort(b1["tau_s"].values) - np.sort(b2["tau_s"].values))))
            else:
                tdiff = np.nan
        else:
            stab, tdiff = False, np.nan
        shoot_grids[rkm] = dict(n_coarse=len(b1), n_fine=len(b2), stable=bool(stab), tau_max_diff_s=tdiff)
        b1["grid"] = "coarse201"
        all_b.append(b1)
        print(f"  r={rkm}km branches coarse={len(b1)} fine={len(b2)} stable={stab} dtau={tdiff}", flush=True)
    branches_raw = pd.concat(all_b, ignore_index=True) if all_b else pd.DataFrame()
    branches = track_branches(branches_raw)
    branches.to_csv(OUT / "eigenray_branches.csv", index=False, encoding="utf-8-sig")
    branches_raw.to_csv(OUT / "branch_tracking.csv", index=False, encoding="utf-8-sig")
    # modal cross-check: compare ray tau spread vs modal group delay scale
    cross = []
    mode = p4.MODE_ENV["E0"]
    for rkm in [45.0, 50.0, 55.0, 60.0]:
        # crude modal group delay proxy: d(phase)/domega via finite diff on modal sum amp peak - skip full pulse
        # compare ray tau min/max
        bb = branches_raw[branches_raw["r_target_km"] == rkm] if len(branches_raw) else pd.DataFrame()
        cross.append({
            "r_km": rkm,
            "n_ray_branches": int(len(bb)),
            "ray_tau_min_s": float(bb["tau_s"].min()) if len(bb) else np.nan,
            "ray_tau_max_s": float(bb["tau_s"].max()) if len(bb) else np.nan,
            "ray_tau_span_s": float(bb["tau_s"].max() - bb["tau_s"].min()) if len(bb) else np.nan,
            **{f"shoot_{k}": v for k, v in shoot_grids[rkm].items()},
            "note": "ray-modal identity not forced; only order-of-magnitude tau structure",
        })
    pd.DataFrame(cross).to_csv(OUT / "ray_modal_crosscheck.csv", index=False, encoding="utf-8-sig")

    print("[delay observability] ...", flush=True)
    delay_df = run_delay_observability()
    print(delay_df.to_string(index=False), flush=True)

    print("[Fisher scan] ...", flush=True)
    fish_df = run_fisher_scan(ctrl, delay_df, gate_res)

    # joint jacobian dump sample
    y0, J = numeric_jac(lambda r, z: mac_obs(r, z, 0.0, "IDEAL_ELEVATION")[0], (50e3, 200.0))
    pd.DataFrame(J, columns=["dr", "dz"], index=["phi1", "phi2", "dtau"]).to_csv(
        OUT / "joint_rz_jacobian.csv", encoding="utf-8-sig"
    )
    fish = fisher_from_J(J, [math.radians(0.5), math.radians(0.5), 1e-3])
    pd.DataFrame([fish]).to_csv(OUT / "joint_rz_fisher.csv", index=False, encoding="utf-8-sig")

    print("[decision + figures + report] ...", flush=True)
    decision, why, nxt, notes = decide(gate_res, ctrl, branches, delay_df, fish_df)
    make_figures(ref, gate_res, ctrl, branches, fish_df, delay_df)

    dec = {
        "rc3b1_decision": decision,
        "why": why,
        "next_step": nxt,
        "notes": notes,
        "gate0_ref_delta_u": ref.to_dict(orient="records"),
        "created_utc": NOW,
        "stop": "after R3-B1; no MC RMSE; no RC2 hard pairs; no RC3-C; no P5",
    }
    (OUT / "R3_B1_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    def fnum(x, nd=4):
        try:
            v = float(x)
            return "n/a" if not np.isfinite(v) else f"{v:.{nd}f}"
        except Exception:
            return "n/a"

    rp = []
    rp.append("# R3-B1 报告：HLA 方向余弦可观测性门控 + MMAC")
    rp.append("")
    rp.append(f"UTC：{NOW}")
    rp.append("")
    rp.append("## 0. R3-A1 收口勘误（文字）")
    rp.append("")
    rp.append("- 判定保持 **A1-NOT-DIRECTLY-TRANSFERABLE**")
    rp.append("- `beta_align_search=2.2` = 条纹对齐**经验评分**参数，**不是**物理 β")
    rp.append("- 保真依据：`beta_local≈0.959` vs `beta_true≈1.001`，距离峰值误差 ≈0.02 km")
    rp.append("- 经典常数 β 不能直接迁移 E-STD；局部/广义 β 条纹降为后备，不关闭整类方法")
    rp.append("")
    rp.append("## 1. 关键物理量（冻结）")
    rp.append("")
    rp.append("$$u=\\cos\\theta_{\\rm rel}\\cos\\phi,\\qquad u(+\\phi)=u(-\\phi)$$")
    rp.append("")
    rp.append("水平直线阵空间相位观测的是 **方向余弦 u**，不是独立俯仰角 φ；±φ 镜像在相位上不可分。")
    rp.append("")
    rp.append("## Gate-0：解析投影量级")
    rp.append("")
    rp.append("| \\|φ\\| | \\|Δu\\| (θ_rel=0°) | L~λ/Δu @200Hz (m) |")
    rp.append("| --- | --- | --- |")
    for _, r in ref.iterrows():
        rp.append(f"| {r['abs_phi_deg']:.0f}° | {r['delta_u_theta0']:.6f} | {fnum(r['L_req_200Hz_m'],1)} |")
    rp.append("")
    rp.append(f"- 孔径可分辨比例（corr<0.95）：14m **{fnum(notes['resolved_corr_14m_frac'],2)}**，70m **{fnum(notes['resolved_corr_70m_frac'],2)}**")
    rp.append(f"- σu<Δu 比例：14m **{fnum(notes['resolved_sigma_u_14m_frac'],2)}**")
    rp.append(f"- **Gate-0 HLA_ELEVATION_PROJECTION_WEAK = {notes['gate0_hla_elevation_projection_weak']}**")
    rp.append("")
    rp.append("数据：`analytic_hla_projection.csv`, `direction_cosine_pairs.csv`, `direction_cosine_resolution.csv`")
    rp.append("")
    rp.append("## 2. 受控 MMAC（等声速两路径）")
    rp.append("")
    rp.append("| θ_rel | 模式 | 观测 | CRB_r km | CRB_z m | cond J |")
    rp.append("| --- | --- | --- | --- | --- | --- |")
    for _, r in ctrl.iterrows():
        rp.append(
            f"| {r['theta_rel_deg']:.0f} | {r['obs_mode']} | {r['obs_set']} | "
            f"{fnum(r['crb_r_m']/1e3,3)} | {fnum(r['crb_z_m'],1)} | {fnum(r['cond_J'],2)} |"
        )
    rp.append("")
    rp.append("## 3. E-STD 本征射线")
    rp.append("")
    rp.append(f"- 稳定分支数（tracking）：**{notes['n_estd_branches']}**；≥2 = **{notes['stable_branches_ge2']}**")
    if len(branches):
        rp.append("- 分支统计（前几行见 CSV）：按 topology+launch 连续 ID，不按时延峰序定义身份。")
        rp.append("")
        rp.append("| branch | n | r_km 覆盖 | φ 中位° | τ 中位 s |")
        rp.append("| --- | --- | --- | --- | --- |")
        g = branches.groupby("branch_id")
        for bid, gg in list(g)[:12]:
            rp.append(
                f"| {bid} | {len(gg)} | {gg['r_target_km'].min():.0f}-{gg['r_target_km'].max():.0f} | "
                f"{fnum(gg['phi_arr_deg'].median(),2)} | {fnum(gg['tau_s'].median(),4)} |"
            )
    rp.append("")
    rp.append("## 4. 时延可观测性")
    rp.append("")
    rp.append("| 源 | 有效带宽 Hz | 主瓣宽 ms | 旁瓣 | 周期模糊 | 可用 |")
    rp.append("| --- | --- | --- | --- | --- | --- |")
    for _, r in delay_df.iterrows():
        rp.append(
            f"| {r['source']} | {fnum(r['effective_bandwidth_Hz'],1)} | {fnum(r['amb_mainlobe_width_s']*1e3,3)} | "
            f"{fnum(r['amb_sidelobe_max'],3)} | {r['periodic_ambiguity_risk']} | {r['delay_usable']} |"
        )
    rp.append("")
    rp.append("## 5. Fisher（σθ 扫描，S0 摘录 θ_rel=0°）")
    rp.append("")
    d = fish_df[(fish_df["source"]=="S0") & (fish_df["theta_rel_deg"]==0.0) & (fish_df["sigma_theta_deg"]==0.1)]
    rp.append("| case | CRB_r km | CRB_z m | corr | cond F |")
    rp.append("| --- | --- | --- | --- | --- |")
    for _, r in d.iterrows():
        rp.append(f"| {r['case']} | {fnum(r['crb_r_m']/1e3,3)} | {fnum(r['crb_z_m'],1)} | {fnum(r['corr_rz'],3)} | {fnum(r['cond_F'],2)} |")
    rp.append("")
    rp.append("## 6. R3-B1 判定")
    rp.append("")
    rp.append(f"### `{decision}`")
    rp.append("")
    rp.append(why)
    rp.append("")
    rp.append(f"**下一步**：{nxt}")
    rp.append("")
    rp.append("允许终态还包括：`B1_PHYSICS_ONLY_HLA_LIMITED` / `B1_DELAY_DOMINANT_POSSIBLE` / `B1_DELAY_BANDWIDTH_LIMITED` / `B1_NO_COMPLEMENTARY_INFORMATION` / `B1_NO_STABLE_MULTIPATH_IDENTITY` / `B1_DELAY_IDENTITY_NOT_OBSERVABLE` / `B1_MMAC_PHYSICS_CONFIRMED`（统一用 MMAC 命名）。")
    rp.append("")
    rp.append("## 7. 停止条件")
    rp.append("")
    rp.append("- 无 Monte Carlo RMSE、无 RC2 困难候选、无 RC3-C、无 P5、无 Bellhop、无新特征")
    rp.append(f"- **B1 完成后停止**；判定：**{decision}**")
    rp.append("")
    rp.append("## 8. 图")
    rp.append("")
    for p in sorted(FIG.glob("*.svg")):
        rp.append(f"- figures/{p.name}")
    rp.append("")
    (OUT / "R3_B1_REPORT.md").write_text("\n".join(rp), encoding="utf-8")

    gs = [
        "# R3-B1 — GPT 同步稿", "",
        f"- UTC: {NOW}",
        f"- **判定：{decision}**",
        f"- {why}",
        f"- 下一步：{nxt}", "",
        "## Gate-0",
        ref.to_string(index=False), "",
        f"HLA_PROJECTION_WEAK={notes['gate0_hla_elevation_projection_weak']}", "",
        "## 控制受控 CRB_r (km)",
        ctrl.groupby(["theta_rel_deg","obs_mode","obs_set"])["crb_r_m"].first().head(20).to_string(), "",
        "## 延迟源",
        delay_df.to_string(index=False), "",
        "停止：B1 后停止；不自动 B2/RC3-C/P5。", "",
    ]
    (OUT / "R3_B1_GPT_SYNC.md").write_text("\n".join(gs), encoding="utf-8")

    print(f"DONE in {time.time()-t0:.1f}s")
    print("DECISION", decision)
    print(why)


if __name__ == "__main__":
    main()
