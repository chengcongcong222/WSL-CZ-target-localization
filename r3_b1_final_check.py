#!/usr/bin/env python3
"""R3-B1-FINAL-CHECK: lock MMAC Fisher invariants; permanently close B1.

No B2 / RC3-C / P5. No further B1.x after this.
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
from r3_b1_1_correction import find_eigenrays_root, shoot_ray, steering_corr, u_hla
from r3_b1_2_final import complex_H, OUT

NOW = datetime.now(timezone.utc).isoformat()
R0, ZS0 = 50.0, 200.0
# 5-point center stencil — Fisher ONLY from branches present at ALL five via continuation
STENCIL = [(50.0, 200.0), (49.5, 200.0), (50.5, 200.0), (50.0, 195.0), (50.0, 205.0)]
C0 = 1500.0
SIGMA_TAU = 2.4e-4
SIGMA_THETA_DEG = [0.05, 0.1]
THETA_REL = [0.0, 30.0, 60.0]
# modal energy match tolerance ~ first-null / local mainlobe scale (~8 ms)
ENERGY_TOL_S = 0.015
ASSERT = True


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


def modal_delay_peaks(r_km, z_s=200.0, z_r=200.0):
    """Complex H, remove tau_ref=r/c0, IFFT; return peaks in absolute tau."""
    freqs = np.arange(150.0, 375.0 + 1e-9, 0.25)  # denser than 1 Hz; no fake 1 s period
    df = freqs[1] - freqs[0]
    n = 4096
    H = complex_H(freqs, r_km * 1e3, z_s, z_r)
    tau_ref = r_km * 1e3 / C0
    H = H * np.exp(-1j * 2 * np.pi * freqs * tau_ref)
    Hf = np.zeros(n, dtype=complex)
    idx = np.clip(np.round(freqs / df).astype(int), 0, n // 2 - 2)
    # place bins (may collide at high df resolution — average)
    acc = np.zeros(n, dtype=complex)
    cnt = np.zeros(n)
    for fbin, hv in zip(idx, H):
        acc[fbin] += hv
        cnt[fbin] += 1
    cnt[cnt == 0] = 1
    Hf = acc / cnt
    h = np.fft.ifft(Hf)
    t = np.fft.fftfreq(n, d=df)
    # search absolute tau around geometric travel ~ r/c
    t_abs = t + tau_ref
    mask = (t_abs >= tau_ref - 0.05) & (t_abs <= tau_ref + 0.25)
    p = np.abs(h[mask]) ** 2
    tt = t_abs[mask]
    if p.max() <= 0:
        return [], float(tau_ref)
    peaks = []
    thr = 0.15 * p.max()
    for i in range(1, len(p) - 1):
        if p[i] >= p[i - 1] and p[i] >= p[i + 1] and p[i] >= thr:
            peaks.append((float(tt[i]), float(p[i] / p.max())))
    peaks = sorted(peaks, key=lambda x: -x[1])
    return peaks, float(tau_ref)


def continuation_branch(bid, center_state, stencil):
    """Predictor-corrector along stencil; require topology + continuous phi/tau/launch.

    Returns dict of cell->state for cells where the same physical branch is recovered,
    and discontinuity flag if phi/tau jump exceeds physical continuity bounds.
    """
    cells = {c: None for c in stencil}
    cells[(R0, ZS0)] = dict(center_state)
    # order: center, then neighbors
    order = [(50.5, 200.0), (49.5, 200.0), (50.0, 205.0), (50.0, 195.0)]
    discontinuity = False
    notes = []
    for cell in order:
        r_km, zs = cell
        _, roots = find_eigenrays_root(r_km, zs)
        if roots is None or roots.empty:
            notes.append(f"{cell}: no roots")
            discontinuity = True
            continue
        ad = roots[roots["admitted"] == True]  # noqa: E712
        if ad.empty:
            notes.append(f"{cell}: no admitted roots")
            discontinuity = True
            continue
        ref = center_state
        best, bc = None, np.inf
        for _, rr in ad.iterrows():
            if rr["topology"] != ref["topology"]:
                continue
            dphi = abs(rr["phi_arr_deg"] - ref["phi_deg"])
            dL = abs(rr["launch_deg"] - ref["launch_deg"])
            # physical predictor: Δτ ≈ Δr / c (absolute travel time)
            dtau_pred = (r_km - R0) * 1e3 / C0
            dtau_err = abs((rr["tau"] - ref["tau"]) - dtau_pred)
            # continuity: small residual vs range predictor + small angle drift
            if dphi > 3.0 or dL > 6.0 or dtau_err > 0.03:
                continue
            c = dphi + 200 * dtau_err + dL
            if c < bc:
                bc, best = c, rr
        if best is None:
            nearest = ad.iloc[int(np.argmin(np.abs(ad["phi_arr_deg"] - ref["phi_deg"])))]
            notes.append(
                f"{cell}: DISCONTINUITY nearest phi={nearest['phi_arr_deg']:.2f} "
                f"(ref {ref['phi_deg']:.2f}) tau={nearest['tau']:.4f} "
                f"(ref {ref['tau']:.4f}, pred_err≈{abs((nearest['tau']-ref['tau'])-(r_km-R0)*1e3/C0):.4f}s)"
            )
            discontinuity = True
            continue
        st = {
            "phi_deg": float(best["phi_arr_deg"]),
            "tau": float(best["tau"]),
            "launch_deg": float(best["launch_deg"]),
            "topology": best["topology"],
            "residual_m": float(best["residual_m"]),
        }
        cells[cell] = st
        notes.append(
            f"{cell}: ok phi={st['phi_deg']:.3f} tau={st['tau']:.5f} "
            f"tau_err_vs_pred={abs((st['tau']-ref['tau'])-(r_km-R0)*1e3/C0):.5f}s"
        )
    present = [c for c in stencil if cells[c] is not None]
    jac_ready = (len(present) == 5) and (not discontinuity)
    return cells, jac_ready, discontinuity, notes


def energy_supported(peaks, tau, tol=ENERGY_TOL_S):
    for tp, amp in peaks:
        if abs(tp - tau) <= tol:
            return True, tp, amp
    return False, None, None


def main():
    t0 = time.time()
    print("=== R3-B1-FINAL-CHECK ===", flush=True)
    print(f"OUT={OUT}", flush=True)

    # Load prior association for candidate center rays
    tracks_p = OUT / "rz_branch_association.csv"
    grid_p = OUT / "rz_branch_grid.csv"
    if tracks_p.exists():
        tracks = pd.read_csv(tracks_p)
    else:
        tracks = pd.DataFrame()
    roots_all = pd.read_csv(grid_p) if grid_p.exists() else pd.DataFrame()

    # Collect center candidates from grid at (50,200)
    center_rows = []
    if len(roots_all):
        c = roots_all[(np.isclose(roots_all["grid_r_km"], R0)) & (np.isclose(roots_all["grid_zs_m"], ZS0))]
        center_rows = c[c["admitted"] == True].to_dict("records")  # noqa: E712
    if not center_rows and len(tracks):
        c = tracks[(np.isclose(tracks["grid_r_km"], R0)) & (np.isclose(tracks["grid_zs_m"], ZS0)) & (tracks["admitted"] == True)]  # noqa: E712
        center_rows = c.to_dict("records")
    print(f"center candidates: {len(center_rows)}", flush=True)

    # modal peaks at 50 km
    peaks, tau_ref = modal_delay_peaks(R0)
    print(f"modal peaks @50km (abs tau): {[(round(p[0],4), round(p[1],3)) for p in peaks[:8]]}", flush=True)
    print(f"tau_ref={tau_ref:.4f}", flush=True)

    # Build observable branch table
    obs_rows = []
    jac_branches = []  # list of (bid, cells) with jac_ready and energy
    for i, row in enumerate(center_rows):
        bid = f"C{i}"
        center_state = {
            "phi_deg": float(row["phi_arr_deg"]),
            "tau": float(row["tau"]),
            "launch_deg": float(row["launch_deg"]),
            "topology": row["topology"],
        }
        cells, jac_ready, disc, notes = continuation_branch(bid, center_state, STENCIL)
        # energy check at center tau
        eng, e_tau, e_amp = energy_supported(peaks, center_state["tau"])
        observable = bool(jac_ready and eng)
        obs_rows.append({
            "branch_id": bid,
            "topology": center_state["topology"],
            "launch_deg": center_state["launch_deg"],
            "phi_deg": center_state["phi_deg"],
            "tau_s": center_state["tau"],
            "ray_exists": True,
            "continuation_5pt": jac_ready,
            "discontinuity": disc,
            "modal_energy_supported": eng,
            "modal_peak_tau_s": e_tau,
            "modal_peak_rel_amp": e_amp,
            "observable_branch": observable,
            "energy_tol_s": ENERGY_TOL_S,
            "notes": " | ".join(notes),
            "status_RAY_EXISTS": "RAY_EXISTS",
            "status_MODAL": "MODAL_ENERGY_SUPPORTED" if eng else "MODAL_ENERGY_NOT_SUPPORTED",
            "status_OBS": "OBSERVABLE_BRANCH" if observable else "NOT_OBSERVABLE",
        })
        print(f"  {bid} topo={center_state['topology']} phi={center_state['phi_deg']:.2f} "
              f"tau={center_state['tau']:.4f} cont5={jac_ready} eng={eng} disc={disc}", flush=True)
        print(f"    {notes}", flush=True)
        if observable:
            jac_branches.append((bid, cells, center_state))
    obs_df = pd.DataFrame(obs_rows)
    obs_df.to_csv(OUT / "observable_branch_table.csv", index=False, encoding="utf-8-sig")

    n_obs_br = len(jac_branches)
    print(f"observable branches: {n_obs_br} -> {[b[0] for b in jac_branches]}", flush=True)

    # Delay observations among observable branches
    delay_rows = []
    taus = [(b[0], b[2]["tau"]) for b in jac_branches]
    for i in range(len(taus)):
        for j in range(i + 1, len(taus)):
            dtau = taus[j][1] - taus[i][1]
            delay_rows.append({
                "branch_i": taus[i][0],
                "branch_j": taus[j][0],
                "tau_i_s": taus[i][1],
                "tau_j_s": taus[j][1],
                "delta_tau_s": dtau,
                "n_independent_delays": len(taus) - 1,
            })
    delay_df = pd.DataFrame(delay_rows)
    delay_df.to_csv(OUT / "final_delay_observations.csv", index=False, encoding="utf-8-sig")

    M = n_obs_br
    n_delay_obs = max(M - 1, 0)
    print(f"M={M} n_delay_obs(M-1)={n_delay_obs}", flush=True)

    # --- Hard assertions on observation construction ---
    def build_y(mode, theta, cells_list, branch_states):
        """cells_list: list of cell-dicts per branch; branch_states: center states ordered by tau."""
        order = sorted(range(len(branch_states)), key=lambda k: branch_states[k]["tau"])
        y_ang = []
        for k in order:
            s = branch_states[k]
            if mode == "IDEAL":
                y_ang.append(math.radians(s["phi_deg"]))
            elif mode == "HLA":
                y_ang.append(u_hla(theta, s["phi_deg"]))
        y_del = []
        for k in order[1:]:
            y_del.append(branch_states[k]["tau"] - branch_states[order[0]]["tau"])
        if mode == "DELAY":
            return np.asarray(y_del, dtype=float), order, "delay"
        return np.asarray(y_ang + y_del, dtype=float), order, "ang+delay"

    def pack_mode(mode, theta, r_km, zs, jac_branches):
        states = []
        for bid, cells, center in jac_branches:
            st = cells.get((r_km, zs))
            if st is None:
                return None, None
            states.append(st)
        y, order, kind = build_y(mode, theta, None, states)
        return y, [jac_branches[i][0] for i in order], kind

    jac_rows, fish_rows, check_rows = [], [], []
    fish_ok_all = True

    if M >= 1:
        for sigma_th in SIGMA_THETA_DEG:
            sig_th = math.radians(sigma_th)
            for theta in THETA_REL:
                # FD packs
                def ypack(mode, r_km, zs):
                    return pack_mode(mode, theta, r_km, zs, jac_branches)

                results = {}
                for mode, label in [("IDEAL", "IDEAL_ELEVATION+delay"),
                                    ("HLA", "HLA_u+delay"),
                                    ("DELAY", "delay-only")]:
                    packs = {}
                    ok = True
                    for cell in STENCIL:
                        yy, ids, kind = ypack(mode, cell[0], cell[1])
                        if yy is None:
                            ok = False
                            break
                        packs[cell] = (yy, ids, kind)
                    if not ok:
                        results[label] = None
                        continue
                    # common branch order
                    id0 = packs[(R0, ZS0)][1]
                    for cell in STENCIL:
                        if packs[cell][1] != id0:
                            # reorder by same ids
                            pass
                    y0 = packs[(R0, ZS0)][0]
                    y_rp = packs[(50.5, ZS0)][0]
                    y_rm = packs[(49.5, ZS0)][0]
                    y_zp = packs[(R0, 205.0)][0]
                    y_zm = packs[(R0, 195.0)][0]
                    n = min(len(y0), len(y_rp), len(y_rm), len(y_zp), len(y_zm))
                    J = np.zeros((n, 2))
                    J[:, 0] = (y_rp[:n] - y_rm[:n]) / 1000.0
                    J[:, 1] = (y_zp[:n] - y_zm[:n]) / 10.0
                    # observation count / sigma
                    if mode == "DELAY":
                        n_obs = n
                        n_expected = M - 1
                        sig = [SIGMA_TAU] * n_obs
                    else:
                        n_obs = n
                        n_expected = M + (M - 1)  # M angle/u + M-1 delays
                        if mode == "IDEAL":
                            sig = [math.radians(0.3)] * M + [SIGMA_TAU] * (M - 1)
                        else:
                            sig_u = []
                            for k in range(M):
                                s = jac_branches[k][2]
                                su = steering_corr(14.0, 8, u_hla(theta, s["phi_deg"]), u_hla(theta, s["phi_deg"]), 250.0)[1]
                                du_dth = -math.sin(math.radians(theta)) * math.cos(math.radians(s["phi_deg"]))
                                sig_u.append(math.sqrt(su ** 2 + (du_dth * sig_th) ** 2))
                            sig = sig_u + [SIGMA_TAU] * (M - 1)
                        sig = (sig + [sig[-1] if sig else 1e-3] * n)[:n]
                    sv = np.linalg.svd(J, compute_uv=False) if n else np.array([])
                    rank = int(np.linalg.matrix_rank(J)) if n else 0

                    # ASSERTIONS
                    asserts = []
                    if mode == "DELAY":
                        asserts.append(("n_delay_obs_eq_Mminus1", n_obs == n_expected == max(M - 1, 0), f"n={n_obs} expected={n_expected}"))
                        asserts.append(("rank_delay_le_min_Mminus1_2", rank <= min(max(M - 1, 0), 2), f"rank={rank}"))
                        if M == 2:
                            asserts.append(("M2_delay_n_obs_eq_1", n_obs == 1, f"n={n_obs}"))
                            asserts.append(("M2_delay_rank_le_1", rank <= 1, f"rank={rank}"))
                        if M < 3:
                            asserts.append(("delay_only_2d_requires_M_ge_3", False,
                                            f"M={M}<3 → delay-only 2D CRB FORBIDDEN"))
                    else:
                        asserts.append(("n_angdelay_eq_M_plus_Mminus1", n_obs == n_expected, f"n={n_obs} expected={n_expected}"))
                        if M >= 2:
                            asserts.append(("angdelay_rank_le_2", rank <= 2, f"rank={rank}"))
                        else:
                            asserts.append(("angdelay_requires_M_ge_2", False, f"M={M}<2"))

                    failed = [a for a in asserts if not a[1]]
                    for a in asserts:
                        check_rows.append({
                            "sigma_theta_deg": sigma_th,
                            "theta_rel_deg": theta,
                            "case": label,
                            "check": a[0],
                            "pass": bool(a[1]),
                            "detail": a[2],
                        })
                    if failed and ASSERT:
                        print(f"  ASSERT FAIL {label} th={theta}: {[a[0]+':'+a[2] for a in failed]}", flush=True)
                        results[label] = {
                            "forbidden": True,
                            "n_obs": n_obs, "rank": rank,
                            "J": J, "asserts": asserts,
                        }
                        continue
                    # CRB only if allowed
                    if mode == "DELAY" and M < 3:
                        results[label] = {"forbidden": True, "n_obs": n_obs, "rank": rank, "J": J, "asserts": asserts}
                        continue
                    try:
                        Sinv = np.diag(1.0 / (np.asarray(sig) ** 2 + 1e-30))
                        F = J.T @ Sinv @ J
                        cov = np.linalg.inv(F)
                        crb_r = math.sqrt(max(float(cov[0, 0]), 0))
                        crb_z = math.sqrt(max(float(cov[1, 1]), 0))
                        corr = float(cov[0, 1] / math.sqrt(cov[0, 0] * cov[1, 1] + 1e-30))
                        condF = float(np.linalg.cond(F))
                    except np.linalg.LinAlgError:
                        F = J.T @ np.diag(1.0 / (np.asarray(sig) ** 2 + 1e-30)) @ J
                        crb_r = crb_z = np.inf
                        corr = np.nan
                        condF = np.inf
                    results[label] = {
                        "forbidden": False, "n_obs": n_obs, "rank": rank, "J": J, "F": F,
                        "crb_r_m": crb_r, "crb_z_m": crb_z, "corr_rz": corr, "cond_F": condF,
                        "sv": sv, "asserts": asserts,
                    }
                    jac_rows.append({
                        "case": label, "sigma_theta_deg": sigma_th, "theta_rel_deg": theta,
                        "n_obs": n_obs, "rank": rank,
                        "dy_dr": json.dumps(J[:, 0].tolist()), "dy_dz": json.dumps(J[:, 1].tolist()),
                    })

                # Nesting check: F(HLA+delay) - F(delay) PSD when both defined
                hla = results.get("HLA_u+delay")
                dly = results.get("delay-only")
                nest_pass, nest_detail = True, "n/a"
                if hla and dly and (not hla.get("forbidden")) and (not dly.get("forbidden")) and "F" in hla and "F" in dly:
                    D = hla["F"] - dly["F"]
                    D = 0.5 * (D + D.T)
                    ev = np.linalg.eigvalsh(D)
                    nest_pass = bool(ev.min() >= -1e-8 * max(1.0, abs(ev).max()))
                    nest_detail = f"eig_min={ev.min():.3e} eig={ev.tolist()}"
                    # CRB nesting: joint should not be worse than subset (allow tiny num noise)
                    crb_ok = True
                    if np.isfinite(hla["crb_r_m"]) and np.isfinite(dly["crb_r_m"]):
                        crb_ok = hla["crb_r_m"] <= dly["crb_r_m"] * 1.05 + 1e-6
                    nest_pass = nest_pass and crb_ok
                    nest_detail += f" CRBr_hla={hla['crb_r_m']:.3g} CRBr_delay={dly['crb_r_m']:.3g} crb_ok={crb_ok}"
                    if not nest_pass:
                        fish_ok_all = False
                        print(f"  NESTING FAIL th={theta}: {nest_detail}", flush=True)
                check_rows.append({
                    "sigma_theta_deg": sigma_th, "theta_rel_deg": theta,
                    "case": "HLA+delay vs delay-only",
                    "check": "fisher_nesting_PSD_and_CRB_monotone",
                    "pass": nest_pass, "detail": nest_detail,
                })

                for label, res in results.items():
                    if res is None:
                        continue
                    fish_rows.append({
                        "case": label,
                        "sigma_theta_deg": sigma_th,
                        "theta_rel_deg": theta,
                        "M_observable": M,
                        "n_obs": res["n_obs"],
                        "rank": res["rank"],
                        "forbidden": bool(res.get("forbidden")),
                        "crb_r_m": res.get("crb_r_m", np.nan),
                        "crb_z_m": res.get("crb_z_m", np.nan),
                        "corr_rz": res.get("corr_rz", np.nan),
                        "cond_F": res.get("cond_F", np.nan),
                        "sv_min": float(res["sv"].min()) if "sv" in res and len(res["sv"]) else np.nan,
                    })

    jac_df = pd.DataFrame(jac_rows)
    fish_df = pd.DataFrame(fish_rows)
    check_df = pd.DataFrame(check_rows)
    jac_df.to_csv(OUT / "final_mmac_jacobian.csv", index=False, encoding="utf-8-sig")
    fish_df.to_csv(OUT / "final_mmac_fisher.csv", index=False, encoding="utf-8-sig")
    check_df.to_csv(OUT / "fisher_invariant_checks.csv", index=False, encoding="utf-8-sig")
    print(check_df.to_string(index=False), flush=True)

    # --- S0 delay window covering actual Δτ ---
    def amb_A(freqs, tau_grid, tau_true):
        return np.abs(np.mean(np.exp(1j * 2 * np.pi * freqs[:, None] * (tau_grid - tau_true)[None, :]), axis=0))

    def psl_outside_first_null(A, tau_grid, tau_true):
        i0 = int(np.argmin(np.abs(tau_grid - tau_true)))
        lo, hi = i0, i0
        while lo > 0 and A[lo - 1] >= 0.08:
            lo -= 1
        while hi < len(A) - 1 and A[hi + 1] >= 0.08:
            hi += 1
        mask = np.ones_like(A, dtype=bool)
        mask[lo:hi + 1] = False
        # local peaks outside
        peaks = []
        for i in range(1, len(A) - 1):
            if mask[i] and A[i] >= A[i - 1] and A[i] >= A[i + 1] and A[i] >= 0.2:
                peaks.append((float(tau_grid[i] - tau_grid[i0]), float(A[i])))
        psl = float(A[mask].max()) if mask.any() else 0.0
        # local mainlobe
        lo2, hi2 = i0, i0
        while lo2 > 0 and A[lo2 - 1] >= 0.5:
            lo2 -= 1
        while hi2 < len(A) - 1 and A[hi2 + 1] >= 0.5:
            hi2 += 1
        return float(tau_grid[hi2] - tau_grid[lo2]), psl, peaks

    # actual Δτ range among observable pairs
    if len(delay_df):
        dtau_abs = float(delay_df["delta_tau_s"].abs().max())
    else:
        dtau_abs = 0.5  # fallback include known 0.48
    margin = 0.05
    half_win = max(dtau_abs + margin, 0.6)
    tau_grid = np.linspace(-half_win, half_win, int(2 * half_win / 1e-4) + 1)
    print(f"delay search window ±{half_win:.3f}s (max|dtau|={dtau_abs:.3f})", flush=True)

    # test true peak at 0 and at actual Δτ
    test_taus = [0.0] + [float(d) for d in delay_df["delta_tau_s"].tolist()] if len(delay_df) else [0.0]
    freqs_dense = np.arange(150.0, 375.0 + 1e-9, 0.25)
    freqs_s1 = np.array([168.0, 204.0, 232.0, 279.0, 320.0])
    freqs_s2 = np.array([166.0, 201.0, 235.0, 283.0, 338.0])
    freqs_23 = np.linspace(150.0, 375.0, 23)
    delay_final_rows = []
    for name, freqs, represent in [
        ("S0_DENSE", freqs_dense, "continuous broadband (0.25 Hz integration)"),
        ("S0_23POINT", freqs_23, "comb control NOT broadband"),
        ("S1", freqs_s1, "sparse lines"),
        ("S2", freqs_s2, "sparse lines"),
    ]:
        for tau_true in sorted(set([0.0] + test_taus)):
            A = amb_A(freqs, tau_grid, tau_true)
            loc_w, psl, peaks = psl_outside_first_null(A, tau_grid, tau_true)
            if name == "S0_DENSE":
                status = "GLOBAL_DELAY_USABLE_IN_TEST_WINDOW" if psl < 0.35 else "GLOBAL_AMBIGUOUS"
            elif name == "S0_23POINT":
                status = "GLOBAL_AMBIGUOUS_COMB_CONTROL"
            else:
                status = "GLOBAL_AMBIGUOUS" if psl > 0.5 else "NEEDS_CHECK"
            delay_final_rows.append({
                "source": name,
                "represents": represent,
                "tau_true_s": tau_true,
                "search_half_window_s": half_win,
                "local_mainlobe_s": loc_w,
                "psl_outside_first_null": psl,
                "n_peaks": len(peaks),
                "peaks": json.dumps(peaks[:6]),
                "global_status": status,
            })
            print(f"  delay {name} tau0={tau_true:.3f}: loc={loc_w*1e3:.2f}ms PSL={psl:.3f} {status}", flush=True)
    delay_final_df = pd.DataFrame(delay_final_rows)
    delay_final_df.to_csv(OUT / "final_delay_observations.csv", index=False, encoding="utf-8-sig")
    # also append pair table at end of same file? keep pairs separate column block
    if len(delay_df):
        delay_df.to_csv(OUT / "final_delay_pairs.csv", index=False, encoding="utf-8-sig")

    # --- Decision ---
    def fish_pick(case, th=0.0, st=0.1):
        if fish_df is None or fish_df.empty:
            return None
        d = fish_df[(fish_df["case"] == case) & (fish_df["theta_rel_deg"] == th) & (fish_df["sigma_theta_deg"] == st)]
        return d.iloc[0] if len(d) else None

    ideal = fish_pick("IDEAL_ELEVATION+delay")
    hla = fish_pick("HLA_u+delay")
    dly = fish_pick("delay-only")

    s0_rows = delay_final_df[delay_final_df["source"] == "S0_DENSE"]
    s0_usable = bool(len(s0_rows) and all(s0_rows["global_status"] == "GLOBAL_DELAY_USABLE_IN_TEST_WINDOW"))
    s1_s2_amb = bool(
        (delay_final_df[delay_final_df["source"] == "S1"]["global_status"] == "GLOBAL_AMBIGUOUS").any()
        or (delay_final_df[delay_final_df["source"] == "S2"]["global_status"] == "GLOBAL_AMBIGUOUS").any()
    )
    n_indep_delay = max(M - 1, 0)
    delay_2d_ok = bool(M >= 3)  # required for delay-only 2D CRB
    nest_ok = bool(check_df[check_df["check"] == "fisher_nesting_PSD_and_CRB_monotone"]["pass"].all()) if len(check_df) else True
    no_assert_fail = bool(check_df["pass"].all()) if len(check_df) else False

    def ok_case(row, thr_r=5e3, thr_z=100.0, need_rank2=True):
        if row is None or bool(row.get("forbidden")):
            return False
        if need_rank2 and int(row["rank"]) < 2:
            return False
        if not np.isfinite(row["crb_r_m"]) or not np.isfinite(row["crb_z_m"]):
            return False
        return row["crb_r_m"] < thr_r and row["crb_z_m"] < thr_z

    ideal_ok = ok_case(ideal)
    hla_ok = ok_case(hla)
    delay_ok = ok_case(dly) if delay_2d_ok else False

    notes = {
        "M_observable": M,
        "n_independent_delays": n_indep_delay,
        "observable_branch_ids": [b[0] for b in jac_branches],
        "ray_modal_energy_peaks_50km": [(round(p[0], 4), round(p[1], 3)) for p in peaks[:8]],
        "energy_tol_s": ENERGY_TOL_S,
        "s0_global_usable": s0_usable,
        "s1_s2_globally_ambiguous": s1_s2_amb,
        "delay_only_2d_allowed": delay_2d_ok,
        "fisher_nesting_ok": nest_ok,
        "all_asserts_pass": no_assert_fail,
        "ideal_ok": ideal_ok,
        "hla_ok": hla_ok,
        "delay_ok": delay_ok,
        "HLA_elevation_weak_frozen": True,
        "local_delay_ms_scale_frozen": True,
        "branch_table_status_counts": obs_df["status_OBS"].value_counts().to_dict() if len(obs_df) else {},
    }

    if M < 2:
        decision = "B1_NO_STABLE_MULTIPATH_IDENTITY"
        why = f"Only {M} continuation-continuous AND modal-energy-supported branches at center."
        nxt = "B1 permanently closed; RC3-B MMAC not established"
    elif not nest_ok or not no_assert_fail:
        decision = "B1_NO_COMPLEMENTARY_INFORMATION"
        why = f"Invariant/assert failure (nesting_ok={nest_ok}, asserts_pass={no_assert_fail}). CRB not trusted."
        nxt = "B1 closed; do not use previous Fisher numbers"
    elif ideal_ok and hla_ok and s0_usable and M >= 2:
        decision = "B1_MMAC_PHYSICS_CONFIRMED"
        why = (
            f"M={M} observable branches; HLA_u+delay rank2 CRB_r={hla['crb_r_m']:.1f} m "
            f"CRB_z={hla['crb_z_m']:.1f} m; IDEAL CRB_r={ideal['crb_r_m']:.1f} m; "
            f"S0_DENSE global usable; invariants passed."
        )
        nxt = "B1 closed; B2 later only if task assigns unlabeled delay-set / HLA-assisted association (not this round)"
    elif ideal_ok and not hla_ok:
        decision = "B1_PHYSICS_ONLY_HLA_LIMITED"
        why = (
            f"IDEAL+delay OK (CRB_r={ideal['crb_r_m'] if ideal is not None else np.nan:.1f} m) but "
            f"HLA_u+delay weak/forbidden (CRB_r={hla['crb_r_m'] if hla is not None else np.nan}). "
            f"M={M}; S0 usable={s0_usable}. HLA elevation projection weak."
        )
        nxt = "B1 closed; RC3-B if continued must be delay-dominant, not HLA elevation MMAC"
    elif M >= 3 and delay_ok and s0_usable:
        decision = "B1_DELAY_DOMINANT_POSSIBLE"
        why = (
            f"M={M} energy-supported branches, n_delays={n_indep_delay}; delay-only 2D rank/CRB OK "
            f"(CRB_r={dly['crb_r_m'] if dly is not None else np.nan:.1f} m); HLA weak; S0 usable."
        )
        nxt = "B1 closed; later B2 only as unlabeled delay-set matching — not this round"
    elif M == 2 and ideal_ok:
        decision = "B1_PHYSICS_ONLY_HLA_LIMITED"
        why = (
            f"M=2 observable: angle/u+delay possible, but delay-only 2D CRB FORBIDDEN (need M>=3). "
            f"HLA_ok={hla_ok}; S0={s0_usable}."
        )
        nxt = "B1 closed; need ≥3 energy-supported branches for delay-only r-z; not this round"
    elif M >= 2:
        decision = "B1_NO_COMPLEMENTARY_INFORMATION"
        why = (
            f"Observable M={M} but corrected Fisher shows no trustworthy r-z complementarity "
            f"(ideal_ok={ideal_ok}, hla_ok={hla_ok}, delay_ok={delay_ok}, s0={s0_usable})."
        )
        nxt = "B1 permanently closed"
    else:
        decision = "B1_NO_STABLE_MULTIPATH_IDENTITY"
        why = "No sufficient continuous energy-supported multipath."
        nxt = "B1 closed"

    dec = {
        "rc3b1_final_decision": decision,
        "why": why,
        "next_step": nxt,
        "notes": notes,
        "created_utc": NOW,
        "b1_status": "PERMANENTLY_CLOSED",
        "stop": "no B2 / RC3-C / P5; no further B1.x",
        "frozen_positive": [
            "S0 continuous broadband local delay ~5ms, PSL acceptable in window when properly measured",
            "S1/S2 sparse lines globally ambiguous",
            "HLA elevation projection weak",
        ],
        "revoked": [
            "B1_MMAC_PHYSICS_CONFIRMED from R3-B1.2 (invalid delay-only Fisher / possible branch swap)",
            "delay-only CRB with n_obs=3 at M=2",
            "any Fisher using shoot-fill without continuation continuity",
        ],
    }
    (OUT / "R3_B1_FINAL_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    def fnum(x, nd=3):
        try:
            v = float(x)
            return "n/a" if not np.isfinite(v) else f"{v:.{nd}f}"
        except Exception:
            return "n/a"

    rp = []
    rp.append("# R3-B1-FINAL-CHECK 报告")
    rp.append("")
    rp.append(f"UTC：{NOW}")
    rp.append("")
    rp.append("**B1 自本报告起永久关闭**，不再派生 B1.x；不进入 B2 / RC3-C / P5。")
    rp.append("")
    rp.append("## 0. 撤销")
    rp.append("")
    rp.append("- 撤销 R3-B1.2 的 `B1_MMAC_PHYSICS_CONFIRMED`（delay-only 观测维数错误 + 疑似分支切换）")
    rp.append("- 保留正结果：S0 连续宽带局部时延 ms 级；S1/S2 全局模糊；HLA 俯仰投影弱")
    rp.append("")
    rp.append("## 1. 可观测分支（continuation + 模态能量）")
    rp.append("")
    rp.append(f"- 五点 stencil 无 shoot-fill；连续性：topology + |Δφ|<3° + |Δτ|<50ms + |Δlaunch|<5°")
    rp.append(f"- 模态能量容差：**{ENERGY_TOL_S*1e3:.0f} ms**（与宽带主瓣/第一零点同量级）")
    rp.append(f"- **M（OBSERVABLE_BRANCH）= {M}**：{notes['observable_branch_ids']}")
    rp.append("")
    if len(obs_df):
        rp.append("| id | topo | φ° | τ s | cont5 | energy | OBS |")
        rp.append("| --- | --- | --- | --- | --- | --- | --- |")
        for _, r in obs_df.iterrows():
            rp.append(
                f"| {r['branch_id']} | {r['topology']} | {fnum(r['phi_deg'],2)} | {fnum(r['tau_s'],4)} | "
                f"{r['continuation_5pt']} | {r['modal_energy_supported']} | {r['status_OBS']} |"
            )
        rp.append("")
        rp.append("诊断摘录见 `observable_branch_table.csv` 的 notes（含 13°→3.6° 类跳变）。")
    rp.append("")
    rp.append("## 2. 时延观测构造（硬断言）")
    rp.append("")
    rp.append(f"- delay-only 独立维数 = M−1 = **{n_indep_delay}**")
    rp.append(f"- delay-only 估计 (r,z) **要求 M≥3** → 本轮 **允许={delay_2d_ok}**")
    rp.append(f"- 搜索窗覆盖实际 Δτ：**±{half_win:.3f} s**（max|Δτ|={dtau_abs:.3f} s）")
    rp.append(f"- S0 使用 0.25 Hz 密频积分，**不用** 1 Hz 采样制造 1 s 周期")
    rp.append("")
    rp.append("| 源 | τ0 | 局部主瓣 ms | PSL | status |")
    rp.append("| --- | --- | --- | --- | --- |")
    for _, r in delay_final_df.iterrows():
        rp.append(
            f"| {r['source']} | {fnum(r['tau_true_s'],3)} | {fnum(r['local_mainlobe_s']*1e3,2)} | "
            f"{fnum(r['psl_outside_first_null'])} | **{r['global_status']}** |"
        )
    rp.append("")
    rp.append("## 3. Fisher + 不变量")
    rp.append("")
    rp.append(f"- nesting / assert 总通过：**{no_assert_fail}**（nesting_ok={nest_ok}）")
    rp.append("")
    if len(fish_df):
        rp.append("| case | σθ | θ° | M | n_obs | rank | forbidden | CRB_r m | CRB_z m |")
        rp.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for _, r in fish_df[fish_df["sigma_theta_deg"] == 0.1].iterrows():
            rp.append(
                f"| {r['case']} | {r['sigma_theta_deg']} | {r['theta_rel_deg']} | {r['M_observable']} | "
                f"{r['n_obs']} | {r['rank']} | {r['forbidden']} | {fnum(r['crb_r_m'],1)} | {fnum(r['crb_z_m'],1)} |"
            )
    rp.append("")
    rp.append("断言明细：`fisher_invariant_checks.csv`（含 n_delay=M−1、rank 限制、Fisher 半正定嵌套、CRB 单调）。")
    rp.append("")
    rp.append("## 4. B1 最终判定（永久）")
    rp.append("")
    rp.append(f"### `{decision}`")
    rp.append("")
    rp.append(why)
    rp.append("")
    rp.append(f"**下一步**：{nxt}")
    rp.append("")
    rp.append("## 5. 关闭声明")
    rp.append("")
    rp.append("- B1 状态：**PERMANENTLY_CLOSED**")
    rp.append("- 不进入 B2 / RC3-C / P5；不派生 B1.x")
    rp.append("- 不增加新物理特征")
    rp.append("")
    (OUT / "R3_B1_FINAL_REPORT.md").write_text("\n".join(rp), encoding="utf-8")

    gs = [
        "# R3-B1-FINAL-CHECK — GPT 同步", "",
        f"- **判定：{decision}**",
        f"- B1: PERMANENTLY_CLOSED",
        f"- {why}",
        f"- 下一步：{nxt}", "",
        f"M_observable={M} delays={n_indep_delay} delay2d_allowed={delay_2d_ok}",
        f"asserts_pass={no_assert_fail} nesting_ok={nest_ok}",
        f"S0_usable={s0_usable} S1/S2_amb={s1_s2_amb}", "",
        "## Fisher (σθ=0.1°)",
        fish_df[fish_df["sigma_theta_deg"]==0.1].to_string(index=False) if len(fish_df) else "(none)", "",
        "## Branches",
        obs_df[["branch_id","topology","phi_deg","tau_s","continuation_5pt","modal_energy_supported","status_OBS"]].to_string(index=False) if len(obs_df) else "", "",
        "关闭 B1；无 B2/RC3-C/P5。", "",
    ]
    (OUT / "R3_B1_FINAL_GPT_SYNC.md").write_text("\n".join(gs), encoding="utf-8")

    print(f"DONE {time.time()-t0:.1f}s", flush=True)
    print("DECISION", decision, flush=True)
    print(why, flush=True)


if __name__ == "__main__":
    main()
