#!/usr/bin/env python3
"""R3-C2.0S: self-adjoint modes + effective-mode convergence + finite-aperture depth signatures.

No full Yang estimator / δ-z / Doppler / f0 / MC / P5.
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
OUT = ROOT / "results" / "R3_C2_Yang_SA_depth"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

H = 5000.0
C0 = 1500.0
FREQS = [201.0, 235.0, 283.0, 338.0]
Z_S = 200.0
Z_R = 200.0
Z_PROBE = [180.0, 190.0, 200.0, 210.0, 220.0]
PAIRS = [(180, 190), (190, 200), (200, 210), (210, 220), (180, 200), (200, 220)]
L_LIST = [0.6e3, 1.2e3, 2.4e3]
NZ_LIST = [2001, 4001, 8001]
THR_LIST = [0.1, 0.2]
C1_CORR = {201.0: 0.579, 235.0: 0.994, 283.0: 1.0, 338.0: 0.067}


def munk_c(z, za=1300.0, c0=1500.0, B=1300.0, eps=0.00737):
    z = np.asarray(z, float)
    eta = 2.0 * (z - za) / B
    return c0 * (1.0 + eps * (np.exp(eta) - 1.0 - eta))


def fem_assemble(z, c, omega):
    """Linear FEM: K φ = γ² M  (for -φ''), Q = ∫(ω²/c²)φv.

    Generalized: (Q − K) φ = k² M φ,  k² = ω²/c²_eff − γ² style.
    BC: φ(0)=0 (Dirichlet eliminate node 0), φ'(H)=0 (natural Neumann).
    """
    z = np.asarray(z, float)
    n = len(z)
    h = float(z[1] - z[0])
    # assemble on all nodes then eliminate 0
    K = np.zeros((n, n))
    M = np.zeros((n, n))
    Q = np.zeros((n, n))
    q = (omega / c) ** 2  # ω²/c²(z) at nodes — use element midpoint for Q if needed
    for e in range(n - 1):
        i, j = e, e + 1
        ke = np.array([[1.0, -1.0], [-1.0, 1.0]]) / h
        me = (h / 6.0) * np.array([[2.0, 1.0], [1.0, 2.0]])
        # Q element: ∫ (ω²/c²) φ_i φ_j ≈ (ω²/c²)_avg * M_e
        qe = 0.5 * (q[i] + q[j]) * me
        for a, ia in enumerate([i, j]):
            for b, ib in enumerate([i, j]):
                K[ia, ib] += ke[a, b]
                M[ia, ib] += me[a, b]
                Q[ia, ib] += qe[a, b]
    # Dirichlet φ0=0: use nodes 1..n-1
    K = K[1:, 1:]
    M = M[1:, 1:]
    Q = Q[1:, 1:]
    # A φ = k² M φ  with A = Q - K  (symmetric)
    A = 0.5 * (Q - K + (Q - K).T)
    M = 0.5 * (M + M.T)
    return A, M, z


def solve_modes_sa(z, c, omega, kr_max, n_keep=None):
    """Symmetric generalized eigen via Cholesky M: A y = k² M y."""
    A, M, z = fem_assemble(z, c, omega)
    # Cholesky M = L L^T
    L = np.linalg.cholesky(M)
    # L^{-1} A L^{-T} y = k² y
    Li = np.linalg.inv(L)
    B = Li @ A @ Li.T
    B = 0.5 * (B + B.T)
    evals, y = np.linalg.eigh(B)
    modes = []
    order = np.argsort(evals)[::-1]
    for j in order:
        k2 = float(evals[j])
        if not np.isfinite(k2) or k2 <= 0:
            continue
        kr = math.sqrt(k2)
        if kr > kr_max + 1e-6:
            continue
        v = Li.T @ y[:, j]  # φ = L^{-T} y
        phi = np.zeros(len(z))
        phi[1:] = v
        # physical L2 normalize
        nrm = math.sqrt(abs(float(np.trapezoid(phi ** 2, z)))) + 1e-30
        phi = phi / nrm
        modes.append({"k_r": kr, "k2": k2, "phi": phi})
    modes = sorted(modes, key=lambda m: -m["k_r"])
    if n_keep:
        modes = modes[:n_keep]
    return modes, z


def inner(f, g, z):
    return float(np.trapezoid(f * g, z))


def shape_corr(p, q, z):
    ip = inner(p, q, z)
    np_ = math.sqrt(abs(inner(p, p, z))) + 1e-30
    nq = math.sqrt(abs(inner(q, q, z))) + 1e-30
    return abs(ip) / (np_ * nq)


def constant_c_unit_final():
    f = 201.0
    omega = 2 * np.pi * f
    rows = []
    ok_all = True
    for n_z in [2001, 4001]:
        z = np.linspace(0.0, H, n_z)
        c = np.full_like(z, C0)
        modes, z = solve_modes_sa(z, c, omega, omega / C0, n_keep=25)
        for n in range(min(20, len(modes))):
            gamma = (n + 0.5) * np.pi / H
            kr_an = math.sqrt(max((omega / C0) ** 2 - gamma ** 2, 0.0))
            phi_an = np.sin(gamma * z)
            phi_an /= math.sqrt(abs(inner(phi_an, phi_an, z))) + 1e-30
            # match by k_r
            j = int(np.argmin([abs(m["k_r"] - kr_an) for m in modes]))
            phi_n = modes[j]["phi"]
            sc = shape_corr(phi_n, phi_an, z)
            err_k = abs(modes[j]["k_r"] - kr_an)
            # BC residual: φ(0), φ'(H)
            h = z[1] - z[0]
            bc = abs(phi_n[0]) + abs((phi_n[-1] - phi_n[-2]) / h)
            ok = err_k < 1e-3 and sc > 0.98 and bc < 1e-3
            ok_all = ok_all and ok
            rows.append({
                "n_z": n_z, "n": n, "gamma_an": gamma,
                "kr_an": kr_an, "kr_num": modes[j]["k_r"], "err_kr": err_k,
                "shape_corr_integral": sc, "bc_residual": bc, "pass": ok,
            })
    return pd.DataFrame(rows), ok_all


def effective_indices(props, z, thr):
    A = np.array([inner(p["phi"], np.ones_like(z) * 0 + 1, z) * 0 for p in props])  # placeholder
    # A_m = φ(zs)*φ(zr) signed
    A = np.array([
        float(np.interp(Z_S, z, p["phi"])) * float(np.interp(Z_R, z, p["phi"]))
        for p in props
    ])
    An = np.abs(A) / (np.max(np.abs(A)) + 1e-30)
    idx = [i for i, a in enumerate(An) if a > thr]
    return idx, A, An


def match_modes(k_ref, phi_ref, z, k_new, phi_new):
    """Greedy min-cost assignment on (Δk, 1-overlap)."""
    used = set()
    pairs = []
    for i, kr in enumerate(k_ref):
        best, best_j = np.inf, -1
        for j, kr2 in enumerate(k_new):
            if j in used:
                continue
            ov = shape_corr(phi_ref[i], phi_new[j], z)
            cost = abs(kr - kr2) / 0.01 + (1.0 - ov) * 10
            if cost < best:
                best, best_j = cost, j
        if best_j >= 0:
            used.add(best_j)
            ov = shape_corr(phi_ref[i], phi_new[best_j], z)
            pairs.append({
                "i_ref": i, "j_new": best_j,
                "k_ref": kr, "k_new": k_new[best_j],
                "delta_k": abs(kr - k_new[best_j]),
                "overlap": ov,
            })
    return pairs


def sa_recover(k_inj, A_signed, L, dr, f):
    """Controlled SA: p(r)=Σ A_m exp(j k_m r); FFT peaks; identity match."""
    n = int(max(256, L / dr))
    r = np.arange(n) * dr
    p = np.zeros(n, dtype=complex)
    for k, a in zip(k_inj, A_signed):
        p += complex(a, 0) * np.exp(1j * k * r)
    P = np.fft.fftshift(np.fft.fft(p * np.hanning(n)))
    kg = np.fft.fftshift(np.fft.fftfreq(n, d=dr)) * 2 * np.pi
    mag = np.abs(P)
    k_nyq = np.pi / dr
    # interior peaks
    peaks = []
    for i in range(2, n - 2):
        if abs(kg[i]) > 0.85 * k_nyq:
            continue
        if mag[i] >= mag[i - 1] and mag[i] >= mag[i + 1] and mag[i] > 0.05 * mag.max():
            peaks.append((float(kg[i]), float(mag[i])))
    # match each injected mode to a peak within Rayleigh
    dkr = 2 * np.pi / L
    used_p = set()
    rows = []
    rec_k = []
    for i, (k, a) in enumerate(zip(k_inj, A_signed)):
        cand = [j for j, (kp, _) in enumerate(peaks) if abs(kp - k) < 0.5 * dkr and j not in used_p]
        if cand:
            j = min(cand, key=lambda jj: abs(peaks[jj][0] - k))
            used_p.add(j)
            rows.append({
                "mode_index": i, "true_k_r": k, "A_signed": a, "A_norm": abs(a),
                "recovered": True, "peak_k": peaks[j][0],
                "delta_k": abs(peaks[j][0] - k),
                "merged_with_modes": 0,
                "status": "RESOLVED",
            })
            rec_k.append(peaks[j][0])
        else:
            # merged if another mode closer to some peak
            rows.append({
                "mode_index": i, "true_k_r": k, "A_signed": a, "A_norm": abs(a),
                "recovered": False, "peak_k": np.nan, "delta_k": np.nan,
                "merged_with_modes": 0, "status": "MISSED",
            })
    for j, (kp, mp) in enumerate(peaks):
        if j not in used_p:
            rows.append({
                "mode_index": np.nan, "true_k_r": np.nan, "A_signed": np.nan,
                "A_norm": np.nan, "recovered": False, "peak_k": kp,
                "delta_k": np.nan, "merged_with_modes": 0, "status": "SPURIOUS",
            })
    # merge stats: peaks that sit between two modes
    for i, (k, a) in enumerate(zip(k_inj, A_signed)):
        near = [j for j, (kp, _) in enumerate(peaks) if abs(kp - k) < dkr]
        if len(near) > 1:
            for r in rows:
                if r.get("mode_index") == i and r["status"] == "RESOLVED":
                    r["merged_with_modes"] = max(0, len(near) - 1)
                    if r["merged_with_modes"] >= 2:
                        r["status"] = "MERGED"
    n_res = sum(1 for r in rows if r["status"] == "RESOLVED")
    n_mer = sum(1 for r in rows if r["status"] == "MERGED")
    n_mis = sum(1 for r in rows if r["status"] == "MISSED")
    n_spu = sum(1 for r in rows if r["status"] == "SPURIOUS")
    n_inj = len(k_inj)
    precision = n_res / max(n_res + n_spu + n_mer, 1)
    recall = (n_res + n_mer) / max(n_inj, 1)
    frac = n_res / max(n_inj, 1)
    return rows, dict(precision=precision, recall=recall, fraction_resolved=frac,
                      n_injected=n_inj, n_resolved=n_res, n_merged=n_mer, n_missed=n_mis, n_spurious=n_spu)


def depth_D(idxs, props, z, zi, zj):
    a = np.array([float(np.interp(zi, z, props[i]["phi"])) for i in idxs])
    b = np.array([float(np.interp(zj, z, props[i]["phi"])) for i in idxs])
    if len(a) < 1:
        return np.nan
    a = a / (np.linalg.norm(a) + 1e-30)
    b = b / (np.linalg.norm(b) + 1e-30)
    return 1.0 - abs(float(np.dot(a, b)))


def main():
    t0 = time.time()
    print("=== R3-C2.0S ===", flush=True)
    # freeze banner
    (OUT / "R3_C2_0S_FROZEN_NOTES.md").write_text(
        "# Frozen / renamed\n\n"
        "- Erratum δ, operator [D2+ω²/c²], k_r≤ω/c_min: FROZEN\n"
        "- Old ModeModel f-invariance: ARTIFACT\n"
        "- Fourier/Rayleigh = pre-gate only\n"
        "- Prior D≈0.9 full-vector: **ORACLE_FULL_MODE_SIGNATURE** (not SA capability)\n",
        encoding="utf-8",
    )

    ut, ok = constant_c_unit_final()
    ut.to_csv(OUT / "constant_c_mode_unit_test_final.csv", index=False, encoding="utf-8-sig")
    print("constant_c_final", ok, flush=True)
    if not ok:
        dec = {"rc3c2_0s_decision": "C2_0S_EFFECTIVE_MODES_NOT_CONVERGED",
               "why": "constant-c self-adjoint FEM unit test failed", "created_utc": NOW}
        (OUT / "R3_C2_0S_DECISION.json").write_text(json.dumps(dec, indent=2), encoding="utf-8")
        return

    # Effective-mode convergence across grids
    conv_rows = []
    set_rows = []
    store = {}  # (f,nz) -> (z, props, A, An)
    for f in FREQS:
        omega = 2 * np.pi * f
        for nz in NZ_LIST:
            z = np.linspace(0.0, H, nz)
            c = munk_c(z)
            cmin = float(np.min(c))
            modes, z = solve_modes_sa(z, c, omega, omega / cmin)
            idx2, A, An = effective_indices(modes, z, 0.2)
            idx1, _, _ = effective_indices(modes, z, 0.1)
            store[(f, nz)] = (z, modes, A, An, idx1, idx2)
            set_rows.append({
                "f_hz": f, "n_z": nz,
                "n_propagating": len(modes),
                "n_eff_0p1": len(idx1),
                "n_eff_0p2": len(idx2),
            })
            print(f"  f={f} nz={nz} nprop={len(modes)} neff2={len(idx2)}", flush=True)
        # match eff 0.2 modes 4001 ↔ 2001 and 4001 ↔ 8001
        z4, m4, A4, An4, _, i4 = store[(f, 4001)]
        for nz_other in [2001, 8001]:
            zo, mo, Ao, Ano, _, io = store[(f, nz_other)]
            # map indices into full props
            k_ref = [m4[i]["k_r"] for i in i4]
            p_ref = [m4[i]["phi"] for i in i4]
            k_new = [mo[j]["k_r"] for j in io]
            p_new = [mo[j]["phi"] for j in io]
            # interpolate z onto 4001 for overlap
            p_new_i = [np.interp(z4, zo, p) for p in p_new]
            pairs = match_modes(k_ref, p_ref, z4, k_new, p_new_i)
            for pr in pairs:
                conv_rows.append({
                    "f_hz": f, "n_z_ref": 4001, "n_z_other": nz_other,
                    "i_ref": pr["i_ref"], "j_new": pr["j_new"],
                    "k_ref": pr["k_ref"], "k_new": pr["k_new"],
                    "delta_k": pr["delta_k"], "overlap": pr["overlap"],
                    "phi180_ref": float(np.interp(180, z4, p_ref[pr["i_ref"]])),
                    "phi180_new": float(np.interp(180, z4, p_new_i[pr["j_new"]])),
                    "phi200_ref": float(np.interp(200, z4, p_ref[pr["i_ref"]])),
                    "phi200_new": float(np.interp(200, z4, p_new_i[pr["j_new"]])),
                    "phi220_ref": float(np.interp(220, z4, p_ref[pr["i_ref"]])),
                    "phi220_new": float(np.interp(220, z4, p_new_i[pr["j_new"]])),
                })
    conv_df = pd.DataFrame(conv_rows)
    set_df = pd.DataFrame(set_rows)
    conv_df.to_csv(OUT / "effective_mode_convergence.csv", index=False, encoding="utf-8-sig")
    set_df.to_csv(OUT / "effective_mode_set_convergence.csv", index=False, encoding="utf-8-sig")

    # convergence quality
    eff_conv = True
    for f in FREQS:
        d2 = conv_df[(conv_df.f_hz == f) & (conv_df.n_z_other == 2001)]
        d8 = conv_df[(conv_df.f_hz == f) & (conv_df.n_z_other == 8001)]
        for d in [d2, d8]:
            if len(d) < 5:
                eff_conv = False
                continue
            if d["overlap"].median() < 0.9 or d["delta_k"].median() > 1e-3:
                eff_conv = False
        n2 = set_df[(set_df.f_hz == f) & (set_df.n_z == 2001)]["n_eff_0p2"].iloc[0]
        n4 = set_df[(set_df.f_hz == f) & (set_df.n_z == 4001)]["n_eff_0p2"].iloc[0]
        n8 = set_df[(set_df.f_hz == f) & (set_df.n_z == 8001)]["n_eff_0p2"].iloc[0]
        if abs(n4 - n2) > 3 or abs(n8 - n4) > 3:
            eff_conv = False
    print("effective_modes_converged", eff_conv, flush=True)

    # Controlled recovery + 3-way depth signatures
    rec_rows = []
    recL_rows = []
    sig_rows = []
    fr_rows = []
    for f in FREQS:
        z, modes, A, An, i1, i2 = store[(f, 4001)]
        # signed A for effective 0.2
        A_signed = A[i2]
        k_eff = [modes[i]["k_r"] for i in i2]
        # rayleigh resolvable at each L
        k_arr = np.array(k_eff)
        order = np.argsort(k_arr)[::-1]
        k_sorted = k_arr[order]
        idx_sorted = [i2[o] for o in order]
        for L in L_LIST:
            dkr = 2 * np.pi / L
            resolvable = []
            for a in range(len(k_sorted)):
                ok = True
                for b in range(len(k_sorted)):
                    if a != b and abs(k_sorted[a] - k_sorted[b]) < dkr:
                        ok = False
                        break
                if ok:
                    resolvable.append(idx_sorted[a])
            # SA recovery
            dr = C0 / f / 4.0
            rows_r, stats = sa_recover(np.array(k_eff), np.array(A_signed, dtype=float), L, dr, f)
            for rr in rows_r:
                rr.update({"f_hz": f, "L_m": L, "dr_name": "lam/4"})
                rec_rows.append(rr)
            recL_rows.append({
                "f_hz": f, "L_m": L,
                "v_r_T600": L / 600.0, "v_r_T1200": L / 1200.0,
                "n_eff_injected": stats["n_injected"],
                "n_rayleigh_resolvable": len(resolvable),
                **stats,
            })
            # recovered mode indices
            rec_idx = [int(r["mode_index"]) for r in rows_r if r["status"] in ("RESOLVED", "MERGED") and not np.isnan(r.get("mode_index", np.nan))]
            for kind, idxs in [
                ("ORACLE_FULL_MODE", i2),
                ("RAYLEIGH_RESOLVABLE_ONLY", resolvable),
                ("ACTUALLY_RECOVERED_ONLY", rec_idx),
            ]:
                for zi, zj in PAIRS:
                    D = depth_D(idxs, modes, z, zi, zj)
                    sig_rows.append({
                        "f_hz": f, "L_m": L, "kind": kind,
                        "z_i": zi, "z_j": zj, "n_modes": len(idxs), "D_abs": D,
                    })
            # aperture recovery already in recL_rows
        # frequency mechanism uses recovered @ L=2.4km
        sigL = [s for s in sig_rows if s["f_hz"] == f and s["L_m"] == 2.4e3 and s["kind"] == "ACTUALLY_RECOVERED_ONLY"]
        d220 = next((s["D_abs"] for s in sigL if s["z_i"] == 200 and s["z_j"] == 220), np.nan)
        d210 = next((s["D_abs"] for s in sigL if s["z_i"] == 200 and s["z_j"] == 210), np.nan)
        fr_rows.append({
            "f_hz": f, "C1_corr_200_220": C1_CORR[f],
            "rec_D_200_220_L2p4": d220, "rec_D_200_210_L2p4": d210,
            "oracle_D_200_220": depth_D(i2, modes, z, 200, 220),
            "oracle_D_200_210": depth_D(i2, modes, z, 200, 210),
        })

    rec_df = pd.DataFrame(rec_rows)
    recL_df = pd.DataFrame(recL_rows)
    sig_df = pd.DataFrame(sig_rows)
    fr_df = pd.DataFrame(fr_rows)
    rec_df.to_csv(OUT / "controlled_effective_mode_recovery.csv", index=False, encoding="utf-8-sig")
    recL_df.to_csv(OUT / "effective_mode_recovery_vs_L.csv", index=False, encoding="utf-8-sig")
    sig_df[sig_df.kind == "ORACLE_FULL_MODE"].to_csv(OUT / "depth_signature_oracle.csv", index=False, encoding="utf-8-sig")
    sig_df[sig_df.kind == "RAYLEIGH_RESOLVABLE_ONLY"].to_csv(OUT / "depth_signature_rayleigh.csv", index=False, encoding="utf-8-sig")
    sig_df[sig_df.kind == "ACTUALLY_RECOVERED_ONLY"].to_csv(OUT / "depth_signature_recovered.csv", index=False, encoding="utf-8-sig")
    fr_df.to_csv(OUT / "frequency_depth_mechanism_observable.csv", index=False, encoding="utf-8-sig")

    # 338 explanation using recovered D
    d338 = float(fr_df.loc[fr_df.f_hz == 338, "rec_D_200_220_L2p4"].iloc[0])
    dmed = float(np.nanmedian(fr_df["rec_D_200_220_L2p4"].to_numpy()))
    if np.isfinite(d338) and d338 > dmed + 0.08:
        c1 = "MODAL_STRUCTURE_SUPPORTED"
    else:
        c1 = "NOT_EXPLAINED"

    # decision
    def recD(f, L):
        s = sig_df[(sig_df.f_hz == f) & (sig_df.L_m == L) & (sig_df.kind == "ACTUALLY_RECOVERED_ONLY")]
        d = s[(s.z_i == 200) & (s.z_j == 220)]
        return float(d["D_abs"].iloc[0]) if len(d) else np.nan

    nrec24 = recL_df[recL_df.L_m == 2.4e3]
    enough = bool((nrec24["n_resolved"] >= 3).any())
    depth_ok = bool(np.nanmax([recD(f, 2.4e3) for f in FREQS]) > 0.08)
    freq_cond = False
    for f in FREQS:
        if recD(f, 2.4e3) > 0.08:
            freq_cond = True
    all_good = all(recD(f, 2.4e3) > 0.08 for f in FREQS)

    if not eff_conv:
        decision = "C2_0S_EFFECTIVE_MODES_NOT_CONVERGED"
        why = "effective mode set / phi(k) not stable 2001↔4001↔8001 after matching"
    elif enough and depth_ok and all_good:
        decision = "C2_0S_EFFECTIVE_MODE_PHYSICS_CONFIRMED"
        why = f"Eff. modes converged; ≥3 recovered at L=2.4km; recovered-mode D_200/220 clear for all f. C1_338={c1}"
    elif freq_cond and not all_good:
        decision = "C2_0S_FREQUENCY_CONDITIONAL"
        why = f"Only some frequencies have recovered-mode depth signature. C1_338={c1}"
    elif not enough:
        decision = "C2_0S_FOURIER_APERTURE_LIMITED"
        why = (
            f"Eff. modes converged but ordinary SA at L≤2.4km recovers few independent modes. "
            f"L=2.4 stats: {nrec24[['f_hz','n_resolved','fraction_resolved']].to_dict(orient='records')}. "
            f"Not a high-res Yang failure."
        )
    else:
        decision = "C2_0S_DEPTH_SIGNATURE_WEAK"
        why = f"Modes recovered but depth signature weak. C1_338={c1}"

    nxt = "stop; no full Yang / δ-z / Doppler / f0 / MC / P5"

    dec = {
        "rc3c2_0s_decision": decision,
        "why": why,
        "next_step": nxt,
        "C1_338": c1,
        "effective_modes_converged": bool(eff_conv),
        "oracle_D_note": "ORACLE_FULL_MODE_SIGNATURE only",
        "recovery_vs_L": recL_df.to_dict(orient="records"),
        "freq": fr_df.to_dict(orient="records"),
        "created_utc": NOW,
    }
    (OUT / "R3_C2_0S_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    rp = []
    rp.append("# R3-C2.0S 报告")
    rp.append("")
    rp.append(f"UTC：{NOW}")
    rp.append("")
    rp.append("冻结：Erratum δ、正确算子、ModeModel 伪影、Fourier=预关；旧 D≈0.9 → **ORACLE_FULL_MODE_SIGNATURE**。")
    rp.append("")
    rp.append("## 有效模态收敛")
    rp.append(f"- **{eff_conv}**（Hungarian/greedy 匹配后 overlap/Δk/集合规模）")
    rp.append(set_df.to_string(index=False))
    rp.append("")
    rp.append("## 有限孔径恢复（signed A_m，身份评估）")
    rp.append(recL_df.to_string(index=False))
    rp.append("")
    rp.append("## 三类深度签名 D(200/220) @ L=2.4 km")
    rp.append("")
    rp.append("| f | ORACLE | RAYLEIGH | RECOVERED | C1 corr |")
    rp.append("| --- | --- | --- | --- | --- |")
    for f in FREQS:
        o = depth_D(store[(f, 4001)][4], store[(f, 4001)][1], store[(f, 4001)][0], 200, 220)
        rp.append(f"| {f} | {o:.3f} | — | {recD(f, 2.4e3)} | {C1_CORR[f]} |")
    rp.append("")
    rp.append(f"C1_338 = **{c1}**")
    rp.append("")
    rp.append(f"## 判定 `{decision}`")
    rp.append("")
    rp.append(why)
    rp.append("")
    rp.append(f"下一步：{nxt}")
    (OUT / "R3_C2_0S_REPORT.md").write_text("\n".join(rp), encoding="utf-8")
    (OUT / "R3_C2_0S_GPT_SYNC.md").write_text(
        f"# R3-C2.0S\n\n**{decision}**\n\n{why}\n\nC1_338={c1}\neff_conv={eff_conv}\n",
        encoding="utf-8",
    )
    print("DECISION", decision)
    print("C1_338", c1)
    print(f"DONE {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
