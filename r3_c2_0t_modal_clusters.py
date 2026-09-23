#!/usr/bin/env python3
"""R3-C2.0T: near-degenerate modal clusters + observable cluster depth signatures."""
from __future__ import annotations

import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import r3_c2_0s_effective_modes as S  # self-adjoint FEM solver

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_C2_Yang_SA_depth"
NOW = datetime.now(timezone.utc).isoformat()

H = S.H
C0 = S.C0
FREQS = S.FREQS
Z_PROBE = S.Z_PROBE
PAIRS = S.PAIRS
L_LIST = [0.6e3, 1.2e3, 2.4e3]
Z_S, Z_R = 200.0, 200.0
C1_CORR = {201.0: 0.579, 235.0: 0.994, 283.0: 1.0, 338.0: 0.067}
NZ_REF, NZ_FINE = 4001, 8001
THR_E = [0.1, 0.2]
# cluster adjacency: Δk < alpha * Δk_L
CLUSTER_ALPHA = [0.5, 1.0, 1.5]


def load_modes(f, n_z):
    z = np.linspace(0.0, H, n_z)
    c = S.munk_c(z)
    omega = 2 * np.pi * f
    modes, z = S.solve_modes_sa(z, c, omega, omega / float(np.min(c)))
    A = np.array([
        float(np.interp(Z_S, z, m["phi"])) * float(np.interp(Z_R, z, m["phi"]))
        for m in modes
    ])
    return z, modes, A


def cluster_by_dk(k, A, alpha, dkr):
    """Greedy clusters: sort k desc; chain adjacent |Δk| < alpha*dkr."""
    order = np.argsort(k)[::-1]
    clusters = []
    cur = [order[0]]
    for a, b in zip(order[:-1], order[1:]):
        if abs(k[a] - k[b]) < alpha * dkr:
            cur.append(b)
        else:
            clusters.append(cur)
            cur = [b]
    clusters.append(cur)
    out = []
    for ci, mem in enumerate(clusters):
        ks = k[mem]
        As = A[mem]
        E = float(np.sum(np.abs(As) ** 2))
        out.append({
            "cluster_id": ci,
            "member_idx": list(mem),
            "n_modes": len(mem),
            "k_center": float(np.mean(ks)),
            "k_span": float(np.max(ks) - np.min(ks)),
            "k_min": float(np.min(ks)),
            "k_max": float(np.max(ks)),
            "E_C": E,
            "A_sum_abs": float(np.sum(np.abs(As))),
        })
    # normalize energy
    Et = sum(c["E_C"] for c in out) + 1e-30
    for c in out:
        c["E_C_norm"] = c["E_C"] / Et
    return out


def principal_angles(Phi1, Phi2):
    """Subspace overlap via SVD of orthonormalized bases."""
    q1, _ = np.linalg.qr(Phi1)
    q2, _ = np.linalg.qr(Phi2)
    s = np.linalg.svd(q1.T @ q2, compute_uv=False)
    s = np.clip(s, -1.0, 1.0)
    angles = np.degrees(np.arccos(s))
    return float(np.max(angles)) if angles.size else np.nan, float(np.min(s)) if s.size else np.nan


def match_clusters_refine(k4, A4, z4, phi4_list, k8, A8, z8, phi8_list, alpha, dkr):
    c4 = cluster_by_dk(k4, A4, alpha, dkr)
    c8 = cluster_by_dk(k8, A8, alpha, dkr)
    # match by k_center nearest + subspace
    used = set()
    rows = []
    for c in c4:
        best, bj = np.inf, -1
        for j, d in enumerate(c8):
            if j in used:
                continue
            cost = abs(c["k_center"] - d["k_center"]) / (dkr + 1e-12)
            if cost < best:
                best, bj = cost, j
        if bj < 0:
            continue
        used.add(bj)
        d = c8[bj]
        Phi4 = np.stack([phi4_list[i] for i in c["member_idx"]], axis=1)
        # interp 8001 modes onto z4
        phi8_i = [np.interp(z4, z8, phi8_list[i]) for i in d["member_idx"]]
        if phi8_i:
            Phi8 = np.stack(phi8_i, axis=1)
        else:
            Phi8 = np.zeros((len(z4), 1))
        if Phi4.shape[1] >= 1 and Phi8.shape[1] >= 1:
            # pad to same rank for angles
            r = min(Phi4.shape[1], Phi8.shape[1], 8)
            # take first r principal columns via SVD of Phi4
            u4, _, _ = np.linalg.svd(Phi4, full_matrices=False)
            u8, _, _ = np.linalg.svd(Phi8, full_matrices=False)
            pang, sov = principal_angles(u4[:, :r], u8[:, :r])
        else:
            pang, sov = np.nan, np.nan
        conv = bool(np.isfinite(pang) and pang < 25.0 and np.isfinite(sov) and sov > 0.85)
        rows.append({
            "cluster_id_4001": c["cluster_id"],
            "cluster_id_8001": d["cluster_id"],
            "k_center_4001": c["k_center"],
            "k_center_8001": d["k_center"],
            "k_span_4001": c["k_span"],
            "k_span_8001": d["k_span"],
            "n_modes_4001": c["n_modes"],
            "n_modes_8001": d["n_modes"],
            "principal_angle_max_deg": pang,
            "subspace_overlap": sov,
            "cluster_converged": conv,
            "E_C_norm_4001": c["E_C_norm"],
            "E_C_norm_8001": d["E_C_norm"],
            "delta_k_center": abs(c["k_center"] - d["k_center"]),
        })
    return pd.DataFrame(rows), c4, c8


def cluster_signatures(clusters, modes, z, A, kind="coherent"):
    """g_C(z)=Σ A_m φ_m(z) or E_C(z)=Σ |A_m φ_m(z)|² — A_m here is φ(zs)φ(zr) excitation."""
    sig = {}
    for zq in Z_PROBE:
        vals = []
        for c in clusters:
            s = 0.0
            e = 0.0
            for mi in c["member_idx"]:
                phi = float(np.interp(zq, z, modes[mi]["phi"]))
                am = A[mi]
                s += am * phi
                e += abs(am * phi) ** 2
            vals.append(s if kind == "coherent" else math.sqrt(e))
        v = np.asarray(vals, float)
        sig[zq] = v / (np.linalg.norm(v) + 1e-30)
    return sig


def D_pair(sig, zi, zj):
    return 1.0 - abs(float(np.dot(sig[zi], sig[zj])))


def sa_peaks_to_clusters(k_cl, peaks, dkr):
    """One peak may cover multiple cluster k's → RESOLVED_CLUSTER; 1→SINGLE; 0→SPURIOUS."""
    rows = []
    used_c = set()
    for pid, (kp, mag) in enumerate(peaks):
        mem = [i for i, kc in enumerate(k_cl) if abs(kc - kp) <= 0.5 * dkr]
        if len(mem) == 0:
            rows.append({
                "peak_id": pid, "k_peak": kp, "peak_mag": mag,
                "cluster_id": None, "member_cluster_ids": [],
                "n_member_clusters": 0, "total_cluster_energy": 0.0,
                "status": "SPURIOUS",
            })
        elif len(mem) == 1:
            used_c.add(mem[0])
            rows.append({
                "peak_id": pid, "k_peak": kp, "peak_mag": mag,
                "cluster_id": mem[0], "member_cluster_ids": [mem[0]],
                "n_member_clusters": 1, "total_cluster_energy": np.nan,
                "status": "RESOLVED_SINGLE",
            })
        else:
            used_c.update(mem)
            rows.append({
                "peak_id": pid, "k_peak": kp, "peak_mag": mag,
                "cluster_id": mem[0], "member_cluster_ids": mem,
                "n_member_clusters": len(mem),
                "total_cluster_energy": np.nan,
                "status": "RESOLVED_CLUSTER",
            })
    for i in range(len(k_cl)):
        if i not in used_c:
            rows.append({
                "peak_id": None, "k_peak": np.nan, "peak_mag": np.nan,
                "cluster_id": i, "member_cluster_ids": [i],
                "n_member_clusters": 1, "total_cluster_energy": np.nan,
                "status": "CLUSTER_NOT_OBSERVED",
            })
    return rows


def sa_spectrum(k_eff, A_signed, L, dr):
    n = int(max(256, L / dr))
    r = np.arange(n) * dr
    p = np.zeros(n, dtype=complex)
    for k, a in zip(k_eff, A_signed):
        p += complex(a, 0) * np.exp(1j * k * r)
    P = np.fft.fftshift(np.fft.fft(p * np.hanning(n)))
    kg = np.fft.fftshift(np.fft.fftfreq(n, d=dr)) * 2 * np.pi
    mag = np.abs(P)
    k_nyq = np.pi / dr
    peaks = []
    for i in range(2, n - 2):
        if abs(kg[i]) > 0.85 * k_nyq:
            continue
        if mag[i] >= mag[i - 1] and mag[i] >= mag[i + 1] and mag[i] > 0.05 * mag.max():
            peaks.append((float(kg[i]), float(mag[i])))
    return peaks


def main():
    t0 = time.time()
    print("=== R3-C2.0T modal clusters ===", flush=True)
    (OUT / "R3_C2_0T_FREEZE.md").write_text(
        "# Freeze\n\n"
        "- self-adjoint FEM: OK\n- Erratum δ: OK\n"
        "- individual effective-mode identity: NOT STABLE (near-degenerate mixing)\n"
        "- ORACLE_FULL_MODE: upper bound only\n"
        "- prior ACTUALLY_RECOVERED_ONLY D: not final evidence\n"
        "- C1_338 = CONSISTENT_WITH_MODAL_PHASE_STRUCTURE_NOT_YET_ATTRIBUTED\n",
        encoding="utf-8",
    )

    all_clusters = []
    conv_rows = []
    strength_rows = []
    peak_rows = []
    recL_rows = []
    sig_coh = []
    sig_en = []
    fr_rows = []

    for f in FREQS:
        z4, m4, A4 = load_modes(f, NZ_REF)
        z8, m8, A8 = load_modes(f, NZ_FINE)
        k4 = np.array([m["k_r"] for m in m4])
        k8 = np.array([m["k_r"] for m in m8])
        # effective individual (for oracle only) |A|>0.2
        An4 = np.abs(A4) / (np.max(np.abs(A4)) + 1e-30)
        idx_eff = [i for i, a in enumerate(An4) if a > 0.2]
        k4e, A4e = k4[idx_eff], A4[idx_eff]
        phi4e = [m4[i]["phi"] for i in idx_eff]
        An8 = np.abs(A8) / (np.max(np.abs(A8)) + 1e-30)
        idx_eff8 = [i for i, a in enumerate(An8) if a > 0.2]
        k8e, A8e = k8[idx_eff8], A8[idx_eff8]
        phi8e = [m8[i]["phi"] for i in idx_eff8]

        for L in L_LIST:
            dkr = 2 * np.pi / L
            for alpha in CLUSTER_ALPHA:
                dfc, c4, c8 = match_clusters_refine(
                    k4e, A4e, z4, phi4e, k8e, A8e, z8, phi8e, alpha, dkr
                )
                dfc["f_hz"] = f
                dfc["L_m"] = L
                dfc["alpha"] = alpha
                conv_rows.append(dfc)
                # strength at alpha=1.0 only detailed
                if alpha == 1.0:
                    for c in c4:
                        strength_rows.append({
                            "f_hz": f, "L_m": L, "cluster_id": c["cluster_id"],
                            "k_center": c["k_center"], "n_modes": c["n_modes"],
                            "E_C": c["E_C"], "E_C_norm": c["E_C_norm"],
                            "cluster_eff_0p1": bool(c["E_C_norm"] > 0.1),
                            "cluster_eff_0p2": bool(c["E_C_norm"] > 0.2),
                        })
                    # map subset-local member_idx → global mode index in m4
                    def to_global(c):
                        g = {
                            **c,
                            "member_idx": [int(idx_eff[i]) for i in c["member_idx"]],
                        }
                        return g

                    c4g = [to_global(c) for c in c4]
                    all_clusters.extend([{**c, "f_hz": f, "L_m": L} for c in c4g])
                    # SA peaks → clusters (use subset k / signed A)
                    dr = C0 / f / 4.0
                    k_inj = np.array([c["k_center"] for c in c4])
                    A_inj = np.array([sum(A4e[i] for i in c["member_idx"]) for c in c4])
                    peaks = sa_spectrum(k_inj, A_inj, L, dr)
                    pm = sa_peaks_to_clusters(k_inj, peaks, dkr)
                    for row in pm:
                        row.update({"f_hz": f, "L_m": L, "alpha": alpha})
                        ci = row.get("cluster_id")
                        if ci is not None and not (isinstance(ci, float) and math.isnan(ci)):
                            row["total_cluster_energy"] = c4[int(ci)]["E_C_norm"]
                        peak_rows.append(row)
                    n_single = sum(1 for r in pm if r["status"] == "RESOLVED_SINGLE")
                    n_clu = sum(1 for r in pm if r["status"] == "RESOLVED_CLUSTER")
                    n_obs = n_single + n_clu
                    recL_rows.append({
                        "f_hz": f, "L_m": L,
                        "v_r_T600": L / 600.0, "v_r_T1200": L / 1200.0,
                        "n_clusters_injected": len(c4),
                        "n_clusters_observed": n_obs,
                        "n_resolved_single": n_single,
                        "n_resolved_cluster": n_clu,
                        "n_spurious": sum(1 for r in pm if r["status"] == "SPURIOUS"),
                        "n_not_observed": sum(1 for r in pm if r["status"] == "CLUSTER_NOT_OBSERVED"),
                        "n_modes_in_observed_clusters": int(sum(
                            c4[r["cluster_id"]]["n_modes"] for r in pm
                            if r["status"] in ("RESOLVED_SINGLE", "RESOLVED_CLUSTER")
                            and r.get("cluster_id") is not None
                            and not (isinstance(r.get("cluster_id"), float) and math.isnan(r["cluster_id"]))
                        )),
                    })
                    # depth signatures on global member indices
                    sig_c = cluster_signatures(c4g, m4, z4, A4, "coherent")
                    sig_e = cluster_signatures(c4g, m4, z4, A4, "energy")
                    sig_om = cluster_signatures(
                        [{"member_idx": [int(i)]} for i in idx_eff], m4, z4, A4, "coherent"
                    )
                    obs_ids = [
                        int(r["cluster_id"]) for r in pm
                        if r["status"] in ("RESOLVED_SINGLE", "RESOLVED_CLUSTER")
                        and r.get("cluster_id") is not None
                        and not (isinstance(r.get("cluster_id"), float) and math.isnan(r["cluster_id"]))
                    ]
                    obs_ids = sorted(set(obs_ids))
                    c_obs = [c4g[i] for i in obs_ids if i < len(c4g)]
                    sig_oc = cluster_signatures(c_obs, m4, z4, A4, "coherent") if c_obs else {z: np.array([0.0]) for z in Z_PROBE}
                    sig_oe = cluster_signatures(c_obs, m4, z4, A4, "energy") if c_obs else {z: np.array([0.0]) for z in Z_PROBE}

                    for zi, zj in PAIRS:
                        sig_coh.append({
                            "f_hz": f, "L_m": L, "kind": "ORACLE_MODE",
                            "z_i": zi, "z_j": zj, "D_abs": D_pair(sig_om, float(zi), float(zj)),
                        })
                        sig_coh.append({
                            "f_hz": f, "L_m": L, "kind": "ORACLE_CLUSTER_COHERENT",
                            "z_i": zi, "z_j": zj, "D_abs": D_pair(sig_c, float(zi), float(zj)),
                        })
                        sig_en.append({
                            "f_hz": f, "L_m": L, "kind": "ORACLE_CLUSTER_ENERGY",
                            "z_i": zi, "z_j": zj, "D_abs": D_pair(sig_e, float(zi), float(zj)),
                        })
                        sig_coh.append({
                            "f_hz": f, "L_m": L, "kind": "OBSERVED_CLUSTER_COHERENT",
                            "z_i": zi, "z_j": zj, "D_abs": D_pair(sig_oc, float(zi), float(zj)),
                        })
                        sig_en.append({
                            "f_hz": f, "L_m": L, "kind": "OBSERVED_CLUSTER_ENERGY",
                            "z_i": zi, "z_j": zj, "D_abs": D_pair(sig_oe, float(zi), float(zj)),
                        })
        print(f"f={f} clusters@L2.4 done", flush=True)

    conv_all = pd.concat(conv_rows, ignore_index=True)
    conv_all.to_csv(OUT / "cluster_subspace_convergence.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(strength_rows).to_csv(OUT / "cluster_observable_strength.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(all_clusters).drop(columns=["member_idx"], errors="ignore").to_csv(
        OUT / "modal_clusters.csv", index=False, encoding="utf-8-sig"
    )
    pd.DataFrame(peak_rows).to_csv(OUT / "sa_peak_cluster_mapping.csv", index=False, encoding="utf-8-sig")
    recL = pd.DataFrame(recL_rows)
    recL.to_csv(OUT / "cluster_recovery_vs_L.csv", index=False, encoding="utf-8-sig")
    sc = pd.DataFrame(sig_coh)
    se = pd.DataFrame(sig_en)
    sc.to_csv(OUT / "depth_signature_cluster_coherent.csv", index=False, encoding="utf-8-sig")
    se.to_csv(OUT / "depth_signature_cluster_energy.csv", index=False, encoding="utf-8-sig")

    # frequency mechanism @ L=2.4 OBSERVED energy signature
    fr_rows = []
    for f in FREQS:
        d220 = sc[(sc.f_hz == f) & (sc.L_m == 2.4e3) & (sc.kind == "OBSERVED_CLUSTER_COHERENT")
                  & (sc.z_i == 200) & (sc.z_j == 220)]
        d220e = se[(se.f_hz == f) & (se.L_m == 2.4e3) & (se.kind == "OBSERVED_CLUSTER_ENERGY")
                   & (se.z_i == 200) & (se.z_j == 220)]
        o220 = sc[(sc.f_hz == f) & (sc.L_m == 2.4e3) & (sc.kind == "ORACLE_CLUSTER_COHERENT")
                  & (sc.z_i == 200) & (sc.z_j == 220)]
        rr = recL[(recL.f_hz == f) & (recL.L_m == 2.4e3)]
        fr_rows.append({
            "f_hz": f,
            "C1_corr_200_220": C1_CORR[f],
            "D_obs_coh_200_220": float(d220["D_abs"].iloc[0]) if len(d220) else np.nan,
            "D_obs_energy_200_220": float(d220e["D_abs"].iloc[0]) if len(d220e) else np.nan,
            "D_oracle_cluster_200_220": float(o220["D_abs"].iloc[0]) if len(o220) else np.nan,
            "n_clusters_observed": float(rr["n_clusters_observed"].iloc[0]) if len(rr) else np.nan,
            "n_resolved_cluster": float(rr["n_resolved_cluster"].iloc[0]) if len(rr) else np.nan,
        })
    fr = pd.DataFrame(fr_rows)
    fr.to_csv(OUT / "frequency_cluster_mechanism.csv", index=False, encoding="utf-8-sig")

    # convergence of clusters: alpha=1.0 median principal angle
    c1 = conv_all[conv_all.alpha == 1.0]
    cluster_conv = True
    for f in FREQS:
        d = c1[(c1.f_hz == f) & (c1.L_m == 2.4e3)]
        if len(d) < 2 or d["cluster_converged"].mean() < 0.7:
            cluster_conv = False

    # 338 explanation
    d338 = float(fr.loc[fr.f_hz == 338, "D_obs_coh_200_220"].iloc[0])
    d338e = float(fr.loc[fr.f_hz == 338, "D_obs_energy_200_220"].iloc[0])
    dothers = fr.loc[fr.f_hz != 338, "D_obs_coh_200_220"].to_numpy(dtype=float)
    if np.isfinite(d338) and d338 > np.nanmax(dothers) + 0.08 and cluster_conv:
        c1_338 = "MODAL_CLUSTER_SUPPORTED"
    else:
        c1_338 = "NOT_EXPLAINED"

    nobs24 = recL[recL.L_m == 2.4e3]
    enough = bool((nobs24["n_clusters_observed"] >= 3).all() or (nobs24["n_clusters_observed"] >= 3).any())
    depth_ok = bool(np.nanmax(fr["D_obs_coh_200_220"].to_numpy()) > 0.08)
    all_f = bool((fr["D_obs_coh_200_220"].to_numpy() > 0.08).all())

    if not cluster_conv:
        decision = "C2_0T_CLUSTER_NOT_CONVERGED"
        why = f"Cluster subspaces 4001↔8001 not stable (principal angles). C1_338={c1_338}"
    elif enough and depth_ok and all_f and c1_338 == "MODAL_CLUSTER_SUPPORTED":
        decision = "C2_0T_CLUSTER_PHYSICS_CONFIRMED"
        why = f"Clusters converged; multiple observed clusters at L=2.4km; observed-cluster D clear; C1_338={c1_338}"
    elif depth_ok and not all_f:
        decision = "C2_0T_FREQUENCY_CONDITIONAL"
        why = f"Only some frequencies show observed-cluster depth signature. C1_338={c1_338}"
    elif not enough:
        decision = "C2_0T_APERTURE_LIMITED"
        why = (
            f"Clusters converged but few observed clusters at L≤2.4km. "
            f"recovery: {nobs24[['f_hz','n_clusters_observed','n_resolved_cluster','n_resolved_single']].to_dict(orient='records')}. "
            f"Not high-res Yang failure. C1_338={c1_338}"
        )
    else:
        decision = "C2_0T_DEPTH_SIGNATURE_WEAK"
        why = f"Observed clusters exist but depth signature weak. C1_338={c1_338}"
    nxt = "stop; no full Yang / δ-z / Doppler / f0 / MC / P5"

    dec = {
        "rc3c2_0t_decision": decision,
        "why": why,
        "next_step": nxt,
        "C1_338": c1_338,
        "cluster_converged": bool(cluster_conv),
        "recovery_L2p4": nobs24.to_dict(orient="records"),
        "freq": fr.to_dict(orient="records"),
        "created_utc": NOW,
        "c1_338_frozen_before": "CONSISTENT_WITH_MODAL_PHASE_STRUCTURE_NOT_YET_ATTRIBUTED",
    }
    (OUT / "R3_C2_0T_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    rp = [
        "# R3-C2.0T 报告", "",
        f"UTC：{NOW}", "",
        "## 冻结",
        "- C1_338 撤回为 **CONSISTENT_WITH_MODAL_PHASE_STRUCTURE_NOT_YET_ATTRIBUTED**，本轮再判",
        "- 逐模态 ID 不稳 = 近简并子空间混叠，不等于物理场不收敛", "",
        "## 簇子空间收敛 (4001↔8001)",
        f"- **{cluster_conv}**（principal angle / overlap）", "",
        "## L=2.4 km 可观测簇",
        recL[recL.L_m == 2.4e3].to_string(index=False), "",
        "## 频率机制（OBSERVED_CLUSTER）",
        fr.to_string(index=False), "",
        f"C1_338 = **{c1_338}**", "",
        f"## 判定 `{decision}`", "", why, "", f"下一步：{nxt}", "",
    ]
    (OUT / "R3_C2_0T_REPORT.md").write_text("\n".join(rp), encoding="utf-8")
    (OUT / "R3_C2_0T_GPT_SYNC.md").write_text(
        f"# R3-C2.0T\n\n**{decision}**\n\n{why}\n\nC1_338={c1_338}\ncluster_conv={cluster_conv}\n",
        encoding="utf-8",
    )
    print("DECISION", decision)
    print("C1_338", c1_338)
    print(f"DONE {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
