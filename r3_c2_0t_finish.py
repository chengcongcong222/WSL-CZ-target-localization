#!/usr/bin/env python3
"""R3-C2.0T finisher: run remaining freqs with incremental save, then aggregate."""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd

import r3_c2_0s_effective_modes as S
from r3_c2_0t_modal_clusters import (
    OUT, NOW, FREQS, L_LIST, CLUSTER_ALPHA, PAIRS, Z_PROBE,
    C1_CORR, NZ_REF, NZ_FINE, Z_S, Z_R, C0,
    load_modes, cluster_by_dk, match_clusters_refine,
    cluster_signatures, D_pair, sa_peaks_to_clusters, sa_spectrum,
)

CACHE = OUT / "_cluster_cache"
CACHE.mkdir(exist_ok=True)


def save_parts(name, rows):
    if rows:
        pd.DataFrame(rows).to_csv(OUT / name, index=False, encoding="utf-8-sig")


def main():
    t0 = time.time()
    conv_rows, strength_rows, peak_rows, recL_rows = [], [], [], []
    sig_coh, sig_en = [], []

    for f in FREQS:
        fp = CACHE / f"modes_{int(f)}.npz"
        if fp.exists():
            z4 = np.load(fp)["z4"]
            # reload via solver still needed for phi lists — use pickle-like json too heavy
        z4, m4, A4 = load_modes(f, NZ_REF)
        z8, m8, A8 = load_modes(f, NZ_FINE)
        k4 = np.array([m["k_r"] for m in m4])
        k8 = np.array([m["k_r"] for m in m8])
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
                dfc, c4, c8 = match_clusters_refine(k4e, A4e, z4, phi4e, k8e, A8e, z8, phi8e, alpha, dkr)
                dfc["f_hz"] = f
                dfc["L_m"] = L
                dfc["alpha"] = alpha
                conv_rows.append(dfc)
                if alpha != 1.0:
                    continue
                for c in c4:
                    strength_rows.append({
                        "f_hz": f, "L_m": L, "cluster_id": c["cluster_id"],
                        "k_center": c["k_center"], "n_modes": c["n_modes"],
                        "E_C": c["E_C"], "E_C_norm": c["E_C_norm"],
                    })
                def to_global(c):
                    return {**c, "member_idx": [int(idx_eff[i]) for i in c["member_idx"]]}
                c4g = [to_global(c) for c in c4]
                dr = C0 / f / 4.0
                k_inj = np.array([c["k_center"] for c in c4])
                A_inj = np.array([sum(A4e[i] for i in c["member_idx"]) for c in c4])
                peaks = sa_spectrum(k_inj, A_inj, L, dr)
                pm = sa_peaks_to_clusters(k_inj, peaks, dkr)
                for row in pm:
                    row.update({"f_hz": f, "L_m": L, "alpha": alpha})
                    peak_rows.append(row)
                n_single = sum(1 for r in pm if r["status"] == "RESOLVED_SINGLE")
                n_clu = sum(1 for r in pm if r["status"] == "RESOLVED_CLUSTER")
                recL_rows.append({
                    "f_hz": f, "L_m": L, "v_r_T600": L / 600.0, "v_r_T1200": L / 1200.0,
                    "n_clusters_injected": len(c4),
                    "n_clusters_observed": n_single + n_clu,
                    "n_resolved_single": n_single,
                    "n_resolved_cluster": n_clu,
                    "n_spurious": sum(1 for r in pm if r["status"] == "SPURIOUS"),
                })
                sig_c = cluster_signatures(c4g, m4, z4, A4, "coherent")
                sig_e = cluster_signatures(c4g, m4, z4, A4, "energy")
                sig_om = cluster_signatures([{"member_idx": [int(i)]} for i in idx_eff], m4, z4, A4, "coherent")
                obs_ids = sorted(set(int(r["cluster_id"]) for r in pm
                                     if r["status"] in ("RESOLVED_SINGLE", "RESOLVED_CLUSTER")
                                     and r.get("cluster_id") is not None
                                     and not (isinstance(r.get("cluster_id"), float) and math.isnan(r["cluster_id"]))))
                c_obs = [c4g[i] for i in obs_ids if i < len(c4g)]
                sig_oc = cluster_signatures(c_obs, m4, z4, A4, "coherent") if c_obs else {z: np.array([0.0]) for z in Z_PROBE}
                sig_oe = cluster_signatures(c_obs, m4, z4, A4, "energy") if c_obs else {z: np.array([0.0]) for z in Z_PROBE}
                for zi, zj in PAIRS:
                    sig_coh.append({"f_hz": f, "L_m": L, "kind": "ORACLE_MODE", "z_i": zi, "z_j": zj,
                                    "D_abs": D_pair(sig_om, float(zi), float(zj))})
                    sig_coh.append({"f_hz": f, "L_m": L, "kind": "ORACLE_CLUSTER_COHERENT", "z_i": zi, "z_j": zj,
                                    "D_abs": D_pair(sig_c, float(zi), float(zj))})
                    sig_en.append({"f_hz": f, "L_m": L, "kind": "ORACLE_CLUSTER_ENERGY", "z_i": zi, "z_j": zj,
                                   "D_abs": D_pair(sig_e, float(zi), float(zj))})
                    sig_coh.append({"f_hz": f, "L_m": L, "kind": "OBSERVED_CLUSTER_COHERENT", "z_i": zi, "z_j": zj,
                                    "D_abs": D_pair(sig_oc, float(zi), float(zj))})
                    sig_en.append({"f_hz": f, "L_m": L, "kind": "OBSERVED_CLUSTER_ENERGY", "z_i": zi, "z_j": zj,
                                   "D_abs": D_pair(sig_oe, float(zi), float(zj))})
        # incremental
        save_parts("cluster_subspace_convergence.csv", pd.concat(conv_rows, ignore_index=True).to_dict("records"))
        save_parts("cluster_observable_strength.csv", strength_rows)
        save_parts("sa_peak_cluster_mapping.csv", peak_rows)
        save_parts("cluster_recovery_vs_L.csv", recL_rows)
        save_parts("depth_signature_cluster_coherent.csv", sig_coh)
        save_parts("depth_signature_cluster_energy.csv", sig_en)
        print(f"saved after f={f} {time.time()-t0:.0f}s", flush=True)

    recL = pd.DataFrame(recL_rows)
    sc = pd.DataFrame(sig_coh)
    se = pd.DataFrame(sig_en)
    conv_all = pd.concat(conv_rows, ignore_index=True)

    fr_rows = []
    for f in FREQS:
        def pick(df, kind, zi=200, zj=220):
            d = df[(df.f_hz == f) & (df.L_m == 2.4e3) & (df.kind == kind) & (df.z_i == zi) & (df.z_j == zj)]
            return float(d["D_abs"].iloc[0]) if len(d) else np.nan
        rr = recL[(recL.f_hz == f) & (recL.L_m == 2.4e3)]
        fr_rows.append({
            "f_hz": f, "C1_corr_200_220": C1_CORR[f],
            "D_obs_coh_200_220": pick(sc, "OBSERVED_CLUSTER_COHERENT"),
            "D_obs_energy_200_220": pick(se, "OBSERVED_CLUSTER_ENERGY"),
            "D_oracle_cluster_200_220": pick(sc, "ORACLE_CLUSTER_COHERENT"),
            "D_oracle_mode_200_220": pick(sc, "ORACLE_MODE"),
            "n_clusters_observed": float(rr["n_clusters_observed"].iloc[0]) if len(rr) else np.nan,
            "n_resolved_cluster": float(rr["n_resolved_cluster"].iloc[0]) if len(rr) else np.nan,
            "n_resolved_single": float(rr["n_resolved_single"].iloc[0]) if len(rr) else np.nan,
        })
    fr = pd.DataFrame(fr_rows)
    fr.to_csv(OUT / "frequency_cluster_mechanism.csv", index=False, encoding="utf-8-sig")

    c1a = conv_all[conv_all.alpha == 1.0]
    cluster_conv = True
    for f in FREQS:
        d = c1a[(c1a.f_hz == f) & (c1a.L_m == 2.4e3)]
        if len(d) < 2 or float(d["cluster_converged"].mean()) < 0.7:
            cluster_conv = False

    d338 = float(fr.loc[fr.f_hz == 338, "D_obs_coh_200_220"].iloc[0])
    dothers = fr.loc[fr.f_hz != 338, "D_obs_coh_200_220"].to_numpy(dtype=float)
    if np.isfinite(d338) and d338 > np.nanmax(dothers) + 0.08 and cluster_conv:
        c1_338 = "MODAL_CLUSTER_SUPPORTED"
    else:
        c1_338 = "NOT_EXPLAINED"

    nobs24 = recL[recL.L_m == 2.4e3]
    enough = bool((nobs24["n_clusters_observed"] >= 3).any())
    depth_ok = bool(np.nanmax(fr["D_obs_coh_200_220"].to_numpy()) > 0.08)
    all_f = bool((fr["D_obs_coh_200_220"].to_numpy() > 0.08).all())

    if not cluster_conv:
        decision = "C2_0T_CLUSTER_NOT_CONVERGED"
        why = f"Cluster subspaces 4001↔8001 not fully stable. C1_338={c1_338}"
    elif enough and depth_ok and all_f and c1_338 == "MODAL_CLUSTER_SUPPORTED":
        decision = "C2_0T_CLUSTER_PHYSICS_CONFIRMED"
        why = f"Clusters converged; observed clusters at L=2.4km with clear D; C1_338={c1_338}"
    elif depth_ok and not all_f:
        decision = "C2_0T_FREQUENCY_CONDITIONAL"
        why = f"Frequency-selective observed-cluster depth signature. C1_338={c1_338}"
    elif not enough:
        decision = "C2_0T_APERTURE_LIMITED"
        why = (
            f"Few observed clusters at L≤2.4km despite converged clusters. "
            f"{nobs24[['f_hz','n_clusters_injected','n_clusters_observed','n_resolved_single','n_resolved_cluster']].to_dict(orient='records')}. "
            f"C1_338={c1_338}"
        )
    else:
        decision = "C2_0T_DEPTH_SIGNATURE_WEAK"
        why = f"Clusters observed but D weak. C1_338={c1_338}"
    nxt = "stop; no full Yang / δ-z / Doppler / f0 / MC / P5"

    dec = {
        "rc3c2_0t_decision": decision, "why": why, "next_step": nxt,
        "C1_338": c1_338, "cluster_converged": bool(cluster_conv),
        "recovery_L2p4": nobs24.to_dict(orient="records"),
        "freq": fr.to_dict(orient="records"),
        "created_utc": NOW,
    }
    (OUT / "R3_C2_0T_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    rp = ["# R3-C2.0T 报告", "", f"UTC：{NOW}", "",
          "C1_338 先冻结为 CONSISTENT_WITH_MODAL_PHASE_STRUCTURE_NOT_YET_ATTRIBUTED。", "",
          f"簇子空间收敛：**{cluster_conv}**", "",
          "L=2.4 km：", nobs24.to_string(index=False), "",
          "频率：", fr.to_string(index=False), "",
          f"C1_338 = **{c1_338}**", "",
          f"## 判定 `{decision}`", "", why, "", f"下一步：{nxt}"]
    (OUT / "R3_C2_0T_REPORT.md").write_text("\n".join(rp), encoding="utf-8")
    (OUT / "R3_C2_0T_GPT_SYNC.md").write_text(
        f"# R3-C2.0T\n\n**{decision}**\n\n{why}\n\nC1_338={c1_338}\n",
        encoding="utf-8",
    )
    print("DECISION", decision)
    print("C1_338", c1_338)
    print(f"DONE {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
