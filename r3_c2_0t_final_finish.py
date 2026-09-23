#!/usr/bin/env python3
"""R3-C2.0T-FINAL finisher: 4001↔8001 energy-weighted, incremental, full decision."""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd

import r3_c2_0s_effective_modes as S
from r3_c2_0t_final import (
    OUT, NOW, FREQS, L_LIST, L_MAIN, PAIRS, Z_R, Z_TRUE, C1_CORR, C0, ALPHA,
    complete_link_clusters, cluster_convergence, sa_field_individual,
    peak_mode_groups, depth_sig, D,
)

NZ_LO, NZ_HI = 4001, 8001


def load(f, n_z):
    z = np.linspace(0.0, S.H, n_z)
    c = S.munk_c(z)
    omega = 2 * np.pi * f
    modes, z = S.solve_modes_sa(z, c, omega, omega / float(np.min(c)))
    A = np.array([
        float(np.interp(Z_TRUE, z, m["phi"])) * float(np.interp(Z_R, z, m["phi"]))
        for m in modes
    ])
    return z, modes, A


def main():
    t0 = time.time()
    cell_rows, conv_rows, spec_rows, map_rows = [], [], [], []
    recL_rows, d_oracle, d_actual, fr_rows, hf_rows = [], [], [], [], []

    for f in FREQS:
        z4, m4, A4 = load(f, NZ_LO)
        z8, m8, A8 = load(f, NZ_HI)
        An4 = np.abs(A4) / (np.max(np.abs(A4)) + 1e-30)
        idx = np.where(An4 > 0.2)[0]
        k4 = np.array([m4[i]["k_r"] for i in idx])
        A4e = A4[idx]
        phi4 = [m4[i]["phi"] for i in idx]
        An8 = np.abs(A8) / (np.max(np.abs(A8)) + 1e-30)
        idx8 = np.where(An8 > 0.2)[0]
        k8 = np.array([m8[i]["k_r"] for i in idx8])
        A8e = A8[idx8]
        phi8 = [m8[i]["phi"] for i in idx8]

        dkr = 2 * np.pi / L_MAIN
        for alpha in [0.5, 1.0]:
            c4 = complete_link_clusters(k4, A4e, alpha, dkr)
            c8 = complete_link_clusters(k8, A8e, alpha, dkr)
            for c in c4:
                cell_rows.append({**{k: v for k, v in c.items() if k != "member_idx"},
                                  "f_hz": f, "alpha": alpha, "grid": NZ_LO})
            dfc = cluster_convergence(c4, c8, k4, k8, z4, z8, phi4, phi8)
            if len(dfc):
                dfc["f_hz"] = f
                dfc["alpha"] = alpha
                conv_rows.append(dfc)
                jmax = int(dfc["E_C_norm_4001"].idxmax())
                hf_rows.append({
                    "f_hz": f, "alpha": alpha,
                    "dominant_cluster_overlap": float(dfc.loc[jmax, "subspace_overlap"]),
                    "dominant_cluster_angle": float(dfc.loc[jmax, "principal_angle_max_deg"]),
                    "dominant_cluster_energy": float(dfc.loc[jmax, "E_C_norm_4001"]),
                    "converged_energy_fraction": float(dfc["converged_energy_fraction"].iloc[0]),
                })
        c4 = complete_link_clusters(k4, A4e, ALPHA, dkr)
        cid = np.zeros(len(k4), dtype=int)
        for c in c4:
            for mi in c["member_idx"]:
                cid[mi] = c["cluster_id"]

        for L in L_LIST:
            dkrL = 2 * np.pi / L
            dr = C0 / f / 4.0
            kg, mag, _ = sa_field_individual(k4, A4e, L, dr)
            groups = peak_mode_groups(k4, A4e, kg, mag, dkrL)
            for g in groups:
                mem = g["member_mode_ids"]
                g = dict(g)
                g["f_hz"] = f
                g["L_m"] = L
                g["cluster_ids"] = json.dumps(sorted(set(int(cid[i]) for i in mem)))
                g["member_mode_ids"] = json.dumps(mem)
                map_rows.append(g)
            obs = [g for g in groups if g["status"] in ("RESOLVED_SINGLE_MODE", "RESOLVED_MODE_GROUP")]
            recL_rows.append({
                "f_hz": f, "L_m": L, "v_r_T600": L / 600.0, "v_r_T1200": L / 1200.0,
                "n_independent_observed_groups": len(obs),
                "n_single": sum(1 for g in groups if g["status"] == "RESOLVED_SINGLE_MODE"),
                "n_group": sum(1 for g in groups if g["status"] == "RESOLVED_MODE_GROUP"),
                "n_spurious": sum(1 for g in groups if g["status"] == "SPURIOUS"),
                "energy_coverage": float(sum(g["member_energy"] for g in obs) / (float(np.sum(A4e ** 2)) + 1e-30)),
            })
            def gl(mem_list):
                return [[int(idx[i]) for i in mem] for mem in mem_list]
            act_mem = []
            for g in obs:
                mem = json.loads(g["member_mode_ids"]) if isinstance(g["member_mode_ids"], str) else g["member_mode_ids"]
                act_mem.append(list(mem))
            sig_o = depth_sig(gl([[i] for i in range(len(k4))]), m4, z4, Z_R, "coherent")
            sig_c = depth_sig(gl([c["member_idx"] for c in c4]), m4, z4, Z_R, "coherent")
            sig_a = depth_sig(gl(act_mem), m4, z4, Z_R, "coherent") if act_mem else {200: np.array([0.0]), 210: np.array([0.0]), 220: np.array([0.0])}
            for zi, zj in PAIRS:
                d_oracle.append({"f_hz": f, "L_m": L, "kind": "ORACLE_FULL_MODE", "z_i": zi, "z_j": zj, "D_abs": D(sig_o, float(zi), float(zj))})
                d_oracle.append({"f_hz": f, "L_m": L, "kind": "ORACLE_RESOLUTION_CELL", "z_i": zi, "z_j": zj, "D_abs": D(sig_c, float(zi), float(zj))})
                d_actual.append({
                    "f_hz": f, "L_m": L, "kind": "ACTUAL_SA_PEAK_GROUP", "z_i": zi, "z_j": zj,
                    "D_abs": D(sig_a, float(zi), float(zj)), "n_groups": len(obs),
                    "scale_flag": "SCALE_UNIDENTIFIABLE_WITH_ONE_GROUP" if len(obs) < 2 else "ok",
                })
        # incremental save
        pd.DataFrame(cell_rows).to_csv(OUT / "resolution_cell_clusters.csv", index=False, encoding="utf-8-sig")
        pd.DataFrame(map_rows).to_csv(OUT / "sa_peak_mode_group_mapping.csv", index=False, encoding="utf-8-sig")
        pd.DataFrame(recL_rows).to_csv(OUT / "sa_group_recovery_vs_L.csv", index=False, encoding="utf-8-sig")
        pd.DataFrame(d_actual).to_csv(OUT / "depth_signature_actual_groups.csv", index=False, encoding="utf-8-sig")
        print(f"saved f={f} {time.time()-t0:.0f}s", flush=True)

    conv_all = pd.concat(conv_rows, ignore_index=True) if conv_rows else pd.DataFrame()
    conv_all.to_csv(OUT / "cluster_energy_convergence_final.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(hf_rows).to_csv(OUT / "high_f_mode_convergence.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(d_oracle).to_csv(OUT / "depth_signature_oracle_final.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(spec_rows).to_csv(OUT / "sa_individual_mode_spectrum.csv", index=False, encoding="utf-8-sig")

    recL = pd.DataFrame(recL_rows)
    dact = pd.DataFrame(d_actual)
    conv_e = {}
    for f in FREQS:
        d = conv_all[(conv_all.f_hz == f) & (conv_all.alpha == ALPHA)] if len(conv_all) else pd.DataFrame()
        conv_e[f] = float(d["converged_energy_fraction"].iloc[0]) if len(d) else 0.0

    fr_rows = []
    for f in FREQS:
        da = dact[(dact.f_hz == f) & (dact.L_m == L_MAIN) & (dact.z_i == 200) & (dact.z_j == 220)]
        rr = recL[(recL.f_hz == f) & (recL.L_m == L_MAIN)]
        hf = [h for h in hf_rows if h["f_hz"] == f and h["alpha"] == ALPHA]
        fr_rows.append({
            "f_hz": f, "C1_corr_200_220": C1_CORR[f],
            "D_actual_200_220": float(da["D_abs"].iloc[0]) if len(da) else np.nan,
            "D_actual_200_210": float(dact[(dact.f_hz == f) & (dact.L_m == L_MAIN) & (dact.z_i == 200) & (dact.z_j == 210)]["D_abs"].iloc[0]) if len(dact) else np.nan,
            "n_obs_groups_L24": int(rr["n_independent_observed_groups"].iloc[0]) if len(rr) else 0,
            "energy_cov_L24": float(rr["energy_coverage"].iloc[0]) if len(rr) else np.nan,
            "conv_energy_frac": conv_e.get(f, 0.0),
            "dom_overlap": hf[0]["dominant_cluster_overlap"] if hf else np.nan,
        })
    fr = pd.DataFrame(fr_rows)
    fr.to_csv(OUT / "frequency_group_mechanism_final.csv", index=False, encoding="utf-8-sig")

    high_f_ok = all(conv_e.get(f, 0) >= 0.9 for f in [283.0, 338.0])
    multi = {f: int(rr) for f, rr in zip(fr.f_hz, fr.n_obs_groups_L24)}
    d220 = dict(zip(fr.f_hz, fr.D_actual_200_220))
    d210 = dict(zip(fr.f_hz, fr.D_actual_200_210))

    if not high_f_ok:
        decision = "C2_0T_HIGH_F_MODE_NOT_CONVERGED"
        why = f"283/338 energy-weighted convergence <0.9: {conv_e}. Use mature normal-mode tools; no more solver R&D."
        nxt = "stop self-built solver"
    elif all(multi.get(f, 0) >= 2 for f in FREQS) and all((d220.get(f, 0) or 0) > 0.08 for f in FREQS):
        decision = "C2_0_FINAL_SA_PHYSICS_CONFIRMED"
        why = f"≥2 groups/f at L=2.4km and clear D. groups={multi} D220={d220}"
        nxt = "eligible for formal Yang estimator"
    elif any(multi.get(f, 0) >= 2 and (d220.get(f, 0) or 0) > 0.08 for f in FREQS):
        decision = "C2_0_FINAL_FREQUENCY_CONDITIONAL"
        why = f"Frequency-conditional. groups={multi} D220={d220} D210={d210}"
        nxt = "Yang only on favorable f if at all"
    elif all(multi.get(f, 0) < 2 for f in FREQS):
        decision = "C2_0_FINAL_APERTURE_LIMITED"
        why = f"<2 independent groups at 2.4km all f. groups={multi}. Close Fourier-SA Yang; go high-res AR/f-k."
        nxt = "close Fourier-SA Yang; next high-res modal (not this round)"
    else:
        decision = "C2_0_FINAL_DEPTH_SIGNATURE_WEAK"
        why = f"Groups exist but D weak. D220={d220} D210={d210}"
        nxt = "close Yang depth"

    d338 = d220.get(338.0, np.nan)
    dref = max([d220.get(f, np.nan) for f in [201.0, 235.0, 283.0]])
    c1_338 = "C1_338_MODAL_GROUP_SUPPORTED" if (
        multi.get(338.0, 0) >= 2 and np.isfinite(d338) and d338 > (dref if np.isfinite(dref) else 0) + 0.08
    ) else "C1_338_NOT_EXPLAINED"

    dec = {
        "rc3c2_0t_final_decision": decision, "why": why, "next_step": nxt,
        "C1_338": c1_338,
        "energy_weighted_convergence": conv_e,
        "n_groups_L24": multi, "D200_220_actual": d220, "D200_210_actual": d210,
        "created_utc": NOW, "no_more_C2_0x": True,
        "note": "high-f compare 4001↔8001 (12001/16001 dense skipped for runtime; same FEM operator)",
    }
    (OUT / "R3_C2_0T_FINAL_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    rp = ["# R3-C2.0T-FINAL", "", f"UTC：{NOW}", "",
          "修正：complete-link、单模态SA、A_m(z)=φ(z)φ(z_r)、能量加权收敛、单组 scale 不可辨识。", "",
          f"收敛：{conv_e}", "",
          recL[recL.L_m == L_MAIN].to_string(index=False), "",
          fr.to_string(index=False), "",
          f"C1_338=**{c1_338}**", "",
          f"## `{decision}`", "", why, "", f"下一步：{nxt}", "",
          "C2.0 最后预关，不再 C2.0x。"]
    (OUT / "R3_C2_0T_FINAL_REPORT.md").write_text("\n".join(rp), encoding="utf-8")
    (OUT / "R3_C2_0T_FINAL_GPT_SYNC.md").write_text(
        f"# R3-C2.0T-FINAL\n\n**{decision}**\n\n{why}\n\nC1_338={c1_338}\n",
        encoding="utf-8",
    )
    print("DECISION", decision)
    print("C1_338", c1_338)
    print(f"DONE {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
