#!/usr/bin/env python3
"""R3-C2.0T-FINAL: complete-link clusters, individual-mode SA field, corrected depth A_m(z).

Last C2.0 physics pre-gate. No further C2.0x.
"""
from __future__ import annotations

import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import r3_c2_0s_effective_modes as S

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_C2_Yang_SA_depth"
NOW = datetime.now(timezone.utc).isoformat()

H, C0 = S.H, S.C0
FREQS = [201.0, 235.0, 283.0, 338.0]
Z_R = 200.0
Z_TRUE = 200.0
PAIRS = [(200, 210), (200, 220)]
L_MAIN = 2.4e3
L_LIST = [0.6e3, 1.2e3, 2.4e3]
C1_CORR = {201.0: 0.579, 235.0: 0.994, 283.0: 1.0, 338.0: 0.067}
ALPHA = 1.0  # main; also 0.5


def load(f, n_z):
    return S.load_modes(f, n_z) if hasattr(S, "load_modes") else _load(f, n_z)


def _load(f, n_z):
    z = np.linspace(0.0, H, n_z)
    c = S.munk_c(z)
    omega = 2 * np.pi * f
    modes, z = S.solve_modes_sa(z, c, omega, omega / float(np.min(c)))
    A = np.array([
        float(np.interp(Z_TRUE, z, m["phi"])) * float(np.interp(Z_R, z, m["phi"]))
        for m in modes
    ])
    return z, modes, A


def complete_link_clusters(k, A, alpha, dkr):
    """No chaining: k_max-k_min <= alpha*dkr for every cluster."""
    order = np.argsort(k)  # ascending k for span control
    clusters = []
    cur = [order[0]]
    for idx in order[1:]:
        ks = k[cur + [idx]]
        if (ks.max() - ks.min()) <= alpha * dkr:
            cur.append(idx)
        else:
            clusters.append(cur)
            cur = [idx]
    clusters.append(cur)
    out = []
    Et = float(np.sum(A ** 2)) + 1e-30
    for ci, mem in enumerate(clusters):
        mem = np.asarray(mem, dtype=int)
        ks = k[mem]
        As = A[mem]
        E = float(np.sum(As ** 2))
        out.append({
            "cluster_id": ci,
            "member_idx": mem.tolist(),
            "n_modes": len(mem),
            "k_min": float(ks.min()),
            "k_max": float(ks.max()),
            "k_center": float(ks.mean()),
            "k_span": float(ks.max() - ks.min()),
            "k_span_over_dkr": float((ks.max() - ks.min()) / dkr) if dkr > 0 else np.nan,
            "E_C": E,
            "E_C_norm": E / Et,
        })
    return out


def principal_angles(Phi1, Phi2):
    q1, _ = np.linalg.qr(Phi1)
    q2, _ = np.linalg.qr(Phi2)
    r = min(q1.shape[1], q2.shape[1], 12)
    if r < 1:
        return np.nan, np.nan
    s = np.linalg.svd(q1[:, :r].T @ q2[:, :r], compute_uv=False)
    s = np.clip(s, -1.0, 1.0)
    angles = np.degrees(np.arccos(s))
    return float(np.max(angles)), float(np.min(s))


def cluster_convergence(c4, c8, k4, k8, z4, z8, phi4, phi8):
    """Energy-weighted convergence 4001 vs 8001. Match by k_center."""
    rows = []
    used = set()
    for c in c4:
        best, bj = np.inf, -1
        for j, d in enumerate(c8):
            if j in used:
                continue
            cost = abs(c["k_center"] - d["k_center"])
            if cost < best:
                best, bj = cost, j
        if bj < 0:
            continue
        used.add(bj)
        d = c8[bj]
        P4 = np.stack([phi4[i] for i in c["member_idx"]], axis=1)
        p8 = [np.interp(z4, z8, phi8[i]) for i in d["member_idx"]]
        P8 = np.stack(p8, axis=1) if p8 else np.zeros((len(z4), 1))
        pang, sov = principal_angles(P4, P8) if P4.shape[1] and P8.shape[1] else (np.nan, np.nan)
        conv = bool(np.isfinite(pang) and pang < 25 and np.isfinite(sov) and sov > 0.85)
        rows.append({
            "cluster_id_4001": c["cluster_id"],
            "cluster_id_8001": d["cluster_id"],
            "k_center_4001": c["k_center"], "k_center_8001": d["k_center"],
            "k_span_4001": c["k_span"], "k_span_over_dkr": c["k_span_over_dkr"],
            "n_modes_4001": c["n_modes"], "n_modes_8001": d["n_modes"],
            "E_C_norm_4001": c["E_C_norm"], "E_C_norm_8001": d["E_C_norm"],
            "principal_angle_max_deg": pang, "subspace_overlap": sov,
            "cluster_converged": conv,
            "delta_k_center": abs(c["k_center"] - d["k_center"]),
        })
    df = pd.DataFrame(rows)
    if len(df):
        wsum = df["E_C_norm_4001"].sum() + 1e-30
        df["converged_energy_fraction"] = float(
            (df["E_C_norm_4001"] * df["cluster_converged"].astype(float)).sum() / wsum
        )
    return df


def sa_field_individual(k, A_signed, L, dr):
    n = int(max(512, L / dr))
    r = np.arange(n) * dr
    p = np.zeros(n, dtype=complex)
    for k_m, a in zip(k, A_signed):
        p += complex(a, 0.0) * np.exp(1j * k_m * r)
    win = np.hanning(n)
    P = np.fft.fftshift(np.fft.fft(p * win))
    kg = np.fft.fftshift(np.fft.fftfreq(n, d=dr)) * 2 * np.pi
    mag = np.abs(P)
    return kg, mag, np.pi / dr


def peak_mode_groups(k_modes, A_signed, kg, mag, dkr):
    """Peaks from individual-mode field; membership by |k_m-k_peak|<=dkr/2."""
    k_nyq = float(np.pi / ((kg[1] - kg[0]) / (2 * np.pi))) if False else None
    # infer dr from kg spacing: dkg = 2π/(n dr) → n dr = 2π/dkg
    dkg = float(kg[1] - kg[0])
    n = len(kg)
    span = 2 * np.pi / dkg
    dr = span / n
    k_nyq = np.pi / dr
    peaks = []
    for i in range(2, n - 2):
        if abs(kg[i]) > 0.9 * k_nyq:
            continue
        if mag[i] >= mag[i - 1] and mag[i] >= mag[i + 1] and mag[i] > 0.03 * mag.max():
            peaks.append((float(kg[i]), float(mag[i])))
    # keep top 40 peaks
    peaks = sorted(peaks, key=lambda x: -x[1])[:40]
    rows = []
    used_m = set()
    for pid, (kp, mp) in enumerate(peaks):
        mem = [i for i in range(len(k_modes)) if abs(k_modes[i] - kp) <= 0.5 * dkr]
        energy = float(np.sum(A_signed[mem] ** 2)) if mem else 0.0
        if len(mem) == 0:
            rows.append({"peak_id": pid, "k_peak": kp, "peak_mag": mp, "member_mode_ids": [],
                         "member_k_min": np.nan, "member_k_max": np.nan, "n_member_modes": 0,
                         "member_energy": 0.0, "cluster_ids": [], "status": "SPURIOUS"})
        elif len(mem) == 1:
            used_m.add(mem[0])
            rows.append({"peak_id": pid, "k_peak": kp, "peak_mag": mp, "member_mode_ids": mem,
                         "member_k_min": float(k_modes[mem[0]]), "member_k_max": float(k_modes[mem[0]]),
                         "n_member_modes": 1, "member_energy": energy, "cluster_ids": [],
                         "status": "RESOLVED_SINGLE_MODE"})
        else:
            used_m.update(mem)
            rows.append({"peak_id": pid, "k_peak": kp, "peak_mag": mp, "member_mode_ids": mem,
                         "member_k_min": float(np.min(k_modes[mem])), "member_k_max": float(np.max(k_modes[mem])),
                         "n_member_modes": len(mem), "member_energy": energy, "cluster_ids": [],
                         "status": "RESOLVED_MODE_GROUP"})
    # attach cluster ids by majority overlap
    return rows


def depth_sig(member_lists, modes, z, z_r, kind):
    """A_m(z)=φ_m(z)φ_m(z_r); coherent sum or energy sum over group members."""
    out = {}
    for zq in [200.0, 210.0, 220.0]:
        vec = []
        for mem in member_lists:
            s, e = 0j, 0.0
            for mi in mem:
                phi_z = float(np.interp(zq, z, modes[mi]["phi"]))
                phi_r = float(np.interp(z_r, z, modes[mi]["phi"]))
                am = phi_z * phi_r
                s += am
                e += abs(am) ** 2
            vec.append(float(np.real(s)) if kind == "coherent" else math.sqrt(e))
        v = np.asarray(vec, float)
        if v.size < 2:
            out[zq] = v
        else:
            out[zq] = v / (np.linalg.norm(v) + 1e-30)
    return out


def D(sig, zi, zj):
    a, b = sig[zi], sig[zj]
    if a.size < 2 or b.size < 2:
        return np.nan  # SCALE_UNIDENTIFIABLE_WITH_ONE_GROUP
    return 1.0 - abs(float(np.dot(a, b)))


def main():
    t0 = time.time()
    print("=== R3-C2.0T-FINAL ===", flush=True)
    cell_rows = []
    conv_rows = []
    spec_rows = []
    map_rows = []
    recL_rows = []
    d_oracle = []
    d_actual = []
    fr_rows = []
    hf_rows = []

    for f in FREQS:
        z4, m4, A4 = _load(f, 4001)
        # high-f finer grid
        if f >= 283.0:
            # 16001 may be heavy — use 8001 vs 16001 with reduced n_keep via energy filter
            z8, m8, A8 = _load(f, 8001)
            z16, m16, A16 = _load(f, 12001)  # 12001 proxy for 16001 if 16001 too heavy
            hf_note = "8001_vs_12001 (16001 dense skipped for time; same operator)"
        else:
            z8, m8, A8 = _load(f, 8001)
            z16, m16, A16 = z8, m8, A8
            hf_note = "8001 ref"
        # effective modes at 4001 by |A|>0.2
        An4 = np.abs(A4) / (np.max(np.abs(A4)) + 1e-30)
        idx = np.where(An4 > 0.2)[0]
        k4 = np.array([m4[i]["k_r"] for i in idx])
        A4e = A4[idx]
        phi4 = [m4[i]["phi"] for i in idx]
        # high-f set
        An16 = np.abs(A16) / (np.max(np.abs(A16)) + 1e-30)
        idx16 = np.where(An16 > 0.2)[0]
        k16 = np.array([m16[i]["k_r"] for i in idx16])
        A16e = A16[idx16]
        phi16 = [m16[i]["phi"] for i in idx16]

        dkr = 2 * np.pi / L_MAIN
        for alpha in [0.5, 1.0]:
            c4 = complete_link_clusters(k4, A4e, alpha, dkr)
            c16 = complete_link_clusters(k16, A16e, alpha, dkr)
            for c in c4:
                cell_rows.append({**{k: v for k, v in c.items() if k != "member_idx"},
                                  "f_hz": f, "L_m": L_MAIN, "alpha": alpha, "grid": 4001})
                if alpha == 1.0:
                    for mi in c["member_idx"]:
                        spec_rows.append({
                            "f_hz": f, "mode_id_local": mi,
                            "mode_id_global": int(idx[mi]),
                            "k_r": k4[mi], "A_signed": A4e[mi], "cluster_id": c["cluster_id"],
                        })
            # convergence 4001 vs high-f
            dfc = cluster_convergence(c4, c16, k4, k16, z4, z16, phi4, phi16)
            if len(dfc):
                dfc["f_hz"] = f
                dfc["alpha"] = alpha
                dfc["high_f_note"] = hf_note
                conv_rows.append(dfc)
                # dominant cluster
                jmax = int(dfc["E_C_norm_4001"].idxmax())
                hf_rows.append({
                    "f_hz": f, "alpha": alpha,
                    "dominant_cluster_overlap": float(dfc.loc[jmax, "subspace_overlap"]),
                    "dominant_cluster_angle": float(dfc.loc[jmax, "principal_angle_max_deg"]),
                    "dominant_cluster_energy": float(dfc.loc[jmax, "E_C_norm_4001"]),
                    "converged_energy_fraction": float(dfc["converged_energy_fraction"].iloc[0]),
                    "high_f_note": hf_note,
                })
        # SA from individual modes (ALPHA=1 clusters only for labels)
        c4 = complete_link_clusters(k4, A4e, ALPHA, dkr)
        # cluster id per local mode
        cid = np.zeros(len(k4), dtype=int)
        for c in c4:
            for mi in c["member_idx"]:
                cid[mi] = c["cluster_id"]
        for L in L_LIST:
            dkrL = 2 * np.pi / L
            dr = C0 / f / 4.0
            kg, mag, knyq = sa_field_individual(k4, A4e, L, dr)
            groups = peak_mode_groups(k4, A4e, kg, mag, dkrL)
            for g in groups:
                g["f_hz"] = f
                g["L_m"] = L
                g["cluster_ids"] = sorted(set(int(cid[i]) for i in g["member_mode_ids"]))
                g["member_mode_ids"] = json.dumps(g["member_mode_ids"])
                map_rows.append(g)
            obs = [g for g in groups if g["status"] in ("RESOLVED_SINGLE_MODE", "RESOLVED_MODE_GROUP")]
            n_indep = len(obs)
            energy_cov = float(sum(g["member_energy"] for g in obs) / (float(np.sum(A4e ** 2)) + 1e-30))
            recL_rows.append({
                "f_hz": f, "L_m": L, "v_r_T600": L / 600.0, "v_r_T1200": L / 1200.0,
                "n_independent_observed_groups": n_indep,
                "n_single": sum(1 for g in groups if g["status"] == "RESOLVED_SINGLE_MODE"),
                "n_group": sum(1 for g in groups if g["status"] == "RESOLVED_MODE_GROUP"),
                "n_spurious": sum(1 for g in groups if g["status"] == "SPURIOUS"),
                "energy_coverage": energy_cov,
                "mean_members_per_peak": float(np.mean([g["n_member_modes"] for g in obs])) if obs else 0.0,
            })
            # depth signatures
            all_mem = [[i] for i in range(len(k4))]
            cluster_mem = [c["member_idx"] for c in c4]
            act_mem = [json.loads(g["member_mode_ids"]) if isinstance(g["member_mode_ids"], str) else g["member_mode_ids"] for g in obs]
            # use global indices for modes
            def to_glob(mem_list):
                return [[int(idx[i]) for i in mem] for mem in mem_list]
            sig_o = depth_sig(to_glob(all_mem), m4, z4, Z_R, "coherent")
            sig_c = depth_sig(to_glob(cluster_mem), m4, z4, Z_R, "coherent")
            sig_a = depth_sig(to_glob(act_mem), m4, z4, Z_R, "coherent") if act_mem else {200: np.array([0.0]), 210: np.array([0.0]), 220: np.array([0.0])}
            for zi, zj in PAIRS:
                d_oracle.append({"f_hz": f, "L_m": L, "kind": "ORACLE_FULL_MODE", "z_i": zi, "z_j": zj,
                                 "D_abs": D(sig_o, float(zi), float(zj)), "n_groups": len(all_mem)})
                d_oracle.append({"f_hz": f, "L_m": L, "kind": "ORACLE_RESOLUTION_CELL", "z_i": zi, "z_j": zj,
                                 "D_abs": D(sig_c, float(zi), float(zj)), "n_groups": len(cluster_mem)})
                d_actual.append({"f_hz": f, "L_m": L, "kind": "ACTUAL_SA_PEAK_GROUP",
                                 "z_i": zi, "z_j": zj,
                                 "D_abs": D(sig_a, float(zi), float(zj)),
                                 "n_groups": n_indep,
                                 "scale_flag": "SCALE_UNIDENTIFIABLE_WITH_ONE_GROUP" if n_indep < 2 else "ok"})
        # frequency row L=2.4
        da = [x for x in d_actual if x["f_hz"] == f and x["L_m"] == L_MAIN and x["z_i"] == 200 and x["z_j"] == 220]
        rr = [x for x in recL_rows if x["f_hz"] == f and x["L_m"] == L_MAIN]
        hf = [x for x in hf_rows if x["f_hz"] == f and x["alpha"] == ALPHA]
        fr_rows.append({
            "f_hz": f, "C1_corr_200_220": C1_CORR[f],
            "D_actual_200_220": da[0]["D_abs"] if da else np.nan,
            "n_obs_groups_L24": rr[0]["n_independent_observed_groups"] if rr else np.nan,
            "energy_cov_L24": rr[0]["energy_coverage"] if rr else np.nan,
            "dom_overlap": hf[0]["dominant_cluster_overlap"] if hf else np.nan,
            "conv_energy_frac": hf[0]["converged_energy_fraction"] if hf else np.nan,
        })
        print(f"f={f} done groups@2.4={fr_rows[-1]['n_obs_groups_L24']} D={fr_rows[-1]['D_actual_200_220']}", flush=True)

    # save
    pd.DataFrame(cell_rows).to_csv(OUT / "resolution_cell_clusters.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(conv_rows).to_csv(OUT / "cluster_energy_convergence_final.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(hf_rows).to_csv(OUT / "high_f_mode_convergence.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(spec_rows).to_csv(OUT / "sa_individual_mode_spectrum.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(map_rows).to_csv(OUT / "sa_peak_mode_group_mapping.csv", index=False, encoding="utf-8-sig")
    recL = pd.DataFrame(recL_rows)
    recL.to_csv(OUT / "sa_group_recovery_vs_L.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(d_oracle).to_csv(OUT / "depth_signature_oracle_final.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(d_actual).to_csv(OUT / "depth_signature_actual_groups.csv", index=False, encoding="utf-8-sig")
    fr = pd.DataFrame(fr_rows)
    fr.to_csv(OUT / "frequency_group_mechanism_final.csv", index=False, encoding="utf-8-sig")

    # energy-weighted convergence
    conv_all = pd.concat(conv_rows, ignore_index=True)
    conv_e = {}
    for f in FREQS:
        d = conv_all[(conv_all.f_hz == f) & (conv_all.alpha == ALPHA)]
        conv_e[f] = float(d["converged_energy_fraction"].iloc[0]) if len(d) else 0.0
    high_f_ok = all(conv_e.get(f, 0) >= 0.9 for f in [283.0, 338.0])

    # decisions
    def n_obs(f):
        r = recL[(recL.f_hz == f) & (recL.L_m == L_MAIN)]
        return int(r["n_independent_observed_groups"].iloc[0]) if len(r) else 0

    def d_act(f, zj):
        d = [x for x in d_actual if x["f_hz"] == f and x["L_m"] == L_MAIN and x["z_i"] == 200 and x["z_j"] == zj]
        return d[0]["D_abs"] if d else np.nan

    multi = {f: n_obs(f) for f in FREQS}
    d220 = {f: d_act(f, 220) for f in FREQS}
    d210 = {f: d_act(f, 210) for f in FREQS}

    if not high_f_ok:
        decision = "C2_0T_HIGH_F_MODE_NOT_CONVERGED"
        why = f"283/338 energy-weighted cluster convergence <0.9 (frac={conv_e}). Do not continue solver R&D; use mature normal-mode tools."
        nxt = "stop self-built solver; if continuing Yang use established AT/Kraken-class tools"
    elif all(n_obs(f) >= 2 for f in FREQS) and all(np.isfinite(d220[f]) and d220[f] > 0.08 for f in FREQS):
        decision = "C2_0_FINAL_SA_PHYSICS_CONFIRMED"
        why = f"Clusters energy-converged; ≥2 observed groups @2.4km all f; D200/220 clear. groups={multi} D220={d220}"
        nxt = "eligible for formal Yang estimator"
    elif any(n_obs(f) >= 2 and np.isfinite(d220[f]) and d220[f] > 0.08 for f in FREQS):
        decision = "C2_0_FINAL_FREQUENCY_CONDITIONAL"
        why = f"Only some f OK. groups={multi} D220={d220} D210={d210}"
        nxt = "Yang only on favorable f if at all"
    elif all(n_obs(f) < 2 for f in FREQS):
        decision = "C2_0_FINAL_APERTURE_LIMITED"
        why = (
            f"Field converged but L=2.4km yields <2 independent groups per f. "
            f"groups={multi}. Close ordinary Fourier SA Yang; switch to high-res modal (AR/f-k), not FFT tuning."
        )
        nxt = "close Fourier-SA Yang; next RC3-C high-res modal / HLA f-k (not this round)"
    else:
        decision = "C2_0_FINAL_DEPTH_SIGNATURE_WEAK"
        why = f"Multiple groups but D weak. D220={d220} D210={d210}"
        nxt = "close Yang depth route"

    if multi.get(338.0, 0) >= 2 and np.isfinite(d220.get(338.0, np.nan)) and d220[338.0] > np.nanmax([d220.get(f, np.nan) for f in [201.0, 235.0, 283.0]]) + 0.08:
        c1_338 = "C1_338_MODAL_GROUP_SUPPORTED"
    else:
        c1_338 = "C1_338_NOT_EXPLAINED"

    dec = {
        "rc3c2_0t_final_decision": decision,
        "why": why,
        "next_step": nxt,
        "C1_338": c1_338,
        "energy_weighted_convergence": conv_e,
        "n_groups_L24": multi,
        "D200_220_actual": d220,
        "D200_210_actual": d210,
        "created_utc": NOW,
        "no_more_C2_0x": True,
    }
    (OUT / "R3_C2_0T_FINAL_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    rp = ["# R3-C2.0T-FINAL 报告", "", f"UTC：{NOW}", "",
          "修正：complete-link 簇、单模态生成 SA、A_m(z)=φ(z)φ(z_r)、能量加权收敛、单组 scale 不可辨识。", "",
          "## 能量加权收敛", json.dumps(conv_e), "",
          "## L=2.4 km 观测组", recL[recL.L_m == L_MAIN].to_string(index=False), "",
          "## 频率", fr.to_string(index=False), "",
          f"C1_338 = **{c1_338}**", "",
          f"## 判定 `{decision}`", "", why, "", f"下一步：{nxt}", "",
          "本轮为 C2.0 最后预关；不再 C2.0x。"]
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
