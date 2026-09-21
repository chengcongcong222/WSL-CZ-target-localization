#!/usr/bin/env python3
"""R3-B1.2 finisher: cache neighbor roots once, compute E-STD Fisher, report, stop."""
from __future__ import annotations

import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from r3_b1_1_correction import find_eigenrays_root, steering_corr, u_hla
from r3_b1_2_final import (
    OUT,
    R0,
    ZS0,
    SIGMA_THETA_DEG,
    THETA_REL,
    decide,
    delay_peaks_final,
    ray_modal_status,
)

NOW = datetime.now(timezone.utc).isoformat()


def main():
    t0 = time.time()
    tracks = pd.read_csv(OUT / "rz_branch_association.csv")
    qual = pd.read_csv(OUT / "rz_branch_quality.csv")
    print("tracks", len(tracks), "qual", len(qual), flush=True)

    # cache all grid states from tracks
    cache = {}  # (bid, r, z) -> state
    for _, row in tracks.iterrows():
        if not bool(row["admitted"]):
            continue
        key = (row["branch_id"], float(row["grid_r_km"]), float(row["grid_zs_m"]))
        cache[key] = {
            "phi_deg": float(row["phi_arr_deg"]),
            "tau": float(row["tau"]),
            "launch_deg": float(row["launch_deg"]),
            "topology": row["topology"],
            "residual_m": float(row["residual_m"]),
        }

    center_bids = [b for b in tracks["branch_id"].unique()
                   if (b, R0, ZS0) in cache]
    print("center bids", center_bids, flush=True)
    fd_pts = [(R0, ZS0), (50.5, ZS0), (49.5, ZS0), (R0, 205.0), (R0, 195.0)]

    def get_state(bid, r_km, zs):
        key = (bid, float(r_km), float(zs))
        if key in cache:
            return cache[key]
        if (bid, R0, ZS0) not in cache:
            return None
        sc = cache[(bid, R0, ZS0)]
        _, roots = find_eigenrays_root(r_km, zs)
        if roots is None or roots.empty:
            return None
        ad = roots[roots["admitted"] == True]  # noqa: E712
        if ad.empty:
            return None
        best, bc = None, np.inf
        for _, rr in ad.iterrows():
            c = abs(rr["launch_deg"] - sc["launch_deg"]) / 3 + abs(rr["phi_arr_deg"] - sc["phi_deg"]) / 3 + abs(rr["tau"] - sc["tau"]) * 1e3 / 80
            if c < bc:
                bc, best = c, rr
        if best is None:
            return None
        st = {
            "phi_deg": float(best["phi_arr_deg"]),
            "tau": float(best["tau"]),
            "launch_deg": float(best["launch_deg"]),
            "topology": best["topology"],
            "residual_m": float(best["residual_m"]),
        }
        cache[key] = st
        print(f"  cached {bid} @ {r_km},{zs}: phi={st['phi_deg']:.3f} tau={st['tau']:.5f} res={st['residual_m']:.3f}", flush=True)
        return st

    # ensure FD coverage
    for bid in center_bids:
        for r_km, zs in fd_pts:
            get_state(bid, r_km, zs)

    def pack(bids, r_km, zs, mode, theta):
        sts = []
        for b in bids:
            s = get_state(b, r_km, zs)
            if s is not None:
                sts.append((b, s))
        if len(sts) < 2:
            return None
        sts = sorted(sts, key=lambda x: x[1]["tau"])
        y = []
        for _, s in sts:
            if mode == "IDEAL":
                y.append(math.radians(s["phi_deg"]))
            else:
                y.append(u_hla(theta, s["phi_deg"]))
        for _, s in sts[1:]:
            y.append(s["tau"] - sts[0][1]["tau"])
        return np.asarray(y), [s for _, s in sts], [b for b, _ in sts]

    rows_j, rows_f = [], []
    bids = center_bids
    for sigma_th in SIGMA_THETA_DEG:
        sig_th = math.radians(sigma_th)
        for theta in THETA_REL:
            for mode, label in [("IDEAL", "IDEAL_ELEVATION+delay"),
                                ("HLA", "HLA_u+delay"),
                                ("DELAY", "delay-only")]:
                packs = [pack(bids, r, z, mode, theta) for r, z in fd_pts]
                if any(p is None for p in packs):
                    print("  skip", label, "theta", theta, "missing pack", flush=True)
                    continue
                common = set(packs[0][2])
                for p in packs[1:]:
                    common &= set(p[2])
                if len(common) < 2:
                    print("  skip", label, "common", common, flush=True)
                    continue
                def restrict(p):
                    y, sts, bs = p
                    keep = sorted([(b, s) for b, s in zip(bs, sts) if b in common], key=lambda x: x[1]["tau"])
                    yy = []
                    for _, s in keep:
                        yy.append(math.radians(s["phi_deg"]) if mode == "IDEAL" else u_hla(theta, s["phi_deg"]))
                    for _, s in keep[1:]:
                        yy.append(s["tau"] - keep[0][1]["tau"])
                    return np.asarray(yy), [s for _, s in keep], [b for b, _ in keep]
                packs = [restrict(p) for p in packs]
                n = min(len(p[0]) for p in packs)
                J = np.zeros((n, 2))
                J[:, 0] = (packs[1][0][:n] - packs[2][0][:n]) / 1000.0
                J[:, 1] = (packs[3][0][:n] - packs[4][0][:n]) / 10.0
                n_br = len(packs[0][1])
                if mode == "IDEAL":
                    sig = [math.radians(0.3)] * n_br + [2.4e-4] * max(n_br - 1, 0)
                elif mode == "HLA":
                    sig = []
                    for s in packs[0][1]:
                        su = steering_corr(14.0, 8, u_hla(theta, s["phi_deg"]), u_hla(theta, s["phi_deg"]), 250.0)[1]
                        du_dth = -math.sin(math.radians(theta)) * math.cos(math.radians(s["phi_deg"]))
                        sig.append(math.sqrt(su ** 2 + (du_dth * sig_th) ** 2))
                    sig += [2.4e-4] * max(n_br - 1, 0)
                else:
                    sig = [2.4e-4] * max(n_br - 1, 0)
                sig = (sig + [sig[-1]] * n)[:n]
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
                    "branch_set": json.dumps(list(common)),
                    "sigma_theta_deg": sigma_th, "theta_rel_deg": theta, "case": label,
                    "n_obs": n, "n_branches": n_br,
                    "dy_dr": json.dumps(J[:, 0].tolist()), "dy_dz": json.dumps(J[:, 1].tolist()),
                    "phi_center_deg": json.dumps([round(s["phi_deg"], 3) for s in packs[0][1]]),
                    "tau_center_s": json.dumps([round(s["tau"], 5) for s in packs[0][1]]),
                })
                rows_f.append({
                    "branch_set": json.dumps(list(common)),
                    "sigma_theta_deg": sigma_th, "theta_rel_deg": theta, "case": label,
                    "n_obs": n, "n_branches": n_br, "rank": rank,
                    "sv": json.dumps(sv.tolist()),
                    "sv_min": float(sv.min()) if sv.size else np.nan,
                    "cond_J": float(sv[0] / sv[-1]) if sv.size and sv[-1] > 1e-30 else np.inf,
                    "crb_r_m": crb_r, "crb_z_m": crb_z, "corr_rz": corr, "cond_F": condF,
                })
                print(f"  Fisher {label} th={theta} sTH={sigma_th}: rank={rank} CRBr={crb_r:.1f} CRBz={crb_z:.1f}", flush=True)

    jac_df, fish_df = pd.DataFrame(rows_j), pd.DataFrame(rows_f)
    jac_df.to_csv(OUT / "estd_rz_jacobian_final.csv", index=False, encoding="utf-8-sig")
    fish_df.to_csv(OUT / "estd_rz_fisher_final.csv", index=False, encoding="utf-8-sig")

    delay_df = delay_peaks_final()
    delay_df.to_csv(OUT / "delay_global_peaks_final.csv", index=False, encoding="utf-8-sig")
    modal_df, modal_status = ray_modal_status()
    modal_df.to_csv(OUT / "ray_modal_delay_crosscheck.csv", index=False, encoding="utf-8-sig")

    decision, why, nxt, notes = decide(fish_df, qual, delay_df, modal_status)
    notes["center_branches"] = center_bids
    notes["fd_branch_sets"] = fish_df["branch_set"].unique().tolist() if len(fish_df) else []

    dec = {
        "rc3b1_2_decision": decision,
        "why": why,
        "next_step": nxt,
        "notes": notes,
        "created_utc": NOW,
        "stop": "after R3-B1.2; no B2/RC3-C/P5",
        "frozen_from_b1_1": [
            "HLA elevation projection weak",
            "local delay mainlobe ms-scale; old 300-400ms revoked",
            "eigenray residual <= 0.1 m",
            "S0_DENSE local delay usable; S1/S2 globally ambiguous",
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
    rp.append("## 0. B1.1 已冻结")
    rp.append("")
    rp.append("- HLA 俯仰投影弱；±φ 不可分")
    rp.append("- 局部时延主瓣 ms 级（旧 300–400 ms **永久撤销**）")
    rp.append("- 射线 residual≤0.1 m；5×5 网格二维关联完成")
    rp.append("")
    rp.append("## 1. 中心二维分支")
    rp.append("")
    rp.append(f"- 中心 (50 km, 200 m) 分支：**{center_bids}**")
    rp.append(f"- jac_ready：{int(qual['jac_ready'].sum()) if len(qual) else 0}")
    if len(qual):
        rp.append("")
        rp.append("| id | topo | cells | jac_ready | φ drift° | τ drift s | max resid m |")
        rp.append("| --- | --- | --- | --- | --- | --- | --- |")
        for _, r in qual.iterrows():
            rp.append(f"| {r['branch_id']} | {r['topology']} | {r['n_cells']} | {r['jac_ready']} | {fnum(r['phi_drift_deg'])} | {fnum(r['tau_drift_s'],5)} | {fnum(r['max_residual_m'])} |")
    rp.append("")
    rp.append("## 2. E-STD 真实 r–z Fisher（B1 核心量）")
    rp.append("")
    if len(fish_df):
        rp.append("| case | σθ° | θ° | rank | CRB_r m | CRB_z m | corr | sv_min | branches |")
        rp.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for _, r in fish_df.iterrows():
            rp.append(
                f"| {r['case']} | {r['sigma_theta_deg']} | {r['theta_rel_deg']} | {r['rank']} | "
                f"{fnum(r['crb_r_m'],1)} | {fnum(r['crb_z_m'],1)} | {fnum(r['corr_rz'])} | "
                f"{fnum(r['sv_min'],6)} | {r['branch_set']} |"
            )
    else:
        rp.append("（Fisher 仍空）")
    rp.append("")
    rp.append("## 3. 时延全局（PSL）")
    rp.append("")
    rp.append("| 源 | 局部主瓣 ms | PSL | global_status |")
    rp.append("| --- | --- | --- | --- |")
    for _, r in delay_df.iterrows():
        rp.append(f"| {r['source']} | {fnum(r['local_mainlobe_width_s']*1e3,2)} | {fnum(r['psl_outside_first_null'])} | **{r['global_status']}** |")
    rp.append("")
    rp.append(f"## 4. ray-modal：`{modal_status}`")
    rp.append("")
    rp.append("## 5. R3-B1 最终判定")
    rp.append("")
    rp.append(f"### `{decision}`")
    rp.append("")
    rp.append(why)
    rp.append("")
    rp.append(f"**下一步**：{nxt}")
    rp.append("")
    rp.append("## 6. 停止")
    rp.append("")
    rp.append("- 不进 B2 / RC3-C / P5；不增加新方法")
    rp.append(f"- **B1.2 完成；判定 `{decision}`**")
    rp.append("")
    (OUT / "R3_B1_2_REPORT.md").write_text("\n".join(rp), encoding="utf-8")

    gs = [
        "# R3-B1.2 — GPT 同步稿", "",
        f"- **判定：{decision}**",
        f"- {why}",
        f"- 下一步：{nxt}", "",
        "## Fisher",
        fish_df.to_string(index=False) if len(fish_df) else "(empty)", "",
        "## Delay",
        delay_df[["source", "local_mainlobe_width_s", "psl_outside_first_null", "global_status"]].to_string(index=False), "",
        f"center branches: {center_bids}", "",
        "停止：无 B2/RC3-C/P5。", "",
    ]
    (OUT / "R3_B1_2_GPT_SYNC.md").write_text("\n".join(gs), encoding="utf-8")

    print(f"DONE {time.time()-t0:.1f}s DECISION={decision}", flush=True)
    print(why, flush=True)


if __name__ == "__main__":
    main()
