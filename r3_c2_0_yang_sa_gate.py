#!/usr/bin/env python3
"""R3-C2.0: Yang SA modal observability pre-gate.

Erratum/full text NOT recovered from network → paper recovery BLOCKED.
Physics pre-gate uses ONLY formulas specified in the task (modal expansion, eta, DFT peaks).
Does NOT invent Yang beamforming equations.
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
OUT = ROOT / "results" / "R3_C2_Yang_SA_depth"
FIG = OUT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

C0 = 1500.0
Z_R = 200.0
FREQS = [201.0, 235.0, 283.0, 338.0]
Z_S_LIST = [180.0, 190.0, 200.0, 210.0, 220.0]
L_LIST = [0.15e3, 0.3e3, 0.6e3, 1.2e3, 2.4e3]
T_FOR_VR = [600.0, 1200.0]
AMPS_THRESH = [0.1, 0.2]
PAIRS = [(180, 190), (190, 200), (200, 210), (210, 220), (180, 200), (200, 220)]
C1_CORR_200_220 = {201.0: 0.579, 235.0: 0.994, 283.0: 1.000, 338.0: 0.067}

MODE = p4.MODE_ENV["E0"]

CONFIG = {
    "package": "R3_C2_Yang_SA_depth / R3_C2_0",
    "created_utc": NOW,
    "anchor": {
        "yang2015": "T.C. Yang, JASA 138(3) 1678-1686 (2015), DOI 10.1121/1.4929748",
        "erratum": "JASA 144(6) 3075 (2018), DOI 10.1121/1.5081712",
        "retrieval": "network fetch of PubMed/AIP blocked (cookies/transport) — erratum content NOT recovered",
    },
    "f0_condition": "KNOWN_F0 / KNOWN_TRUE_RANGE_INCREMENT (physics ceiling only)",
    "aperture_note": "L = |r(t_end)-r(t_start)| radial range change, NOT total track length",
    "eta": "eta_m = L |k_r,m+1-k_r,m| / (2π); Δk_SA ≈ 2π/L",
    "stop": ["no full Yang estimator", "no Doppler error", "no f0 error", "no MC", "no P5"],
}


def mode_table(f, z_s, z_r=Z_R, amp_thr=0.1):
    ms = MODE.modes(float(f))
    kr = [float(k) for k, _ in ms]
    phis = [float(np.interp(z_s, MODE.z, phi)) for _, phi in ms]
    phir = [float(np.interp(z_r, MODE.z, phi)) for _, phi in ms]
    A = np.array([a * b for a, b in zip(phis, phir)])
    amax = float(np.max(np.abs(A))) if A.size else 1.0
    An = np.abs(A) / (amax + 1e-30)
    dk = [abs(kr[i + 1] - kr[i]) for i in range(len(kr) - 1)] + [np.nan]
    rows = []
    for i in range(len(kr)):
        rows.append({
            "f_hz": f, "z_s_m": z_s,
            "mode_id": i,
            "k_r": kr[i],
            "adj_delta_k_r": dk[i],
            "phi_m_zs": phis[i],
            "phi_m_zr": phir[i],
            "A_m": float(A[i]),
            "A_norm": float(An[i]),
            "effective_A0p1": bool(An[i] > 0.1),
            "effective_A0p2": bool(An[i] > 0.2),
        })
    return rows, kr, An


def eta_modes(L, kr, An, thr=0.1):
    """Adjacent effective-mode separability index."""
    idx = [i for i, a in enumerate(An) if a > thr]
    if len(idx) < 2:
        return [], np.nan
    etas = []
    for a, b in zip(idx[:-1], idx[1:]):
        dk = abs(kr[b] - kr[a])
        eta = L * dk / (2 * np.pi)
        etas.append({"i": a, "j": b, "dk": dk, "eta": float(eta), "resolvable": bool(eta > 1.0)})
    return etas, float(np.median([e["eta"] for e in etas]))


def controlled_sa_recovery(f, z_s, L, dr, thr=0.1):
    """Generate p(r)=Σ A_m exp(j k_rm r) on [0,L]; DFT peak vs predicted k_r."""
    rows, kr, An = mode_table(f, z_s, amp_thr=thr)
    n = int(max(32, L / dr))
    r = np.arange(n) * dr
    # excitation uses A_norm as complex amplitude (sign from phi products)
    A = np.array([rw["A_m"] for rw in rows])
    p = np.zeros(n, dtype=complex)
    for i, k in enumerate(kr):
        p += A[i] * np.exp(1j * k * r)
    # wavenumber spectrum via FFT along r
    P = np.fft.fftshift(np.fft.fft(p))
    kgrid = np.fft.fftshift(np.fft.fftfreq(n, d=dr)) * 2 * np.pi
    mag = np.abs(P)
    # peaks
    peaks = []
    for i in range(1, n - 1):
        if mag[i] >= mag[i - 1] and mag[i] >= mag[i + 1] and mag[i] > 0.15 * mag.max():
            peaks.append((float(kgrid[i]), float(mag[i])))
    peaks = sorted(peaks, key=lambda x: -x[1])[:8]
    # match predicted effective k
    rec = []
    for i, (k, m) in enumerate(peaks):
        # nearest predicted
        pred = min(kr, key=lambda kp: abs(kp - k))
        rec.append({
            "f_hz": f, "z_s_m": z_s, "L_m": L, "dr_m": dr,
            "recovered_k": k, "pred_k_nearest": float(pred),
            "peak_err": abs(k - pred),
            "peak_mag": m,
            "is_eff_mode": bool(min(abs(kr[j] - pred) for j in range(len(kr)) if An[j] > thr) < 1e-6),
        })
    eta_list, med_eta = eta_modes(L, kr, An, thr)
    n_res = sum(1 for e in eta_list if e["resolvable"])
    return rec, eta_list, med_eta, n_res


def depth_signature(z_s_list, f, thr=0.2):
    """Normalized mode-depth vector a(z)=[phi_m(z_s)] on effective modes at z_r."""
    rows, kr, An = mode_table(f, Z_R, amp_thr=thr)  # excitation at z_r for listing modes
    # use modes effective at z_s=200 as reference set
    rows2, _, An2 = mode_table(f, 200.0, amp_thr=thr)
    idx = [i for i, a in enumerate(An2) if a > thr]
    vecs = {}
    for z in z_s_list:
        phis = [float(np.interp(z, MODE.z, MODE.modes(f)[i][1])) for i in idx]
        v = np.asarray(phis, float)
        v = v / (np.linalg.norm(v) + 1e-30)
        vecs[z] = v
    out = []
    for zi, zj in PAIRS:
        a, b = vecs[float(zi)], vecs[float(zj)]
        corr = float(np.dot(a, b))
        D = 1.0 - abs(corr)
        out.append({
            "f_hz": f, "z_i": zi, "z_j": zj,
            "n_modes": len(idx),
            "corr": corr, "D_abs": D,
            "mode_ids": json.dumps(idx),
        })
    return out, idx, vecs


def main():
    t0 = time.time()
    print("=== R3-C2.0 Yang SA modal observability pre-gate ===", flush=True)
    (OUT / "R3_C2_0_CONFIG.json").write_text(json.dumps(CONFIG, ensure_ascii=False, indent=2), encoding="utf-8")

    # --- Paper recovery ---
    am = [
        "# ANCHOR — Yang 2015 / Erratum 2018",
        "",
        f"UTC：{NOW}",
        "",
        "## 文献",
        "",
        "- T.C. Yang, *Source depth estimation based on synthetic aperture beamforming for a moving source*, JASA 138(3), 1678–1686 (2015). DOI 10.1121/1.4929748",
        "- Erratum: JASA 144(6), 3075 (2018). DOI 10.1121/1.5081712",
        "",
        "## 检索状态",
        "",
        "- **PubMed / AIP / DOI 本轮网络访问失败**（cookie/transport）",
        "- **Erratum 具体修正公式：NOT_RECOVERED**",
        "- 因此 **不得编写/猜测 Yang 正式波束形成公式细节**",
        "",
        "## 任务书已给出的物理量（本轮唯一使用）",
        "",
        "- 模态场：\\(p(r,f)=\\sum_m A_m e^{jk_{rm}r}\\)，\\(A_m\\propto\\phi_m(z_s)\\phi_m(z_r)\\)",
        "- 合成孔径：\\(L=|r(t_{end})-r(t_{start})|\\)（**径向相对距离变化**）",
        "- 分辨尺度：\\(\\Delta k_{SA}\\approx 2\\pi/L\\)，\\(\\eta_m=L|k_{r,m+1}-k_{r,m}|/(2\\pi)\\)",
        "- 可控恢复：对 \\(p(r)\\) 做空间谱（FFT）比对预测 \\(k_r\\) 峰",
        "",
        "## 任务书未覆盖、本轮禁止臆造",
        "",
        "- Yang 正文 steering 相位约定、SA beam 具体核、Doppler–Δr 精确式、Erratum 符号修正",
        "- PLL 实测细节",
        "",
        "## f0 条件",
        "",
        "- **KNOWN_F0 / KNOWN_TRUE_RANGE_INCREMENT**（理论上限）",
        "",
    ]
    (OUT / "ANCHOR_YANG2015_METHOD.md").write_text("\n".join(am), encoding="utf-8")
    (OUT / "YANG2015_ERRATUM_NOTE.md").write_text(
        "\n".join([
            "# YANG2015 Erratum",
            "",
            f"UTC：{NOW}",
            "",
            "**状态：NOT_RECOVERED**",
            "",
            "DOI 10.1121/1.5081712 (JASA 144(6), 3075, 2018) 本轮无法从 PubMed/AIP 抓取正文。",
            "",
            "按任务规程：**不得自行猜测被修正的公式/符号**。",
            "",
            "影响：任何依赖 Yang 原文精确相位/波束核的实现保持 **C2_0_PAPER_RECOVERY_BLOCKED**，",
            "直至 Erratum 内容恢复。本轮仅输出任务书给定的模态/孔径可观测性物理量。",
            "",
        ]),
        encoding="utf-8",
    )

    # --- Mode tables ---
    mode_rows = []
    for f in FREQS:
        for z_s in Z_S_LIST:
            rows, kr, An = mode_table(f, z_s)
            mode_rows.extend(rows)
    pd.DataFrame(mode_rows).to_csv(OUT / "estd_mode_table.csv", index=False, encoding="utf-8-sig")

    # --- Aperture / eta ---
    ap_rows = []
    for f in FREQS:
        rows, kr, An = mode_table(f, 200.0)
        for thr in AMPS_THRESH:
            for L in L_LIST:
                eta_list, med = eta_modes(L, kr, An, thr)
                n_res = sum(1 for e in eta_list if e["resolvable"]) if eta_list else 0
                v_need = {T: L / T for T in T_FOR_VR}
                ap_rows.append({
                    "f_hz": f, "amp_thr": thr, "L_m": L, "L_km": L / 1e3,
                    "delta_k_SA": 2 * np.pi / L,
                    "median_eta": med,
                    "n_adjacent_pairs": len(eta_list),
                    "n_resolvable_pairs": n_res,
                    "v_r_needed_T600": v_need[600.0],
                    "v_r_needed_T1200": v_need[1200.0],
                    "note": "L is radial Δr not total track",
                })
    ap_df = pd.DataFrame(ap_rows)
    ap_df.to_csv(OUT / "sa_aperture_resolution.csv", index=False, encoding="utf-8-sig")

    # --- Controlled SA recovery ---
    rec_rows = []
    for f in FREQS:
        for L in [0.6e3, 1.2e3, 2.4e3]:
            for dr_rel, dname in [(8.0, "lam/8"), (4.0, "lam/4"), (2.0, "lam/2")]:
                lam = C0 / f
                dr = lam / dr_rel
                rec, _, med_eta, n_res = controlled_sa_recovery(f, 200.0, L, dr, thr=0.1)
                for rr in rec:
                    rr["dr_name"] = dname
                    rr["med_eta"] = med_eta
                    rr["n_resolvable_pairs"] = n_res
                    rec_rows.append(rr)
    rec_df = pd.DataFrame(rec_rows)
    rec_df.to_csv(OUT / "controlled_mode_recovery.csv", index=False, encoding="utf-8-sig")

    # --- Depth mode signatures ---
    dep_rows = []
    for f in FREQS:
        out, idx, vecs = depth_signature(Z_S_LIST, f, thr=0.2)
        dep_rows.extend(out)
        print(f"  depth sig {f}Hz modes={idx}", flush=True)
        for row in out:
            print(f"    {row['z_i']}/{row['z_j']}: D={row['D_abs']:.4f}", flush=True)
    dep_df = pd.DataFrame(dep_rows)
    dep_df.to_csv(OUT / "depth_mode_signature.csv", index=False, encoding="utf-8-sig")

    # --- Frequency mechanism vs C1 ---
    freq_rows = []
    for f in FREQS:
        rows, kr, An = mode_table(f, 200.0)
        n_eff = int(np.sum(An > 0.1))
        n_eff2 = int(np.sum(An > 0.2))
        eta_list, med = eta_modes(2.4e3, kr, An, 0.1)
        n_res = sum(1 for e in eta_list if e["resolvable"]) if eta_list else 0
        drow = dep_df[(dep_df.f_hz == f) & (dep_df.z_i == 200) & (dep_df.z_j == 220)]
        D220 = float(drow["D_abs"].iloc[0]) if len(drow) else np.nan
        d210 = dep_df[(dep_df.f_hz == f) & (dep_df.z_i == 200) & (dep_df.z_j == 210)]
        D210 = float(d210["D_abs"].iloc[0]) if len(d210) else np.nan
        freq_rows.append({
            "f_hz": f,
            "C1_corr_200_220": C1_CORR_200_220[f],
            "C1_D_abs_implied": 1.0 - abs(C1_CORR_200_220[f]),
            "n_eff_modes_A0p1": n_eff,
            "n_eff_modes_A0p2": n_eff2,
            "n_resolvable_eta_L2p4km": n_res,
            "median_eta_L2p4km": med,
            "mode_D_200_220": D220,
            "mode_D_200_210": D210,
        })
    freq_df = pd.DataFrame(freq_rows)
    freq_df.to_csv(OUT / "frequency_depth_mechanism.csv", index=False, encoding="utf-8-sig")
    print(freq_df.to_string(index=False), flush=True)

    # --- Figures (simple SVG) ---
    def write_svg(path, title, xs, ys_list, xlab, ylab, names):
        w, h = 760, 360
        p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}">',
             f'<rect width="{w}" height="{h}" fill="#f7f4ef"/>',
             f'<text x="380" y="28" text-anchor="middle" font-family="sans-serif" font-size="14" font-weight="600">{title}</text>',
             f'<rect x="60" y="50" width="620" height="250" fill="#fff" stroke="#ccc"/>']
        all_y = [y for ys in ys_list for y in ys if y is not None and np.isfinite(y)]
        ymin, ymax = (min(all_y + [0]), max(all_y + [1]))
        xmin, xmax = min(xs), max(xs) if xs else (0, 1)
        if xmax - xmin < 1e-12:
            xmax = xmin + 1
        cols = ["#b45309", "#0f766e", "#1d4ed8", "#8a8a8a"]
        for i, ys in enumerate(ys_list):
            pts = []
            for x, y in zip(xs, ys):
                if y is None or not np.isfinite(y):
                    continue
                X = 60 + (x - xmin) / (xmax - xmin) * 620
                Y = 50 + 250 - (y - ymin) / (ymax - ymin + 1e-15) * 250
                pts.append(f"{X:.1f},{Y:.1f}")
            if len(pts) >= 2:
                p.append(f'<polyline fill="none" stroke="{cols[i%4]}" stroke-width="2" points="{" ".join(pts)}"/>')
            p.append(f'<rect x="520" y="{60+i*16}" width="12" height="3" fill="{cols[i%4]}"/>')
            p.append(f'<text x="540" y="{64+i*16}" font-size="11" font-family="sans-serif">{names[i]}</text>')
        p.append(f'<text x="380" y="330" text-anchor="middle" font-size="12" font-family="sans-serif">{xlab}</text>')
        p.append(f'<text x="30" y="180" font-size="12" font-family="sans-serif" transform="rotate(-90 30 180)">{ylab}</text></svg>')
        Path(path).write_text("\n".join(p), encoding="utf-8")

    # fig1 eta vs L for 338Hz
    sub = ap_df[(ap_df.f_hz == 338.0) & (ap_df.amp_thr == 0.1)]
    write_svg(FIG / "fig1_mode_spacing_vs_aperture.svg",
              "fig1 η vs L (338Hz, thr=0.1)", sub["L_km"].tolist(),
              [sub["median_eta"].tolist(), sub["n_resolvable_pairs"].tolist()],
              "L (km)", "eta / n_res", ["median eta", "n resolvable"])
    # fig2 peak err
    sub = rec_df[(rec_df.f_hz == 338.0) & (rec_df.dr_name == "lam/4")]
    if len(sub):
        write_svg(FIG / "fig2_controlled_mode_recovery.svg",
                  "fig2 recovered k peak error (338Hz, λ/4)",
                  list(range(len(sub))), [sub["peak_err"].tolist()],
                  "peak idx", "|Δk|", ["peak_err"])
    # fig3 depth signature D
    write_svg(FIG / "fig3_depth_mode_signature.svg",
              "fig3 mode-depth D_abs (200/220)",
              FREQS, [[float(dep_df[(dep_df.f_hz==f)&(dep_df.z_i==200)&(dep_df.z_j==220)]["D_abs"].iloc[0]) for f in FREQS]],
              "f Hz", "D_abs", ["200/220"])
    # fig4 C1 corr vs mode D
    write_svg(FIG / "fig4_frequency_comparison.svg",
              "fig4 C1 corr vs mode D_abs 200/220",
              FREQS,
              [[1-abs(C1_CORR_200_220[f]) for f in FREQS],
               [float(dep_df[(dep_df.f_hz==f)&(dep_df.z_i==200)&(dep_df.z_j==220)]["D_abs"].iloc[0]) for f in FREQS]],
              "f Hz", "D", ["C1 1-|corr|", "mode vector D"])

    # --- Decision ---
    # Physics summary at L=2.4km thr=0.1
    phys = []
    for f in FREQS:
        rows, kr, An = mode_table(f, 200.0)
        eta_list, med = eta_modes(2.4e3, kr, An, 0.1)
        n_res = sum(1 for e in eta_list if e["resolvable"]) if eta_list else 0
        d220 = float(dep_df[(dep_df.f_hz==f)&(dep_df.z_i==200)&(dep_df.z_j==220)]["D_abs"].iloc[0])
        d210 = float(dep_df[(dep_df.f_hz==f)&(dep_df.z_i==200)&(dep_df.z_j==210)]["D_abs"].iloc[0])
        phys.append(dict(f=f, n_res=n_res, med_eta=med, D220=d220, D210=d210))
    any_res = any(p["n_res"] >= 2 for p in phys)
    depth_sensitive = any(p["D220"] > 0.08 for p in phys)
    freq_cond = any(p["n_res"] >= 1 and p["D220"] > 0.08 for p in phys) and not all(p["D220"] > 0.08 for p in phys)

    # HARD GATE: erratum not recovered
    paper_blocked = True
    if paper_blocked:
        decision = "C2_0_PAPER_RECOVERY_BLOCKED"
        why = (
            "Yang 2015 full text + 2018 Erratum (DOI 10.1121/1.5081712) NOT recovered "
            "(PubMed/AIP/DOI network blocked). Per protocol, no Yang beamforming formulas invented. "
            "Task-specified modal/aperture physics pre-gate computed below as supporting numbers only."
        )
        nxt = "recover Yang2015+Erratum text before any SA depth estimator; then re-open C2 gate"
    elif not any_res:
        decision = "C2_0_APERTURE_LIMITED"
        why = f"Even at L=2.4 km, effective modes not separated. phys={phys}"
        nxt = "Zhu/Yang intensity-SA not enough; consider longer radial Δr or other RC3-C methods"
    elif depth_sensitive and any(p["n_res"] >= 2 for p in phys):
        decision = "C2_0_SA_MODE_PHYSICS_CONFIRMED"
        why = f"Multiple effective modes resolvable and depth-sensitive. phys={phys}"
        nxt = "eligible for full Yang estimator later (needs erratum)"
    elif freq_cond:
        decision = "C2_0_FREQUENCY_CONDITIONAL"
        why = f"Only some frequencies (e.g. 338 Hz) have mode separability + depth signature. phys={phys}"
        nxt = "frequency-selective SA depth only; still needs Yang+erratum for estimator"
    else:
        decision = "C2_0_DEPTH_SIGNATURE_WEAK"
        why = f"Modes separable but depth vectors similar for 180/200/220. phys={phys}"
        nxt = "close Yang-SA for fine depth near 200 m unless new physics"

    notes = {
        "erratum_recovered": False,
        "paper_recovery": "BLOCKED",
        "physics_pre_gate": phys,
        "any_resolvable_modes_L2p4km": any_res,
        "depth_signature_sensitive_200_220": depth_sensitive,
        "C1_338Hz_explained": (
            "338 Hz has stronger mode-vector D_200/220 / more effective modes than 235/283 — "
            "consistent with frequency-selective modal depth information (see frequency_depth_mechanism.csv)"
        ) if any(p["f"] == 338.0 and p["D220"] > 0.15 for p in phys) else "not fully explained — keep as numerical cue only",
        "aperture_is_radial": True,
    }
    dec = {
        "rc3c2_0_decision": decision,
        "why": why,
        "next_step": nxt,
        "notes": notes,
        "created_utc": NOW,
        "C1_frozen": "Zhu-QF CLOSED; Zhu-QM C1_QM_MULTIFREQ_ONLY_CONDITIONAL; HLA no depth DOF",
        "stop": "no full Yang / Doppler / f0 error / MC / P5",
    }
    (OUT / "R3_C2_0_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    rp = []
    rp.append("# R3-C2.0 报告：Yang 合成孔径模态可观测性预关口")
    rp.append("")
    rp.append(f"UTC：{NOW}")
    rp.append("")
    rp.append("## 0. C1 冻结")
    rp.append("")
    rp.append("- Zhu-QF：CLOSED；Zhu-QM：`C1_QM_MULTIFREQ_ONLY_CONDITIONAL`（20 m 候选信号 ≠ 定深精度）")
    rp.append("- HLA 无额外深度自由度；不再 C1.x")
    rp.append("")
    rp.append("## 1. 文献恢复")
    rp.append("")
    rp.append("- Yang 2015 DOI 10.1121/1.4929748；Erratum 2018 DOI 10.1121/1.5081712")
    rp.append("- **网络未取回 Erratum 内容 → `C2_0_PAPER_RECOVERY_BLOCKED`**，不臆造 Yang 公式")
    rp.append("- f0 条件：KNOWN_F0 / KNOWN_TRUE_RANGE_INCREMENT")
    rp.append("")
    rp.append("## 2. 孔径定义")
    rp.append("")
    rp.append("**L = 径向相对距离变化 |Δr|**，不是航程。")
    rp.append("")
    rp.append("| L km | v_r(T=600) m/s | v_r(T=1200) m/s |")
    rp.append("| --- | --- | --- |")
    for L in L_LIST:
        rp.append(f"| {L/1e3:.2f} | {L/600:.4f} | {L/1200:.4f} |")
    rp.append("")
    rp.append("## 3. 模态与 η（z_s=200, thr=0.1, L=2.4 km）")
    rp.append("")
    rp.append("| f Hz | n_eff | median η | n 可分对 | D(200/220) | D(200/210) |")
    rp.append("| --- | --- | --- | --- | --- | --- |")
    for f in FREQS:
        d220 = float(dep_df[(dep_df.f_hz==f)&(dep_df.z_i==200)&(dep_df.z_j==220)]["D_abs"].iloc[0])
        d210 = float(dep_df[(dep_df.f_hz==f)&(dep_df.z_i==200)&(dep_df.z_j==210)]["D_abs"].iloc[0])
        rows, kr, An = mode_table(f, 200.0)
        el, med = eta_modes(2.4e3, kr, An, 0.1)
        nres = sum(1 for e in el if e["resolvable"]) if el else 0
        rp.append(f"| {f:.0f} | {int(np.sum(An>0.1))} | {med:.3f} | {nres} | {d220:.4f} | {d210:.4f} |")
    rp.append("")
    rp.append("## 4. C1 的 338 Hz 现象")
    rp.append("")
    rp.append(freq_df.to_string(index=False))
    rp.append("")
    rp.append(notes["C1_338Hz_explained"])
    rp.append("")
    rp.append("## 5. 判定")
    rp.append("")
    rp.append(f"### `{decision}`")
    rp.append("")
    rp.append(why)
    rp.append("")
    rp.append(f"**下一步**：{nxt}")
    rp.append("")
    rp.append("## 6. 停止")
    rp.append("")
    rp.append("- 不做完整 Yang 深度估计 / Doppler / f0 误差 / MC / P5")
    rp.append("")
    (OUT / "R3_C2_0_REPORT.md").write_text("\n".join(rp), encoding="utf-8")
    (OUT / "R3_C2_0_GPT_SYNC.md").write_text(
        f"# R3-C2.0\n\n**{decision}**\n\n{why}\n\nnext: {nxt}\n\nerratum NOT_RECOVERED\n",
        encoding="utf-8",
    )

    print("DECISION", decision)
    print(why)
    print(f"DONE {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
