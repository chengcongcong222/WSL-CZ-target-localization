#!/usr/bin/env python3
"""R3-C1 finisher: stricter depth-info decision + figures from existing CSVs."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent / "results" / "R3_C1_depth_motion"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()


def _esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def write_line_svg(path, title, series, xlab, ylab, note=""):
    w, h = 760, 360
    ml, mr, mt, mb = 70, 140, 48, 50
    pw, ph = w - ml - mr, h - mt - mb
    # x may be numeric or categorical labels
    all_x = []
    for s in series:
        all_x.extend(s["x"])
    numeric = all(isinstance(x, (int, float)) and np.isfinite(x) for x in all_x) if all_x else True
    ys = []
    for s in series:
        ys.extend([y for y in s["y"] if y is not None and np.isfinite(y)])
    if not ys:
        ys = [0, 1]
    y0, y1 = min(ys + [0]), max(ys + [0.1])
    cols = ["#b45309", "#0f766e", "#1d4ed8", "#8a8a8a", "#9f1239"]
    p = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
         f'<rect width="{w}" height="{h}" fill="#f7f4ef"/>',
         f'<text x="{w/2}" y="26" text-anchor="middle" font-family="-apple-system,PingFang SC,Microsoft YaHei,sans-serif" font-size="14" font-weight="600">{_esc(title)}</text>',
         f'<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="#fff" stroke="#ccc"/>']
    if numeric and all_x:
        x0, x1 = min(all_x), max(all_x)
        if x1 - x0 < 1e-12:
            x1 = x0 + 1

        def sx(x):
            return ml + (x - x0) / (x1 - x0) * pw
        def sy(y):
            return mt + ph - (y - y0) / (y1 - y0 + 1e-15) * ph
        for i, s in enumerate(series):
            pts = [f"{sx(x):.1f},{sy(y):.1f}" for x, y in zip(s["x"], s["y"]) if y is not None and np.isfinite(y)]
            if len(pts) >= 2:
                p.append(f'<polyline fill="none" stroke="{cols[i%len(cols)]}" stroke-width="1.8" points="{" ".join(pts)}"/>')
            p.append(f'<rect x="{w-mr+8}" y="{mt+8+i*16}" width="12" height="3" fill="{cols[i%len(cols)]}"/>')
            p.append(f'<text x="{w-mr+24}" y="{mt+12+i*16}" font-size="10" font-family="sans-serif">{_esc(s["name"])}</text>')
        p.append(f'<text x="{ml+pw/2}" y="{h-16}" text-anchor="middle" font-size="11" font-family="sans-serif">{_esc(xlab)}</text>')
    else:
        # categorical: bar-like points
        cats = list(dict.fromkeys(str(x) for x in all_x))
        for i, s in enumerate(series):
            for j, xv in enumerate(s["x"]):
                y = s["y"][j]
                if y is None or not np.isfinite(y):
                    continue
                X = ml + (j + 0.5) * pw / max(len(cats), 1) + i * 8
                Y = mt + ph - (y - y0) / (y1 - y0 + 1e-15) * ph
                p.append(f'<circle cx="{X:.1f}" cy="{Y:.1f}" r="4" fill="{cols[i%len(cols)]}"/>')
            p.append(f'<rect x="{w-mr+8}" y="{mt+8+i*16}" width="12" height="3" fill="{cols[i%len(cols)]}"/>')
            p.append(f'<text x="{w-mr+24}" y="{mt+12+i*16}" font-size="10" font-family="sans-serif">{_esc(s["name"])}</text>')
        for j, c in enumerate(cats):
            X = ml + (j + 0.5) * pw / max(len(cats), 1)
            p.append(f'<text x="{X:.1f}" y="{mt+ph+16}" text-anchor="middle" font-size="9" font-family="sans-serif">{_esc(c[:12])}</text>')
    p.append(f'<text x="{ml}" y="{h-4}" font-size="10" fill="#777" font-family="sans-serif">{_esc(note)} {_esc(ylab)}</text></svg>')
    path.write_text("\n".join(p), encoding="utf-8")


def main():
    paper = pd.read_csv(OUT / "paper_reproduction.csv")
    amb = pd.read_csv(OUT / "estd_depth_ambiguity.csv")
    pair = pd.read_csv(OUT / "depth_pair_distinguishability.csv")
    track = pd.read_csv(OUT / "track_length_depth_boundary.csv")
    cmpdf = pd.read_csv(OUT / "single_vs_hla_vs_motion.csv")

    def fnum(x, nd=2):
        try:
            v = float(x)
            return "n/a" if not np.isfinite(v) else f"{v:.{nd}f}"
        except Exception:
            return "n/a"

    # Paper repro: peak on true z is necessary but NOT sufficient if FWHM~full and PSL~1
    def paper_quality(df):
        out = {}
        for m in df["method"].unique():
            d = df[df["method"] == m]
            out[m] = {
                "mean_abs_err": float(d["err_z_m"].abs().mean()),
                "mean_fwhm": float(d["fwhm_m"].mean()),
                "mean_psl": float(d["psl"].mean()),
                "frac_err_lt_15m": float((d["err_z_m"].abs() < 15).mean()),
            }
        return out

    pq = paper_quality(paper)
    # pass if any method has localization: mean FWHM < 40 m OR (err small AND psl < 0.6)
    paper_pass = any(
        v["mean_fwhm"] < 40.0 or (v["frac_err_lt_15m"] > 0.7 and v["mean_psl"] < 0.6)
        for v in pq.values()
    )
    # even if peak lands on truth, flat J is not depth localization
    paper_peak_only = any(v["frac_err_lt_15m"] > 0.7 for v in pq.values())

    d = amb[(amb["source_level"] == "C-S2_235") & (amb["T_s"] == 600.0)]
    notes = {
        "paper_quality": pq,
        "paper_pass_strict": paper_pass,
        "paper_peak_on_truth_but_flat": paper_peak_only,
        "C-S0": "METHOD_NOT_APPLICABLE_WITHOUT_TRACKABLE_LINE",
    }
    if len(d):
        notes["cz_mean_abs_err_m"] = float(d["err_z_m"].abs().mean())
        notes["cz_mean_fwhm_m"] = float(d["fwhm_m"].mean())
        notes["cz_mean_psl"] = float(d["psl"].mean())
        notes["cz_frac_fwhm_lt_40m"] = float((d["fwhm_m"] < 40).mean())
        notes["cz_true_in_peak"] = float((d["err_z_m"].abs() < 25).mean())
    if len(pair):
        notes["pair_mean_D"] = float(pair["D_1minus_corr"].mean())
        notes["pair_max_D"] = float(pair["D_1minus_corr"].max())
        notes["pair_resolvable_frac"] = float(pair["distinguishable_D_gt_0p05"].mean())

    # CZ depth information: not "err=0 self-match", but narrow peak + distinct pair structure
    cz_localized = bool(notes.get("cz_frac_fwhm_lt_40m", 0) >= 0.4 and notes.get("cz_mean_psl", 1) < 0.7)
    cz_pairs = bool(notes.get("pair_resolvable_frac", 0) >= 0.4 and notes.get("pair_mean_D", 0) > 0.05)
    info_present = cz_localized or (notes.get("cz_mean_fwhm_m", 999) < 60 and cz_pairs)

    if not paper_pass:
        decision = "C1_PAPER_REPRO_FAIL"
        why = (
            f"Strict paper-condition localization failed: methods show peak-at-truth but "
            f"flat ambiguity (FWHM~full, PSL~1) → QF/matched-structure not yet a valid "
            f"depth localizer in this implementation. quality={pq}"
        )
        nxt = "fix Zhu Q_F / intensity-structure implementation before CZ transfer claims"
    elif info_present and cz_pairs:
        decision = "C1_PAPER_REPRO_PASS_CZ_DEPTH_CONFIRMED"
        why = (
            f"Paper repro OK; CZ depth peaks reasonably localized "
            f"(FWHM mean={notes.get('cz_mean_fwhm_m')} m, frac<40m={notes.get('cz_frac_fwhm_lt_40m')}, "
            f"PSL={notes.get('cz_mean_psl')}) and pair D resolvable "
            f"(frac={notes.get('pair_resolvable_frac')}, mean D={notes.get('pair_mean_D')})."
        )
        nxt = "eligible for C2 later (not this round)"
    elif paper_pass and info_present:
        decision = "C1_PAPER_REPRO_PASS_CZ_CONDITIONAL"
        why = (
            f"Paper OK; CZ some depth response but pair distinguishability weak "
            f"(pair_frac={notes.get('pair_resolvable_frac')}, FWHM={notes.get('cz_mean_fwhm_m')})."
        )
        nxt = "C2 only on favorable windows if at all"
    elif paper_pass:
        decision = "C1_PAPER_REPRO_PASS_CZ_NOT_TRANSFERABLE"
        why = (
            f"Paper-method logic accepted at algorithm level, but E-STD first-CZ modal track "
            f"intensity does not yield narrow unique depth peaks "
            f"(mean FWHM={notes.get('cz_mean_fwhm_m')} m, PSL={notes.get('cz_mean_psl')}, "
            f"pair mean D={notes.get('pair_mean_D')}, frac={notes.get('pair_resolvable_frac')}). "
            f"Self-match err≈0 must not be read as depth accuracy."
        )
        nxt = "next RC3-C candidate in order: Yang 2015 SA beamforming → HLA modal/k-spectrum → Emmetière 2019"
    else:
        decision = "C1_PAPER_REPRO_FAIL"
        why = "implementation not validated"
        nxt = "fix paper method"

    # figures
    p1 = paper[(paper["z_true_m"] == 55.0)]
    series = []
    for m in p1["method"].unique():
        dd = p1[p1["method"] == m].sort_values("r1_m")
        series.append({"name": f"{m[:18]} FWHM", "x": dd["r1_m"].tolist(), "y": dd["fwhm_m"].tolist()})
        series.append({"name": f"{m[:18]} |err|", "x": dd["r1_m"].tolist(), "y": dd["err_z_m"].abs().tolist()})
    write_line_svg(FIG / "fig1_paper_repro.svg", "图1  论文条件：FWHM/|err| vs 轨迹终点 (z=55m)",
                   series, "r1 (m)", "m", note="峰值贴真值≠深度可定位（看FWHM/PSL）")
    if len(d):
        write_line_svg(FIG / "fig2_estd_depth_ambiguity.svg",
                       "图2  E-STD 深度模糊 (C-S2 235Hz,T=600s)",
                       [{"name": "FWHM", "x": d["z_true_m"].tolist(), "y": d["fwhm_m"].tolist()},
                        {"name": "PSL", "x": d["z_true_m"].tolist(), "y": d["psl"].tolist()}],
                       "z_true m", "m / ratio", note="对照P4.5 z UNRESOLVED")
    if len(pair):
        write_line_svg(FIG / "fig3_depth_pair_distinguishability.svg",
                       "图3  深度对 D=1−corr (235Hz,T=600s)",
                       [{"name": "D", "x": list(range(len(pair))), "y": pair["D_1minus_corr"].tolist()}],
                       "pair idx", "D", note="D低→轨迹声强结构几乎重合")
    if len(track):
        d4 = track[track["z_true_m"] == 200].sort_values("track_len_m")
        if len(d4):
            write_line_svg(FIG / "fig4_track_length_boundary.svg",
                           "图4  轨迹长度 U*T vs FWHM (z=200m,235Hz)",
                           [{"name": "FWHM", "x": d4["track_len_m"].tolist(), "y": d4["fwhm_m"].tolist()}],
                           "track m", "m", note="")
    if len(cmpdf):
        write_line_svg(FIG / "fig5_static_vs_motion.svg",
                       "图5  单通道运动 / HLA运动 / HLA静态",
                       [{"name": r["config"], "x": [float(i)], "y": [r["fwhm_m"]]} for i, r in cmpdf.iterrows()],
                       "config index", "FWHM m", note="")

    dec = {
        "rc3c1_decision": decision,
        "why": why,
        "next_step": nxt,
        "notes": notes,
        "created_utc": NOW,
        "stop": "after R3-C1; no C2/P5/RC3-B",
        "b1_status": "PERMANENTLY_CLOSED B1_NO_STABLE_MULTIPATH_IDENTITY",
    }
    (OUT / "R3_C1_DECISION.json").write_text(json.dumps(dec, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    rp = []
    rp.append("# R3-C1 报告：运动/声强结构深度（Zhu 2023）")
    rp.append("")
    rp.append(f"UTC：{NOW}")
    rp.append("")
    rp.append("B1 已关闭：当前 CZ 仅 1 条能量可观测多径，HLA-MMAC 不成立。")
    rp.append("")
    rp.append("## 0. 锚定")
    rp.append("")
    rp.append("- **Zhu et al. 2023** JASA-EL：深海 HLA + 移动源；improved Fourier / matched intensity structure")
    rp.append("- 观测为**轨迹声强结构**，非 P4.5 CZ 轮廓")
    rp.append("- **C-S0**：`METHOD_NOT_APPLICABLE_WITHOUT_TRACKABLE_LINE`")
    rp.append("")
    rp.append("## 1. 论文条件复现（严格口径）")
    rp.append("")
    rp.append("峰值落在 true z **不等于**深度可定位；必须同时看 **FWHM / PSL**。")
    rp.append("")
    rp.append("| method | mean\\|err\\| m | mean FWHM m | mean PSL | frac\\|err\\|<15m |")
    rp.append("| --- | --- | --- | --- | --- |")
    for m, v in pq.items():
        rp.append(f"| {m} | {fnum(v['mean_abs_err'])} | {fnum(v['mean_fwhm'])} | {fnum(v['mean_psl'])} | {fnum(v['frac_err_lt_15m'])} |")
    rp.append("")
    rp.append(f"- 严格 paper_pass = **{paper_pass}**（当前实现下 ambiguity 几乎铺满深度窗 → **复现未形成有效深度定位**）")
    rp.append("")
    rp.append("## 2. E-STD 第一 CZ（模态 I' + Zhu 匹配）")
    rp.append("")
    rp.append(f"- C-S2 235 Hz, T=600 s：mean FWHM=**{fnum(notes.get('cz_mean_fwhm_m'))} m**，mean PSL=**{fnum(notes.get('cz_mean_psl'))}**，frac FWHM<40 m=**{fnum(notes.get('cz_frac_fwhm_lt_40m'))}**")
    rp.append(f"- 深度对 mean D=**{fnum(notes.get('pair_mean_D'))}**，可分比例(D>0.05)=**{fnum(notes.get('pair_resolvable_frac'))}**")
    rp.append(f"- 自匹配 err≈0 **不得**解读为深度 RMSE；结构几乎重合时 corr 峰在 true z 但无分辨力")
    rp.append("")
    rp.append("### 与 P4.5 对照")
    rp.append("")
    rp.append("- P4.5：z 弱/UNRESOLVED")
    rp.append(f"- R3-C1：在 matched 环境、无噪、CW 轨迹声强下，CZ 深度峰仍 **宽 / 旁瓣高** → 相对 P4.5 **未形成可用的独立深度约束**")
    rp.append("")
    rp.append("## 3. R3-C1 判定")
    rp.append("")
    rp.append(f"### `{decision}`")
    rp.append("")
    rp.append(why)
    rp.append("")
    rp.append(f"**下一步**：{nxt}")
    rp.append("")
    rp.append("## 4. 停止")
    rp.append("")
    rp.append("- 不进 C2 / P5 / RC3-B")
    rp.append("- 若迁移失败，RC3-C 下一候选按序：Yang 2015 → HLA 模态/波数谱 → Emmetière 2019")
    rp.append(f"- **R3-C1 完成后停止**")
    rp.append("")
    (OUT / "R3_C1_REPORT.md").write_text("\n".join(rp), encoding="utf-8")

    gs = [
        "# R3-C1 — GPT 同步", "",
        f"- **判定：{decision}**",
        f"- {why}",
        f"- 下一步：{nxt}", "",
        f"paper_pass_strict={paper_pass} quality={json.dumps(pq, ensure_ascii=False)}", "",
        f"CZ FWHM={notes.get('cz_mean_fwhm_m')} PSL={notes.get('cz_mean_psl')} pairD={notes.get('pair_mean_D')} pair_frac={notes.get('pair_resolvable_frac')}", "",
        "C-S0: METHOD_NOT_APPLICABLE_WITHOUT_TRACKABLE_LINE", "",
        "停止：无 C2/P5/RC3-B。", "",
    ]
    (OUT / "R3_C1_GPT_SYNC.md").write_text("\n".join(gs), encoding="utf-8")
    print("DECISION", decision)
    print(why)


if __name__ == "__main__":
    main()
