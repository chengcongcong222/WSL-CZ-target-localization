#!/usr/bin/env python3
"""Finish P4.5: run missing scenarios, aggregate, figures, report."""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd

import p4_g4_performance_boundary as p4
import p4_5_parameter_accuracy as p45

OUT = p45.OUT
FIG = p45.FIG
NOW = p45.NOW


def ensure():
    p45.ensure_tables()
    gkin, gfull, n_k, n_f = p45.build_grid()
    return gkin, gfull, n_k, n_f


def main():
    t0 = time.time()
    gkin, gfull, n_k, n_f = ensure()
    print(f"grid kin={n_k} full={n_f}", flush=True)
    (OUT / "P4_5_CONFIG.json").write_text(
        json.dumps(p45.CONFIG, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    trial_path = OUT / "parameter_accuracy_trials.csv"
    existing = pd.read_csv(trial_path) if trial_path.exists() else pd.DataFrame()
    have = set(existing["scenario"].unique()) if len(existing) else set()
    print("have scenarios", have, "n", len(existing), flush=True)

    cache_bear, cache_feat, cache_sig = {}, {}, {}
    parts = [existing] if len(existing) else []
    for key in ["B-L", "B-M", "B-U"]:
        if key in have and len(existing[existing["scenario"] == key]) >= p45.N_TRUTH * p45.N_NOISE:
            print(f"skip {key} already done", flush=True)
            continue
        print(f"[MC] {key} ...", flush=True)
        df_t, pri = p45.run_mc_one_scenario(key, gkin, gfull, n_k, n_f, cache_bear, cache_feat, cache_sig)
        parts.append(df_t)
        trial_df = pd.concat([x for x in parts if len(x)], ignore_index=True)
        trial_df.to_csv(trial_path, index=False, encoding="utf-8-sig")
        print(f"  saved total={len(trial_df)}", flush=True)

    trial_df = pd.read_csv(trial_path)
    print("final trials", len(trial_df), trial_df.groupby(["scenario", "layer"]).size(), flush=True)
    pri = p45.prior_stats()

    print("aggregate ...", flush=True)
    acc = p45.aggregate_accuracy(trial_df, pri)
    acc.to_csv(OUT / "parameter_accuracy_baseline_vs_method.csv", index=False, encoding="utf-8-sig")

    cov_rows = []
    for _, r in acc.iterrows():
        cov_rows.append({
            "scenario": r["scenario"], "parameter": r["parameter"],
            "rc2_coverage": r.get("rc2_coverage"), "rc23_coverage": r.get("rc23_coverage"),
            "rc2_ci95_width": r.get("rc2_ci95_width"), "rc23_ci95_width": r.get("rc23_ci95_width"),
            "rc2_status": r.get("rc2_status"), "rc23_status": r.get("rc23_status"),
            "status": r.get("status"), "rmse_improvement": r.get("rmse_improvement"),
        })
    pd.DataFrame(cov_rows).to_csv(OUT / "parameter_coverage.csv", index=False, encoding="utf-8-sig")

    print("sensitivity B-M ...", flush=True)
    df1, df2, hw = p45.rc3_sensitivity(gfull, "B-M")
    df1.to_csv(OUT / "rc3_1d_sensitivity.csv", index=False, encoding="utf-8-sig")
    df2.to_csv(OUT / "rc3_rv_sensitivity.csv", index=False, encoding="utf-8-sig")
    df2[df2["slice"] == "rz"].to_csv(OUT / "rc3_rz_sensitivity.csv", index=False, encoding="utf-8-sig")
    print("half-widths", hw, flush=True)

    print("figures ...", flush=True)
    p45.write_figA(FIG / "figA_rmse_baseline_vs_method.svg", acc)
    p45.write_figB(FIG / "figB_rmse_improvement.svg", acc)
    p45.write_figC(FIG / "figC_rc3_parameter_sensitivity.svg", df1, df2, hw)

    def fnum(x, nd=3):
        try:
            if x is None:
                return "n/a"
            v = float(x)
            if not np.isfinite(v):
                return "n/a"
            return f"{v:.{nd}f}"
        except Exception:
            return "n/a"

    md = []
    md.append("# P4.5 报告：五维参数统一估计精度基准")
    md.append("")
    md.append(f"UTC：{NOW}  ·  目录：`results/P4_5_parameter_accuracy/`")
    md.append("")
    md.append("本阶段**不**开发新方法；把成果从“候选排除率”转换为 **先验 → RC2 → RC2+RC3** 的参数 RMSE/MAE/覆盖率。")
    md.append("")
    md.append("## 0. 实验设计")
    md.append("")
    md.append(f"- 每场景 **N_truth={p45.N_TRUTH}** × **N_noise={p45.N_NOISE}** = {p45.N_TRUTH*p45.N_NOISE} 次估计")
    md.append("- 真值连续 **off-grid** 采样；M1/M2 共用真值、方位噪声、候选空间与优化框架")
    md.append("- M0=先验中心参考；M1=仅方位；M2=方位+CZ（时间轮廓特征，未知 S(f) 消去）")
    md.append("- 点估计：粗网格 argmin → 局部细化；CI：J≤J_min+Δχ²")
    md.append("- 覆盖率不足 → INVALID/OVERCONFIDENT，不把小 RMSE 称高精度")
    md.append("")
    md.append("## 1. 主表 `parameter_accuracy_baseline_vs_method.csv`")
    md.append("")
    md.append("| 场景 | 参数 | 先验std | RC2 RMSE | RC2覆盖 | RC2+RC3 RMSE | RC2+RC3覆盖 | 改善 | status |")
    md.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for _, r in acc.iterrows():
        md.append(
            f"| {r['scenario']} | {r['parameter']} ({r['unit']}) | {fnum(r['prior_std'])} | "
            f"{fnum(r.get('rc2_rmse'))} | {fnum(r.get('rc2_coverage'),2)} | "
            f"{fnum(r.get('rc23_rmse'))} | {fnum(r.get('rc23_coverage'),2)} | "
            f"{fnum(r.get('rmse_improvement'),3)} | {r.get('status','')} |"
        )
    md.append("")
    md.append("## 2. 论文式表述（B-M 主结果）")
    md.append("")
    for pname in ["r", "v", "theta", "psi", "z"]:
        row = acc[(acc["scenario"] == "B-M") & (acc["parameter"] == pname)]
        if not len(row):
            continue
        row = row.iloc[0]
        st = str(row.get("status", ""))
        if "UNRESOLVED" in st:
            md.append(f"- **{pname}**：RC2 结构性弱/不可辨识；CZ 后 RMSE={fnum(row.get('rc23_rmse'))} {row['unit']}，95%宽={fnum(row.get('rc23_ci95_width'))}，覆盖={fnum(row.get('rc23_coverage'),2)}。")
        elif "INVALID" in st or "OVERCONFIDENT" in st:
            md.append(f"- **{pname}**：RMSE={fnum(row.get('rc23_rmse'))} 但覆盖={fnum(row.get('rc23_coverage'),2)} → **{st}**，不称高精度。")
        else:
            md.append(
                f"- **{pname}**：RC2 RMSE={fnum(row.get('rc2_rmse'))} {row['unit']} → "
                f"CZ 后 {fnum(row.get('rc23_rmse'))} {row['unit']}；"
                f"降低 {fnum(100*(row.get('rmse_improvement') or 0),1)}%；覆盖={fnum(row.get('rc23_coverage'),2)}。"
            )
    md.append("")
    md.append("## 3. CZ 直接敏感度（B-M 切片）")
    md.append("")
    md.append("| 参数 | Δ@+9.21 | 解读 |")
    md.append("| --- | --- | --- |")
    interp = {
        "r": "CZ 对距离直接敏感" if hw.get("r", 99) < 5 else "距离敏感度有限",
        "v": "单独弱，多来自 r–v 耦合" if hw.get("v", 0) > 0.3 else "对速度亦有直接敏感度",
        "z": "剖面平，深度仍难独立约束" if hw.get("z", 0) > 25 else "深度有一定敏感度",
        "psi": "平，航向主要来自 RC2" if hw.get("psi", 0) > 5 else "航向有 CZ 信息",
        "theta": "θ 主要由方位提供" if hw.get("theta", 0) > 2 else "θ 有 CZ 信息",
    }
    for dim in ["r", "v", "z", "psi", "theta"]:
        md.append(f"| {dim} | {fnum(hw.get(dim),4)} | {interp.get(dim,'')} |")
    md.append("")
    md.append("**结论口径**：CZ 直接观测贡献主要来自 **距离 r**；速度改善更多来自距离轨迹耦合；深度/航向改善有限。**不得**用“候选排除率”反推参数精度。")
    md.append("")
    md.append("## 4. 图与文件")
    md.append("")
    md.append("- figures/figA_rmse_baseline_vs_method.svg")
    md.append("- figures/figB_rmse_improvement.svg")
    md.append("- figures/figC_rc3_parameter_sensitivity.svg")
    for pth in sorted(OUT.rglob("*")):
        if pth.is_file() and pth.suffix != ".svg":
            md.append(f"- results/P4_5_parameter_accuracy/{pth.relative_to(OUT)}")
    md.append("")
    md.append("**P4.5 完成后停止；未进入 P5。**")
    (OUT / "P4_5_REPORT.md").write_text("\n".join(md), encoding="utf-8")

    gs = ["# P4.5 — GPT 同步稿", "", f"- UTC: {NOW}",
          f"- trials: {len(trial_df)}",
          "- 主表: parameter_accuracy_baseline_vs_method.csv", "",
          "## B-M 摘录", "",
          "| 参数 | RC2 RMSE | RC2+RC3 RMSE | 改善 | 覆盖M2 | status |",
          "| --- | --- | --- | --- | --- | --- |"]
    for _, r in acc[acc["scenario"] == "B-M"].iterrows():
        gs.append(
            f"| {r['parameter']} | {fnum(r.get('rc2_rmse'))} | {fnum(r.get('rc23_rmse'))} | "
            f"{fnum(r.get('rmse_improvement'),3)} | {fnum(r.get('rc23_coverage'),2)} | {r.get('status','')} |"
        )
    gs += ["", "## CZ 半高宽", "", json.dumps(hw, ensure_ascii=False), "",
           "> P4.5 完成后停止。", ""]
    (OUT / "P4_5_GPT_SYNC.md").write_text("\n".join(gs), encoding="utf-8")

    print(f"DONE in {time.time()-t0:.1f}s", flush=True)
    print(acc.to_string(index=False)[:2500], flush=True)


if __name__ == "__main__":
    main()
