#!/usr/bin/env python3
"""Post-process P3 summaries after mechanism-layer bugfix."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent / "results" / "P3_RC3_increment"
sum_df = pd.read_csv(OUT / "rc2_vs_rc3_summary.csv")
selected = pd.read_csv(OUT / "selected_hard_pairs.csv")


def mean_rej(df, source, layer, mech=None):
    d = df[(df["source"] == source) & (df["layer"] == layer) & (df["snr_db"].astype(str) == "noiseless")]
    if mech:
        d = d[d["mechanism"] == mech]
    if not len(d):
        return float("nan")
    return float(pd.to_numeric(d["reject_rate_wrong"], errors="coerce").mean())


rows = []
for mech in ["A", "B", "C"]:
    for src in ["S0", "S1", "S2"]:
        rows.append({
            "mechanism": mech,
            "source": src,
            "mean_reject_RC2": mean_rej(sum_df, src, "RC2", mech),
            "mean_reject_RC2CZ": mean_rej(sum_df, src, "RC2+CZ", mech),
            "mean_reject_RC2Dop": mean_rej(sum_df, src, "RC2+Doppler", mech),
            "mean_reject_combined": mean_rej(sum_df, src, "RC2+CZ+Doppler", mech),
        })
mech_df = pd.DataFrame(rows)
mech_df.to_csv(OUT / "mechanism_summary.csv", index=False, encoding="utf-8-sig")
print(mech_df.to_string(index=False))

layers = ["RC2", "RC2+CZ", "RC2+Doppler", "RC2+CZ+Doppler"]
rej = {src: {ly: mean_rej(sum_df, src, ly) for ly in layers} for src in ["S0", "S1", "S2"]}
print("REJ", json.dumps(rej, indent=2))

r0c2 = rej["S0"]["RC2"]
r0cz = rej["S0"]["RC2+CZ"]
cz_gain = r0cz - r0c2
line_vals = []
for s in ["S1", "S2"]:
    for ly in ["RC2+Doppler", "RC2+CZ+Doppler"]:
        v = rej[s][ly]
        if np.isfinite(v):
            line_vals.append(v)
line_max = max(line_vals) if line_vals else float("nan")

if np.isfinite(cz_gain) and cz_gain >= 0.15:
    g3 = "PROPAGATION_INCREMENT_CONFIRMED"
    why = (
        f"CZ increased wrong-candidate rejection on S0 hard pairs by "
        f"delta={cz_gain:.3f} (RC2={r0c2:.3f} -> RC2+CZ={r0cz:.3f}). "
        f"Doppler control branch is weak on constant-v_rad collinear pairs "
        f"(absorbed by unknown source frequency bias)."
    )
elif np.isfinite(line_max) and line_max >= 0.2:
    g3 = "LINE_INCREMENT_ONLY"
    why = f"CZ gain weak (delta={cz_gain}); line/Doppler max={line_max:.3f}."
else:
    g3 = "NO_USEFUL_INCREMENT_IN_TESTED_MODEL"
    why = f"Under frozen E-STD model: CZ delta={cz_gain}, line max={line_max}."

c0 = rej["S0"]
c1 = rej["S1"]
c2 = rej["S2"]
client = {
    "S0": (
        f"无稳定谱线时：RC2错误候选排除率约 {c0['RC2']:.0%}；"
        f"加入会聚区传播后约 {c0['RC2+CZ']:.0%}；"
        f"Doppler支路对恒定径向速度差不可用（被未知源频偏置吸收）。"
    ),
    "S1": (
        f"典型机械谱线时：RC2约 {c1['RC2']:.0%}；CZ约 {c1['RC2+CZ']:.0%}；"
        f"bearing+Doppler约 {c1['RC2+Doppler']:.0%}；联合约 {c1['RC2+CZ+Doppler']:.0%}。"
        f"（线谱频点少，CZ轮廓自由度低于S0连续谱。）"
    ),
    "S2": (
        f"理想稳定多谱线时：RC2约 {c2['RC2']:.0%}；CZ约 {c2['RC2+CZ']:.0%}；"
        f"bearing+Doppler约 {c2['RC2+Doppler']:.0%}；联合约 {c2['RC2+CZ+Doppler']:.0%}。"
        f"（S2为上界工况，非真实UUV谱。）"
    ),
}

meta = {
    "g3_result": g3,
    "g3_why": why,
    "rej_rc2_s0": r0c2,
    "rej_cz_s0": r0cz,
    "cz_gain_s0": cz_gain,
    "rej_matrix": rej,
    "line_gain_max": line_max,
    "client_lines": client,
    "n_selected_pairs": int(len(selected)),
    "notes": [
        "RC2-only rejection is not always 0: some hard pairs retain small bearing RMSE.",
        "CZ layer: multi-freq amplitude profile with unknown S(f) LS elimination + time-normalized contrast.",
        "Doppler: constant v_rad differences act as unknown f0 bias; only time-varying Doppler is informative.",
        "KRAKEN not installed on host; modal theory model recorded in CONFIG (no Bellhop audit).",
        "Mode-count absolute amplitude is not fully converged; range-profile CONTRAST is the stable discrimination quantity.",
        "P2 CRLB rank-deficient = UNDEFINED/UNBOUNDED; maneuver wording downgraded; P2 not re-run.",
    ],
}
(OUT / "g3_decision.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
print("G3", g3)
print(json.dumps(client, ensure_ascii=False, indent=2))

# Patch reports with corrected numbers
rep = OUT / "P3_RC3_REPORT.md"
if rep.exists():
    txt = rep.read_text(encoding="utf-8")
    # rewrite key result paragraphs by appending a correction section if needed
    corr = []
    corr.append("\n## 5b. 汇总修正（分层统计）\n")
    corr.append("机制 × 源 × 层的正确平均错误排除率：\n")
    corr.append("| 机制 | 源 | RC2 | RC2+CZ | +Doppler | combined |")
    corr.append("| --- | --- | --- | --- | --- | --- |")
    for _, row in mech_df.iterrows():
        def f(v):
            try:
                if v is None or (isinstance(v, float) and not np.isfinite(v)):
                    return "n/a"
                return f"{float(v):.3f}"
            except Exception:
                return "n/a"
        corr.append(
            f"| {row['mechanism']} | {row['source']} | {f(row['mean_reject_RC2'])} | {f(row['mean_reject_RC2CZ'])} | "
            f"{f(row['mean_reject_RC2Dop'])} | {f(row['mean_reject_combined'])} |"
        )
    corr.append("")
    corr.append(f"**G3 = `{g3}`**\n")
    corr.append(why + "\n")
    corr.append("### Q3 甲方可理解表述（修正后）\n")
    for k, v in client.items():
        corr.append(f"- **{v}**")
    corr.append("")
    corr.append("### 解释要点\n")
    corr.append("1. **CZ 主增量在 S0 连续谱**：多频幅度轮廓 + 未知 S(f) 消元后，不同 r(t) 轨迹可分。")
    corr.append("2. **S1/S2 线谱频点少**，CZ 轮廓自由度低于连续谱，排除率低于 S0；Doppler 对恒定 Δv_rad 的 r–v 对几乎无效（未知 f0 偏置可吸收常数频移）。")
    corr.append("3. **Q2**：传播信息本质上随距离剖面变化；若 ψ/θ0 差异不改变 r(t) 与声场轮廓，CZ 对该类歧义**无明显独立增量**——属合法理论结果。")
    corr.append("4. **Q4**：本模型下增量主要落在**距离维**及由距离轨迹区分的速度补偿对；深度固定不报增量。")
    if "## 5b. 汇总修正" not in txt:
        rep.write_text(txt + "\n" + "\n".join(corr) + "\n", encoding="utf-8")

# Patch GPT sync
gpath = OUT / "P3_RC3_GPT_SYNC.md"
gpath.write_text(
    "\n".join([
        "# P3 RC3 — GPT 同步稿",
        "",
        f"- G3 判定: **{g3}**",
        f"- 原因: {why}",
        f"- 入选困难候选: {len(selected)} 对（A/B/C）",
        "",
        "## P2 口径（已改，未重跑）",
        "",
        "- 秩亏 CRLB(r) = **UNDEFINED/UNBOUNDED**（禁止 0.000 对外）",
        "- 小机动：改变局部几何，**不足以形成实用距离约束**",
        "- RC2 角色：不是远距离测距，而是约束方位/航向并给出 r–v 困难候选",
        "",
        "## P3 主结果（分层修正后）",
        "",
        f"- {client['S0']}",
        f"- {client['S1']}",
        f"- {client['S2']}",
        f"- CZ 增量 Δ(S0) = {cz_gain:.3f}（RC2 {r0c2:.3f} → CZ {r0cz:.3f}）",
        "",
        "## 机制摘要",
        "",
        mech_df.to_string(index=False),
        "",
        "## 下一轮（GPT）",
        "",
        "> 判定 G3 结论是否成立，冻结研究内容三的阶段性表述；**不**自动进入 P4。",
        "",
    ]),
    encoding="utf-8",
)
print("reports patched")
print("files", sorted(p.name for p in OUT.iterdir()))
