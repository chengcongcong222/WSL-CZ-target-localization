#!/usr/bin/env python3
"""R3-C2.2A: Yang 2015 formula evidence recovery only. No simulation, no estimator."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "R3_C2_Yang_SA_depth" / "R3_C2_2A"
OUT.mkdir(parents=True, exist_ok=True)
NOW = datetime.now(timezone.utc).isoformat()

PRIMARY = (
    "T. C. Yang, Source depth estimation based on synthetic aperture beamfoming "
    "for a moving source, JASA 138(3), 1678-1686 (2015). DOI 10.1121/1.4929748. "
    "PMID 26428805."
)
ERRATUM = (
    "Erratum: JASA 144(6), 3075 (2018). DOI 10.1121/1.5081712. PMID 30599696."
)

RETRIEVAL = {
    "doi_main": "10.1121/1.4929748",
    "doi_erratum": "10.1121/1.5081712",
    "pmid_main": "26428805",
    "pmid_erratum": "30599696",
    "openalex_main": "closed; has_fulltext=false; has_content.pdf=false; any_repository_has_fulltext=false",
    "unpaywall_main": "is_oa=false; oa_status=closed; has_repository_copy=false; oa_locations=[]",
    "unpaywall_erratum": "is_oa=false; oa_status=closed; has_repository_copy=false",
    "semantic_scholar": "openAccessPdf empty/CLOSED; abstract elided by publisher",
    "crossref_pdf_main": "https://pubs.aip.org/asa/jasa/article-pdf/138/3/1678/15317235/1678_1_online.pdf -> HTTP 403",
    "crossref_pdf_erratum": "https://pubs.aip.org/asa/jasa/article-pdf/144/6/3075/14737039/3075_1_online.pdf -> HTTP 403",
    "aip_landing": "HTTP 403",
    "pubmed_html": "cookie wall",
    "europepmc": "abstract only; hasPDF=N; inPMC=N; Subscription required",
    "core": "totalHits=1 metadata only; fullText Not available for public API users",
    "researchgate_282528816": "HTTP 403",
    "hindawi_liang_7824671_pdf": "HTTP 403 (secondary literature)",
    "wiley_liang_pdf": "HTTP 403 (secondary literature)",
    "local_user_pdfs": "none matching yang/4929748/5081712/JASA",
    "sci_hub": "NOT USED (unauthorized redistribution)",
    "abstract_main_recovered": True,
    "abstract_erratum_recovered": False,
    "eq_1_to_9_fulltext": "NOT_OBTAINED",
    "erratum_formula_text_this_session": "NOT_REOBTAINED_FROM_PDF; forms carried from prior recovery note + protocol/user restatement",
}

# Abstract is evidence for framework words only, never for equation forms.
ABSTRACT_CLAIMS = [
    "single hydrophone; moving CW source",
    "synthetic aperture created by source motion",
    "signal at each range steered by range-dependent phase relative to starting point",
    "range increment (aperture) from Doppler shift estimated from data, knowing original signal frequency",
    "depth estimated from beam output, assuming mode depth functions from nominal sound speed and bottom profile",
    "PLL removes fast random phase; leaves deterministic range-dependent phase due to range change",
    "if vertical array covers depth span of interest, beam output can estimate depth without acoustic environment knowledge",
]

# Eq.(5) original printed form is NOT recovered; only Erratum-corrected form is documented.
EQ_ROWS = [
    {
        "equation_id": "Eq.(1)",
        "original_formula": "NOT_RECOVERED",
        "symbol_definition": "NOT_RECOVERED",
        "physical_meaning": "synthetic aperture beamforming integral/sum; steering phase (Erratum cites this eq for δ synchronization)",
        "input_required": "NOT_RECOVERED",
        "output": "NOT_RECOVERED",
        "used_for_depth_estimation": "UNKNOWN (likely SA beam construction, not depth map itself)",
        "primary_source_page": "JASA 138, 1678-1686 — full text not obtained",
        "confidence": "NOT_RECOVERED",
        "key_lock": "A. SA integral/sum + steering phase",
    },
    {
        "equation_id": "Eq.(2)",
        "original_formula": "NOT_RECOVERED",
        "symbol_definition": "NOT_RECOVERED",
        "physical_meaning": "NOT_RECOVERED",
        "input_required": "NOT_RECOVERED",
        "output": "NOT_RECOVERED",
        "used_for_depth_estimation": "NOT_RECOVERED",
        "primary_source_page": "JASA 138, 1678-1686 — full text not obtained",
        "confidence": "NOT_RECOVERED",
        "key_lock": "",
    },
    {
        "equation_id": "Eq.(3)",
        "original_formula": "NOT_RECOVERED",
        "symbol_definition": "NOT_RECOVERED",
        "physical_meaning": "NOT_RECOVERED",
        "input_required": "NOT_RECOVERED",
        "output": "NOT_RECOVERED",
        "used_for_depth_estimation": "NOT_RECOVERED",
        "primary_source_page": "JASA 138, 1678-1686 — full text not obtained",
        "confidence": "NOT_RECOVERED",
        "key_lock": "",
    },
    {
        "equation_id": "Eq.(4)",
        "original_formula": "NOT_RECOVERED",
        "symbol_definition": "NOT_RECOVERED",
        "physical_meaning": "NOT_RECOVERED",
        "input_required": "NOT_RECOVERED",
        "output": "NOT_RECOVERED",
        "used_for_depth_estimation": "NOT_RECOVERED",
        "primary_source_page": "JASA 138, 1678-1686 — full text not obtained",
        "confidence": "NOT_RECOVERED",
        "key_lock": "",
    },
    {
        "equation_id": "Eq.(5)",
        "original_formula": (
            "PRINTED FORM NOT_RECOVERED; "
            "Erratum-corrected form: g(k_m, z_r) = b_m * phi_m(z_r) * exp(i * k_m * delta)"
        ),
        "symbol_definition": "b_m contents NOT_RECOVERED; phi_m(z_r) mode depth function at receiver; k_m modal wavenumber; delta initial range offset (Erratum)",
        "physical_meaning": (
            "modal spectral coefficient g at wavenumber peak k_r=k_m; "
            "Erratum: missing exp(i k_m delta) in original relative to true data range synchronization"
        ),
        "input_required": "g spectrum or beam output at k_m; phi_m(z_r); k_m; delta",
        "output": "complex modal coefficient g(k_m, z_r)",
        "used_for_depth_estimation": "YES — feeds depth ambiguity / matching (exact map via Eq.(6) NOT_RECOVERED)",
        "primary_source_page": "Erratum JASA 144, 3075 (2018) for corrected form; original Eq.(5) print not obtained",
        "confidence": "ERRATUM_FORM_DOCUMENTED; ORIGINAL_PRINT NOT_RECOVERED; b_m PARTIAL/UNKNOWN",
        "key_lock": "B. why g=b_m phi_m(z_r) at k_r=k_m; Erratum phase factor",
    },
    {
        "equation_id": "Eq.(6)",
        "original_formula": "NOT_RECOVERED",
        "symbol_definition": "NOT_RECOVERED",
        "physical_meaning": "source-depth ambiguity / estimator (Erratum says δ phase alters this function)",
        "input_required": "NOT_RECOVERED",
        "output": "NOT_RECOVERED",
        "used_for_depth_estimation": "YES (per abstract/erratum context) — exact form NOT_RECOVERED",
        "primary_source_page": "JASA 138 + Erratum cite — full text not obtained",
        "confidence": "NOT_RECOVERED",
        "key_lock": "C. concrete source-depth ambiguity / estimator",
    },
    {
        "equation_id": "Eq.(7)",
        "original_formula": "NOT_RECOVERED",
        "symbol_definition": "NOT_RECOVERED",
        "physical_meaning": "NOT_RECOVERED",
        "input_required": "NOT_RECOVERED",
        "output": "NOT_RECOVERED",
        "used_for_depth_estimation": "NOT_RECOVERED",
        "primary_source_page": "JASA 138, 1678-1686 — full text not obtained",
        "confidence": "NOT_RECOVERED",
        "key_lock": "",
    },
    {
        "equation_id": "Eq.(8)",
        "original_formula": "NOT_RECOVERED",
        "symbol_definition": "NOT_RECOVERED",
        "physical_meaning": "NOT_RECOVERED",
        "input_required": "NOT_RECOVERED",
        "output": "NOT_RECOVERED",
        "used_for_depth_estimation": "NOT_RECOVERED",
        "primary_source_page": "JASA 138, 1678-1686 — full text not obtained",
        "confidence": "NOT_RECOVERED",
        "key_lock": "",
    },
    {
        "equation_id": "Eq.(9)",
        "original_formula": "NOT_RECOVERED",
        "symbol_definition": "NOT_RECOVERED",
        "physical_meaning": "related to Eq.(6); array/data condition (Erratum cites this eq)",
        "input_required": "NOT_RECOVERED",
        "output": "NOT_RECOVERED",
        "used_for_depth_estimation": "UNKNOWN — possibly VLA direct beam-depth path from abstract",
        "primary_source_page": "JASA 138 + Erratum cite — full text not obtained",
        "confidence": "NOT_RECOVERED",
        "key_lock": "D. relation to Eq.(6); receiver-array / data condition",
    },
]

SYM_ROWS = [
    (
        "k_m / k_{r,m}",
        "modal horizontal wavenumber",
        "Erratum + abstract framework",
        "NAME_ONLY",
        "spectral peak location of g; NOT a recovered equation definition",
    ),
    (
        "phi_m(z)",
        "mode depth function from nominal environment",
        "abstract: mode depth functions from nominal sound speed and bottom profile",
        "NAME_ONLY",
        "which normalization (unit peak / flux / etc.) NOT_RECOVERED",
    ),
    (
        "b_m",
        "modal coefficient multiplying phi_m(z_r) in g(k_m,z_r)",
        "Erratum form g = b_m phi_m(z_r) exp(i k_m delta)",
        "SYMBOL_PRESENT_FACTORS_UNKNOWN",
        "whether includes phi_m(z_s), 1/sqrt(k_m), source level, attenuation, range spreading, other mode factors: NOT_RECOVERED",
    ),
    (
        "delta δ",
        "initial range offset: δ = r1 − r1_data = r2 − r2_data",
        "Erratum 2018 (prior recovery note + protocol restatement)",
        "RECOVERED_FROM_ERRATUM",
        "enters g as exp(i k_m δ); changes depth ambiguity, not wavenumber spectral density; δ=0 ↔ original simulation/data processing synchronization; joint vs staged search with z: NOT_RECOVERED",
    ),
    (
        "z_r",
        "receiver depth",
        "framework",
        "NAME_ONLY",
        "",
    ),
    (
        "z_s / source depth",
        "estimated parameter",
        "framework",
        "NAME_ONLY",
        "",
    ),
    (
        "f0",
        "original CW frequency used with Doppler for aperture increment",
        "abstract",
        "ABSTRACT_ONLY",
        "exact Doppler–Δr formula NOT_RECOVERED",
    ),
    (
        "Δr(t) / range increment",
        "synthetic aperture sample coordinate (radial range change)",
        "abstract: range increment determined by Doppler",
        "ABSTRACT_ONLY",
        "exact integration variable in Eq.(1) NOT_RECOVERED",
    ),
    (
        "PLL",
        "phase locked loop removing fast random phase",
        "abstract",
        "ABSTRACT_ONLY",
        "input/output equations NOT_RECOVERED",
    ),
    (
        "g(k_m, z_r)",
        "beam/spectral output at modal wavenumber for assumed or true z_r",
        "Erratum-corrected form only",
        "ERRATUM_FORM",
        "whether g is complex beam output, spectral amplitude, or normalized modal coefficient: NOT_RECOVERED",
    ),
    (
        "beam output",
        "observable used for depth estimation",
        "abstract",
        "ABSTRACT_ONLY",
        "exact mathematical object compared to phi_m(z) NOT_RECOVERED",
    ),
    (
        "source level / unknown amplitude",
        "nuisance amplitude of CW source",
        "not mentioned in abstract",
        "NOT_RECOVERED",
        "normalization / relative amplitudes / analytic cancel / correlation ambiguity: NOT_RECOVERED",
    ),
]


def write_equation_table() -> None:
    pd.DataFrame(EQ_ROWS).to_csv(
        OUT / "YANG2015_EQUATION_TABLE.csv", index=False, encoding="utf-8-sig"
    )


def write_symbol_table() -> None:
    rows = [
        {
            "symbol": s,
            "definition": d,
            "evidence": e,
            "status": st,
            "open_questions": oq,
        }
        for s, d, e, st, oq in SYM_ROWS
    ]
    pd.DataFrame(rows).to_csv(
        OUT / "YANG2015_SYMBOL_TABLE.csv", index=False, encoding="utf-8-sig"
    )


def write_anchor() -> None:
    lines = [
        "# YANG2015_FULL_METHOD_ANCHOR",
        "",
        f"UTC: {NOW}",
        "",
        "## 三态分离（冻结）",
        "",
        "- KRAKEN 数据基线：**VALIDATED** (`C2_1_PARSER_VALIDATED`)",
        "- Yang 原方法公式恢复：**PARTIAL** (`YANG_FORMULA_RECOVERY_PARTIAL`)",
        "- Yang 在 E-STD 是否可用：**UNDECIDED** (`YANG_ROUTE_UNDECIDED`)",
        "",
        "## 文献",
        "",
        f"- 主文：{PRIMARY}",
        f"- Erratum：{ERRATUM}",
        "",
        "## 检索状态（本轮）",
        "",
        "- OpenAlex / Unpaywall：**closed access**，无机构库全文",
        "- Crossref 给出官方 PDF 直链，AIP/Silverchair：**HTTP 403**",
        "- EuropePMC：仅摘要，`hasPDF=N`",
        "- CORE：仅元数据",
        "- ResearchGate 282528816：**HTTP 403**",
        "- 次级 Liang et al. 2018 PDF：**HTTP 403**（本轮不做次级公式替代）",
        "- 本地用户目录：无相关 PDF",
        "- Sci-Hub：**未使用**",
        "- **Eq.(1)–Eq.(9) 主文全文：NOT_OBTAINED**",
        "",
        "## 仅摘要级框架（禁止当作公式）",
        "",
        *[f"- {x}" for x in ABSTRACT_CLAIMS],
        "",
        "## Erratum 级已锁定（非主文原式）",
        "",
        "- δ = r1 − r1_data = r2 − r2_data",
        "- g(k_m, z_r) = b_m φ_m(z_r) exp(i k_m δ)",
        "- exp(i k_m δ) 不改变 wavenumber spectral density，**会**改变 source-depth ambiguity",
        "- 原文模拟/数据处理对应 δ = 0 同步条件",
        "- “不需要源距离”表述不严格成立",
        "",
        "## 禁止事项",
        "",
        "- 禁止用摘要/次级文献重构 Eq.(1)–(9)",
        "- 禁止自造 Yang estimator / depth score",
        "- 禁止把 A_m = φ_m(z_s)φ_m(z_r) 当作 Yang 正式公式",
        "- 禁止输出 Yang 性能判定（APERTURE_LIMITED / DEPTH_AMBIGUOUS 等）",
        "",
    ]
    (OUT / "YANG2015_FULL_METHOD_ANCHOR.md").write_text("\n".join(lines), encoding="utf-8")


def write_erratum_mapping() -> None:
    lines = [
        "# YANG2015_ERRATUM_MAPPING",
        "",
        f"UTC: {NOW}",
        "",
        "Erratum 引用主文 **Eq.(1), Eq.(5), Eq.(6), Eq.(9)**。",
        "",
        "| 位置 | Erratum/协议要求回答 | 本轮证据 | 状态 |",
        "| --- | --- | --- | --- |",
        "| Eq.(1) | SA 如何积分/求和；steering phase | 仅摘要：range-dependent phase relative to starting point | **NOT_RECOVERED** |",
        "| Eq.(5) | 为何 k_r=k_m 时 g=b_m φ_m(z_r)；Erratum 为何加 exp(i k_m δ) | Erratum 形式已记；原推导与原印刷式 | **FORM_FROM_ERRATUM; DERIVATION NOT_RECOVERED** |",
        "| Eq.(6) | 具体 source-depth ambiguity / estimator | 仅知 δ 相位会改变该函数 | **NOT_RECOVERED** |",
        "| Eq.(9) | 与 Eq.(6) 关系；接收阵/数据条件 | 仅摘要有 VLA 直接估深路径 | **NOT_RECOVERED** |",
        "",
        "## δ 关系",
        "",
        "- 定义：δ = r1 − r1_data = r2 − r2_data（初始距离偏移，两时刻同步差）",
        "- 进入位置：g(k_m, z_r) 的相位因子 exp(i k_m δ)（Erratum）",
        "- δ = 0：原文 simulation / data processing 隐含同步",
        "- δ 未知时是 (z, δ) 联合搜索还是先 δ 后 z：**NOT_RECOVERED**（不得推测）",
        "",
        "## 证据来源分级",
        "",
        "- Erratum 形式：来自先前会话 `YANG2015_ERRATUM_NOTE.md`（标 RECOVERED）+ 本任务书/用户复述；本轮未能重下 PDF 复核字节级原文",
        "- 主文 Eq 印刷式：**全部 NOT_RECOVERED**",
        "",
    ]
    (OUT / "YANG2015_ERRATUM_MAPPING.md").write_text("\n".join(lines), encoding="utf-8")


def write_observable_chain() -> None:
    lines = [
        "# YANG2015_OBSERVABLE_CHAIN",
        "",
        f"UTC: {NOW}",
        "",
        "## 必须区分的 observable（主文未锁定前均不得实现）",
        "",
        "| 候选对象 | 本轮状态 |",
        "| --- | --- |",
        "| complex beam output | **NOT_RECOVERED**（摘要出现 beam output 一词） |",
        "| spectral amplitude | **NOT_RECOVERED** |",
        "| spectral energy | **NOT_RECOVERED** |",
        "| normalized modal coefficient | **NOT_RECOVERED** |",
        "| modal shading coefficient | **NOT_RECOVERED** |",
        "",
        "摘要仅说明 depth 从 **beam output** 估计，并用 nominal 环境的 **mode depth functions**。",
        "深度匹配到底比较哪一个量：**NOT_RECOVERED**。",
        "",
        "## b_m 分解（决定 KRAKEN mode table 如何变成 Yang observable）",
        "",
        "必须明确 b_m 是否包含：",
        "",
        "- φ_m(z_s)：**NOT_RECOVERED**",
        "- 1/sqrt(k_m)：**NOT_RECOVERED**",
        "- source level：**NOT_RECOVERED**",
        "- attenuation：**NOT_RECOVERED**",
        "- range spreading：**NOT_RECOVERED**",
        "- other mode factors：**NOT_RECOVERED**",
        "",
        "**禁止**继续把 `A_m = φ_m(z_s) φ_m(z_r)` 写成 Yang 正式公式。",
        "",
        "## source-level nuisance",
        "",
        "- 是否归一化：**NOT_RECOVERED**",
        "- 是否只使用模态相对幅度：**NOT_RECOVERED**",
        "- 是否存在比例因子解析消除：**NOT_RECOVERED**",
        "- 是否通过相关/内积形成 ambiguity：**NOT_RECOVERED**",
        "",
        "## Doppler / PLL（只恢复用途，不仿真）",
        "",
        "| 项目 | 摘要级 | 公式级 |",
        "| --- | --- | --- |",
        "| Doppler → range increment Δr | 由数据估计 Doppler，已知原始频率 f0 | **NOT_RECOVERED** |",
        "| PLL 输入 | 含随机+确定性相位的实数据相位 | **NOT_RECOVERED** |",
        "| PLL 输出 | 去掉快变随机分量后的确定性 range-dependent 相位 | **NOT_RECOVERED** |",
        "",
        "## 次级文献",
        "",
        "- Liang et al. 2018, Match-Mode Autoregressive Method..., DOI 10.1155/2018/7824671：仅可交叉核验对 Yang 的描述与 SA 孔径限制；**本轮 PDF 403，未做公式对照**",
        "- 不得用 AR/Hankel 公式冒充 Yang 2015 原公式",
        "",
    ]
    (OUT / "YANG2015_OBSERVABLE_CHAIN.md").write_text("\n".join(lines), encoding="utf-8")


def write_freeze() -> None:
    lines = [
        "# C2_1 freeze + audit-code notes (no KRAKEN rerun)",
        "",
        f"UTC: {NOW}",
        "",
        "## 冻结",
        "",
        "- `C2_1_PARSER_VALIDATED`：`.mod` 读取链可靠（用户独立二进制复核一致）",
        "- zgrid 150:2:250 为 KRAKEN 真实深度节点，可合法计算 φ_m(z)",
        "- self-built FEM 保持 **NOT_ADMISSIBLE**",
        "- 283 Hz：`M_mod=902`, `M_prt=903`，保留 `MODE_COUNT_DISCREPANCY`；后续以 `.mod` 结构为准",
        "",
        "## 审计代码漏洞（本轮不改数据、不重跑 KRAKEN；仅记录）",
        "",
        "1. 非 235 Hz `depth_ok` 含 `(... or True)`，使该项退化",
        "2. 跨 z_s 一致性只比了前两个运行，未覆盖 z_s=220（用户三组两两独立复核 corr=1、Δk=0）",
        "3. zgrid case 的 `k_checksum_ok` / `phi_ok` 有硬编码 True 成分",
        "",
        "原则：**DECISION 中的 True 本身不是证据，必须追溯 True 的产生方式。**",
        "",
        "本轮按 R3-C2.2A 范围**不**改 parser gate 代码；留给后续独立审计补丁。",
        "",
    ]
    (OUT / "C2_1_FREEZE_AND_AUDIT_NOTES.md").write_text("\n".join(lines), encoding="utf-8")


def write_report_and_decision() -> None:
    decision = "C2_2A_PRIMARY_FORMULA_NOT_RECOVERED"
    why = (
        "Yang 2015 Eq.(1)–Eq.(9) 全文未取得（closed access；官方 PDF/落地页 403；"
        "OpenAlex/Unpaywall 无仓库全文；CORE/仅元数据）。"
        "摘要仅提供方法框架。Erratum 提供 g(k_m,z_r)=b_m φ_m(z_r) exp(i k_m δ) 的修正形式及 δ 定义，"
        "并指出影响 Eq.(1)(5)(6)(9)，但不恢复主文全部原式。"
        "b_m 分解、深度 ambiguity Eq.(6)/(9)、SA 求和核、source-level 处理均为 NOT_RECOVERED。"
        "未重构任何 estimator。"
    )

    dec = {
        "stage": "R3-C2.2A",
        "rc3c2_2a_decision": decision,
        "why": why,
        "three_state": {
            "kraken_data_baseline": "C2_1_PARSER_VALIDATED",
            "yang_formula_recovery": "YANG_FORMULA_RECOVERY_PARTIAL",
            "yang_route_in_E_STD": "YANG_ROUTE_UNDECIDED",
        },
        "allowed_terminal_only": [
            "C2_2A_FORMULA_RECOVERED",
            "C2_2A_PRIMARY_FORMULA_NOT_RECOVERED",
        ],
        "not_done": [
            "new Yang simulation",
            "APERTURE_LIMITED / DEPTH_AMBIGUOUS verdicts",
            "close Yang candidate",
            "delta scan",
            "Doppler error sim",
            "Monte Carlo",
            "AR methods",
            "P5",
            "KRAKEN new scenarios",
            "z_hat/FWHM/PSL",
            "self-built estimator from abstract/secondary",
        ],
        "equation_recovery_summary": {
            "Eq.(1)": "NOT_RECOVERED",
            "Eq.(2)": "NOT_RECOVERED",
            "Eq.(3)": "NOT_RECOVERED",
            "Eq.(4)": "NOT_RECOVERED",
            "Eq.(5)": "ERRATUM_CORRECTED_FORM_ONLY",
            "Eq.(6)": "NOT_RECOVERED",
            "Eq.(7)": "NOT_RECOVERED",
            "Eq.(8)": "NOT_RECOVERED",
            "Eq.(9)": "NOT_RECOVERED",
        },
        "key_locks": {
            "A_eq1_sa_sum_steering": "NOT_RECOVERED",
            "B_eq5_g_form": "ERRATUM_FORM_ONLY",
            "C_eq6_depth_ambiguity": "NOT_RECOVERED",
            "D_eq9_relation": "NOT_RECOVERED",
            "b_m_factorization": "NOT_RECOVERED",
            "observable_object": "NOT_RECOVERED",
            "source_level_nuisance": "NOT_RECOVERED",
            "delta_search_strategy": "NOT_RECOVERED",
        },
        "created_utc": NOW,
        "retrieval": RETRIEVAL,
    }
    (OUT / "R3_C2_2A_DECISION.json").write_text(
        json.dumps(dec, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    rp = f"""# R3-C2.2A 报告：Yang 2015 原方法公式证据恢复

UTC：{NOW}

## 0. 冻结当前基础

- KRAKEN 数据基线：**VALIDATED**（`C2_1_PARSER_VALIDATED`）
- Yang 原方法公式恢复：**PARTIAL**
- Yang 在 E-STD 是否可用：**UNDECIDED**
- 283：`M_mod=902` / `M_prt=903` 保留 `MODE_COUNT_DISCREPANCY`
- self-built FEM：**NOT_ADMISSIBLE**

详见 `C2_1_FREEZE_AND_AUDIT_NOTES.md`。

## 1–2. 原文获取与逐公式证据

主文：{PRIMARY}

Erratum：{ERRATUM}

检索结果：**Eq.(1)–Eq.(9) 全文 NOT_OBTAINED**。证据表见 `YANG2015_EQUATION_TABLE.csv`。

仅 Erratum 给出的 g 修正形式可写入表；其余一律 **NOT_RECOVERED**。

## 3. 四个关键位置

| 锁 | 问题 | 状态 |
| --- | --- | --- |
| A | Eq.(1) SA 积分/求和与 steering phase | **NOT_RECOVERED** |
| B | Eq.(5) g=b_m φ_m(z_r) 及 Erratum 相位 | **仅 Erratum 形式** |
| C | Eq.(6) depth ambiguity / estimator | **NOT_RECOVERED** |
| D | Eq.(9) 与 Eq.(6) 关系及阵列条件 | **NOT_RECOVERED** |

完整前**禁止**写 Yang estimator。

## 4–7. b_m / observable / source-level / δ

见 `YANG2015_SYMBOL_TABLE.csv` 与 `YANG2015_OBSERVABLE_CHAIN.md`。

结论：b_m 分解、比较对象、源级处理、δ 搜索策略均为 **NOT_RECOVERED**。
**禁止**把 `A_m=φ_m(z_s)φ_m(z_r)` 或 `sqrt(member_energy)` 自制 score 写成 Yang 公式。

## 8. Doppler / PLL

仅摘要级用途已记录；公式 **NOT_RECOVERED**。本轮不仿真误差。

## 9. 次级文献

Liang et al. 2018 等只可交叉核验描述与 SA 孔径限制；本轮 PDF 未取得，**未**用其次级公式替代 Yang 原式。

## 10. 输出清单

- `YANG2015_FULL_METHOD_ANCHOR.md`
- `YANG2015_EQUATION_TABLE.csv`
- `YANG2015_SYMBOL_TABLE.csv`
- `YANG2015_ERRATUM_MAPPING.md`
- `YANG2015_OBSERVABLE_CHAIN.md`
- `C2_1_FREEZE_AND_AUDIT_NOTES.md`
- `R3_C2_2A_REPORT.md`
- `R3_C2_2A_DECISION.json`
- `R3_C2_2A_GPT_SYNC.md`

## 11. 判定

### `{decision}`

{why}

按用户指令：**宁可停在此状态，也不按摘要/次级文章补 estimator。**

下一阶段（需另行授权）：取得 2015 主文 PDF 后做 **忠实实现 → 先复现论文条件 → 再迁移 E-STD**。

## 12. 停止声明

- 不运行 KRAKEN 新场景
- 不生成 z_hat / FWHM / PSL
- 不做 δ 扫描 / Doppler 误差 / MC
- 不进入 AR
- 不进 P5
- 不评价 APERTURE_LIMITED / DEPTH_AMBIGUOUS
- 不关闭 Yang 候选
"""
    (OUT / "R3_C2_2A_REPORT.md").write_text(rp, encoding="utf-8")

    sync = f"""# R3-C2.2A GPT SYNC

UTC：{NOW}

## 三态

- KRAKEN 数据基线：VALIDATED
- Yang 原方法公式恢复：PARTIAL
- Yang 在 E-STD 是否可用：UNDECIDED

## 判定

**{decision}**

{why}

## 关键

- Eq.(1)-(9) 主文全文未取得
- Eq.(5) 仅有 Erratum 修正形式 g=b_m φ_m(z_r) exp(i k_m δ)
- b_m / observable / source-level / δ 搜索策略 NOT_RECOVERED
- 未写任何 Yang estimator
- 完成即停止

## 下一步（未执行）

取得主文全文后再做忠实实现；在此之前不得性能判定、不得关闭 Yang。
"""
    (OUT / "R3_C2_2A_GPT_SYNC.md").write_text(sync, encoding="utf-8")
    print("DECISION", decision)


def main() -> None:
    write_equation_table()
    write_symbol_table()
    write_anchor()
    write_erratum_mapping()
    write_observable_chain()
    write_freeze()
    write_report_and_decision()


if __name__ == "__main__":
    main()
