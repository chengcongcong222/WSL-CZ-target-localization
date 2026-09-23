# YANG2015_FULL_METHOD_ANCHOR

UTC: 2026-09-23T08:46:33.049990+00:00

## 三态分离（冻结）

- KRAKEN 数据基线：**VALIDATED** (`C2_1_PARSER_VALIDATED`)
- Yang 原方法公式恢复：**PARTIAL** (`YANG_FORMULA_RECOVERY_PARTIAL`)
- Yang 在 E-STD 是否可用：**UNDECIDED** (`YANG_ROUTE_UNDECIDED`)

## 文献

- 主文：T. C. Yang, Source depth estimation based on synthetic aperture beamfoming for a moving source, JASA 138(3), 1678-1686 (2015). DOI 10.1121/1.4929748. PMID 26428805.
- Erratum：Erratum: JASA 144(6), 3075 (2018). DOI 10.1121/1.5081712. PMID 30599696.

## 检索状态（本轮）

- OpenAlex / Unpaywall：**closed access**，无机构库全文
- Crossref 给出官方 PDF 直链，AIP/Silverchair：**HTTP 403**
- EuropePMC：仅摘要，`hasPDF=N`
- CORE：仅元数据
- ResearchGate 282528816：**HTTP 403**
- 次级 Liang et al. 2018 PDF：**HTTP 403**（本轮不做次级公式替代）
- 本地用户目录：无相关 PDF
- Sci-Hub：**未使用**
- **Eq.(1)–Eq.(9) 主文全文：NOT_OBTAINED**

## 仅摘要级框架（禁止当作公式）

- single hydrophone; moving CW source
- synthetic aperture created by source motion
- signal at each range steered by range-dependent phase relative to starting point
- range increment (aperture) from Doppler shift estimated from data, knowing original signal frequency
- depth estimated from beam output, assuming mode depth functions from nominal sound speed and bottom profile
- PLL removes fast random phase; leaves deterministic range-dependent phase due to range change
- if vertical array covers depth span of interest, beam output can estimate depth without acoustic environment knowledge

## Erratum 级已锁定（非主文原式）

- δ = r1 − r1_data = r2 − r2_data
- g(k_m, z_r) = b_m φ_m(z_r) exp(i k_m δ)
- exp(i k_m δ) 不改变 wavenumber spectral density，**会**改变 source-depth ambiguity
- 原文模拟/数据处理对应 δ = 0 同步条件
- “不需要源距离”表述不严格成立

## 禁止事项

- 禁止用摘要/次级文献重构 Eq.(1)–(9)
- 禁止自造 Yang estimator / depth score
- 禁止把 A_m = φ_m(z_s)φ_m(z_r) 当作 Yang 正式公式
- 禁止输出 Yang 性能判定（APERTURE_LIMITED / DEPTH_AMBIGUOUS 等）
