# R3-C2.2A 报告：Yang 2015 原方法公式证据恢复

UTC：2026-09-23T08:46:33.049990+00:00

## 0. 冻结当前基础

- KRAKEN 数据基线：**VALIDATED**（`C2_1_PARSER_VALIDATED`）
- Yang 原方法公式恢复：**PARTIAL**
- Yang 在 E-STD 是否可用：**UNDECIDED**
- 283：`M_mod=902` / `M_prt=903` 保留 `MODE_COUNT_DISCREPANCY`
- self-built FEM：**NOT_ADMISSIBLE**

详见 `C2_1_FREEZE_AND_AUDIT_NOTES.md`。

## 1–2. 原文获取与逐公式证据

主文：T. C. Yang, Source depth estimation based on synthetic aperture beamfoming for a moving source, JASA 138(3), 1678-1686 (2015). DOI 10.1121/1.4929748. PMID 26428805.

Erratum：Erratum: JASA 144(6), 3075 (2018). DOI 10.1121/1.5081712. PMID 30599696.

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

### `C2_2A_PRIMARY_FORMULA_NOT_RECOVERED`

Yang 2015 Eq.(1)–Eq.(9) 全文未取得（closed access；官方 PDF/落地页 403；OpenAlex/Unpaywall 无仓库全文；CORE/仅元数据）。摘要仅提供方法框架。Erratum 提供 g(k_m,z_r)=b_m φ_m(z_r) exp(i k_m δ) 的修正形式及 δ 定义，并指出影响 Eq.(1)(5)(6)(9)，但不恢复主文全部原式。b_m 分解、深度 ambiguity Eq.(6)/(9)、SA 求和核、source-level 处理均为 NOT_RECOVERED。未重构任何 estimator。

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
