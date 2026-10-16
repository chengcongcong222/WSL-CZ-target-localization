# 最终交付清单及GPT整体验收标准

## 文件树（在项目仓库 `results/ENV_CAP_FULL_CAMPAIGN/`）

```
E0_BASELINE_AUDIT.md
E0_FREEZE.json
E1_ENV_SOURCE_INDEX.csv
E1_CZ_ENVIRONMENTS.csv
E1_CZ_GATE.json
E1_CZ_PHYSICS_EVIDENCE.md
E1_CZ_RANGE_SHAPE_DIAGNOSTICS.csv
E1_CZ_NEGATIVE_EXAMPLES.csv
E2_PAPER_METHOD_IDENTITY.md
E2_YANG_ERRATUM_FORMULA_AUDIT.md
E2_YANG_MMAR_COMPARISON.csv
E2_MODE_DENSITY_AND_APERTURE.csv
E2_METHOD_GATE.json
E3_PERFORMANCE_BY_CASE.csv
E3_OBSERVABLE_EXTRACTION.csv
E3_NUISANCE_PROFILED_INFORMATION.csv
E3_ENV_SENSITIVITY.csv
E3_GATE.json
E4_RESOURCE_INVERSION.csv
E4_NEEDS_FROM_CLIENT.csv
E4_PARETO_RESOURCE_FEASIBILITY.csv
E4_GATE.json
ENV_CAP_EVIDENCE_LEDGER.csv
ENV_CAP_FINAL_DECISION.json
ENV_CAP_FINAL_SCIENTIFIC_REPORT.md
GPT_SYNC_ENV_CAP_FINAL.md
APPLICATION_USABLE_CLAIMS.md
APPLICATION_FORBIDDEN_CLAIMS.md
REVIEW_ENTRYPOINT.md
RUN_STATE.json
SESSION_HANDOFF.md
INDEPENDENT_AUDIT_SUMMARY.json
FILE_HASHES_SHA256.csv
STOP.json
figures/ (CZ、模式密度、敏感性、资源Pareto、证据层级)
scripts/ (全部新增科学源与独立冷复核脚本)
raw_local_archive_index.csv (D盘大文件哈希索引)
```

没有计算的阶段仍提供同名表头/`NOT_EVALUATED`状态和原因，不能填伪数或无意义的空图。失败、拒判、计算超限的原始数据、trace和决策同等重要。

## 终稿必须回答的10个问题

1. 什么环境具有**实际数值证据的第一CZ**？首焦区的距离区间、源收深/海深与SSP条件是什么？是否有数值/模型不稳定？
2. 哪些环境只是浅海/非CZ论文控制、无CZ或来源不完整？为什么没有被纳入主比较？
3. Yang论文和勘误各哪些公式、源频/径向运动、delta前提已经复现？哪些受控oracle没有在现实提取层成立？
4. CZ名义环境中理想深度信息与可提取模态身份之间差多少？是模态密度、孔径、噪声、f0还是SSP失配限制？
5. Liang/MMAR带来什么新增（若有）？其论文复现完整性怎样？AR峰值与Hankel幅度是否严格区分？
6. 不同海洋环境/SSP会改变哪些参数信息弱方向、CZ位置、深度候选和不辨识结构？哪些只是物理上限？
7. 报告每项条件精度及不利/失败分母，特别是左右舷错误分支、源频未知、深度区间覆盖和连续目标状态的困难案例。
8. 预先设定的若干示例深度、距离、速度等级对应哪些必要条件、**已测条件可行**配置和未决资源？需要哪些SSP/海底/导航/时频/硬件条件？
9. 相对单HLA、新增垂向阵/非共线横向阵/辅助节点各增加了什么**独立观测**，并付出何种成本？不要拿T3-W1有符号方位做左右舷已解决证据。
10. 已知负结果限制到哪里？继续研究还是到此停止？哪些图和结论可以安全进入项目申请、哪些绝对不能写？

## 核心科研图最低集合

`F01_ENVIRONMENT_SSP_AND_CZ`：SSP、射线折返、TL首CZ位置/宽及环境标签；
`F02_CZ_GATE_AND_NEGATIVES`：数值稳定与未准入环境，明确未评估区域；
`F03_MODAL_SPACING_APERTURE`：模式波数间隔及有效径向孔径Rayleigh诊断；
`F04_YANG_MMAR_EVIDENCE_LAYERS`：理想信息、模态提取和性能递进；
`F05_SENSITIVITY_BOUNDARIES`：正负实验的环境/孔径/SNR/delta敏感性，有数据才画；
`F06_RESOURCE_PRECISION_PARETO`：不同参数等级与资源成本的条件可行域/未知域（若只有必要界则画“仅必要界”，不杜撰达标曲线）；
`F07_APPLICATION_ROUTE`：RC0首CZ场景准入→RC1可辨识→RC2提取与多分支→RC3精度和资源反推。

图号可以不变，单张可含子图，标题/图例必须标证据层级/单位/参数口径。正式生成科研图PDF或申请DOCX时，Codex应按相关格式技能或专业标准做正确字体/边距检查（本任务包不直接生成那些文件）。

## 独立审计

独立审计程序不能直接导入主分析的研究判定函数；必须从记录的原始数据/给定数学公式重算：环境数据身份、CZ ray/TL关键数值、方法基准、主要信息投影与QR/SVD一致性、分位数和拒判、资源Pareto非支配性、逐图数值来源、输入与结果哈希。留出数据与开发状态分开；无新数据的结果不允许宣称外部重复验证。

审计结果记录`check, error, tolerance, PASS, evidence_path, scope`。某阶段没有科学准入时，只读源身份审计PASS不等于科学数值认证PASS。

## 研究负责人接受条件

- **不是必须性能PASS**。失败也可收口，但必须有完整的原始证据/复现和边界说明。
- **不准漏写硬条件**：第一CZ、向量/信号获取、未知源和环境处理、左右舷、协同硬件成本、数值模型真实性。
- **不准一图改变主观结论**：正例与负例、理想与非理想、名义与敏感性并列。
- 物理CZ未认证时，不允许交付“深海首CZ经验深度精度”；仍可交付有价值的资料缺口和理论资源需求。
- 结题须有完整Git SHA、分支名、与远端一致、`STOP=true`、最终报告摘要、源数据索引、预算消费及下一步建议。R4历史进度保持0%，不要人为加信用。
