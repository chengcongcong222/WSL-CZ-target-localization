# YANG2015_ERRATUM_MAPPING

UTC: 2026-09-23T08:46:33.049990+00:00

Erratum 引用主文 **Eq.(1), Eq.(5), Eq.(6), Eq.(9)**。

| 位置 | Erratum/协议要求回答 | 本轮证据 | 状态 |
| --- | --- | --- | --- |
| Eq.(1) | SA 如何积分/求和；steering phase | 仅摘要：range-dependent phase relative to starting point | **NOT_RECOVERED** |
| Eq.(5) | 为何 k_r=k_m 时 g=b_m φ_m(z_r)；Erratum 为何加 exp(i k_m δ) | Erratum 形式已记；原推导与原印刷式 | **FORM_FROM_ERRATUM; DERIVATION NOT_RECOVERED** |
| Eq.(6) | 具体 source-depth ambiguity / estimator | 仅知 δ 相位会改变该函数 | **NOT_RECOVERED** |
| Eq.(9) | 与 Eq.(6) 关系；接收阵/数据条件 | 仅摘要有 VLA 直接估深路径 | **NOT_RECOVERED** |

## δ 关系

- 定义：δ = r1 − r1_data = r2 − r2_data（初始距离偏移，两时刻同步差）
- 进入位置：g(k_m, z_r) 的相位因子 exp(i k_m δ)（Erratum）
- δ = 0：原文 simulation / data processing 隐含同步
- δ 未知时是 (z, δ) 联合搜索还是先 δ 后 z：**NOT_RECOVERED**（不得推测）

## 证据来源分级

- Erratum 形式：来自先前会话 `YANG2015_ERRATUM_NOTE.md`（标 RECOVERED）+ 本任务书/用户复述；本轮未能重下 PDF 复核字节级原文
- 主文 Eq 印刷式：**全部 NOT_RECOVERED**
