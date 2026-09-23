# YANG2015_PAPER_REPRO_CONFIG

UTC: 2026-09-23T11:30:59.914058+00:00

## DIRECTLY_SPECIFIED_IN_YANG2015（论文模拟）

| 参数 | 值 | 出处 |
| --- | --- | --- |
| source depth z_s | 4 m（shallow）, 50 m（deep） | Sec.II |
| receiver depth z_r | 18 m, 70 m | Sec.II |
| frequency | 350 Hz | Sec.II |
| speed | 5 kn | Sec.II |
| motion | moving **away** from receiver | Sec.II |
| initial range r1 | 5010 m（模拟生成用；**接收端算法不把 r1 当已知输入**） | Sec.II |
| sampling | 1 s | Sec.II |
| synthetic range span | **约 4990 m** | Sec.II |
| VLA demo | 20 elements, dz=4 m, covering 2-80 m | Sec.IV B |
| SWellEx96 real data | 127 Hz zs=9 m; 130 Hz zs=54 m; z_r=bottom phone; ~35 min | Sec.III |

## REFERENCED_TO_REF7

Ref.7 = T. C. Yang, *Data-based matched-mode source localization for a moving source*, JASA **135**, 1218-1230 (2014).

| 项目 | 状态 |
| --- | --- |
| downward refractive SSP（具体剖面） | **NEEDS_REF7** |
| water depth | **NEEDS_REF7** |
| bottom properties | **NEEDS_REF7** |
| mode attenuation / model details | **NEEDS_REF7** |
| 其它传播配置 | **NEEDS_REF7** |
| Eq.(1) 在 Ref.7 中称 generalized Hankel transform | 已记录 |

**PAPER_REPRO_ENV = `PAPER_REPRO_NEEDS_REF7`**

## 模态峰识别（必须分层）

| 标签 | 定义 | 论文依据 |
| --- | --- | --- |
| **ORACLE_MODE_ID** | 用模型 true k_m 标注谱峰 | Fig.1 "+" true modal wavenumbers |
| **PRACTICAL_MODE_ID** | mode spacing + 相邻 spacing 增长率 + 模拟谱辅助；含 educated guess | Sec.III 实测流程 |

**禁止混为一谈。** 忠实复现论文仿真用 ORACLE_MODE_ID。

## 孔径条件（Sec.IV C）

- 分辨 n/m 阶模态：L > lambda_nm = 2 pi / |k_n-k_m|
- 示例：lambda_12 约 2.3 km（**论文环境示例**）
- 1 km 短孔径 -> 论文建议高分辨算法（Ref.12 AR 等）

**2.3 km 不得直接迁移为 E-STD 固定阈值**；E-STD 须用 KRAKEN k_m 重算。

## PAPER_REPRODUCTION_BASELINE

- **delta = 0**（Erratum 修正 carry-forward；无 Erratum PDF 本轮重核）
- 第一轮忠实复现 **不得**把 2.4 km E-STD 场景塞进论文条件
- 论文孔径 **约 4.99 km** 优先
