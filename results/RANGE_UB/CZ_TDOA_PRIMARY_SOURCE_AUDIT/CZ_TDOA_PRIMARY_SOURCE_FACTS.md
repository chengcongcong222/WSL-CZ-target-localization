# CZ_TDOA_PRIMARY_SOURCE_FACTS

UTC: 2026-09-29

**状态：`CZ_TDOA_PRIMARY_SOURCE_NOT_RECOVERED`（全文 PDF 未获取）**

以下事实来源于**期刊官网官方摘要页**（aipub.cn，2026-09-29 访问），已用 SHA256 锁定。
标注 `[ABSTRACT]` 表示摘要明确支持；标注 `[NEEDS_PDF]` 表示需全文核验。

---

## 基本信息

| 项 | 值 | 来源 |
|---|---|---|
| 作者 | 徐嘉璘, 郭良浩 | [ABSTRACT] |
| 单位 | 中国科学院声学研究所 | [ABSTRACT] |
| 期刊 | 应用声学, 2024, 43(2): 237-251 | [ABSTRACT] |
| DOI | 10.11684/j.issn.1000-310X.2024.02.001 | [ABSTRACT] |
| 投稿/修订 | 2022-11-15 / 2024-02-26 | [ABSTRACT] |

---

## 方法核心

| 项 | 值 | 置信 |
|---|---|---|
| 理论基础 | 虚源理论推导多途到达时延结构 | [ABSTRACT] |
| 核心量 | 多途时延差（multipath time delay difference） | [ABSTRACT] |
| 关键曲线 | 时延差曲线随**接收深度**变化的近似表达式 | [ABSTRACT] |
| 测距原理 | 时延差曲线与声源水平距离的关系 | [ABSTRACT] |
| 目标声源 | **近海面**声源 | [ABSTRACT] |
| 信号类型 | **爆炸声源**（exploding/exploding sound source） | [ABSTRACT] |
| 定位类型 | 被动定位（关键词）；但使用爆炸声源=主动/合作式 | [ABSTRACT] |

## 实验验证

| 项 | 值 | 置信 |
|---|---|---|
| 海试 | 南海实验 | [ABSTRACT] |
| 收发距离范围 | **42 km ~ 52 km** | [ABSTRACT] |
| 距离估计误差 | **0.6% ~ 6.1%** | [ABSTRACT] |
| 信号 | 会聚区爆炸声源 | [ABSTRACT] |
| 比较基准 | 实际爆炸距离（known explosion range） | [ABSTRACT] |

## 关键词揭示的条件

| 关键词 | 含义 | 迁移影响 |
|---|---|---|
| 深海会聚区 | 与 E-STD 一致 | 兼容 |
| 多途时延结构 | 核心 observable | 需脉冲/宽带分离多途 |
| 被动定位 | 被动接收 | 但源是爆炸声源 |
| **大接收深度** | 接收深度范围大，可能需要垂直孔径 | **HLA 不兼容** |

---

## 待全文核验 [NEEDS_PDF]

以下必须从原文逐项恢复，当前摘要未提供：

- 声源深度（近海面具体多少 m？）
- 接收深度范围（"大接收深度"具体 0–3000 m？）
- 接收设备类型（单水听器 / VLA / 垂直阵？）
- 阵元数、孔径、深度跨度
- 海深、SSP
- 第一会聚区距离范围
- 使用哪个/哪些会聚区
- 信号带宽、脉冲宽度
- 多途数（2 途？3 途？）
- 各多途物理定义
- 时延差提取方式（互相关？匹配？）
- 距离匹配公式
- 是否同时估深
- 仿真条件细节
- 海试详细配置

---

## 声明核验

| 声明 | 状态 | 依据 |
|---|---|---|
| 42–52 km | **CONFIRMED_BY_PRIMARY_ABSTRACT** | 官方摘要原文 |
| 0.6%–6.1% | **CONFIRMED_BY_PRIMARY_ABSTRACT** | 官方摘要原文 |
| 基于虚源理论 | **CONFIRMED_BY_PRIMARY_ABSTRACT** | 官方摘要原文 |
| 需要接收深度变化曲线 | **CONFIRMED_BY_PRIMARY_ABSTRACT** | 摘要"时延差曲线随接收深度变化" |
| 爆炸声源 | **CONFIRMED_BY_PRIMARY_ABSTRACT** | 官方摘要原文 |
