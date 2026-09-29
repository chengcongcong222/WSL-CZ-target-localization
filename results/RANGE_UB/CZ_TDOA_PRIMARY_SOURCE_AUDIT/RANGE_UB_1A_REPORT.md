# RANGE_UB_1A_REPORT — CZ Multipath TDOA Primary-Source Lock & Transferability Audit

UTC: 2026-09-29

基线：`ad91ce0bd11bf3a511bdf99a84fc38f1a53b5329`

## 判定

### `CZ_TDOA_REQUIRES_ADDITIONAL_APERTURE_OR_SIGNAL`

论文方法与 E-STD 存在**三项 BLOCKING 不兼容**：

1. **REQUIRES_BROADBAND_OR_IMPULSIVE_SIGNAL** — 论文用爆炸声源（脉冲）分离多途时延峰；E-STD 是被动窄带线谱。
2. **REQUIRES_VERTICAL_APERTURE** — 论文核心是"时延差曲线随接收深度变化"；E-STD 是近水平 HLA。
3. **TDOA_OBSERVABLE_NOT_AVAILABLE** — 窄带连续信号无到达时刻，无法提取 TDOA。

---

## 原始文献锁定

| 项 | 值 |
|---|---|
| 作者 | 徐嘉璘, 郭良浩（中科院声学所） |
| 题目 | 深海会聚区声场多途时延差分析及声源距离估计 |
| 期刊 | 应用声学, 2024, 43(2): 237-251 |
| DOI | 10.11684/j.issn.1000-310X.2024.02.001 |
| **PDF 状态** | **`CZ_TDOA_PRIMARY_SOURCE_NOT_RECOVERED`**（反爬拦截） |
| 摘要状态 | 官方摘要已获取并锁定（SHA256: 0bad8fb5…） |

---

## 已确认声明（摘要原文）

| 声明 | 状态 |
|---|---|
| 42–52 km 收发距离 | **CONFIRMED** |
| 0.6%–6.1% 估计误差 | **CONFIRMED** |
| 虚源理论 | **CONFIRMED** |
| 爆炸声源 | **CONFIRMED** |
| 大接收深度 | **CONFIRMED** |
| 时延差随接收深度变化 | **CONFIRMED** |

这些数字**不能**作为 E-STD 性能依据——条件不匹配。

---

## 迁移审计（BLOCKING）

| 论文条件 | E-STD | 后果 |
|---|---|---|
| 爆炸声源（脉冲宽带） | 被动窄带 4 线 | `REQUIRES_BROADBAND_OR_IMPULSIVE_SIGNAL` |
| 大接收深度（垂直采样） | HLA z≈200 m | `REQUIRES_VERTICAL_APERTURE` |
| 多途时延峰互差 | 窄带无到达时刻 | `TDOA_OBSERVABLE_NOT_AVAILABLE` |

### WARNING

- 源深：近海面 vs UUV 200 m
- 距离：42–52 km 验证 vs 50–60 km 使用（50–52 重叠）
- 几何：可能双基（爆炸声源）vs 单站被动

---

## 对主线的意义

1. **CZ-TDOA 不能作为当前 HLA 的增强模块**。
2. **可作为"理想声学距离信息上界"**（ORACLE_UPPER_BOUND）：若有多途时延 + 垂直孔径，距离可达 0.6%–6.1% 精度。
3. **与主线互补**：主线（RC2+RC3+转向）处理 [r,θ,v,ψ]；TDOA 上界说明距离信息的物理极限在哪里。

---

## 下一步

- 若用户提供完整 PDF：进入 `RANGE-UB-1B` 做公式恢复 + ORACLE 仿真。
- 若不提供：TDOA 只作 `ORACLE_UPPER_BOUND_ONLY` 引用，不进入算法链。

## 未做

TDOA 仿真、人造时延、oracle 定位、RC2/RC3 集成、噪声扫描、平台机动、P5。
