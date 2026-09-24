# SPATIAL_SAMPLING_ALIAS_AUDIT

UTC: 2026-09-24T09:51:45.113589+00:00

## 推导（只做解析，不评 MMAR 成败）

运动距离样本：r_i = r_0 + i Δr

模态相位因子（Eq.2 空间部分）：y_i ∝ exp[−j k_m r_i]

⇒ 相邻样本：**y_i ∝ exp[−j i k_m Δr]**

AR 归一化角频率与物理水平波数：

$$
\omega_s \equiv \pm k_m\,\Delta r \pmod{2\pi}
$$

存在 **符号** 与 **分支** 选择；须由可控复指数单元测试 + 论文谱横轴（Fig.3/4 约 **k=1.32–1.50 rad/m**）共同确定。

## 历史 Yang/Ref.7 采样 Δr≈2.5 m 解析检查

| k (rad/m) | kΔr | (kΔr)/π | 空间折叠？ |
| --- | --- | --- | --- |
| 1.32 | 3.300 | 1.050 | YES |
| 1.50 | 3.750 | 1.194 | YES |

**k_max Δr = 3.750 > π** ⇒ 按空间 Nyquist 会发生 **波数折叠/混叠**（不能用 Δr=2.5 m 无折叠覆盖 1.5 rad/m）。

## 无折叠 control

要求 k_max Δr < π ⇒ Δr < π/k_max ≈ **2.0944 m**

预注册 control：**Δr = 0.9 π / k_max ≈ 1.8850 m**（表中 NO_FOLD_CONTROL）。

## 边界

本审计 **不** 评价 MMAR 成败；只冻结实现解释与采样约束。
`PAPER_RANGE_SAMPLE_INTERVAL_NOT_EXPLICIT` 仍成立；Δt=1 s 仅 `REF7_CONSISTENT_SAMPLING_ASSUMPTION`。
