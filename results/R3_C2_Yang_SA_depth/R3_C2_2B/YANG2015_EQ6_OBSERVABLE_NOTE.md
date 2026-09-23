# YANG2015_EQ6_OBSERVABLE_NOTE

UTC: 2026-09-23T11:30:59.914058+00:00

## Eq.(6) 印刷形式

$$
D(z)=\left|\sum_{m=1}^{M} \phi_m(z)\,\frac{g(k_m,z_r)}{\bar{\phi}_m(z_r)}\right|^2
$$

正则化逆（**视觉锁定**）：

$$
\bar{\phi}_m^{-1}(z_r)=\frac{\phi_m(z_r)}{\phi_m^2(z_r)+\Delta^2}
$$

故

$$
\frac{g(k_m,z_r)}{\bar{\phi}_m(z_r)} = g(k_m,z_r)\cdot\bar{\phi}_m^{-1}(z_r)
= g(k_m,z_r)\,\frac{\phi_m(z_r)}{\phi_m^2(z_r)+\Delta^2}
$$

## observable 到底是什么

| 问题 | 锁定答案 |
| --- | --- |
| 输入是复值还是幅值？ | **复值谱峰** $g(k_m,z_r)$ 进入求和后再取模平方（Bartlett） |
| 是否只用 spectrum amplitudes？ | 附录原文：*Eq.(6) uses only the spectrum amplitudes, not the corresponding mode wavenumbers* |
| $k_m$ 数值是否进深度 score？ | **NO**。$k_m$ 只用于 **mode identification**（峰<->模态编号） |
| 印刷记号 $g(k_m,z_r)$ | 表示在模态波数峰处的谱值；实现取该峰处的复谱 |

**禁止** `sqrt(member_energy)` 或任何自造 score。

## 模态编号错误的影响

论文 Sec.III：

- 编号 $n$ 与 $m$ 错配 -> 深度误差量级 **$H/n - H/m$**（$H$ 为水深）
- 实测：用模拟谱 + 相邻模态波数差及其增长率辅助；作者承认 **educated guess**
- 可通过微调 residual Doppler 平移谱横轴使数据谱与模型谱对齐

## 正则化 Delta

| 项 | 值 |
| --- | --- |
| 符号 | $\Delta$ |
| 精确定义 | $\bar{\phi}_m^{-1}=\phi_m/(\phi_m^2+\Delta^2)$ |
| 论文推荐量级 | *on the order of one tenth of the maximum value of the mode depth functions* |
| 是否唯一数值 | **否** -> 状态 **`EMPIRICAL_NOT_UNIQUE`** |
| 作用 | 防止 $\phi_m(z_r)$ 零交叉处逆模态函数爆炸 |

**不得**把 0.1 写成论文固定参数。

## 深度分布归一化（Sec.II 原文）

normalized depth distribution = Eq.(6) 的 $D(z)$ **除以搜索深度网格上 $D(z)$ 的总和**。

这是论文展示用归一化分布，**不是**另造 correlation score。
