# CZ envelope source/gain nuisance 草案

本轮仅符号模型与未来审计规则。真实 UUV 漂移幅度、时尺度和噪声/接收增益没有冻结；不以旧 model stress 的 dB 数值冒充真实源统计。

## N0 — STATIC_PER_LINE_OFFSET

`a_f(t)+g_f(t)=A_fw`，每 line、每旧窗口内为未知常数，跨 line 值任意。C_w 消除 A_fw；不需要 common source level、known spectrum 或跨频 calibration。E1/E2 为线性 projection，E3 的逐 line normalized intensity 也对常数严格不变。允许窗口之间常数不同，不在 M0 中制造一个不存在的跨窗 level step。

N0 是第一轮 matched-model mechanism 的基本 control，不是“真实 UUV 已足够稳定”。A1 存档的 acoustic feature 没有证明实测声噪声分布；bearing 含噪不等于 acoustic level 已有实测噪声。

## N1 — SLOW_BOUNDED_DRIFT

合并源级与不能独立区分的 beam gain 漂移，避免对 a_f 和 g_f 作任意重复拟合。建议最小模型 `d_fw(t)=B_fw b_w1(t)`，`|B_fw|<=D_fw`；D 是预运行由负责人依据合同/诊断目标选择的参数，不给真实 UUV 数值。本轮标为 `DRIFT_BOUND_NOT_FROZEN`。可讨论 common drift 与 line-specific drift 两个预注册结构，但不能看结果后改成有利的结构。

在 E1 中 linear drift 恰进入一阶系数；若该自由度无界，一阶传播 trend 被吸收，二阶系数仍可作为有限假设方向。bounded drift 保留多少信息必须未来验证，不能把 D 随 range 候选调小。若把独立 quadratic drift 也设成无界，E1 的两个系数可全部吸收：这是结构限制，不是优化器失败。

E2 的 drift 通过同一 block means 映射；E3 在 level 域加 d 后再指数/归一化，非线性 nuisance 不能在最终 feature 上偷做线性加法。模型/观测端处理保持一致。未来审计必须记录 drift family、bounds、每候选拟合 nuisance、是否贴边与 residual 不确定性。

gain 的未知常数在 N0；时变 gain 必须有 metadata 支持或进入 N1。任意未知 gain variation 与 N2 一样危险。“慢漂”与“慢 CZ envelope”可能同尺度；不能先 detrend 再假定只删了源变化。

## N2 — UNCONSTRAINED_LIMIT

对于任意观测和候选 propagation，令 `a_f(t)+g_f(t)=L_f(t)-P_f(R(t;x),z_s,z_r,e)`，即可在允许的无限 nuisance 下解释它（噪声也可并入）。因此仅靠 level shape 没有独立 range 约束：**`NON_IDENTIFIABLE_NUISANCE_LIMIT`**。这是纯符号推论，不是新数值结果，不要求运行一个大规模 N2 主实验。

输出这项 limit 可以防止把任意 flexible source model 当 robustness；不能从 N2 的人为无限自由度推断真实 UUV 物理不可辨识。

## 其他 nuisance 与 future support 处理

初始水平状态保留 `(r,theta,v,psi)`；z 为继承 finite profile nuisance，匹配环境 E-STD、receiver depth/导航使用旧 proxy 元数据。首次 N0 matched 设计不发展 depth/SSP，也不把评估 truth z、truth motion 固定给 estimator。旧 profile label 与 effective stored modal depth 的 mapping 须归档；连续深度全域覆盖不在该 matched finite-profile claim 内。

若未来加环境/阵深误差，应事先定义集合与 score 中的 union/minimization，不可以从 truth 选择 matched template。本轮不生成误差样本或 source-drift observations。N1 和 acoustic noise 的生成/验证均要新授权、来源标签和 pre-run freeze。

## feature acceptance band 与禁止误用的旧 J

support compatibility 的 `tau_F` 与 outcome 的 range-width thresholds 是两类不同条件。推荐 `tau_F = feature-noise bound + validated numeric error allowance`（或由已知 covariance 建立预注册 region）；未知项目声噪声不能默认为真实零噪声。N0 的旧 acoustic proxy 可以明确声明 synthetic noiseless-level control，但非真实接收 validation。

线性 E1/E2 的 feature uncertainty 可由固定线性映射对 noise/gain bounds 或 covariance 传播，E3 需包含 nonlinear normalization 的传播与误差界。数值误差使用 neutral fixtures/self-consistency 依据，不能读取 truth recovery 后定容差。当前 `tau_F / acoustic_noise_model` **NOT_FROZEN**，不能由 `J<0.001 dB`、`J<=Jmin+0.5` 或新运行的 Jmin 自行替代；旧 A1 Gate 一律不改。

N1 兼容 score 在预先规定 D 内 profile，并保留全部兼容状态。所有 nuisance 估计只能来自 observation 与先验合同；truth 仅交 evaluator。pre-run 未确定 nuisance/score tolerance 时不得执行。
