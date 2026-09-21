# P3 RC3 — GPT 同步稿

- G3 判定: **PROPAGATION_INCREMENT_CONFIRMED**
- 原因: CZ increased wrong-candidate rejection on S0 hard pairs by delta=0.967 (RC2=0.000 -> RC2+CZ=0.967). Doppler control branch is weak on constant-v_rad collinear pairs (absorbed by unknown source frequency bias).
- 入选困难候选: 30 对（A/B/C）

## P2 口径（已改，未重跑）

- 秩亏 CRLB(r) = **UNDEFINED/UNBOUNDED**（禁止 0.000 对外）
- 小机动：改变局部几何，**不足以形成实用距离约束**
- RC2 角色：不是远距离测距，而是约束方位/航向并给出 r–v 困难候选

## P3 主结果（分层修正后）

- 无稳定谱线时：RC2错误候选排除率约 0%；加入会聚区传播后约 97%；Doppler支路对恒定径向速度差不可用（被未知源频偏置吸收）。
- 典型机械谱线时：RC2约 0%；CZ约 23%；bearing+Doppler约 0%；联合约 20%。（线谱频点少，CZ轮廓自由度低于S0连续谱。）
- 理想稳定多谱线时：RC2约 0%；CZ约 97%；bearing+Doppler约 0%；联合约 97%。（S2为上界工况，非真实UUV谱。）
- CZ 增量 Δ(S0) = 0.967（RC2 0.000 → CZ 0.967）

## 机制摘要

mechanism source  mean_reject_RC2  mean_reject_RC2CZ  mean_reject_RC2Dop  mean_reject_combined
        A     S0              0.0                0.9                 NaN                   0.9
        A     S1              0.0                0.2                 0.0                   0.2
        A     S2              0.0                0.9                 0.0                   0.9
        B     S0              0.0                1.0                 NaN                   1.0
        B     S1              0.0                0.2                 0.0                   0.2
        B     S2              0.0                1.0                 0.0                   1.0
        C     S0              0.0                1.0                 NaN                   1.0
        C     S1              0.0                0.3                 0.0                   0.2
        C     S2              0.0                1.0                 0.0                   1.0

## 下一轮（GPT）

> 判定 G3 结论是否成立，冻结研究内容三的阶段性表述；**不**自动进入 P4。
