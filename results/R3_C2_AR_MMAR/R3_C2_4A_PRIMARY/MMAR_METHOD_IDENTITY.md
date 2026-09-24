# MMAR_METHOD_IDENTITY

UTC: 2026-09-24T09:03:08.658004+00:00

方法链（主文确认）：

**HLA beamforming (Eq1–4) → AR 估 k̂ (Eq17–19) → generalized Hankel 在 k̂ 处取 g(k̂) (Eq7) → ordered-subset mode matching (Eq24–27) → matched-mode D(z) (Eq20–23)**

| 模块 | 来源 | 角色 |
| --- | --- | --- |
| HLA beamforming | Liang2018 | 相对 Yang 单听器 SAB **新增** |
| moving-range sequence | 同 Yang SAB 类 | 运动合成距离 |
| AR wavenumber | Liang2018 核心 | **只定 k̂** |
| Hankel amplitudes | Yang/Ref 传统 | 在 k̂ 处取幅度 |
| ordered mode matching | Liang2018 Eq24–27 | **不是 nearest / Hungarian** |
| matched-mode depth | Liang2018 Eq20–23 | 矩阵形式，非 Yang Eq6 原文 |

**不得称为 Yang2015 的 AR 版。**
AR peak height **不得**作 modal amplitude。
