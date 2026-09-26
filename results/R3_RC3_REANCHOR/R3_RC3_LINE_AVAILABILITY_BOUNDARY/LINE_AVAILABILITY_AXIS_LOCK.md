# LINE_AVAILABILITY_AXIS_LOCK

UTC: 2026-09-26T03:58:23.823471+00:00

`S1_BRIDGE_LINE_AVAILABILITY_ONLY` / `TRACKABLE_LINE_AVAILABILITY_AXIS`

- 环境：E0_NOMINAL（无 SSP mismatch / 无 IID TL error / 无频漂）
- F4 = {201,235,283,338} Hz；全部 **15** 个非空子集冻结运行
- 每频独立去均值 TL 形状；多频拼接
- 双侧 profile z∈150:5:250
- 主量：ΔJ = J_alt* − J_true*
- 子集 PASS：fraction_tested(ΔJ>0)≥0.8 且 median ΔJ>0（路线 Gate，非统计检验）
- 仍假定：所选线谱持续存在、中心频率稳定、幅度变化不污染 TL

**不是** 完整 S1；不含 168/204/232/279/320 Hz 迁移。
