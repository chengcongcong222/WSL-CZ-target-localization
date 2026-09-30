# 1E_B_INVALIDATION_AUDIT

UTC: 2026-09-30T07:49:13.619437+00:00

旧1E-B判定降级为 SUPERSEDED_INVALID_FREQUENCY_MODEL_MANIFEST。

原因:
1. 所有非标称漂移频率 model_exists=False
2. manifest 失败未停止执行
3. q!=0 时间片被静默置零
4. D>0 结果无效

保留旧输出作为失败审计。
