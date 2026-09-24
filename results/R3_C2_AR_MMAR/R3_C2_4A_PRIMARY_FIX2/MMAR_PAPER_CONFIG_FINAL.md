# MMAR_PAPER_CONFIG_FINAL

UTC: 2026-09-24T09:59:41.308845+00:00
PRIMARY SHA256: `9364914d34f5408f80065a31cd7fbcb5195c21047c97b2eb118348191b861834`
唯一配置入口（4B 不得从 superseded V1/V2 继承参数）

| 项 | 值 | 状态 |
| --- | --- | --- |
| frequency | 350 Hz | PRINTED_PRIMARY |
| HLA N_elem | 2L+1=**11**（L=5） | PRINTED_PRIMARY |
| d | **λ/2**；参考声速 **主文未明确** | `PAPER_LAMBDA_REFERENCE_NOT_EXPLICIT` |
| HLA depth | 70 m | PRINTED_PRIMARY |
| z_s | 4 m / 50 m | PRINTED_PRIMARY |
| r0 | 5010 m（算法当未知） | PRINTED_PRIMARY |
| v | 2.5 m/s | PRINTED_PRIMARY |
| motion | 沿 beam 方向远离 | PRINTED_PRIMARY |
| look angle | θ̂=θ=θ | PRINTED_PRIMARY（ORACLE_LOOK_DIRECTION） |
| depth during aperture | fixed | PRINTED_PRIMARY |
| range span | **1990 m**（insufficient） / **4990 m**（sufficient） | PRINTED_PRIMARY |
| SNR | 20 / 5 / −5 dB | PRINTED_PRIMARY |
| field | KRAKEN | PRINTED_PRIMARY |
| bottom | referred to Yang2015 Ref.[7] | REFERENCED |
| MC | C0=500；correct \|ẑ−z\|≤5 m；90% CI Eq.30 | PRINTED_PRIMARY |
| δ | **0** = `ORACLE_OFFSET_ALIGNMENT` | 项目约定 |
| Δ regularizer | ~½ max\|φ\| | PRINTED_PRIMARY（Eq.23） |

## 非 Liang 主文明示

| 项 | 说明 |
| --- | --- |
| **Δr / Δt** | 主文本节**未给**时间/距离采样间隔 |
| 若用 Δt=1 s | 只能标 **`REF7_CONSISTENT_SAMPLING_ASSUMPTION`** |
| 空间无折叠 | Δr < π/k_max；见 `SPATIAL_SAMPLING_ALIAS_AUDIT.md` |
| AR p | 见 `AR_ORDER_PREREGISTRATION.md`（两分支） |
