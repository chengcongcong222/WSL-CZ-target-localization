# DZ_INPUT_CONTRACT

UTC: 2026-09-25T03:59:05.755902+00:00

未来 D(z) 输入（本轮只锁契约，不计算 b、D(z)）：

| 字段 | 来源 |
| --- | --- |
| matched theoretical mode IDs | Eq.(24) / nearest unique |
| k̂_m | AR p=7 printed branch |
| g(k̂_m,z_r) | **Eq.(7) Hankel**（复数） |
| φ_m(z_r), φ_m(z) | KRAKEN |
| δ | **0 = ORACLE_OFFSET_ALIGNMENT** |
| Δ regularizer | Liang ~½ max\|φ\| |

**禁止** AR peak height → amplitude。
**禁止**本轮计算 b=(Φ+U)^-1g 或 D(z)。

状态：机制验证，**非** Liang Fig.5/6 strict reproduction。
