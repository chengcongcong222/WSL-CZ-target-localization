# 283 Hz mode count audit

UTC: 2026-09-23T08:12:42.907979+00:00

- `.mod` M=902（record 5） vs `.prt` M=903 → **MODE_COUNT_DISCREPANCY**
- **不**用 prt 覆盖 mod；**不**静默删除 283
- 下游结构以 `.mod` 的 M=902 为准，差异待 AT/KRAKEN 源码或官方 reader 解释
