# R3-C2.1-FIX

UTC：2026-09-23T07:54:09.501452+00:00

`.mod` 按 AT 布局解析：尾部完整 `k[M]`，M 与 `.prt` 一致，k 与 prt checksum 一致。

旧 `C2_1_YANG_APERTURE_LIMITED` → **SUPERSEDED_PENDING_PARSER_FIX**

## M / 峰组
{"201.0": 2.0, "235.0": 3.0, "283.0": 0.0, "338.0": 3.0}

## `C2_1_YANG_FORMULA_BLOCKED`

Parser OK; multiple peak groups appear (n_peaks={201.0: 2.0, 235.0: 3.0, 283.0: 0.0, 338.0: 3.0}) but YANG_FORMULA_RECOVERY_PARTIAL — cannot close Yang on self-made depth score alone.

下一步：recover full Yang depth estimator equations before ADMISSION/DEPTH_AMBIGUOUS