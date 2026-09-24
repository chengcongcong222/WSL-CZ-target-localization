# SNR_EQ28_AUDIT

UTC: 2026-09-24T16:02:27.087144+00:00

Liang Eq.(28): SNR = 10 log10(Ps/Pn)|_{r=r0}, r0=5010 m；**array gain 已含在 SNR 中**。

修正：**Ps = |B(r0)|²**（r=r0 处 beam output），禁止 whole-aperture mean。

噪声：Pn = Ps / 10^(SNR/10)；生成后反算 `SNR_AT_R0_DB` 验证。
不按谱结果改噪声。
