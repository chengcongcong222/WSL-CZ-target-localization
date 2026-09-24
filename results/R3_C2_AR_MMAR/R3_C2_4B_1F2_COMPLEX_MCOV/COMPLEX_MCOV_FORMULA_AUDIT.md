# COMPLEX_MCOV_FORMULA_AUDIT

UTC: 2026-09-24T11:37:52.613113+00:00

`COMPLEX_BACKWARD_CONJUGATION_REQUIRED`

```text
Af[n,:] = [y[n-1], ..., y[n-p]]
bf[n]   = -y[n]
Ab[n,:] = [conj(y[n-p+1]), ..., conj(y[n])]
bb[n]   = -conj(y[n-p])
```

对应 e_b*[n] = y*[n-p] + Σ a[k] y*[n-p+k]（backward 用 a* 后整体共轭）。

p=1 恒等式：y=e^{-jknΔr} ⇒ a1=−e^{-jkΔr}；ω_peak=wrap(−kΔr)。

Gate：**P1_COMPLEX_EXP_IDENTITY_VALIDATED**
