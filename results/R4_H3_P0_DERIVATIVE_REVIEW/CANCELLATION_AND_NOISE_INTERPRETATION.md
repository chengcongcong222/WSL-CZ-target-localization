# Cancellation and noise weighting interpretation

| Resource | R1/D1 effective information range | Absolute-floor / relative information range (same sigma) |
|---|---:|---:|
| B0 | 1.0282977 to 1.0709162 | 0.048363206 to 1.9467194 |
| B1 | 1.034704 to 1.071432 | 0.077275854 to 2.2055683 |
| B2 | 1.0281877 to 1.0723646 | 0.049533546 to 1.9406946 |

These noise models use the same nominal reference, but different covariance weights. Their information ratio describes frozen ideal noise assumptions; it is not measured SNR or receiver feasibility.

| Method | Max derivative modal cancellation | Max pressure-stencil cancellation |
|---|---:|---:|
| D0.5 | 1131.79 | 720.60427 |
| D1 | 784.95918 | 248.79883 |
| D2 | 388.324 | 56.62162 |
| R1 | 1025.6915 | 971.07385 |
| R2 | 1054.2957 | 484.81151 |

Lowest nominal amplitude: 3.6220821e-06, H12_M/B0/n80001, [time,element,frequency] index [2, 6, 12] (zero based). Per-record snapshot and frequency order are retained in the CSV. All seven source-depth pressures and five derivatives are saved in each STENCILS NPZ; source depths and receiver positions remain frozen.

Large cancellation factors explain why derivative serialization and finite differences deserve separate attention. They do not bound the true derivative error. R1 adds approximately2.82% to 7.24% effective information relative to D1; therefore the choice of stencil has a measurable effect even though the two high-order forms are internally close.

The observed near-h^2  pattern and smaller high-order discrepancy support a substantial truncation component. They do not isolate it from modal serialization, source-knot smoothness or remaining propagation errors. A rigorous error budget has not been established.
