# Auxiliary measurement requirement confirmation

Gate2A maxima are SINGLE FACTOR ONLY. Joint requirements below come from Gate2B simultaneous error models; the columns cannot be interchanged. Client values remain UNKNOWN.

| Quantity | Single-factor tested max | Jointly validated requirement | Client value | Status |
|---|---|---|---|---|
| T5 Anchor A: baseline km | 5.0 | 5.0 at tested lambda=0.5 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T5 Anchor A: per-node random bearing sigma deg | 0.05 | 0.05 at tested lambda=0.5 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T5 Anchor A: common bias magnitude deg | 0.1 | 0.05 at tested lambda=0.5 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T5 Anchor A: differential HALF-bias magnitude deg | 0.05 | 0.025 at tested lambda=0.5 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T5 Anchor A: per-node per-axis position sigma m | 50 | 25.0 at tested lambda=0.5 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T5 Anchor A: known actual deployment abs(beta) deg | 20 | 10.0 at tested lambda=0.5 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T5 Anchor B: baseline km | 7.0 | 7.0 at tested lambda=0.75 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T5 Anchor B: per-node random bearing sigma deg | 0.075 | 0.075 at tested lambda=0.75 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T5 Anchor B: common bias magnitude deg | 0.1 | 0.075 at tested lambda=0.75 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T5 Anchor B: differential HALF-bias magnitude deg | 0.05 | 0.0375 at tested lambda=0.75 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T5 Anchor B: per-node per-axis position sigma m | 50 | 37.5 at tested lambda=0.75 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T5 Anchor B: known actual deployment abs(beta) deg | 20 | 15.0 at tested lambda=0.75 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T10 Anchor A: baseline km | 5.0 | 5.0 at tested lambda=0.75 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T10 Anchor A: per-node random bearing sigma deg | 0.05 | 0.05 at tested lambda=0.75 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T10 Anchor A: common bias magnitude deg | 0.1 | 0.075 at tested lambda=0.75 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T10 Anchor A: differential HALF-bias magnitude deg | 0.1 | 0.075 at tested lambda=0.75 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T10 Anchor A: per-node per-axis position sigma m | 100 | 75.0 at tested lambda=0.75 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T10 Anchor A: known actual deployment abs(beta) deg | 20 | 15.0 at tested lambda=0.75 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T10 Anchor B: baseline km | 7.0 | 7.0 at tested lambda=0.75 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T10 Anchor B: per-node random bearing sigma deg | 0.075 | 0.075 at tested lambda=0.75 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T10 Anchor B: common bias magnitude deg | 0.1 | 0.075 at tested lambda=0.75 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T10 Anchor B: differential HALF-bias magnitude deg | 0.1 | 0.075 at tested lambda=0.75 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T10 Anchor B: per-node per-axis position sigma m | 200 | 150.0 at tested lambda=0.75 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |
| T10 Anchor B: known actual deployment abs(beta) deg | 20 | 15.0 at tested lambda=0.75 | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |

Each joint vector must be considered as a complete vector with its anchor random sigma. Signed corners tested; no continuous interior guarantee. Differential HALF-bias is (b2-b1)/2, total differential twice this value. Position sigma is Cartesian per-axis1sigma, not radial RMS. Deployment changes true geometry while navigation perturbs estimator positions. Gate2A20deg is a tested grid edge, not a proven tolerance limit.

Primary half point requirements:

- Anchor A: T5 half point PASS; common magnitude0.05deg, half-diff0.025deg, position25m/axis, deployment10deg, frozen random sigma0.05deg.
- Anchor B: T5 half point PASS; common magnitude0.05deg, half-diff0.025deg, position25m/axis, deployment10deg, frozen random sigma0.075deg.

All synchronization, target association/front-back resolution and array/attitude calibration requirements require client evidence. Synchronization/association NOT_NUMERICALLY_VALIDATED, no timing tolerance is claimed. Separate hardware measured random and systematic errors, covariance, node/aperture configuration, directed-bearing accuracy statistics, navigation position covariance, heading calibration, deployable baseline, position updates and clock provenance are required.

No actual equipment admission; R4-A1-NEW NOT_OPENED; R4=0%. STOP numerical architecture work and obtain client hardware confirmation.
