# E1-G0 nuisance lock
N0 each log F unknown, receiver reference known, source drift zero.
N1 N0 + unknown shared fractional source drift.
N2 N1 + independent node fractional constant rho and linear k; primary NO_CALIBRATION.
N3 N1 + independent additive Hz B,K, fractional errors fixed zero: different model.
N4_BOTH N2 + additive B,K; N4_LINE N2 + per-line drift dl; broader controls separately retained. Shared d/per-line dl redundancy removed only by SVD gauge quotient.
N5 N2 + independent receiver rho sigma5e-6 and kappa sigma5e-9/s calibration priors. Hypothetical requirement, NOT existing hardware evidence. No source-frequency prior.
Each model full registered covariance and source-time convention from design. Fractional errors dimensionless; Hz errors never substituted for fractional clock uncertainty.
N0 is an ideal-reference upper control, not proposed prior-free result.
SUN-MEASUREMENT uses locked Sun Eq4/unknown F, no UKF, PAPER_ADAPTED_MEASUREMENT_INFORMATION_CONTROL. Source-state random-walk process covariance is omitted in this local static information control and explicitly NOT full original dynamic CRLB. SUN complete MFB-AUKF ORIGINAL_ESTIMATOR_NOT_FULLY_LOCKED.
All rank thresholds, covariance assumptions, calibration requirements and corner values fixed before scientific evaluation. No numerical result used to select model.
