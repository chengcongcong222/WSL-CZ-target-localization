# Preregistered route decision

Primary Anchor A; B remains secondary. All scientific/application completion credit remains zero. All empirical quantiles use nearest rank over 500 samples including failed runs as INF. PROJECT gates unchanged: range/speed <=10%, bearing <=1 deg, heading <=5 deg. E1 route also requires every case failure <=1%. STRONG speed/CRLB <=5% is diagnostic only.

Priority I: E1 A has all 12 speed PROJECT cases, all other three PROJECT metrics and failure gates -> GAUSSIAN_MLE_SPEED_ROUTE_WORTH_FRESH_VALIDATION; recommend FRESH_VALIDATION only, never run it automatically.

If any A FULL saved FIM lacks full rank -> CRLB_NUMERICAL_RANK_NOT_CLOSED, STOP.

Priority II: E1 is not 12/12 speed PROJECT, every case |E1 speed_P95-E2 speed_P95| <=0.02 in relative fraction (2 percentage points), and at least one A case max saved FULL speed_196 equivalent >0.10 -> CURRENT_1200S_GEOMETRY_INFORMATION_LIMIT_CONFIRMED_BY_CRLB. Recommend OBSERVATION_DESIGN / R4_A1_NEW_SPEED_OBSERVATION_DESIGN_GATE only. Accurate meaning: local Gaussian 1.96-sigma equivalent exceeds target for at least one tested primary geometry; not a universal empirical P95 lower bound.

Priority III: all 12 A case max FULL speed_196 <=0.10 but E1 not all speed PROJECT -> ESTIMATOR_OR_NONLINEAR_EFFICIENCY_GAP_REMAINS. Recommend ESTIMATOR_RESEARCH only; execution stops for research-lead decision.

Other: E1/E2 case P95 maximum gap >0.02 -> OPTIMIZER_OR_INITIALIZER_INEFFICIENCY_PRESENT, recommend ESTIMATOR_RESEARCH only. Uncovered pattern -> INFORMATION_EFFICIENCY_AUDIT_INCONCLUSIVE, STOP. E1/E2 paired cost/state/speed gaps additionally reported, no post-result threshold changes.

No new scenarios, filters, Doppler, range-rate, TDOA, acoustics, depth or automatic next stage. Commit A must be pushed before tests involving FIM, any E1/E2 re-solve or CRLB. Commit B follows full saved-scene reconstruction and is pushed, then STOP. Remote equality and seal verify must complete. Independent audit may require investigation of numerical integrity, but frozen scientific policy cannot be retuned after execution.
