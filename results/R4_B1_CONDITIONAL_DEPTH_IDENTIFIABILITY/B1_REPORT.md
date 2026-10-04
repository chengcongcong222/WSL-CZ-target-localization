# R4 B1 conditional depth profile identifiability

Decision: **B1_DEPTH_IDENTIFIABLE_ONLY_UNDER_TIGHT_HORIZONTAL_CONDITIONING**. Execution is valid; research-lead audit is pending. Confirmed R4 progress remains 0%; proposed B1/overall progress is 0%.

The accepted historical six-case definition is S7 (201/235/283 Hz) x 180/200/220 m plus S8 (201/235/283/338 Hz) x the same depths. They share one horizontal scenario and noiseless matched acoustic observation rule; they are not six independent scenes or Monte Carlo replicates. The three-frequency primary route is S7; S8 retains the historical four-line companion without converting it to a new three-line score. All six are included in every all-case Gate.

Design was committed and pushed first at 0f93b162b15cd89340a6a2fba9ac26411c0c5fb2, parent 20846a9e1ebb10518762587635856c71763da790. Cases were selected by the complete accepted final_nominal_candidate_set, never by depth success. All twelve final survivor rows (psi=4 and psi=5) are separately cold reconstructed in B1_HISTORICAL_SCORE_REPLAY.csv. The old artifact stores profiled minima, not complete J(z) or a standalone acoustic observation array: B1_OBSERVATIONS.csv explicitly reconstructs the accepted generation rule after design push. No new signal/noise assumption is introduced.

H0 uses exact generation/reference horizontal state as an authorized controlled input. H1 uses all eight legal +/- inherited grid steps (1 km,0.5 degree,0.2 m/s,1 degree), one axis at a time. This describes the historically inherited local grid panel, not a promise that these errors are realistic post-acquisition errors or a demonstrated continuous tolerance radius. Score is the accepted per-frequency/per-window demeaned RMS over 121 times, without bearing or Q2 representation. Receiver label is 200 m; nearest stored source/receiver modal depths and all 21 depth labels are saved, with no interpolation.

| Case | H0 z-hat | Second margin (dB) | Neighbor margin (dB) | Curvature (dB) | Interior minima | W50 (m) | H1 <=10 m | H1 rank<=3 | Worst H1 (m) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| S7_Z180 | 180.0 | 1.11621 | 2.00563 | 4.05247 | 7 | 0.0 | 5/8 | 6/8 | 30.0 |
| S7_Z200 | 200.0 | 0.982204 | 0.982204 | 3.78467 | 7 | 5.0 | 5/8 | 4/8 | 50.0 |
| S7_Z220 | 220.0 | 0.853971 | 1.7071 | 3.4558 | 6 | 0.0 | 7/8 | 4/8 | 40.0 |
| S8_Z180 | 180.0 | 1.20415 | 1.91685 | 4.04349 | 7 | 0.0 | 5/8 | 4/8 | 25.0 |
| S8_Z200 | 200.0 | 0.9953 | 0.9953 | 3.49833 | 7 | 5.0 | 5/8 | 4/8 | 50.0 |
| S8_Z220 | 220.0 | 0.794481 | 1.56346 | 3.15885 | 6 | 0.0 | 5/8 | 4/8 | 35.0 |


| H1 axis | Profiles | Exact | <=5 m | <=10 m | Rank<=3 | Worst error (m) |
|---|---:|---:|---:|---:|---:|---:|
| r_km | 12 | 0 | 4 | 4 | 2 | 50.0 |
| theta_deg | 12 | 12 | 12 | 12 | 12 | 0.0 |
| v_mps | 12 | 0 | 1 | 4 | 0 | 35.0 |
| psi_deg | 12 | 12 | 12 | 12 | 12 | 0.0 |


H0 unique true-depth minima: 6/6. Positive second/neighbor margins: 6/6 and 6/6. Positive curvature: 6/6. H1 exact/within5/within10/rank<=3: 24/29/32/26 of 48. Worst H1 error: 50.0 m. Other strict local minima and boundary minima are preserved; a unique global truth minimum does not imply one local valley. Full profiles/metrics are exported without dropping any row.

W50 is the span between terminal labels of the contiguous component J_norm<=0.5 containing the deterministic global minimum. It is not a confidence interval; a width of 0 m means one retained depth label, not zero uncertainty. No spline, parabola, continuous depth or off-grid z is used. The 5 m nuisance labels may map to nearby stored physical samples; numerical mapping is retained explicitly per frequency.

Numerical tolerance was computed by the preregistered rule max(128*eps*max(1,max_abs_J),10*max_repeat_or_independent_residual), not from a depth margin. Result: 4.30873114965e-10 dB. Cold repeat maximum: 0 dB; independent score/observation maximum: 4.30873114965e-11 dB. Independent parser uses struct, geometry uses scalar math, propagation uses ordered compensated mode sums, and score uses math.fsum. This supplies a numerical implementation check of the same matched physical model, not a second physical model validation. 11 unit tests pass; runtime validation has 8505 passing checks, including 7241 historical byte hashes. Saved metric/Gate reconstruction is separately recorded.

The CZ certified global acquisition branch is CLOSED: WINDOW_SHAPE_Q2 + N0 + whole-cell conservative partition over 45-60 km cold-start 4D support did not establish closure under registered budgets. Its mathematical lower-bound mechanism remains historical evidence. No C4/F4, grid/tau/enclosure fix, budget increase or 21-case CZ development is run. A1 remains BLOCKED / NOT_COMPLETED, 0/15%.

Maximum claim is ORACLE/EXACT-HORIZONTAL-CONDITIONED DEPTH MECHANISM in synthetic matched E0, on-grid source labels, with H1 outcomes explicitly bounded by the fixed panel. A CONDITIONAL result supports the exact-conditioned mechanism and identifies horizontal/depth coupling; it does not establish an untested tighter continuous neighborhood. It gives no cold-start joint inversion, end-to-end depth performance, environmental robustness, real ocean/UUV observation or final accuracy guarantee. B2 is NOT_OPENED; A2/A3/A4/B3/B4/C and P5 are NOT_OPENED. Stop after execution commit/push and wait independent audit.

See [all saved profile curves](B1_DEPTH_PROFILE_CURVES.svg), [profiles](B1_DEPTH_PROFILES.csv), [metrics](B1_PROFILE_METRICS.csv), [decision](B1_DECISION.json), [validation](B1_VALIDATION.json).
