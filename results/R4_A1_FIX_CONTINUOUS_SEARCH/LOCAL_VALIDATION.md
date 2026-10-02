# Local execution validation

This is local consistency validation, not the research lead's independent scientific audit. A1-FIX remains blocked; R4 remains 0%.

- 18 tests pass via `python -m unittest discover -s tests -v`. New repair tests independently construct continuous bearing observations, disable truth error/neighborhood helpers during estimation, compare theta profiling against scalar minimization, check geometric bounds, and ensure DE theta-offset diversity and continuous likelihood feasibility.
- 338 checks pass via `python r4_a1_fix_audit.py`. All 36 cases reconstruct bearing costs, exact shared-depth TRIPLE scores, final errors, survivor counts and RC2 hypotheses. The largest exact-score reconstruction discrepancy is 4.312e-11 dB (tolerance 1e-7 dB). Direct bearing formulas and independent RMS reductions are used; modal propagation shares the audited forward function. All freeze hashes and selected frozen legacy challenge bytes agree.
- All 27 frozen A1 failure cases replay, and the R3 385-node cloud plus all TRIPLE scores/survivors at three depths agree with frozen evidence. No old result file is rewritten.
- Continuous forward self-match: noiseless cost 0 within 1e-20 rad-squared tolerance; maximum exact profiled J 4.786e-14 dB within 1e-10 dB tolerance.
- Pressure interpolation: 256 fixed random ranges per frequency, all 21 nuisance labels; maximum TL discrepancy 1.117e-5 dB within the frozen 0.001 dB tolerance. Final reported candidates use exact modal scoring.
- 315/316 local optimization runs report convergence. All 42 global DE islands reach the fixed iteration budget rather than report convergence. A high local convergence fraction therefore does not establish global basin recovery. DE nfev counts vectorized objective batches; counts/upper bounds are explicit in CONVERGENCE_AUDIT.csv.
- Both scientific figures were visually inspected. They identify the nominal-only, fixed small-panel scope and distinguish regression from held-out confirmation. They do not claim a sensor tolerance.

The frozen method's eight local basins per island and four exact-polish basins are caps. Population collapse or deduplication may yield fewer starts. MODE_INVENTORY.csv distinguishes the certified single noiseless bearing solution from finite noisy candidate clusters; candidate clustering is not an exhaustive physical-mode count.

Initial implementation V1 had a degenerate global theta-offset coordinate. Its 22 completed cases, frozen panel/configuration and executor snapshot are preserved under INITIAL_IMPLEMENTATION, and excluded from V2 confirmation. The correction changes initialization diversity without changing search budgets or likelihood/acoustic thresholds. A fresh six-case deterministic holdout (seed 2026100203) was frozen before any V2 holdout observations. It contains six noiseless controls and twelve nominal 0.1-degree realizations. No extra trials were added to make the Gate pass.

CONTINUOUS_ALIAS_DIAGNOSTIC.csv checks exported states for separation at the frozen coordinate resolution and joint matched-observation tolerances of 1e-8 degree bearing RMS and 1e-6 dB acoustic RMS. No distinct joint exact match is demonstrated. This finite test does not prove global uniqueness under noisy observations.
