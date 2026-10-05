# Paired saved-realization interventions

Every original6000A+6000Bscene produces exactlyD0-D4. No performance-based subset or changedsampling. All scenes reconstructed with the frozen historicalgenerator and originalNPZuniform/normalarrays; never use any RNGcall.

D0 CURRENT_BASELINE:copyparentrunCSV byte-for-byte as identitycontrol, reuse saved estimates/errors, no optimizer rerun. New absolute speederror computed post-evaluation; original columns retained identically.
D1 ORACLE_BIAS_REMOVED:subtracttruegeneratedb1=common-HALFdiff,b2=common+HALFdiff from observedbearings; same noisyestimatedpositions/randombearingnoise; exactfrozen4stateestimator.
D2 ORACLE_NAVIGATION:truepositions tosame4stateestimator; originalbiased/randombearings unchanged.
D3 ORACLE_BIAS_AND_NAV:subtractgeneratedbias,use truepositions; keep exactsame randombearingnoise, actualbetadeployment and4stateestimator.
D1-D3=ORACLE_DIAGNOSTIC_ONLY; not deployablealgorithms, not applicationmetrics. r/theta estimates reference MAIN positionavailableto respectiveestimator; truthreferencealwaysoriginaltrueMAIN0. This referencechange in navigationoracle is explicit. Speed compares globaltargetvelocity inallbranches.

D4 BIAS_AWARE_6STATE:originalbiasedbearings,originalnoisypositions. UniqueobservationonlyE1initializer,common/diffboth0,one boundedunregularized sixstateTRF. Nooldrunwarmstart,truthbias,truthstateor multistart. Bounds designA common.05deg/HALFdiff.025deg,B.075/.0375deg. Onlycandidatealgorithm.

Every failed initialization/solve/nonfiniteestimate counts infinite errorinall4metricsandabsolute speed; no outlierremoval,replacementor redraw. Per-case500unconditional nearest-rankmedian/P90/P95/P99,failure rates. BiaserrP95 failureINF; fractionatbound uses validfitsdenominator explicit; no hiding solverfailures.
