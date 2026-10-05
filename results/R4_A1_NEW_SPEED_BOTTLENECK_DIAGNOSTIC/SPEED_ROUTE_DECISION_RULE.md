# Frozen development decision rules

OnlyAnchorA drivesroute;Bsecondary. Speedthreshold remains10%,strongspeed5%. All12case-levelP95 required; nearest rankceil(.95*500)-1. OtherPROJECTlimits:range10%,bearing1deg,heading5deg;D4solverfailurerate<=1%INEACHcase. Pooledmetricsneverrescue failingcase.

D1all12speed10% => STATIC_BEARING_BIAS_IS_SUFFICIENT_TO_EXPLAIN_PROJECT_SPEED_FAILURE withinfrozenpanel,notuniquecauseproof.
D1notall12 and D2orD3all12 => NAVIGATION_ERROR_MATERIALLY_LIMITS_SPEED. 'Materialclosure' frozenasall12, notpost-hocsubjectivecutoff.
D3anycase>10% => randombearing/frozengeometry/current4stateestimatorinsufficientinthiscontrol. Formalroute label SPEED_INFORMATION_INSUFFICIENT_AT_CURRENT_1200S_APPLICATION_SCENARIO is scope-limited; not universalno-information theorem.

Routepriority: ifD4all12fourPROJECTmetrics plus<=1%failure/case => BIAS_AWARE_SPEED_ROUTE_WORTH_FRESH_VALIDATION; nextrecommended R4_A1_NEW_BIAS_AWARE_FRESH_VALIDATION,neverautomatic. Else ifD3speedfailsanycase => SPEED_INFORMATION_INSUFFICIENT_AT_CURRENT_1200S_APPLICATION_SCENARIO;STOP. Else ifD1all12speedpasses => BIAS_INFORMATION_PRESENT_BUT_ESTIMATOR_NOT_YET_ESTABLISHED;STOP/noFIX2. Otherwise SPEED_ROUTE_NOT_ESTABLISHED_UNDER_CURRENT_SAVED_SCENE_DIAGNOSTIC;STOP. D4all12speed5% separately reports BIAS_AWARE_ROUTE_MEETS_STRONG_SPEED_ON_DEVELOPMENT_DATA, no scientificcredit.

No newperformanceclaim,randomseed,truthcase,noisedraw orMCrealization. R4-A1=0%;R4=0%;depthnotopened. No newbaseline,RMS,time,turn,acoustics ordepth. After completebothanchors,D0-D4,independentvalidation,CommitBpushSTOPforlead audit. Freshrequiresfutureauthorization,fullyfrozenestimator,newseedsandfresh/off-gridconfirmationpanel. PassingfuturePROJECT/STRONGalone canrestore15%progress, not thisdevelopmentdiagnostic.
