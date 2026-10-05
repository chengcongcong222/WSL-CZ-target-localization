# AUX Gate2B frozen joint static error budget

Parent Gate2A independently accepted: SINGLE_FACTOR ONLY. This is the last geometry-only numerical Gate before client hardware confirmation.

Exactly two frozen anchors: A5km/random0.05deg, B7km/random0.075deg; random terms do not scale. T5 vectors both(0.10deg common,0.05deg differential HALF,50m position/axis,20deg actual deployment). T10 A(0.10,0.10,100,20), B(0.10,0.10,200,20). Each vector scales simultaneously by lambda[0,.25,.5,.75,1]. Every positive lambda includes8 signed common/diff/beta combinations,11 integer ranges50..60km and both mirrors. Zero uses one sign identity.1452cells per family,2904total,30000trials each.

Primary: at least one T5 anchor at lambda=.50 has worst P95<=5% => AUX1_HALF_T5_JOINT_BUDGET_ESTABLISHED and route AUX1_JOINT_STATIC_REQUIREMENT_ESTABLISHED_CLIENT_CONFIRMATION_REQUIRED. Otherwise primary NOT_ESTABLISHED; if independent T10 family half point passes at least one anchor => AUX1_JOINT_SUPPORTS_10PCT_NOT_5PCT, else AUX1_STATIC_GEOMETRY_TOO_FRAGILE_FOR_PROJECT_TARGET. Complete both ladders regardless primary results. No post-hoc anchor/threshold/seed changes.

All failures infinite; nearest-rank unconditional quantiles. Report P95/P99, failure rates and absolute2D position P95; joint aggregate worst over signs/ranges/mirrors. Range output uses estimated MAIN reference, compared with true MAIN range; global position separately compared with truth. True deployment generates bearings; noisy positions enter estimator. b1=common-half_diff, b2=common+half_diff; half_diff is not full difference.

New seed2026100504 PCG64 six-column draws frozen before MC. Independent slope reconstruction, full grid regeneration, per-cell quantiles/failures, mirror controls, old Gate2A exact shared-draw quantile reproduction and preregistered two-sample DKW zero-control criterion required. The policy JSON specifies familywise alpha=.001 over44 unique controls. Failure invalidates execution. Shared draws and reflected mirrors are not independent repeat samples.

Tested corners and discrete ranges only: no continuous hyperbox or interior guarantee. Maximum tested lambda is an observed discrete frontier, not a tolerance boundary. Tested beta20deg is a grid-edge observation, not a maximum physical tolerance. Actual hardware UNKNOWN; time/association NOT_NUMERICALLY_VALIDATED. No propagation/depth/SSP/TDOA/tracker/motion/time-skew/association-error experiment. R4=0%; A1-NEW NOT_OPENED.

Two commits: freeze joint static error budget (push verify before MC), evaluate joint static error budget (push verify thenSTOP). Regardless outcome STOP numerical architecture work; client confirmation or explicit research-lead route decision required. No automatic Gate2C. Preserve old Gate0/1/2A artifacts; updated client table in this stage supersedes operational presentation without editing frozen originals.
