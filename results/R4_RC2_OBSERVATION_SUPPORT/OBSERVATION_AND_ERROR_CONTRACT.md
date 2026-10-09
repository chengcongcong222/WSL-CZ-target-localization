# Observation and error contract

U0 uses all21 saved0.1deg P/Q bearing arrays,121epochs. Archived relative_TL is not used. Exact platform motion and zero navigation/system bias are inherited, not extended to physical hardware. Physical domain r45..60km,theta-5..5deg,v1..3m/s,heading-15..15deg is a restrictive inherited prior.

Preferred H01/H06/H12 raw bearings are absent. Dual-node archives contain random draws and truths but not received bearings/reported positions. Reconstructing via scene(truth,...) is forbidden here. U is input-incomplete. No same-scene single/dual quantitative comparison or engineering coverage claim is possible.

A conservative reflection union retains unresolved left/right sectors even though the stored generator produced signed bearings. It never adds information or removes signed-compatible states. Gaussian coverage is conditional on the archived121-epoch iid vector law, exact navigation, constant motion and physical domain. Cross-case dependence is allowed by a union bound. Real correlated noise/navigation or array sign resolution is not certified.
