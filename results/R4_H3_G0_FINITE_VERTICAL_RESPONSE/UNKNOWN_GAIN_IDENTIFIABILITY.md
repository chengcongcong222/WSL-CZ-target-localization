# Registered nuisance models and depth identifiability

C0: unknown arbitrary complex source at every frequency/snapshot; array responses known. Mechanism upper bound only.
C1a: source as C0; fixed per-element complex gain shared across every frequency and snapshot.
C1b (PRIMARY): source as C0; per-element/per-frequency complex gain shared across 0/600/1200 s. This retains only temporal changes that cannot be absorbed by fixed frequency-dependent calibration. Gains at different times are NOT fitted separately.
C1g_FREE_MODAL_GAIN: completely free mode-dependent complex gain. At regular nonzero source coefficients, dphi_m=phi_m*(dphi_m/phi_m), hence depth Jacobian lies in modal-gain tangent space. This is a mandatory ZERO-depth control, not an admitted useful-information model. It does not estimate or expose oracle mode amplitudes to a practical receiver.
C1g_PHYSICALLY_CONSTRAINED_GROUP_GAIN: NOT_ADMITTED_NO_INDEPENDENT_GROUP_CONTRACT. The current cache/reference set contains no separately justified smoothness/coherence/group-gain constraint for CZ unknown propagation. No group count or low-dimensional gain is chosen to generate a positive result. Thus no classification B is available in this execution.
C2: independent unknown complex response for each element/frequency/snapshot; its tangent is the full raw observation space and all target information is zero.
C1e: SSP, source/receiver depth error and pose nuisance NOT_EVALUATED with these fixed caches. All positive results remain conditional on exact environment and geometry.

Free modal gains are more flexible than C1b channel calibration. A C1b conditional positive cannot refute their analytic zero information. Any application claim requires an independently justified propagation/calibration restriction first.
