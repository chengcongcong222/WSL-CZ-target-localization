# AUX1 measurement model and physical geometry

MAIN retains single towed HLA. AUX adds one mobile passive bearing-producing node. No auxiliary HLA/VLA, active sonar, cooperative beacon or known source waveform is assumed. A plain omnidirectional hydrophone does not by itself provide bearing; a real bearing chain remains a future requirement.

At common time, theta_i=atan2(y_T-y_i,x_T-x_i)+n_i, i=1,2. Exact known node positions; directed bearings associated with the same target; independent unbiased equal-sigma Gaussian errors. Front/back, target association, clock/latency, attitude, positioning, bias/correlation and detection are not validated.

MAIN=(0,0), target=(R,0). Target-centered smaller LoB angle alpha; baseline B=norm(p2-p1). Auxiliary target distance d satisfies d^2-2R cos(alpha)d+R^2-B^2=0. Both d=R cos(alpha)+/-sqrt(B^2-R^2 sin(alpha)^2) roots and both reflected node placements are registered; positive real d required. Feasible iff B>=R sin(alpha) in this grid. B,R,alpha are not independent physical knobs. Impossible cells are retained; no invented placements or replacement geometry. No numerical geometry/score/MC evaluated before freeze/push.

AUX p2=(R-d cos(alpha),side*d sin(alpha)), true directed AUX bearing=-side*alpha, MAIN bearing=0. Estimator receives only positions and observed directions to intersect forward rays. R truth sets the frozen synthetic scene and evaluation only; no truth clipping/range constraint/initializer in estimator.

H_i=[-q_iy,q_ix]/norm(q_i)^2; F=H^T Sigma^-1 H; local covariance F^-1. H/F condition numbers and range/cross-range sigma describe local linearized lower-bound/approximation, not actual-estimator guarantee. Two nonparallel bearings observe instantaneous horizontal position; no instantaneous z,v,psi information. Temporal paired positions could support motion, but no tracker or motion accuracy here.

1920 requested cells. 20000 Gaussian draws per feasible cell after design commit/push only, PCG64seed2026100501. Common random numbers across cells, sign reversal for mirrors: mirrors are symmetry controls, not independent trials. No acoustic level/propagation/CZ/depth/SSP/optimizer/5D estimator.
