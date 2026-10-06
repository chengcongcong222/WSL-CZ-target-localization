# E1 deterministic algebra policy
Registered before checking. These are dimensionless toy identities, not acoustic scenarios, propagation, signal synthesis, Monte Carlo, estimator performance, or a coverage certificate.

| ID | Fixed inputs | Expected identity |
|---|---|---|
| A | F=[2,3], D=[0.9,0.8], a=1.1 | F'=aF, D'=D/a preserves all products |
| B | t=[-1,0,1], source drift column t, radial acceleration column -t | joint rank 1 |
| C | two-node block diagonal [1,t], affine columns per node | nuisance projection removes affine motion |
| C2 | same nuisance, curvature [1,-2,1] on node 1 | curvature survives projection |
| D | line loading [1,2,3] | shared Doppler/shared fractional drift rank 1 |
| D2 | additive Hz bias [1,1,1] versus line loading [1,2,3] | combined rank 2; distinct clock models |
| E | two-node common source column [1,1], differential [-1,1] | differential survives common source, not free per-node nuisance |
| F | same-column motion/reference, reference prior precision 1 | acoustic effective info 0, with prior 0.5 |
| SUN-W | n=7, kappa=-3 | printed sum kappa*(n+1)/(n+kappa)=-6, not 1 |
| SUN-F | T=10, constant velocity component 2 | printed next velocity 20 differs from constant velocity 2 |
| SUN-H | p=[3,4], v=[-1,2], F=200, c=1500 | derivative of Eq.4 wrt position is generally nonzero |

Projection uses P=I-N pinv(N), rank relative tolerance 1e-10. Projection norm <=1e-12 counts zero. Expected nonzero norms >1e-6. These checks validate structural statements only. Full target dynamics/bearing rank remains uncomputed this turn. SUN checks diagnose printed-formula inconsistencies; they do not establish an author-approved correction.
