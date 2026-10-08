# Analytic pressure and state chain
Each mode contributes T_m=w_m sqrt(2π/(k_m ρ)) exp(-i k_m ρ-iπ/4).
Retain all cached modes and complex wavenumbers. ∂ρT_m=T_m(-i k_m-1/(2ρ)).
For the target p(t)=[r0 cos θ0+vx t,r0 sin θ0+vy t], the exact horizontal
chain matrix is [[cos θ0,-r0 sin θ0,t,0],[sin θ0,r0 cos θ0,0,t]].
Each element displacement divided by its range contracts with this matrix.
The normalized log-complex-field derivative divides ∂p by nominal pressure;
state scales remain [50000 m,1 rad,2 m/s,2 m/s,200 m].
Depth alone uses the already cached ±1 m and ±0.5 m outputs, never interpolation.
The independent bearing tangent rotates the source around each node center:
δp=[-center_y,center_x]δβ. No angular finite difference enters the primary chain.

Statistics are unchanged from the pilot: nominal fixed raw Gaussian covariance,
shared original-frequency noise, source scalars eliminated separately per node/time/frequency,
fixed complex element gains jointly shared across all frequencies and snapshots (C1),
and saturated response C2. No covariance derivative contributes information.
The linear CA chart exactly represents the rank-one spatial normalization locally.
P0/P1 and P0_DUAL/P2 comparisons use identical underlying noise/nuisance assumptions.
Finite candidates reuse saved pressures and the same joint fixed-gain log-chart feasible
fit; it is not a certificate for the global nonlinear gain optimum or global localization.
Frozen synthetic evaluation states define local points, never estimator priors.
