# Whole-cell enclosure: scope and proof

This checkpoint establishes a conservative mathematical enclosure for the specified fixed-input modal level model. It establishes neither useful contraction at the frozen terminal resolution nor development coverage/convergence. No P/Q numerical level or feature was loaded. Large cells can give zero mismatch lower bounds and enormous upper bounds; such cells must remain unresolved.

## Exact inputs and arithmetic

Python-FLINT 0.9.0 Arb/Acb with 160-bit midpoint precision evaluates real/complex balls with outward error radii. Binary modal coefficients, k components, the inherited binary64 platform positions, candidate endpoint numbers and observation coefficients are converted from their exact integer ratios. Candidate trigonometric degree conversion uses a certified real pi. Reference reconstruction uses 256 bits. Do not execute this backend concurrently in threads: the precision context is global. Package and source identities are in the freeze manifest.

The backend guarantees enclosures for elementary ball operations; it does not guarantee tight bounds. See the primary [Arb API](https://python-flint.readthedocs.io/en/latest/arb.html), [Acb API](https://python-flint.readthedocs.io/en/latest/acb.html), and [precision and ball semantics](https://python-flint.readthedocs.io/en/latest/general.html). The verified official wheel digest and installation fallback are recorded in calibration; the existing pip vendor import was broken, so the pinned wheel was hash-checked and extracted on D.

## Trajectory and bearing

For each time t, a Cartesian box is propagated from all four state intervals:

`dx = 1000*r*cos(theta) + v*t*cos(psi) - xp(t)`

`dy = 1000*r*sin(theta) + v*t*sin(psi) - yp(t)`.

The full inherited box has certified dx>0 and R=sqrt(dx^2+dy^2)>0. Thus atan(dy/dx) encloses the principal bearing without an atan branch crossing. Dependency inflation can widen intervals but cannot remove states.

RC2 cost remains the SUM of squared wrapped radian residuals over 121 samples, not mean cost or degree residuals. If a residual interval is certified inside [-pi,pi], its absolute lower and upper endpoints produce valid squared bounds. If wrapping is uncertain, the contribution is conservatively [0,pi^2]. The observation-side cutoff is supplied unchanged; for the 21 future cases it must be read from the frozen manifest, never fitted from envelope results. This PRE-RUN does not recompute any case cutoff.

## Finite profile and levels

Each inherited label z=150,155,...,250 maps independently, per frequency, to the nearest stored modal source-depth sample. Receiver depth maps from 200 m. The 63 label-frequency mappings are exported; duplicated effective depths remain explicitly enumerated. There is no continuous depth estimation.

For each cell and depth:

`p_f(t) = sum_m phi_m(z)*phi_m(200)*sqrt(2*pi/(Re(k_m)*R(t))) * exp(-alpha_m*R(t) - i*Re(k_m)*R(t) - i*pi/4)`

with alpha=-Im(k). The parser verifies positive real wavenumber and finite modal inputs. Acb enclosures contain every possible p. If a<=|p|<=b, level bounds are constructed as

`20*log10(max(a,p_floor)) <= L <= 20*log10(max(b,p_floor))`,

with the inherited binary64 p_floor=1e-30. **The positive endpoints are logged before their union.** A broad amplitude ball may straddle zero; logging that broad ball would incorrectly create a numerical failure. No mode truncation, interpolation, sample-based remainder, phase observation or source spectrum input is introduced. Shared modal factors inside a single all-depth batch are common expressions, not a persistent unmetered cache.

## Certified Q2 projection and mismatch

For each n=61 or 60, let x_i=2i-(n-1), q_i=x_i^2-mean(x^2). Define b1=x/sqrt(sum x^2), b2=q/sqrt(sum q^2). Exact symmetry gives sum x=sum q=sum x*q=0, so these are orthonormal under the inherited discrete Euclidean inner product and orthogonal to the constant. Numerator sums are exact rationals; square roots and division are enclosed.

Explicit per-line/window centering followed by each b dot product encloses the ideal coefficients. Any constant per-line/window N0 offset cancels algebraically. The predicted 12 intervals enclose the image of the whole joint state cell. They can overestimate that image by discarding coefficient correlations, which is conservative.

For observed binary64 coefficient y_j treated as an exact point and interval I_j,

`LB_z^2 = sum_j distance(y_j,I_j)^2 / 6`

is no greater than any same-depth feature mismatch in the cell. Using the largest absolute endpoint instead yields UB_z. Whole-profile bounds are `min_z LB_z` and `min_z UB_z`; all 21 labels must be completed before an acoustic rejection is admitted. The minimum upper bound is valid because one fixed depth with a uniform upper bound suffices over the entire cell. Bearing-only rejection may precede acoustic profiling and records that branch explicitly.

Binary64 API and original forward reconstruction differ slightly from this specified ideal-real modal model. tau_F=1e-7 dB is the preregistered fixture-derived numerical compatibility allowance, not a proof of a uniform error bound for all binary64 operations and not an acoustic noise model. The whole-cell certificate itself bounds the specified ideal model exactly; future reference retention must still be independently evaluated. Fixed-corner calibration and enclosure containment tests are regressions supplementing the analytic argument, never rejection certificates.

## Classification and coverage

Only a certain Arb comparison `bearing_LB > frozen_cutoff`, or completed `feature_LB > tau_F`, can reject a whole cell. A sample, a local solver or an ambiguous comparison cannot reject it. A uniform compatible upper bound yields RETAINED_COMPATIBLE. Everything else is RETAINED_UNRESOLVED and is split only under the frozen rule/budget. BOUNDARY_CENSORED is an orthogonal flag and does not delete support.

The queue initially contains the complete domain. Splitting replaces one closed box by two closed children whose union is the parent. Each terminal box is either certified rejected or retained. On a cap or invalid bound, the current box and every queued box are retained. Induction on this partition proves that no unprocessed region disappears. A cap, terminal loose bound or invalid bound gives coverage_closed=false; invalid arithmetic also gives execution_valid=false. Empty support is a future scientific failure even when its exclusions were certified.

Each raw coarse and raw fine run has a new engine, cache, queue and counter object and starts at the complete domain. No cumulative retained set enters the Gate. All six operation categories and aggregate hard limits are charged before the relevant work. Atomic all-profile admission prevents partial work from bypassing a cap. Charged counters represent admitted model/bound calls, including calls that later fail; attempts and denied transactions are separately logged. An invalid calculation cannot certify coverage.

The 60/120 s runtime envelopes are calibration targets used to derive strict operation caps. They are **not hard wall-clock guarantees**; an individual admitted cell can exceed its measured fixture cost. There is no claim that the 14/29-cell budgets can resolve the frozen grids. Loose intervals and an incomplete run produce broader support and a failed future Gate, never a point-grid substitution or additional budget.
