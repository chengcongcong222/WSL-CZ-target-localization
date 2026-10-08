# Structural single-HLA rank
At a snapshot the element ranges are functions solely of node-center range ρ_t
and instantaneous bearing β_t. In the frozen horizontally range-dependent
environment the pressure also depends on the shared source depth z.
Consequently J_horizontal=J_ρ D_ρ + J_β D_β. Profiling independent β_t at
three snapshots removes the second term, leaving at most three center-range
directions and one depth direction: rank(P1) <= 4.
Let D_ρ be the 3×4 scaled horizontal center-range Jacobian. Its signed
3×3 cofactors construct h in ker(D_ρ), independently of the acoustic field.
The full state null vector is n=[h,0]. Whitening, CA transformation and
further source/gain nuisance profiling cannot restore this lost direction.
Both depth stencils have n_z=0 and therefore must annihilate this same n.
Measure the unmodified projected analytic Jacobian response to n BEFORE any
quotient/projection. Keep the original 1e-10 singular-rank threshold.
Coordinate estimates with a nonzero null-space component have infinite uncertainty,
not zero variance. Profiles that remain estimable use orthogonal nuisance removal.
P2 has six center ranges; its geometry may have horizontal rank four. With depth,
rank five is allowed but measured rather than imposed.
The independent audit constructs n via the geometric SVD and pressure gradients
via a separate per-mode scalar sum and Cartesian contraction.
