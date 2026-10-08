# Plane-wave null: derivation and pre-registered numerical checks

Let P_m(f)=S(f)A(f;r,z) exp(i k(f)u·d_m), with common u independent of range/depth in this control. Then a_m=b exp(i[k(f+)−k(f−)]u·d_m), where b=S(f+)S*(f−)A(f+;r,z)A*(f−;r,z).
Consequently C_A[m,n]=exp(iΔk u·(d_m−d_n))/M. The nonzero common spectrum, common propagation amplitude and common propagation phase cancel exactly. Range/depth cannot create independent information at fixed bearing. Array bearing can change over a moving trajectory; such bearing geometry is not an independent propagation observable.

A fixed frequency-independent element gain g_m contributes |g_m|² to a_m: its phase cancels, its amplitude does not. Unknown amplitude gains must be profiled across all registered frequencies/times in C1. Arbitrary element/frequency gains can absorb any snapshot field; C2 is explicitly a saturated **snapshot-local** response absorption control. A frequency-dependent response constrained constant over a moving trajectory is a different model and is not claimed to eliminate all temporal information.

The numerical checks use the frozen 8×2m array, all registered pairs, different r,z and nonflat complex spectra. Maximum entrywise CA discrepancy must be <1e-10. Single-element normalized CA must equal1. Arbitrary response absorption discrepancy must be <1e-12. Any failure stops all physical interpretation. Numerical result is appended after the design commit; no scientific checks run before freezing.
