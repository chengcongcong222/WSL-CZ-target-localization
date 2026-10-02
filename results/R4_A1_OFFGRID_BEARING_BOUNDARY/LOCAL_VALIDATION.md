# R4-A1 local validation

This verifies local numerical consistency. The research lead's independent scientific Gate audit remains pending, and R4 progress remains 0%.

- `python -m unittest discover -s tests -v`: 14 tests pass, including the nine existing R3 audit tests and five new R4 boundary tests. The latter cover circular angles, all nearest-node ties, Manhattan adjacency without edge wrapping, direct residual agreement and observation immutability, and nuisance-depth mapping/RMS scoring.
- `python r4_a1_result_audit.py`: 162 checks pass. For all 27 pilot cases, independent geometry and direct circular residual sums across the full stored grid reproduce exact RC2 node IDs. Saved catalogs then reproduce top1, survivor IDs, errors, floors, excesses and survivor envelopes. Additional checks verify continuous-truth observations, exact truth relative TL, three recomputed modal-score samples, all RC2 summary fractions, design hashes and source hashes. Results are in RESULT_RECONSTRUCTION_AUDIT.csv.
- R3 replay: exact 385-node control cloud and all 1155 TRIPLE candidate scores at three nuisance-truth depths agree with frozen evidence; survivor IDs agree exactly. No frozen R3 source or result file changed.
- Both scientific figures were visually inspected. They explicitly distinguish RC2-only statistics from the 27-case nominal end-to-end pilot. No complete multi-sigma accuracy boundary is claimed.

The modal-score sample recomputation shares the forward function with the executor; the bearing geometry/residual and final metric reconstruction use independent formulas. API separation checks and artifact reconstruction do not constitute a formal proof that every possible future implementation is free of truth leakage.

Input-source hashes were captured after execution. DESIGN_FREEZE_MANIFEST.json records the pre-run design freeze. R4 evidence bytes are preserved through local Git attributes; legacy CSV input validation uses the additionally recorded LF-normalized hash for checkout portability.
