# Audit implementation correction

Frozen design SHA1eb535fca52077e628777e7a9ab63e218f0d551f retained; original r4_aux_gate0_audit.py byte-identical. A precommit feasible-cell-count edit accidentally changed two SHA256 identifiers to SHA240. Original audit failed at the first hash call. An initial correction fixed only the function name; it passed numerical checks then failed at the policy-hash field lookup. Neither attempt reached figure/report/finalization writes. Primary geometry execution was not repeated.

Corrected audit source differs in exactly two SHA identifier tokens: hashlib.sha240 -> hashlib.sha256 and policy_sha240 -> policy_sha256. All numerical formulas, scenes, draws, metrics, tolerances, Gates and presentation logic unchanged. Correction JSON binds both source hashes and failed audit attempts. Use r4_aux_gate0_audit_corrected.py --finalize/--verify. This is an audit-infrastructure correction, no scientific credit.
