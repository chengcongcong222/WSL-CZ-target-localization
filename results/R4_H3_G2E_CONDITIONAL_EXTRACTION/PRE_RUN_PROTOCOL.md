# Frozen execution protocol

1. Push the design commit and verify remote/main equality.
2. Run r4_h3_g2e_audit.py pre: independent synthetic controls only.
3. If any precontrol fails: classify IMPLEMENTATION_INVALID and STOP without executing cells.
4. Run r4_h3_g2e.py exactly once. RUN_ONCE marker prohibits silent reruns.
5. Run independent cold reconstruction r4_h3_g2e_audit.py audit, then frozen reporting.
6. Push execution commit, independently verify remote HEAD and STOP.

The execution and cold-audit cap is 3600 seconds; delivered storage cap is 1,000,000,000 bytes.
No extra optimizer, solver, frequencies, time-domain simulation or depth localization is authorized.
If a frozen interface fails, preserve partial outputs and classify it; do not patch and rerun.
