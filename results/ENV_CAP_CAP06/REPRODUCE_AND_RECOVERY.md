# Reproduce / recover CAP06

Scientific files and thresholds are locked by CODE_FREEZE.json at Commit A000ac28f43b8689c45a22a5e1ce72263be7b70d3. Read RUN_STATE.json and every units/*_SCIENCE_STARTED.json before any action. A started unit cannot retry. Do not run development_tests.py again as a scientific result; its saved P0 controls and failed first fixture are already retained.

Runtime: Python3.10.11; versions and paths in P0_PRECHECK.json. Dependencies are loaded read-only from the existing D runtime. Set PYTHONDONTWRITEBYTECODE=1, PYTHONUTF8=1, OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1. All writes go to D:\ProjectStorage\WSL-CZ\ENV_CAP_CAP06. Dot-sourcing the historic runtime environment changes cwd, so Set-Location to this worktree afterwards.

Large public receiver samples, private labels, extracted products and P3 boxes are on D and indexed by SHA256. They are intentionally not duplicated in Git. Private labels are allowed only for generation and post-extraction evaluation. The receiver/support processes contain read guards. Main receiver input is physically the same file for S0/S2. This is an in-process coding separation, not a security proof against hostile software.

Formal phases use code/run_phase.py with P1/P2/P3 once each; the wrapper records code/input/config identities, wall time and process-tree memory and enforces the lower frozen budget. Independent audits do not import the scientific calculation functions. Reporting programs only read saved outputs and do not regenerate signals or optimize states. Their post-freeze creation changes no science threshold or result selection.

CP1, CP2 and FINAL are automatically committed and pushed to research/env-cap-cap06. No merge to main. Final SHA cannot be self-embedded in the same commit; remote verification is returned to the user and recorded outside the committed artifact after push. Historical R4 remains0%.
