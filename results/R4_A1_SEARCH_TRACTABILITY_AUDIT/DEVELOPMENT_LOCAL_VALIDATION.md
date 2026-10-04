# Development local validation

Single frozen experiment: 126 raw runs / 2646 branches. No rerun and no algorithm repair.

Mechanical reconstruction: 14479 passed / 0 failed. See DEVELOPMENT_RECONSTRUCTION_AUDIT.csv and DEVELOPMENT_SCORE_RECONSTRUCTION.csv.

Pre-release tests: 74 passed after creating the missing runtime_tmp parent directory. The earlier setup attempt had 67 passed / 7 WinError 3 setup errors before any noisy run; no frozen algorithm was changed. See DEVELOPMENT_ENVIRONMENT_NOTES.json. Post-run unchanged tests:

........................................................................ [ 97%]
..                                                                       [100%]
74 passed in 4.45s


Frozen METHOD, EXECUTION and accepted CHECKPOINT hashes revalidated after execution. New evidence receives its own manifest; historical evidence remains unchanged. No fresh experiment was released.
