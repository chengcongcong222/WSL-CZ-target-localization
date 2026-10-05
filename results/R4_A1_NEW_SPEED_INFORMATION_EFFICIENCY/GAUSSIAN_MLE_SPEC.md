# Gaussian L2-MLE frozen specification

DEVELOPMENT ONLY. E0 reads accepted D3 results without solving; every result/state/error field must be identical. E1 uses the same saved D3 true nodes, bias-removed bearings, sigma and 121 epochs / 242 observations. Only loss changes from soft_l1 to linear. State [x0,y0,vx,vy]; scales [50000 m,50000 m,2 m/s,2 m/s]. Original ray-intersection + five IRLS initializer iterations, original analytic wrapped bearing residual/Jacobian. scipy least_squares TRF, max_nfev=100, ftol/xtol/gtol=1e-10. f_scale=1.5 is retained but mathematically inactive under linear loss. No bounds, priors, rejection, multistart or truth access in E1 API.

E2 uses the identical L2 solver with exact true Cartesian initial state, explicitly ORACLE_INITIALIZATION_DIAGNOSTIC_ONLY. No perturbation, extra start or claim. All A and B scenes run once (6000 each). Failures remain INF in all Gate quantiles. E0 never re-solves. Signed speed bias/median/sample std ddof=1/RMSE are reported; finite-only statistics are descriptive, never a way to discard failed runs. eta_v is defined only with zero failures and finite positive variance, using mean conditional saved-geometry local CRLB variance divided by signed empirical variance.

No new noise, truth, MC, time, geometry, turn, relative motion, acoustics or depth. Any subsequent stage requires new explicit authorization.
