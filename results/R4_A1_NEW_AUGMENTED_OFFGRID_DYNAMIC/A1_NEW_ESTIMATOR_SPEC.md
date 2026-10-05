# Observation-only estimator freeze

API estimate(times,node_positions,bearings,sigma_deg); no truth/state/range/speed/heading/error argument. Separate generator/evaluator permittedtruth; standaloneestimator imports onlyNumPy/SciPy. No truth warm starts, target-dependent formation or historic state-grid seeds. Runtime TypeError truthargument rejection and AST audit required.

E1 directed-ray intersection determinanttol1e-12. Onlyparallel/nonfinite/behind-sensor intersections invalid; require>=4valid times withpositive time span. No truth/error/threshold based rejection. Fit Cartesianposition on[1,t/1200] ordinary LS thenexactly5IRLSiterations. Pointweights=min(1,1.5median(norm residual)/norm residual), scale floor10m; all valid points retained. Initial x0/y0/vx/vy observationonly.

E2 optimize[x0,y0,vx,vy] in globalmetric coordinates scaled[50000m,50000m,2m/s,2m/s]. All242wrapped angularresiduals / randombearing sigma, includingepoch pairswithinvalidE1. scipy least_squares trf, analyticJacobian, soft_l1loss, f_scale1.5, no parameter bounds, max_nfev100;ftol/xtol/gtol1e-10. One frozen localrefinement from E1; not a global uniquenesscertificate. Nonconverged/nonfinitesolver invalid ->all4trutherrors infinity. No dropped/replacedrealizations.

Naive local covariance pinv(JTJ,rcond1e-12)*residualvariance withdf242-4, scaled tophysicalunits; conditionnumber ofscaledJacobian, residualRMS, SE saved. Diagnostics only; not calibrated uncertainties under nav/systematicmisspecification. No covariance-basedperformanceGate.

Independent audit: all12000savedscenes cold reconstructed, tangent slope E1 andIRLS independently rebuilt, full error andsummary recomputation; predefinedfirst/lastrealizationpercase/anchor48scenes rerunwithindependent3-point finite-differenceJacobiansameobjective/settings. Agreement tolerances1m Cartesianinitialposition and.001m/s Cartesianvelocity; failure invalidates executionaudit. No post-freeze numericalretuning. Tests5 cover noiselessoffgrid, truth API, analyticJacobian, degenerateparallelfailure, navigationinputuse.
