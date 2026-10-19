# Independent algebra and interval controls

For y_m=s g_m exp(i k x_m cos(alpha-beta))+n_m, multiplying s by a complex scalar is profiled by least squares. A common phase/amplitude cannot create range or source-depth sensitivity in a single plane-wave snapshot: d(normalized spatial vector)/dr = d/dz =0. cos(alpha)=cos(-alpha), so no true side can enter extraction.

Known-angle calibration compares channel/reference ratios at +/-theta: their phase difference is2k(x_m-x_0)sin(theta)sin(beta). Fixed channel response cancels. This is identifiable only with the extra known-angle calibration transmissions; target data alone do not provide that reference.

For log channel powers, D=[-1,I] cancels common source level. D D^T=I+11^T: reference noise induces off-diagonal covariance. Free per-time differential gain spans the entire contrast space; projected Fisher information is zero. G1 and G2 remove fixed and linear components; projection cannot increase information. Toy source/depth labels do not enter the extractor.

For a state box, the target-minus-reported-node position interval is affine in(x0,y0,vx,vy), expanded by navigation/time bounds. All x values are positive in the frozen physical domain. atan2 extrema lie on rectangle corners. Taking the complete heading/calibration interval and cosine image yields an enclosing axial response. Boxes are rejected only if this entire image is disjoint from every candidate interval at an epoch. Missing observations reject nothing; unresolved boxes remain. This guarantees algorithmic outer inclusion conditional on input sets, not the statistical validity of all sets.

Development tests and first failed fixture are preserved. Attempt1 expected a detection during a deliberately absent target; the positive-control fixture was corrected before freeze, without modifying any formal seed or observation.
