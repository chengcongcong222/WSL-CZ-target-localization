# Random and systematic measurement error budget

Anchor sigma_random is per-node independent zero-mean Gaussian directed-bearing RANDOM error only: A0.05deg, B0.075deg. It is not total effective system RMS. If component independence and zero means are substantiated, sigma_random^2=sigma_DOA^2+sigma_heading_random^2+sigma_cal_random^2. Otherwise covariance terms must be included; no numerical component allocation is made here.

Fixed/systematic residual b1,b2 remains separate, using common=(b1+b2)/2 and differential=(b2-b1)/2. Differential is a half-difference; a0.05deg axis cap corresponds to0.10deg full inter-node difference. Source components need not physically be independent merely because parametrized this way. Do not square a fixed bias into an independent Gaussian variance and claim that it was tested.

Node navigation error is independent2D Gaussian per-axis1-sigma in metres for each node. Real correlated navigation, attitude bias, lever arm, time skew or target association remain unverified. Deployment beta is known actual geometry, not node-position estimation error or bearing bias.

Each family is single factor. Bias axis tolerance requires the other bias coordinate zero; combined bias map is separate. Axis maxima may not be simultaneously satisfied. Navigation/deployment limits are each with other nonideal errors zero. Joint engineering point requires separately authorized Gate2B pre-run freeze; no joint safety statement here. T5/T10 scopes are discrete11-range finite-MC only, no continuous tolerance theorem. No acoustic/depth/SSP/TDOA/tracker.
