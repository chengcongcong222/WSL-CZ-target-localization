# E1-G0 Jacobian and covariance lock
Coordinates follow existing code: +x initial nominal look direction, atan2(y,x). x=[p0x,p0y,vx,vy]. S=diag(50000,50000,2,2). Target p=p0+t*v. Main speed2m/s, turn15deg after600s; at exactly600s use pre-turn velocity. AUX adds fixed [0,mirror*5000] and shares translation. Deployment beta=0 nominal control; not 500 saved random placements.

e=(p-node)/r; w=v-vnode; q=e.w; q_p=(w-e*q)/r; q_x=[q_p,t*q_p+e].
bearing_x=[-e_y/r,e_x/r,t*(-e_y/r),t*e_x/r].
u=t-r/c-delta_t; u_x=-[e,t*e]/c.
mu=exp(a_l+(d+dl_l)*u/1200+rho_j+k_j*t/1200)*(1-q/c)+B_j+K_j*t/1200.
mu_x=(mu-B-K*t/1200)*((d+dl)/1200*u_x-q_x/(c*(1-q/c))).
eta derivatives: a_l,rho_j -> Doppler component; d,dl_l -> component*u/1200; k_j -> component*t/1200; B_j ->1; K_j ->t/1200. No unknown emitted f supplied as calibration.

Navigation derivatives at each node/epoch are minus position derivatives. R_local=diag(sigma_b²,sigma_f²,...)+25² J_nav J_nav.T. Same-epoch same-node bearing and frequency share nav errors; cross-lines also correlated. White instrument errors are ideal independent G0 inputs, not extracted same-waveform features. Static uniform common/half-differential bearing offsets integrated with variance bound²/3 as Gaussian moment-matched rank2 covariance. Not a uniform-distribution exact Fisher bound. No navigation velocity errors or target process prior added.

Whiten block covariance by Cholesky then rank2 inverse-square-root for static bearing covariance. Profile with orthonormal nuisance left basis U: R=(I-UU.T)A_x, I_eff=R.T R. Independently build Schur using full-rank SVD-selected nuisance coordinate basis. N5 appends receiver calibration prior rows before profiling. Target inversion via scaled R SVD; speed gradient S*[0,0,vx/v,vy/v]. Null-space overlap -> UNBOUNDED, never finite pseudoinverse-zero. Weak vector is scaled four-state null/weak combination, physical loading=S*vector.
