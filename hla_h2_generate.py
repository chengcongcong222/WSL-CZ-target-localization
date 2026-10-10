"""Generator only. Hidden truths/innovations never returned to an estimator."""
import numpy as np,pandas as pd
import hla_h2_core as c
SOURCE_PANEL=c.ROOT/'results/R4_SINGLE_HLA_SOURCE_PROFILED_HORIZONTAL_PILOT/OFFGRID_TRUTH_EVALUATION_ONLY.csv'
def generate():
    if (c.OUT/'OBSERVATIONS.npz').exists():raise RuntimeError('Observations already generated: no rerun')
    panel=pd.read_csv(SOURCE_PANEL);h=panel[list(c.AXES)].to_numpy(float)
    beta,q,rho=c.geometry(h)
    checks=[]
    for g,state in enumerate(h):
        th,psi=np.radians(state[[1,3]]);vt=state[2]*np.array([np.cos(psi),np.sin(psi)])
        p0=1000*state[0]*np.array([np.cos(th),np.sin(th)])
        for j in range(6):
            times=np.array([200*j,200*(j+1)],float)
            vo=2*np.array([1.,0.]) if j<3 else 2*np.array([np.cos(np.radians(15)),np.sin(np.radians(15))])
            rel=p0+times[:,None]*vt-c.platform(times)
            rate=rel@(vt-vo)/np.linalg.norm(rel,axis=1)
            checks.append(dict(geometry=g,window=j+1,endpoint_instantaneous_rate_start=rate[0],endpoint_instantaneous_rate_end=rate[1],within_window_sign_change=bool(rate[0]*rate[1]<0),near_zero_net_rate=bool(abs(q[g,j])<.02),net_rate=q[g,j],retained=True))
    pd.DataFrame(checks).to_csv(c.OUT/'WINDOW_SIGN_AND_NEAR_ZERO_DIAGNOSTICS.csv',index=False,float_format='%.17g')
    eb=np.empty((8,32,121));eq=np.empty((8,32,6))
    for g in range(8):
        for rep in range(32):
            eb[g,rep]=np.random.default_rng(np.random.SeedSequence([2026101002,g,rep,0])).standard_normal(121)
            eq[g,rep]=np.random.default_rng(np.random.SeedSequence([2026101002,g,rep,1])).standard_normal(6)
    observed_b=c.wrap(beta[:,None,:]+c.SIG_B*eb)
    signed=q[:,None,None,:]+np.array(c.SIGMAS)[None,None,:,None]*eq[:,:,None,:]
    unsigned=abs(signed)
    biased=abs(q[:,None,:]+.10+.05*eq)
    np.savez_compressed(c.OUT/'OBSERVATIONS.npz',beta=observed_b,unsigned=unsigned,unsigned_biased=biased,sigmas=np.array(c.SIGMAS))
    np.savez_compressed(c.OUT/'SIGNED_CONTROL_OBSERVATIONS.npz',beta=observed_b,signed=signed,sigmas=np.array(c.SIGMAS))
    np.savez_compressed(c.OUT/'GENERATOR_TRUTH_AND_INNOVATIONS.npz',horizontal_truth=h,beta_truth=beta,qbar_truth=q,range_truth=rho,depth_metadata=panel.z_label.to_numpy(),bias_truth=np.array(.10),bearing_innovation=eb,radial_innovation=eq)
    c.write_json(c.OUT/'OBSERVATION_GENERATION_MANIFEST.json',dict(source_panel_sha256=c.sha(SOURCE_PANEL),seed=2026101002,PRNG='NumPy default_rng PCG64 / SeedSequence [seed,geometry_zero_based,replicate_zero_based,stream_zero_based]',numpy_version=np.__version__,base_innovations=256,shared_across_methods_and_scales=True,sha256={p:c.sha(c.OUT/p) for p in ['OBSERVATIONS.npz','SIGNED_CONTROL_OBSERVATIONS.npz','GENERATOR_TRUTH_AND_INNOVATIONS.npz']},signed_control_is_separate=True,source_depth='METADATA_ONLY_NOT_USED',INJECTED_FEATURE_ONLY=True))
if __name__=='__main__':generate()
