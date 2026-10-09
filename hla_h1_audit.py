"""Independent execution audit: verify each saved original or reprofiled depth, not the minimum of an original-depth record."""
import json,time
import numpy as np,pandas as pd
from scipy.stats import norm
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import hla_h1_core as c
def csv(name,df):df.to_csv(c.OUT/name,index=False,float_format='%.17g')
def wilson(p,n):
    z=norm.ppf(.975);d=1+z*z/n
    center=(p+z*z/(2*n))/d;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return center-half,center+half
def report():
    manifest=json.loads((c.OUT/'EXECUTION_MANIFEST.json').read_text())
    f=pd.read_csv(c.OUT/'CONFIGURATION_RESULTS.csv') if (c.OUT/'CONFIGURATION_RESULTS.csv').exists() else pd.DataFrame()
    agg=[]
    for key,group in f.groupby(['method','geometry','sigma','source']):
        method,g,sig,src=key;n=len(group);ret=group.truth_retained.mean();lo,hi=wilson(ret,n)
        r=dict(method=method,geometry=g,sigma=sig,source=src,n=n,truth_retention=ret,truth_retention_Wilson95_low=lo,truth_retention_Wilson95_high=hi,joint_truth_retention=group.truth_joint_retained.mean(),empty_rate=group['empty'].mean(),median_count=group['count'].median())
        for ax in c.AXES:
            for col in ('span','point_error','worst_error'):
                v=group[f'{col}_{ax}'];r[f'median_{col}_{ax}']=v.median();r[f'max_{col}_{ax}']=v.max()
        agg.append(r)
    a=pd.DataFrame(agg)
    csv('PRIMARY_METRICS.csv',a)
    csv('COVERAGE_AND_SET_WIDTHS.csv',f[[x for x in f.columns if x.startswith(('min_','max_','span_','tags_','range_')) or x in ['method','geometry','sigma','source','replicate','truth_retained','truth_joint_retained','empty','count']]])
    csv('FULL_SET_WORST_ERRORS.csv',f[[x for x in f.columns if 'worst' in x or x in ['method','geometry','sigma','source','replicate','count','range_branches','range_branch_gap_km']]])
    primary=[]
    for src in ('S1','S2'):
        for g in [f'G{i:02}' for i in range(1,7)]:
            panel=a[(a.geometry==g)&(a.sigma==.25)&(a.source==src)].set_index('method')
            if not set(('M1','M0','BEARING'))<=set(panel.index):continue
            m1,m0,b=panel.loc['M1'],panel.loc['M0'],panel.loc['BEARING']
            shrink=[1-m1[f'median_span_{ax}']/b[f'median_span_{ax}'] if b[f'median_span_{ax}']>0 else np.nan for ax in ('r_km','v_mps','psi_deg')]
            gains=[1-m1[f'median_span_{ax}']/m0[f'median_span_{ax}'] if m0[f'median_span_{ax}']>0 else np.nan for ax in ('r_km','v_mps','psi_deg')]
            A=bool(m1.truth_retention>=.9 and any(v>=.3 for v in shrink))
            B=bool(m1.truth_retention-m0.truth_retention>=.1 or (min(m1.truth_retention,m0.truth_retention)>=.9 and any(v>=.2 for v in gains)))
            primary.append(dict(source=src,geometry=g,truth_retention_M1=m1.truth_retention,truth_retention_M0=m0.truth_retention,A=A,B=B,**dict(zip(['range_shrink_vs_bearing','speed_shrink_vs_bearing','heading_shrink_vs_bearing'],shrink))))
    p=pd.DataFrame(primary);csv('PRIMARY_SIGNAL_BY_GEOMETRY.csv',p)
    csv('SOURCE_MODEL_TRADEOFF.csv',a)
    off=pd.read_csv(c.OUT/'OFFGRID_ESTIMATES.csv') if (c.OUT/'OFFGRID_ESTIMATES.csv').exists() else pd.DataFrame()
    if len(off):
        truths=pd.read_csv(c.OUT/'OFFGRID_TRUTH_EVALUATION_ONLY.csv');export=np.load(c.OUT/'OFFGRID_EXPORTED_CANDIDATES.npz')
        for ix,row in off.iterrows():
            true=truths.loc[int(row.panel),list(c.AXES)].to_numpy(float)
            e=c.errors(row[list(c.AXES)].to_numpy(float),true)
            candidates=export[f'{row.method}_{int(row.config)}'];accepted=candidates[candidates[:,4]<=c.CUT[row.method]]
            for j,ax in enumerate(c.AXES):
                off.loc[ix,f'point_error_{ax}']=e[j]
                off.loc[ix,f'exported_span_{ax}']=np.ptp(accepted[:,j]) if len(accepted) else np.nan
                off.loc[ix,f'exported_worst_error_{ax}']=c.errors(accepted[:,:4],true)[:,j].max() if len(accepted) else np.nan
            off.loc[ix,'point_relative_r']=e[0]/true[0];off.loc[ix,'point_relative_v']=e[2]/true[2]
        csv('OFFGRID_RESULTS.csv',off)
    else:csv('OFFGRID_RESULTS.csv',pd.DataFrame(columns=['method','status']))
    # Independent audit uses explicit double centering, not the main projection helper.
    checks=[]
    m=c.Model()
    states=np.load(c.OUT/'GRID_STATES.npy')
    def independent_sse(level,observation,method):
        e=level-observation;total=np.zeros(e.shape[:-2])
        for sl in (slice(0,61),slice(61,121)):
            q=e[...,sl];q=q-q.mean(axis=-1,keepdims=True)
            if method=='M1':q=q-q.mean(axis=-2,keepdims=True)
            total+=np.sum(q*q,axis=(-2,-1))
        return total
    samples=[]
    for g in range(6):
        path=c.OUT/f'observations/G{g+1:02}.npz'
        if not path.exists():continue
        obs=np.load(path);accepted=np.load(c.OUT/f'accepted/G{g+1:02}.npz')
        # Every configuration's analytical source/true residual independently rebuilt.
        truelevel=m.levels(c.GEOMS[g],True)[0,int((c.ZS[g]-150)/5)]
        for ci,(sig,src,rep) in enumerate(obs['configs']):
            for method in c.METHODS:
                yy=obs['levels'][ci]
                if method=='ORACLE':yy=yy-c.sources()[int(src)]
                bcost=np.square(c.wrap(c.geometry(c.GEOMS[g])[0][0]-obs['bearings'][int(rep)])/np.radians(.1)).sum()
                ts=bcost if method=='BEARING' else bcost+independent_sse(truelevel,yy,'M1' if method=='M1' else 'M0')/sig**2
                row=f[(f.geometry==f'G{g+1:02}')&(f.config==ci)&(f.method==method)].iloc[0]
                if method!='BEARING':checks.append(dict(check='all_true_scores',error=abs(ts-row.truth_joint_score),PASS=abs(ts-row.truth_joint_score)<1e-6))
        # Fixed source/sigma/repeat entries, all retained candidate states directly rechecked.
        chosen=[0,32,64,96,128,160,192,224,256,288,320,352]
        for ci in chosen:
            sig,src,rep=obs['configs'][ci]
            for method in c.METHODS:
                ids=accepted[f'{method}_{ci}_ids'];saved=accepted[f'{method}_{ci}_scores']
                for start in range(0,len(ids),16):
                    ii=ids[start:start+16];bc=np.square(c.wrap(c.geometry(states[ii])[0]-obs['bearings'][int(rep)])/np.radians(.1)).sum(axis=1)
                    if method=='BEARING':sc=bc
                    else:
                        yy=obs['levels'][ci]-(c.sources()[int(src)] if method=='ORACLE' else 0)
                        sc=bc+independent_sse(m.levels(states[ii],True),yy,'M1' if method=='M1' else 'M0').min(axis=1)/sig**2
                    # Empirical interpolation check remains separate from acceptance near threshold.
                    for v,w in zip(sc,saved[start:start+16]):
                        err=abs(v-w);checks.append(dict(check='direct_accepted_score',error=float(err),PASS=bool(v<=c.CUT[method]+1e-8)))
    if len(off):
        # Independent Cartesian trajectory geometry already preflight-certified.
        data=np.load(c.OUT/'observations/OFFGRID.npz');export=np.load(c.OUT/'OFFGRID_EXPORTED_CANDIDATES.npz')
        for key in export.files:
            method,ci=key.rsplit('_',1);ci=int(ci);arr=export[key]
            for start in range(0,len(arr),16):
                s=arr[start:start+16,:4];bc=np.square(c.wrap(c.geometry(s)[0]-data['bearings'][ci])/np.radians(.1)).sum(axis=1)
                if method=='BEARING':score=bc
                else:
                    costs=independent_sse(m.levels(s,True),data['levels'][ci],method)/.25**2
                    labels=arr[start:start+16,5]
                    zi=np.rint((labels-150)/5).astype(int)
                    score=bc+costs[np.arange(len(s)),zi]
                for v,w in zip(score,arr[start:start+16,4]):
                    checks.append(dict(check='all_offgrid_terminals',error=float(abs(v-w)),PASS=bool(abs(v-w)<1e-6)))
    c.write_json(c.OUT/'VALIDATION.json',dict(PASS=all(v['PASS'] for v in checks),checks=len(checks),FAIL=sum(not v['PASS'] for v in checks),max_score_difference=max((v['error'] for v in checks),default=0),scope='All main true residuals, fixed registered configurations all accepted states, all offgrid terminal states; preflight independent modal/Cartesian/QR controls. No global interpolation error certificate.',counts=pd.Series([v['check'] for v in checks]).value_counts().to_dict()))
    csv('COLD_REBUILD_CHECKS.csv',pd.DataFrame(checks))
    valid=all(v['PASS'] for v in checks)
    counts={src:dict(A=int(p[p.source==src].A.sum()),B=int(p[p.source==src].B.sum()),A_and_B=int((p[p.source==src].A&p[p.source==src].B).sum())) for src in ('S1','S2')} if len(p) else {}
    gain=bool(counts and all(v['A']>=4 and v['A_and_B']>=4 for v in counts.values()))
    status='HLA_H1_SOURCE_PROFILED_HORIZONTAL_GAIN_SUPPORTED_CONDITIONALLY' if gain else 'HLA_H1_SOURCE_INVARIANCE_VERIFIED_BUT_CONTRACTION_LIMITED'
    if manifest['status']!='COMPLETE':status='HLA_H1_SEARCH_OR_NUMERICAL_LIMITED'
    if not valid:status='HLA_H1_IMPLEMENTATION_OR_INPUT_INVALID'
    mismatch=a[(a.method=='M1')&(a.source=='S3')]
    limits=[]
    if len(mismatch) and (mismatch.truth_retention<.9).any():limits.append('HLA_H1_GAIN_LIMITED_BY_NONCOMMON_SOURCE_VARIATION')
    decision=dict(stage='HLA-H1',classification=status,parallel_limits=limits,signal=counts,execution=manifest['status'],main_rows=len(f),offgrid_rows=len(off),R4_percent=0,scope='EXPLORATORY_FEATURE_LEVEL_NOMINAL_MATCHED_SINGLE_HLA',no_continuous_coverage=True,depth_success='NOT_EVALUATED',real_signal_extraction='NOT_EVALUATED',source_accuracy='CONDITION_COMMON_SOURCE_ONLY',STOP=True)
    c.write_json(c.OUT/'DECISION.json',decision)
    # Three required standalone scientific plots.
    colors={'BEARING':'#555555','M0':'#cc6677','M1':'#0072b2','ORACLE':'#009e73'}
    fig,axs=plt.subplots(2,2,figsize=(11,7))
    sub=a[(a.sigma==.25)&(a.source=='S2')]
    for ax,tag in zip(axs.flat,c.AXES):
        for method in c.METHODS:
            v=sub[sub.method==method].sort_values('geometry')
            ax.plot(v.geometry,v[f'median_span_{tag}'],marker='o',label=method,color=colors[method])
        ax.set_title(tag+' median finite-set span');ax.grid(alpha=.2)
    axs.flat[0].legend();fig.suptitle('S2 / 0.25 dB: empty sets excluded from width; zero means one grid label');fig.tight_layout();fig.savefig(c.OUT/'FOUR_PARAMETER_SET_WIDTHS.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(12,4))
    for ax,tag in zip(axs,('r_km','v_mps','psi_deg')):
        for method in c.METHODS:
            v=sub[sub.method==method];ax.scatter(v[f'median_span_{tag}'],v.truth_retention,label=method,color=colors[method])
        ax.set_xlabel(tag+' median span');ax.set_ylabel('Empirical true horizontal retention');ax.set_ylim(-.05,1.05);ax.grid(alpha=.2)
    axs[0].legend();fig.tight_layout();fig.savefig(c.OUT/'COVERAGE_WIDTH_TRADEOFF.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4))
    for method in c.METHODS:
        v=a[(a.sigma==.25)&(a.method==method)].groupby('source').truth_retention.mean()
        ax.plot(v.index,v.values,marker='o',label=method,color=colors[method])
    ax.set_ylabel('Mean of six empirical scene retentions');ax.set_ylim(-.05,1.05);ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(c.OUT/'SOURCE_CONDITION_METHOD_GAIN.png',dpi=160);plt.close(fig)
    c.write_json(c.OUT/'OUTPUT_MANIFEST.json',{'sha256':{str(p.relative_to(c.OUT)):c.sha(p) for p in c.OUT.rglob('*') if p.is_file() and p.name not in ['OUTPUT_MANIFEST.json','GPT_SYNC.md']}})
    print(json.dumps(decision,indent=2))
if __name__=='__main__':
    with threadpool_limits(limits=4):report()
