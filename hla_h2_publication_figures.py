"""Publication display only: common INF category for each axis.
All scalar values and frozen experiment code remain unchanged.
Preserve the original generated plots for comparison.
"""
import shutil,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import hla_h2_core as c
LABELS=('range_relative','bearing_deg','speed_relative','heading_deg')
TITLES=('Range error (%)','Bearing error (deg)','Speed error (%)','Heading error (deg)')
SCALES=(100,1,100,1)
REFS=(10,1,10,5)
def footer(fig,stats,methods,sigmas):
    lines=['A=accepted R=abstain N=not evaluated; denominator=32 per geometry']
    for m in methods:
        for sigma in sigmas:
            group=stats[(stats.method==m)&(stats.sigma==sigma)].sort_values('geometry')
            if len(group):lines.append(m+f' {sigma:g}: '+', '.join(f"G{int(r.geometry)+1:02} {int(r.accepted)}/{int(r.abstain)}/{int(r.planned-r.evaluated)}" for _,r in group.iterrows()))
    fig.text(.02,.01,'\n'.join(lines),fontsize=7,family='monospace',va='bottom')
    fig.tight_layout(rect=[0,.14 if len(lines)<5 else .18,1,.96])
def panel(ax,curves,reference=None):
    vals=np.concatenate([np.asarray(y,float) for x,y,kw in curves])
    finite=vals[np.isfinite(vals)]
    top=max([float(v) for v in finite]+([reference] if reference is not None else [])+[.01])
    infinity=top*1.2
    hasinf=False
    for x,y,kw in curves:
        x=np.asarray(x,float);y=np.asarray(y,float);bad=~np.isfinite(y)
        plotted=y.copy();plotted[bad]=np.nan
        ax.plot(x,plotted,**kw)
        if bad.any():
            hasinf=True;ax.scatter(x[bad],np.full(bad.sum(),infinity),marker='^',color=kw['color'])
    if reference is not None:ax.axhline(reference,color='black',ls='--',label='reference')
    if hasinf:
        ax.set_ylim(-top*.04,infinity+top*.07)
        ticks=[v for v in ax.get_yticks() if 0<=v<=top*1.04]
        ax.set_yticks(ticks+[infinity]);ax.set_yticklabels([f'{v:g}' for v in ticks]+['INF'])
        ax.axhspan(top*1.08,infinity+top*.07,color='gray',alpha=.12)
        ax.text(.01,.93,'INF category (not a finite error)',transform=ax.transAxes,fontsize=8)
    else:ax.set_ylim(bottom=-top*.04)
    ax.grid(alpha=.2)
def render():
    stats=pd.read_csv(c.OUT/'METRICS_BY_GEOMETRY.csv')
    names=['PRECISION_VS_FOUR_PARAMETER_ERRORS.png','SIGNED_UNSIGNED_BOUNDARY.png','RADIAL_ZERO_POINT_BOUNDARY.png']
    originals=c.OUT/'ORIGINAL_FROZEN_REPORT_FIGURES';originals.mkdir(exist_ok=False)
    for name in names:shutil.copyfile(c.OUT/name,originals/name)
    fig,axs=plt.subplots(2,2,figsize=(13,10))
    for i,(ax,label,scale,title) in enumerate(zip(axs.flat,LABELS,SCALES,TITLES)):
        curves=[]
        for g in range(8):
            group=stats[(stats.method=='U')&(stats.geometry==g)].sort_values('sigma')
            curves.append((group.sigma,group[label+'_unconditional_P95']*scale,dict(label=f'G{g+1:02}',color=plt.cm.tab10(g),marker='o')))
        panel(ax,curves,REFS[i]);ax.set_title(title+' - exploratory unconditional P95');ax.set_xlabel('sigma_q (m/s)')
    axs.flat[0].legend(ncol=3,fontsize=8);fig.suptitle('INJECTED FEATURE ONLY; all 8 geometries; abstain included; no precision interpolation')
    footer(fig,stats,['U'],c.SIGMAS);fig.savefig(c.OUT/names[0],dpi=160);plt.close(fig)
    for methods,colors,title,name in [
        (['S','U'],['#009e73','#0072b2'],'Sign boundary; sigma_q=0.05 m/s; calibrated zero',names[1]),
        (['U','U-mismatch','UB-matched'],['#0072b2','#cc6677','#e69f00'],'Zero point: clean / +0.10 ignored / +0.10 fitted in [-0.20,+0.20]',names[2])]:
        fig,axs=plt.subplots(2,2,figsize=(13,10))
        for ax,label,scale,pt in zip(axs.flat,LABELS,SCALES,TITLES):
            curves=[]
            for method,col in zip(methods,colors):
                group=stats[(stats.method==method)&(stats.sigma==.05)].sort_values('geometry')
                curves.append((group.geometry+1,group[label+'_unconditional_median']*scale,dict(label=method,color=col,marker='o')))
            panel(ax,curves);ax.set_title(pt+' - unconditional median');ax.set_xlabel('all 8 development geometries')
        axs.flat[0].legend();fig.suptitle('INJECTED FEATURE ONLY; '+title)
        footer(fig,stats,methods,[.05]);fig.savefig(c.OUT/name,dpi=160);plt.close(fig)
    c.write_json(c.OUT/'FIGURE_DISPLAY_PROVENANCE.json',dict(scope='Display correction only: per-axis common separate INF category replaces per-curve near-zero placement; rejection denominators unchanged',scientific_values_changed=False,frozen_estimator_or_report_code_changed=False,original_figures_preserved=True,rendered_from='METRICS_BY_GEOMETRY.csv',input_sha256=c.sha(c.OUT/'METRICS_BY_GEOMETRY.csv')))
if __name__=='__main__':render()
