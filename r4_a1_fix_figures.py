"""Scientific figures from saved repair-Gate case data only."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import r4_a1_fix_continuous as fix


def figures():
    data=pd.read_csv(fix.OUT/'CASE_LEVEL_RESULTS.csv')
    nominal=data[data.sigma_deg==.1].reset_index(drop=True)
    x=np.arange(len(nominal)); labels=[s.replace('_s0.1_seed',' / ') for s in nominal.case_id]
    colors=['#176b8e' if g=='REGRESSION' else '#c55a11' for g in nominal.group]
    fig,axes=plt.subplots(3,1,figsize=(12,9),sharex=True,layout='constrained')
    axes[0].scatter(x,nominal.Jmin,c=colors)
    axes[0].axhline(.001,color='black',linestyle='--',label='Frozen matched-basin diagnostic: 0.001 dB')
    axes[0].set_yscale('log'); axes[0].set_ylabel('Exact profiled acoustic J (dB)'); axes[0].legend(fontsize=8)
    axes[1].plot(x,nominal.top1_rel_r*100,'o',label='Top1 range')
    axes[1].plot(x,nominal.survivor_worst_rel_r*100,'x',label='Survivor worst range')
    axes[1].set_ylabel('Range relative error (%)'); axes[1].legend(fontsize=8)
    axes[2].plot(x,nominal.top1_abs_psi_deg,'o',label='Top1 heading')
    axes[2].plot(x,nominal.survivor_worst_abs_psi_deg,'x',label='Survivor worst heading')
    axes[2].set_ylabel('Circular heading error (deg)'); axes[2].legend(fontsize=8)
    axes[2].set_xticks(x,labels,rotation=60,ha='right',fontsize=8)
    for ax in axes:
        ax.axvline(8.5,color='gray',alpha=.5); ax.grid(alpha=.2)
    fig.suptitle('R4-A1-FIX V2: fixed-budget continuous search at nominal 0.1 deg\n9 regression cases (blue) + 12 preregistered holdout cases (orange); no sensor-tolerance claim',fontsize=11)
    fig.savefig(fix.OUT/'CONTINUOUS_SEARCH_BASIN_AND_ENVELOPES.png',dpi=160)
    fig.savefig(fix.OUT/'CONTINUOUS_SEARCH_BASIN_AND_ENVELOPES.pdf'); plt.close(fig)
    old=pd.read_csv(fix.legacy.OUT/'CASE_LEVEL_RESULTS.csv'); legacy=old[old.seed==410001].set_index('panel_id')
    dev=nominal[nominal.group=='REGRESSION'].set_index('panel_id')
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for ax,metric,label,scale in zip(axes,['top1_rel_r','top1_abs_psi_deg'],['Range relative error (%)','Heading error (deg)'],[100,1]):
        xx=np.arange(len(dev))
        ax.bar(xx-.18,legacy.loc[dev.index,metric]*scale,width=.36,label='Frozen-grid legacy')
        ax.bar(xx+.18,dev[metric]*scale,width=.36,label='Continuous V2')
        ax.set_xticks(xx,dev.index); ax.set_ylabel(label); ax.legend(fontsize=8); ax.grid(axis='y',alpha=.2)
    fig.suptitle('Paired regression seed 410001 only; improvements do not constitute Gate PASS',fontsize=11)
    fig.savefig(fix.OUT/'PAIRED_LEGACY_CONTINUOUS_ERRORS.png',dpi=160)
    fig.savefig(fix.OUT/'PAIRED_LEGACY_CONTINUOUS_ERRORS.pdf'); plt.close(fig)


if __name__=='__main__':figures()
