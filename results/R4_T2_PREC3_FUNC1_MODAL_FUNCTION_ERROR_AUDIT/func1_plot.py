"""Post-execution figures, using only already sealed FUNC1 outputs."""
from pathlib import Path
import csv,json,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parent
def rows(n):
 with (OUT/n).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
def save(fig,name):
 fig.savefig(OUT/(name+'.png'),dpi=200);fig.savefig(OUT/(name+'.pdf'));plt.close(fig)
samples=np.load(OUT/'FUNC1_DENSE_SAMPLES.npz',allow_pickle=False);z=samples['z_m'];phase=rows('FUNC1_NORMALIZATION_AND_PHASE_AUDIT.csv');cross=rows('FUNC1_FE_FD_CROSS_GRID.csv')
fig,axes=plt.subplots(1,3,figsize=(15,4),layout='constrained')
for ax,mode in zip(axes[:2],[1,481]):
 for N,color in [(80000,'#2563eb'),(160000,'#0f766e')]:
  f=samples[f'FE_m{mode}_N{N}'];g=samples[f'FD_m{mode}_N{N}'];sig=int(next(r for r in phase if r['mode']==str(mode) and r['N']==str(N) and r['grid']=='S3')['sign']);ax.plot(z/1000,abs(f-sig*g),color=color,label=f'N={N}',lw=.6)
 ax.set_title(f'Mode{mode}: paired FE-FD');ax.set_xlabel('Depth (km)');ax.set_ylabel('Absolute modal-function difference');ax.legend(fontsize=8)
sig=int(next(r for r in cross if r['comparison']=='SAME_METHOD_CROSS_GRID' and r['method']=='FE' and r['mode']=='961')['sign'])
axes[2].plot(z/1000,abs(samples['FE_m961_N80000']-sig*samples['FE_m961_N160000']),color='#9333ea',lw=.6);axes[2].set_title('Mode961: FE own cross-grid only');axes[2].set_xlabel('Depth (km)');axes[2].set_ylabel('Absolute modal-function difference')
for ax in axes:ax.grid(alpha=.2)
save(fig,'FUNC1_MODAL_FUNCTION_DIFFERENCES')
d=rows('FUNC1_DENSE_SAMPLING_COMPARISON.csv');labels=['SSP_UNWEIGHTED','SSP_TRAPZ','S1','S2','S3'];pretty=['20m discrete','20m trapz','1m','0.5m','0.25m']
fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained')
for ax,mode in zip(axes,[1,481]):
 for N,color in [(80000,'#2563eb'),(160000,'#0f766e')]:
  rr={r['grid']:float(r['e_sign'])*100 for r in d if r['mode']==str(mode) and r['N']==str(N)}
  ax.plot(range(5),[rr[n] for n in labels],'o-',color=color,label=f'N={N}')
 ax.set_xticks(range(5),pretty,rotation=20);ax.set_title(f'Mode{mode}: evaluation grid');ax.set_ylabel('Sign-only discrepancy (%)');ax.grid(alpha=.2);ax.legend(fontsize=8)
 if mode==481:ax.axhline(.2,color='#dc2626',ls='--',label='Historical 0.2% (unchanged)');ax.legend(fontsize=8)
save(fig,'FUNC1_SAMPLING_COMPARISON')
b=rows('FUNC1_DEPTH_RESOLVED_ERROR.csv')
fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained')
for ax,mode in zip(axes,[1,481]):
 for N,color,shift in [(80000,'#2563eb',-.035),(160000,'#0f766e',.035)]:
  rr=[r for r in b if r['mode']==str(mode) and r['N']==str(N)];ax.bar([(float(r['z_low_m'])+100)/1000+shift for r in rr],[100*float(r['fraction_of_total_squared_error']) for r in rr],width=.07,color=color,label=f'N={N}')
 ax.set_title(f'Mode{mode}:200m error bands');ax.set_xlabel('Depth (km)');ax.set_ylabel('Share of total squared difference (%)');ax.legend(fontsize=8);ax.grid(axis='y',alpha=.2)
save(fig,'FUNC1_DEPTH_ERROR_DISTRIBUTION')
c=rows('FUNC1_LOCAL_MODAL_COUPLING.csv');resources=['A0_0','A_P4M','B1_V4','B2_DUP'];colors=['#334155','#2563eb','#0f766e','#9333ea'];probes=['P01','P02','P03','P04','P05']
fig,axes=plt.subplots(2,2,figsize=(12,7),layout='constrained')
for column,N in enumerate([80000,160000]):
 for res,color in zip(resources,colors):
  rr={r['probe']:r for r in c if r['mode']=='481' and r['N']==str(N) and r['resource']==res}
  axes[0,column].plot(range(5),[float(rr[p]['raw_relative_difference'])*100 for p in probes],'o-',color=color,label=res)
  axes[1,column].plot(range(5),[float(rr[p]['array_shape_difference'])*100 for p in probes],'o-',color=color,label=res)
 axes[0,column].set_title(f'Mode481, N={N}: raw Q');axes[1,column].set_title(f'Mode481, N={N}: Q shape diagnostic')
 for ax in axes[:,column]:ax.set_xticks(range(5),probes);ax.grid(alpha=.2);ax.legend(fontsize=7)
axes[0,0].set_ylabel('Uncalibrated coupling discrepancy (%)');axes[1,0].set_ylabel('Normalized single-mode Q shape (%)')
fig.suptitle('Q=phi(source)phi(receiver); no radial kernel or full-mode field. Mode1: all 40 couplings flagged weak.')
save(fig,'FUNC1_LOCAL_COUPLING')
print('4PNG+4PDF figures from completed saved data; no new solver or scientific analysis')
