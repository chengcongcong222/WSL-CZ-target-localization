from pathlib import Path
import csv,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parent
def rows(p):
 with p.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
a=rows(OUT/'PREC3_REAL_K_INDEPENDENT_VALIDATION.csv');r=rows(OUT/'PREC3_INDEPENDENT_COMPLEX_AND_PHI.csv')
fig,ax=plt.subplots(1,2,figsize=(12,4.4),layout='constrained')
for f,color in [(201,'#2563eb'),(235,'#0d9488'),(283,'#9333ea')]:
 rr=[x for x in a if int(x['freq_hz'])==f];xx=[i for i,v in enumerate(a) if int(v['freq_hz'])==f]
 ax[0].scatter(xx,[float(v['GALERKIN_vs_PREC2']) for v in rr],color=color,label=f'{f}Hz: Galerkin vs FD')
ax[0].plot(range(15),[float(v['GALERKIN_R4_gap']) for v in a],'x',color='#334155',label='Galerkin R4 gap')
ax[0].axhline(1e-4/60000,color='#dc2626',ls='--',label='Registered single-mode budget')
ax[0].set_yscale('symlog',linthresh=1e-14);ax[0].set_xlabel('Frozen mode index (see CSV)');ax[0].set_ylabel('Wavenumber difference (1/m)');ax[0].set_title('A: finite-subset real-k consistency');ax[0].legend(fontsize=7)
for N,color in [(80000,'#2563eb'),(160000,'#0d9488')]:
 rr=[v for v in r if int(v['N'])==N];ax[1].plot([int(v['mode']) for v in rr],[100*float(v['modal_function_FE_vs_FD_relative']) for v in rr],'o-',color=color,label=f'N={N}')
ax[1].axhline(.2,color='#dc2626',ls='--',label='Frozen modal-function threshold')
ax[1].set_yscale('log');ax[1].set_xlabel('201Hz main matched mode');ax[1].set_ylabel('Relative function difference (%)');ax[1].set_title('B: Galerkin / FD functions on251SSP nodes');ax[1].legend(fontsize=8)
for x in ax:x.grid(True,alpha=.18)
fig.savefig(OUT/'PREC3_NUMERICAL_DIAGNOSTICS.png',dpi=200);fig.savefig(OUT/'PREC3_NUMERICAL_DIAGNOSTICS.pdf')
plt.close(fig);print('Figure generated exclusively from sealed completed results; C not plotted')
