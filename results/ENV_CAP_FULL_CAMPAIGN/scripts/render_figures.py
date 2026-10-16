"""Publication artifacts from stored results only; no new scientific execution."""
from env_provider import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import textwrap
from matplotlib.patches import FancyBboxPatch
F=O/'figures';F.mkdir(exist_ok=True);plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'savefig.bbox':'tight'})
meta=[]
def save(fig,name,sources,level,note):
 fig.get_layout_engine().set(rect=(0,.10,1,.90));fig.text(.01,.012,textwrap.fill(note,140),fontsize=8,color='#6b3742',va='bottom');fig.savefig(F/(name+'.png'),dpi=170);fig.savefig(F/(name+'.pdf'));plt.close(fig)
 meta.append({'figure':name,'source_paths':';'.join(sources),'source_SHA256':';'.join(sha(O/p) if (O/p).exists() else sha(R/p) for p in sources),'evidence_level':level,'scope':note})
cs=list(csv.DictReader((O/'E1_CZ_ENVIRONMENTS.csv').open()));small=[c for c in cs if c['env_id'].startswith('SCS17')]
fig,ax=plt.subplots(1,3,figsize=(15,5),layout='constrained')
for env in ['SCS17_JW','SCS17_JA','SCS17_ANNUAL']:
 z,c=loadprofile(O/'inputs'/f'{env}.csv');ax[0].plot(c,z,label=env)
ax[0].invert_yaxis();ax[0].set(xlabel='Sound speed (m/s)',ylabel='Depth (m)',title='WOA23 conditional SSP / one source family');ax[0].legend(fontsize=8)
D=W/'E1/SCS17_JW_z200'
for i,(a,ns,nb,x) in enumerate(readrays(D/'R1.ray')):
 if i%8:continue
 hits=np.where((x[:-1,1]>200)&(x[1:,1]<=200))[0]
 if not len(hits):continue
 j=hits[0]
 if x[:j+2,1].max()>=3908:continue
 ax[1].plot(x[:j+2,0]/1000,x[:j+2,1],lw=.8)
ax[1].set(xlabel='Range (km)',ylabel='Depth (m)',title='First deep-turning ray segments / no bottom contact',xlim=(0,60));ax[1].invert_yaxis()
for lev in ['I0','I1','I2']:
 a=np.load(D/(lev+'_ENERGY_CURVE.npz'));ax[2].plot(a['r_km'],a['TL_dB'],label=lev,lw=1)
ax[2].set(xlabel='Range (km)',ylabel='Band-energy TL (dB)',title='JW z=200m, receiver175..225m',xlim=(35,65),ylim=(115,60));ax[2].legend()
save(fig,'F01_ENVIRONMENT_SSP_AND_CZ',['inputs/SCS17_JW.csv','inputs/SCS17_JA.csv','inputs/SCS17_ANNUAL.csv','E1_CZ_ENVIRONMENTS.csv'],'PROVISIONAL_RAY_ENERGY_DIAGNOSTIC','Global E1 admission WITHHELD: raw provisional examples; no real-ocean or coherent-pressure accuracy claim.')
fig,ax=plt.subplots(1,2,figsize=(13,6),layout='constrained')
labels=[c['case_id'] for c in cs];vals=[1 if c['CZ_class']=='CZ2_RAY_TL_STABLE' else (0 if c['CZ_class']=='CZ_NOT_PRESENT_IN_TESTED_DOMAIN' else -1) for c in cs]
colors=['#c79b42' if v==1 else ('#8d969f' if v==0 else '#a94546') for v in vals];ax[0].barh(labels,vals,color=colors);ax[0].invert_yaxis();ax[0].set(xlim=(-1.2,1.2),xticks=[-1,0,1],xticklabels=['Uncertified','No return','Provisional'],title='Main output labels are NOT final admission')
fails=json.loads((O/'E1_INDEPENDENT_COLD_AUDIT.json').read_text());a={s:sum(not x['PASS'] and x['check'].startswith(s) for x in fails['rows']) for s in labels};ax[1].barh(labels,[a[s] for s in labels],color='#a94546');ax[1].invert_yaxis();ax[1].set(xlabel='Failed independent checks (count)',title='12 FAIL /274 cold checks; formal CZ2 count=0')
save(fig,'F02_CZ_GATE_AND_NEGATIVES',['E1_CZ_ENVIRONMENTS.csv','E1_INDEPENDENT_COLD_AUDIT.json','E1_FINAL_ADJUDICATION.json'],'AUDIT_AND_NEGATIVE_EVIDENCE','Do not interpret an unclosed numeric gate as absence of physical CZ. Three constant-speed controls have no refracted return.')
fig,ax=plt.subplots(1,2,figsize=(12,4.5),layout='constrained');mod=list(csv.DictReader((O/'E2_MODE_DENSITY_AND_APERTURE.csv').open()));sh=[a for a in mod if a['env_id']=='YANG_SHALLOW_NON_CZ_H88'];L=[float(a['radial_aperture_m']) for a in sh];n=[int(a['rayleigh_unique_all_modes']) for a in sh];ax[0].semilogx(L,n,'o-');ax[0].set(xlabel='Relative radial aperture Lr (m)',ylabel='Rayleigh-isolated modes (all23)',title='Archived shallow modes / NOT a CZ transfer',ylim=(-1,25));ax[0].grid(alpha=.2)
x=np.logspace(-5,-2,200);ax[1].loglog(x,2*np.pi/x/1000);ax[1].set(xlabel='Adjacent wavenumber gap (rad/m)',ylabel='Ordinary Fourier aperture scale (km)',title='Lr=2pi / gap: mechanism-specific scale');ax[1].grid(alpha=.2)
save(fig,'F03_MODAL_SPACING_APERTURE',['E2_MODE_DENSITY_AND_APERTURE.csv','E4_ANALYTIC_DESIGN_SCALES.csv'],'HISTORICAL_MODAL_DIAGNOSTIC_AND_L0','No new CZ mode density. Radial synthetic aperture is not physical HLA length; Fourier scale is not a universal high-resolution limit.')
fig,ax=plt.subplots(figsize=(13,4),layout='constrained');ax.axis('off')
texts=[('Paper control','Shallow H=88m\nCW / known radial coordinate\ntrue mode IDs / delta=0\n12 stored-score checks match','#c2d8cf'),('CZ forward','New ray + incoherent energy\nGlobal cold gate unclosed\nComplete coherent modes\nNOT ADMITTED','#ead5ba'),('Observable extraction','Yang Fourier / MMAR AR k\nHankel complex amplitudes\nUnknown f0, delta, gains\nNOT EVALUATED','#e6c8c8'),('Depth / full support','Unknown horizontal branches\nheld-out estimator / rejection\ncoverage and wrong-side risk\nNOT EVALUATED','#e6c8c8')]
for i,(title,body,col) in enumerate(texts):
 x=.01+i*.25;ax.add_patch(FancyBboxPatch((x,.22),.225,.58,boxstyle='round,pad=.01',facecolor=col,edgecolor='#555'));ax.text(x+.112,.73,title,ha='center',weight='bold');ax.text(x+.112,.49,body,ha='center',va='center',linespacing=1.6,fontsize=9)
 if i<3:ax.annotate('',xy=(x+.25,.5),xytext=(x+.225,.5),arrowprops={'arrowstyle':'->'})
save(fig,'F04_YANG_MMAR_EVIDENCE_LAYERS',['E2_FINAL_ADJUDICATION.json','E2_METHOD_GATE.json','E3_GATE.json'],'EVIDENCE_SCOPE_MAP','E2 total controls25/26PASS; int32 mirror control failed. Stored paper matches do not repair scientific admission.')
fig,ax=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
for env in ['SCS17_JW','SCS17_JA','SCS17_ANNUAL']:
 a=[c for c in small if c['env_id']==env];z=[float(c['source_z_m']) for c in a];ax[0].plot(z,[float(c['peak_km']) for c in a],'o-',label=env);ax[1].plot(z,[float(c['numerical_P95_dB']) for c in a],'o-',label=env)
ax[0].set(xlabel='Source depth; receiver band centered there (m)',ylabel='Provisional peak range (km)',title='Coupled source/receiver-band and SSP changes');ax[1].set(xlabel='Source / receiver-band center (m)',ylabel='Largest refinement P95 difference (dB)',title='Descriptive E1 numeric sensitivity only');ax[0].legend(fontsize=8);ax[1].legend(fontsize=8)
save(fig,'F05_SENSITIVITY_BOUNDARIES',['E3_ENV_SENSITIVITY.csv','E1_FINAL_ADJUDICATION.json'],'PROVISIONAL_ENVIRONMENT_DIAGNOSTIC','No Yang/MMAR performance or SNR curves were generated. All current-CZ estimator cells remain NOT_EVALUATED.')
fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained');q=list(csv.DictReader((O/'E4_PARETO_RESOURCE_FEASIBILITY.csv').open()));fields=[('range','Range worst-scene P95 (%)'),('speed','Speed worst-scene P95 (m/s)'),('bearing','Bearing worst-scene P95 (deg)'),('heading','Heading worst-scene P95 (deg)')]
for ax,(param,label) in zip(axs.ravel(),fields):
 for quality,style in [('Q1','o-'),('Q2','s--')]:
  a=[v for v in q if v['config'].startswith(quality)];ax.plot([float(v['node_separation_km']) for v in a],[float(v[param+'_worst_case_P95']) for v in a],style,label=quality)
 ax.set(xlabel='Node separation (km)',ylabel=label,title=param);ax.legend();ax.grid(alpha=.2)
save(fig,'F06_RESOURCE_PRECISION_PARETO',['E4_PARETO_RESOURCE_FEASIBILITY.csv','results/R4_T3_W1_WIDE_BASELINE_DEMO/FOUR_PARAMETER_DEMO_TABLE.csv'],'HISTORICAL_CONDITIONAL_KINEMATIC_MC_ONLY','Historical T3 point-estimator tradeoffs, NOT new CZ acoustic feasibility or minimal equipment. All8 finite resource points nondominated.')
fig,ax=plt.subplots(figsize=(13,5),layout='constrained');ax.axis('off')
steps=[('Environment','Traceable SSP / bathymetry\nfirst refracted CZ + numeric audit\ncurrent gate: WITHHELD'),('Independent information','Source / clock / calibration\nprofiled weak directions\nL0 controls; E2 implementation FAIL'),('Extracted observables','Actual coherent receive data\nno true f0, mode IDs or trajectory\ncurrent status: NOT EVALUATED'),('Support and resources','Retain competing branches\nheld-out precision / refusal\ncurrent evidence: L0 + historicT3')]
for i,(title,body) in enumerate(steps):
 x=.02+i*.245;ax.add_patch(FancyBboxPatch((x,.35),.22,.45,boxstyle='round,pad=.01',facecolor='#dfe5eb',edgecolor='#6d7781'));ax.text(x+.11,.73,title,ha='center',fontsize=11,weight='bold');ax.text(x+.11,.54,body,ha='center',va='center',fontsize=9,linespacing=1.5)
 if i<3:ax.annotate('',xy=(x+.245,.58),xytext=(x+.22,.58),arrowprops={'arrowstyle':'->'})
ax.text(.5,.16,'Application: study the environment–observable–resource boundary; do not promise established five-state accuracy.',ha='center',fontsize=12)
save(fig,'F07_APPLICATION_ROUTE',['ENV_CAP_FINAL_DECISION.json','E3_GATE.json','E4_GATE.json'],'RESEARCH_ROUTE_WITH_EXPLICIT_UNCLOSED_GATES','No automatic next experiment. Old T2/FUNC2 STOP and R4=0 remain unchanged.')
table(O/'FIGURE_DATA_PROVENANCE.csv',meta)
print('Generated7 standalone PNG/PDF figure pairs from stored results only.')
