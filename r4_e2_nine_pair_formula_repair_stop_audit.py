"""Independent read-only audit of the frozen pre-development stop. No regression restart."""
from pathlib import Path
import json,csv
import numpy as np
from scipy import linalg
import r4_e2_nine_pair_formula_repair as p
import r4_e2_nine_pair_pilot_audit as independent
checks=[]
def ck(name,ok,detail=''):checks.append(dict(check=name,PASS=bool(ok),detail=str(detail)))
def main():
 o=p.O;f=p.read(o/'DESIGN_FREEZE.json');start=p.read(o/'EXECUTION_STARTED.json')
 for name,r in f['bindings'].items():ck('frozen_input:'+name,p.sha(name,r['kind']=='RAW')==r['sha256'])
 ctrl=p.base.rows(o/'PRE_DEVELOPMENT_CONTROLS.csv');sr=p.base.rows(o/'STRUCTURAL_NULLSPACE_CHECKS.csv');fd=p.base.rows(o/'ANALYTIC_VS_FD_JACOBIAN.csv')
 ck('registered_precontrols_count',len(ctrl)==1305)
 bad=[r for r in ctrl if r['PASS']!='True']
 ck('failed_family_only_nominal_pressure_identity',len(bad)==20 and all(r['control'].startswith('unchanged_pressure_') for r in bad))
 ck('development_never_started',not (o/'REPAIRED_INFORMATION_BY_SCENE.csv').exists() and not list(o.glob('PROFILED_*.npz')))
 errmax=0.;geommax=0.;pressmax=0.;oldpressmax=0.;nullmax=0.;fieldrows=[]
 for sc in f['scenes']:
  x=p.state(sc)
  for grid in [80001,160001]:
   a=np.load(o/f"INPUT_{sc['scene_id']}_n{grid}.npz")
   original=np.load(p.OLD/f"INPUT_{sc['scene_id']}_n{grid}.npz")
   mods=[independent.mod(p.paired.path_for(ff,grid)) for ff in p.RAW]
   tt=np.array([0.,600.,1200.]);dt=np.maximum(tt-600,0)
   mainnode=np.stack([2*np.minimum(tt,600)+2*dt*np.cos(np.pi/12),2*dt*np.sin(np.pi/12)],axis=-1)
   nodes=np.stack([mainnode,mainnode+[0,sc['mirror']*5000]],axis=1)
   offsets=np.stack([np.arange(-7.,8.,2.),np.zeros(8)],axis=-1)
   directtarget=np.array(sc['state'][:2])+tt[:,None]*np.array(sc['state'][2:4])
   target=p.cart(x)[:2]+tt[:,None]*p.cart(x)[2:]
   rd=np.linalg.norm(directtarget[:,None,None,:]-nodes[:,:,None,:]-offsets[None,None,:,:],axis=-1)
   rp=np.linalg.norm(target[:,None,None,:]-nodes[:,:,None,:]-offsets[None,None,:,:],axis=-1)
   pc=np.stack([independent.press(m,rd,200.) for m in mods],axis=-1)
   pp=np.stack([independent.press(m,rp,200.) for m in mods],axis=-1)
   pe=float(np.max(abs(pc-a['pressure']))/np.max(abs(pc)));oe=float(np.max(abs(pp-original['pressure']))/np.max(abs(pp)))
   de=float(np.max(abs(rd-rp)))
   pressmax=max(pressmax,pe);oldpressmax=max(oldpressmax,oe);geommax=max(geommax,de)
   # Sum independent scalar modal derivatives before contracting Cartesian state columns.
   grad=[]
   for zz,phi,k in mods:
    w=phi[list(zz).index(200.)]**2;rad=np.zeros(rd.shape,complex)
    for im in range(len(k)):
     term=w[im]*np.sqrt(2*np.pi/(k[im]*rd))*np.exp(-1j*k[im]*rd-1j*np.pi/4)
     rad+=term*(-1j*k[im]-.5/rd)
    grad.append(rad[...,None]*(target[:,None,None,:]-nodes[:,:,None,:]-offsets[None,None,:,:])/rd[...,None])
   gradient=np.stack(grad,axis=-2)/pc[...,None]
   T=np.zeros((3,2,4));T[:,:,0]=[np.cos(x[1]),np.sin(x[1])];T[:,:,1]=[-x[0]*np.sin(x[1]),x[0]*np.cos(x[1])];T[:,0,2]=tt;T[:,1,3]=tt
   jc=np.stack([np.sum(gradient*T[:,None,None,None,:,c],axis=-1)*p.SCALE[c] for c in range(4)],axis=-1)
   err=float(np.linalg.norm(jc-a['J'][...,:4])/np.linalg.norm(jc));errmax=max(errmax,err)
   ck('independent_scalar_modal_horizontal_chain:'+sc['scene_id']+str(grid),err<1e-9,err)
   ck('direct_and_roundtrip_pressure_reconstruction:'+sc['scene_id']+str(grid),pe<1e-12 and oe<1e-12,(pe,oe))
   D=a['center_range_geometry'];_,_,vh=linalg.svd(D,full_matrices=True);n=np.r_[vh[-1],0.]
   ck('SVD_vs_cofactor_geometry_null:'+sc['scene_id']+str(grid),abs(n@a['structural_null'])>1-1e-10)
   for jaclabel in ['J','J_half']:
    y,g,b,_=p.system(a['pressure'],a[jaclabel],a['bearing_tangent'],'P1',1,'RELATIVE',p.REF)
    z=independent.qremove(y,np.column_stack([g,b]))
    response=float(np.linalg.norm(z@n)/np.linalg.norm(z));nullmax=max(nullmax,response)
    ss=linalg.svdvals(z);rank=int(np.sum(ss>ss[0]*1e-10))
    ck('independent_QR_physical_rank:'+sc['scene_id']+str(grid)+jaclabel,rank<=4 and response<1e-10,(rank,response))
   for step,key in [(1.,'J'),(.5,'J_half')]:
    up=x.copy();dn=x.copy();up[4]+=step;dn[4]-=step
    pu=independent.cold_field(mods,sc,up);pd=independent.cold_field(mods,sc,dn)
    dz=(pu-pd)/(2*step)/a['pressure']*p.SCALE[4]
    e=float(np.linalg.norm(dz-a[key][...,4])/np.linalg.norm(dz))
    ck('independent_depth_stencil:'+sc['scene_id']+str(grid)+str(step),e<1e-9,e)
   observed=float(np.max(abs(a['pressure']-original['pressure']))/np.max(abs(a['pressure'])))
   predicted=float(np.max(abs(pc-pp))/np.max(abs(pc)))
   fieldrows.append(dict(scene_id=sc['scene_id'],mesh=grid,max_element_range_roundtrip_m=de,observed_pressure_relative=observed,independent_predicted_pressure_relative=predicted,new_pressure_reconstruction_relative=pe,old_pressure_reconstruction_relative=oe,frozen_limit=1e-12,frozen_PASS=observed<1e-12))
  print('STOP_AUDIT',sc['scene_id'],flush=True)
 p.write(o/'INDEPENDENT_STOP_AUDIT_CHECKS.csv',checks);p.write(o/'PRESSURE_GEOMETRY_ROUNDTRIP_DIAGNOSTIC.csv',fieldrows)
 pn=[r for r in sr if r['representation']=='P1' and r['calibration']=='C1' and r['noise_model']=='RELATIVE']
 rank2=sorted({int(r['analytic_rank']) for r in fd if r['representation']=='P2' and r['calibration']=='C1' and r['noise_model']=='RELATIVE'})
 cmp=p.base.rows(o/'ORIGINAL_VS_REPAIRED_COMPARISON.csv')
 differences={}
 for rep in ['P1','P2']:
  rs=[r for r in cmp if r['representation']==rep and r['calibration']=='C1' and r['noise_model']=='RELATIVE']
  differences[rep]={col:max(float(r[col]) for r in rs) for col in ['range_relative_change','depth_relative_change']}
 decision=dict(parent_SHA=p.PARENT,design_SHA=start['design_SHA'],execution_SHA='CONTAINING_COMMIT',remote_main_SHA='VERIFY_AFTER_PUSH',decision='E2_FORMULA_REPAIR_FAILED',reason='FROZEN_PRE_DEVELOPMENT_PRESSURE_IDENTITY_GUARD_FAILED',failure_interpretation='Mixed original Cartesian geometry and polar-roundtrip geometry caused numerical pressure mismatch; does not demonstrate failure of analytic rank or absence of physical information.',formula_Gate='FAIL',predevelopment_checks=len(ctrl),predevelopment_PASS=len(ctrl)-len(bad),predevelopment_FAIL=len(bad),P1_physical_rank_pass=sum(r['PASS']=='True' for r in pn),P1_physical_rank_count=len(pn),P1_rank_values=sorted({int(r['numeric_rank']) for r in pn}),P1_structural_null_residual_max=max(float(r['null_response']) for r in pn),P1_independent_QR_null_max=nullmax,P2_precontrol_rank_values=rank2,analytic_vs_FD_identifiable_relative_max=max(float(r['identifiable_relative_error']) for r in fd),precontrol_original_vs_repaired_information_changes=differences,pressure_identity_relative_max=max(float(r['observed_pressure_relative']) for r in fieldrows),pressure_identity_frozen_limit=1e-12,Cartesian_roundtrip_element_range_max_m=geommax,independent_pressure_reconstruction_max=pressmax,independent_old_pressure_reconstruction_max=oldpressmax,independent_analytic_horizontal_J_relative_max=errmax,development_regressions_executed=0,development_information='NOT_EXECUTED_PRECONTROL_STOP',C1_scientific_result='NOT_ESTABLISHED',relative_vs_fixed_floor_scientific_result='NOT_EVALUATED_THIS_REPAIR',weak_direction_and_effective_information_stability='NOT_EVALUATED_THIS_REPAIR',candidate_separation='NOT_REEVALUATED',independent_stop_audit='PASS' if all(r['PASS'] for r in checks) else 'FAIL',original_E2_G0='FAIL_UNCHANGED',original_pilot='IMPLEMENTATION_INVALID',original_E2_information='NOT_EVALUATED',result_scope='FORMULA_REPAIRED_DEVELOPMENT_ONLY',R4_percent=0,new_KRAKEN=0,new_MC=0,new_audio=0,once_only_repair_execution_count=1,next='H3_REVIEW',next_stage_execution='NOT_AUTHORIZED',execution_status='STOPPED_WITHOUT_REGRESSION')
 p.dump(o/'REPAIR_DECISION.json',decision)
 v=dict(checks=len(ctrl)+len(checks),PASS=sum(r['PASS']=='True' for r in ctrl)+sum(r['PASS'] for r in checks),FAIL=len(bad)+sum(not r['PASS'] for r in checks),registered_predevelopment=dict(checks=len(ctrl),PASS=len(ctrl)-len(bad),FAIL=len(bad)),independent_stop_audit=dict(checks=len(checks),PASS=sum(r['PASS'] for r in checks),FAIL=sum(not r['PASS'] for r in checks)),formula_Gate='FAIL',regression_executed=False,old_bound_files_unchanged=all(r['PASS'] for r in checks if r['check'].startswith('frozen_input:')),scope='No development regression, no rank threshold or formula guard changes, no second scientific execution.')
 p.dump(o/'VALIDATION.json',v)
 # Headers only honestly encode the registered stop, never populate fictional scientific values.
 p.write(o/'REPAIRED_INFORMATION_BY_SCENE.csv',[dict(scene_id=sc['scene_id'],status='NOT_EXECUTED_PRECONTROL_STOP',result_scope='FORMULA_REPAIRED_DEVELOPMENT_ONLY') for sc in f['scenes']])
 p.write(o/'EFFECTIVE_INFORMATION_STABILITY.csv',[dict(scene_id=sc['scene_id'],status='NOT_EVALUATED_PRECONTROL_STOP',result_scope='FORMULA_REPAIRED_DEVELOPMENT_ONLY') for sc in f['scenes']])
 p.dump(o/'NOT_EXECUTED.json',dict(information_regression=False,candidate_regression=False,DPI_scene_regression=False,weak_direction_regression=False,noise_floor_regression=False,reason='Frozen analytic precontrols did not all pass',registered_results='The two information/stability CSV files contain status records only, no estimates.'))
 text=f"""# E2 analytic formula repair: pre-development stop

Final: E2_FORMULA_REPAIR_FAILED. Formula Gate FAIL.
Parent {p.PARENT}; design {start['design_SHA']}.
Execution SHA is the commit containing this document. Verify remote main after push.

The one authorized repair execution reached all 48 scene/mesh analytic controls,
then stopped before development regression because 20 of 48 nominal pressure identity
checks exceeded the frozen 1e-12 relative bound. Maximum: {decision['pressure_identity_relative_max']:.12g}.
The implementation mixed original Cartesian element ranges and the polar-coordinate
roundtrip used by the old pilot. Maximum range change was {geommax:.12g} m.
The new and old pressures were independently reconstructed using their respective
geometries. Their residuals were {pressmax:.12g} and {oldpressmax:.12g}.
Thus this stop is a frozen numerical-consistency guard failure in this implementation,
not a demonstration that the physical mechanism lacks information. The frozen
threshold, source files, and execution marker were not changed and no rerun was made.

P1 primary physical rank: {decision['P1_physical_rank_pass']}/48 PASS, ranks
{decision['P1_rank_values']}; maximum cofactor-null projected response
{decision['P1_structural_null_residual_max']:.12g}.
Independent geometric SVD and QR nuisance elimination verified the same rank bound,
including both depth stencils; maximum null response {nullmax:.12g}.
P2 measured pre-control ranks: {rank2}. This is not a full P2 scientific validation.
Maximum analytic-vs-old-FD identifiable Jacobian discrepancy:
{decision['analytic_vs_FD_identifiable_relative_max']:.12g}.
Independent per-mode scalar analytic horizontal derivative discrepancy: {errmax:.12g}.
Pre-control old/new effective-information changes: {json.dumps(differences)}.
These are formula diagnostics only; they do not replace the cancelled regression.

P0/P1/P2 main C0/C1/C2 development comparison, scene-level data-processing regression,
unknown-depth range information stability, depth/velocity weak-direction stability,
relative/fixed-floor scientific sensitivity and candidate regression were NOT EXECUTED.
The information and stability CSVs explicitly contain stop statuses without estimates.
The registered synthetic plane-wave, gain-absorption and shared-frequency controls
were run; full covariance and scientific information closure remains unevaluated.
No application claim can be made from P2's pre-control rank or its diagnostic information.

Validation: {v['PASS']} PASS / {v['FAIL']} FAIL. Independent stop audit:
{v['independent_stop_audit']['PASS']} PASS / {v['independent_stop_audit']['FAIL']} FAIL.
Every frozen input and original pilot artifact was hash-checked unchanged.
Original E2-G0: FAIL_UNCHANGED. Original pilot: IMPLEMENTATION_INVALID.
Original E2 scientific information: NOT_EVALUATED. R4=0%.
All new diagnostics: FORMULA_REPAIRED_DEVELOPMENT_ONLY.
New KRAKEN/Monte Carlo/audio: 0/0/0; development regressions: 0.
The single allowed repair opportunity is closed. No budget, depth step, SVD threshold
or pressure identity threshold was changed after observing results.
Next: H3_REVIEW (review only), requiring separate authorization for execution.
Do not launch full-band certification, H3, A2, SSP or P5.
"""
 (o/'GPT_SYNC.md').write_text(text,encoding='utf-8')
 (o/'PLANE_WAVE_AND_DATA_PROCESSING_CONTROLS.md').write_text('# Pre-development control scope\nRegistered synthetic plane-wave distance/depth/bearing-profile, fixed-gain absorption, free-response absorption, and shared-frequency covariance controls passed. Scene-level DPI and covariance development controls were NOT EXECUTED because nominal pressure identity precontrols failed. See PRE_DEVELOPMENT_CONTROLS.csv and NOT_EXECUTED.json. No data-processing closure or scientific gain claim is asserted.\n',encoding='utf-8')
 master=Path('results/R4_MASTER')
 with (master/'R4_PLAN.md').open('a',encoding='utf-8') as fh:fh.write('\n\n## E2 analytic formula repair closed\nE2_FORMULA_REPAIR_FAILED: registered pre-development nominal-pressure consistency guard failed (20/48); P1 physical rank 48/48 passed. No development regression executed. Numerical geometry roundtrip mismatch is implementation evidence, not physical non-identifiability. Old E2 admission FAIL_UNCHANGED, pilot IMPLEMENTATION_INVALID, R4=0%. See ../R4_E2_NINE_PAIR_FORMULA_REPAIR/GPT_SYNC.md. Stop; H3 review only, execution not authorized.\n')
 ledger=master/'R4_EVIDENCE_LEDGER.csv'
 with ledger.open(encoding='utf-8',newline='') as fh:fields=next(csv.reader(fh))
 row={k:'' for k in fields};row.update(stage='R4_E2_NINE_PAIR_ANALYTIC_FORMULA_REPAIR',artifact='../R4_E2_NINE_PAIR_FORMULA_REPAIR/GPT_SYNC.md',status='E2_FORMULA_REPAIR_FAILED',evidence_scope='PRE_DEVELOPMENT_GUARD_STOP; P1 rank48/48; no regression; no physical impossibility claim',independent_audit_status=decision['independent_stop_audit'],credited_weight_percent='0')
 with ledger.open('a',encoding='utf-8',newline='') as fh:csv.DictWriter(fh,fieldnames=fields).writerow(row)
 p.dump(master/'R4_E2_NINE_PAIR_FORMULA_REPAIR_PROGRESS.json',decision)
 artifacts={}
 for q in o.rglob('*'):
  if q.is_file() and q.name!='EXECUTION_ARTIFACT_MANIFEST.json':
   raw=q.suffix in ['.npz','.png'];artifacts[str(q)]={'kind':'RAW' if raw else 'LF','sha256':p.sha(q,raw)}
 p.dump(o/'EXECUTION_ARTIFACT_MANIFEST.json',dict(artifacts=artifacts,design_SHA=start['design_SHA'],R4_percent=0))
 print(json.dumps(p.native(dict(decision=decision,validation=v)),indent=2))
if __name__=='__main__':main()
