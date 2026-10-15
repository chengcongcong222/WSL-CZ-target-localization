"""One T3-W1 main: verify nominal information then generate paired public archives
and fit once per configuration. No acoustic solver or historical-result mutation.
"""
from pathlib import Path
import json,csv,hashlib,math,time,traceback,os
import numpy as np
import t3_w1_generate as gen
import t3_w1_information as inf
from t3_w1_model import estimate
from t3_w1_statistics import csvread,csvwrite,evaluate,summarize,decision
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent;LOCAL=Path('D:/ProjectStorage/WSL-CZ/T3_W1')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def meta(p):return {'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)}

def main():
    cfg=read(OUT/'T3_W1_DESIGN_FREEZE.json');panel=csvread(ROOT/cfg['truth_panel']);geom=[];infos=[];numerics=[];matrices=[];known_matrices=[];covariances=[];effective=[];variance_contributions=[];public=[];private=[];paired_hash={};start=time.monotonic()
    t,mainp=gen.platform()
    for qname,quality in cfg['qualities'].items():
        for bname,scheme in cfg['geometries'].items():
            for tr in panel:
                state=gen.cartesian(tr);target=state[:2]+t[:,None]*state[2:];s6=np.r_[state,0.,0.]
                for side in [-1,1]:
                    _,P=gen.node_geometry(scheme,side);distance=np.linalg.norm(target-P[:,1],axis=1);main_distance=np.linalg.norm(target-P[:,0],axis=1)
                    key={'quality':qname,'geometry':bname,'case_id':tr['case_id'],'side':side}
                    geom.append({**key,'aux_target_min_m':float(min(distance)),'aux_target_max_m':float(max(distance)),'main_target_min_m':float(min(main_distance)),'main_target_max_m':float(max(main_distance)),'node_separation_m':float(np.hypot(scheme['forward_km'],scheme['lateral_km'])*1000),'aux_outside45to60km_epochs':int(np.sum((distance<45000)|(distance>60000))),'main_outside45to60km_epochs':int(np.sum((main_distance<45000)|(main_distance>60000))),'epochs':121,'forward_km':scheme['forward_km'],'lateral_km':scheme['lateral_km'],'resource_scope':'two_nodes;fixed_global_offset;deployability_UNKNOWN;inrange_not_detection'})
                    angle,var=inf.finite_difference_check(s6,t,P,quality);numerics.append({**key,'max_scaled_mean_J_error':angle,'max_relative_variance_J_error':var,'mean_tolerance':1e-7,'variance_tolerance':1e-6,'PASS':angle<=1e-7 and var<=1e-6})
                    if not numerics[-1]['PASS']:
                        csvwrite(OUT/'JACOBIAN_CROSSCHECK.csv',numerics);save(OUT/'T3_W1_DECISION.json',{'status':'IMPLEMENTATION_INVALID','stage':'JACOBIAN_INFORMATION','MC_configurations':0,'R4_percent':0,'STOP':True});return
                    row,F,Fs,E,C=inf.analyze(s6,t,P,quality,False);r2,F2,_,E2,C2=inf.analyze(s6,t,P,quality,True);index=len(matrices);infos.extend([{**key,'matrix_index':index,**row},{**key,'matrix_index':index,**r2}]);matrices.append(F);known_matrices.append(F2);effective.append(E);covariances.append(C if C is not None else np.full((6,6),np.nan))
                    H,dv,A=inf.design(s6,t,P,quality);meanF=A[:242].T@A[:242]/inf.SCALE[:,None]/inf.SCALE[None,:];variance_contributions.append(F-meanF)
    csvwrite(OUT/'GEOMETRY_AND_RESOURCE_AUDIT.csv',geom);csvwrite(OUT/'JACOBIAN_CROSSCHECK.csv',numerics);csvwrite(OUT/'LOCAL_INFORMATION.csv',infos)
    np.savez_compressed(OUT/'INFORMATION_MATRICES.npz',F_data=np.array(matrices),F_bias_known=np.array(known_matrices),F_variance_term=np.array(variance_contributions),effective_velocity_scaled=np.array(effective),C_data_or_nan_when_unidentifiable=np.array(covariances),parameter_scale=inf.SCALE)
    save(OUT/'INFORMATION_SEAL.json',{'Jacobian_controls':len(numerics),'information_rows':len(infos),'matrices':len(matrices),'files':[meta(OUT/name) for name in ['GEOMETRY_AND_RESOURCE_AUDIT.csv','JACOBIAN_CROSSCHECK.csv','LOCAL_INFORMATION.csv','INFORMATION_MATRICES.npz']],'uniform_hard_prior_not_added_to_Fisher':True})
    print('INFORMATION_COMPLETE',len(numerics),len(infos),flush=True)
    # Generate and seal every public block before fitting any test record.
    (LOCAL/'private').mkdir(exist_ok=True);(LOCAL/'public').mkdir(exist_ok=True)
    for g,tr in enumerate(panel):
        rng=np.random.Generator(np.random.PCG64(np.random.SeedSequence([cfg['seed'],g])));unif=rng.random((500,2));normal=rng.standard_normal((500,121,6));sides=np.where(np.arange(500)%2==0,-1,1)
        priv=LOCAL/'private'/f"{tr['case_id']}_innovations.npz";np.savez_compressed(priv,uniforms=unif,normals=normal,sides=sides);private.append({'case_id':tr['case_id'],**meta(priv)})
        for qname,quality in cfg['qualities'].items():
            for bname,scheme in cfg['geometries'].items():
                t,P,y,bias=gen.observations(tr,scheme,quality,normal,unif,sides)
                p=LOCAL/'public'/f"{qname}_{bname}_{tr['case_id']}.npz";np.savez_compressed(p,times=t,reported_positions=P,bearings=y,sides=sides,replicate_ids=np.arange(500))
                mainhash=hashlib.sha256(t.tobytes()+np.ascontiguousarray(P[:,:,0]).tobytes()+np.ascontiguousarray(y[:,:,0]).tobytes()).hexdigest();basekey=qname,tr['case_id']
                if basekey in paired_hash:assert mainhash==paired_hash[basekey],'MAIN_INPUT_PAIR_IDENTITY'
                else:paired_hash[basekey]=mainhash
                public.append({'config':qname+'_'+bname,'quality':qname,'geometry':bname,'case_id':tr['case_id'],'records':500,'mirror_minus':250,'mirror_plus':250,'MAIN_byte_hash':mainhash,'innovation_sha256':sha(priv),**meta(p)})
        print('ARCHIVED',tr['case_id'],flush=True)
    csvwrite(OUT/'PUBLIC_BLOCK_IDENTITIES.csv',public);csvwrite(OUT/'PRIVATE_INNOVATION_IDENTITIES.csv',private);save(OUT/'OBSERVATION_SEAL.json',{'public_blocks':len(public),'private_blocks':len(private),'paired_base_innovations':6000,'method_configurations':48000,'all_inputs_generated_before_any_fit':True,'public_index_sha256':sha(OUT/'PUBLIC_BLOCK_IDENTITIES.csv'),'private_index_sha256':sha(OUT/'PRIVATE_INNOVATION_IDENTITIES.csv'),'seed':cfg['seed']})
    fields=['config','quality','geometry','case_id','replicate','draw_index','side','accepted','failure_reason','r_hat_km','theta_hat_deg','v_hat_mps','psi_hat_deg','common_hat_deg','half_diff_hat_deg','solver_status','nfev','cost','optimality','initializer_fallback','initializer_intersections','initial_r_km','initial_theta_deg','initial_v_mps','initial_psi_deg','initial_common_deg','initial_diff_deg','state_at_bound','bias_at_bound','common_at_bound','diff_at_bound','rank_augmented_likelihood_jac','smallest_largest_jac_ratio','mean_residual_RMS_standard']
    runpath=OUT/'ESTIMATOR_OUTPUTS.csv';completed=0
    with runpath.open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,lineterminator='\n');writer.writeheader()
        # Frozen order Q1/Q2, B0/B1/B2/B3, H01..H12, rep0..499.
        for qname,quality in cfg['qualities'].items():
            for bname in cfg['geometries']:
                for g,tr in enumerate(panel):
                    block=next(r for r in public if r['config']==qname+'_'+bname and r['case_id']==tr['case_id']);archive=np.load(block['path']);t=archive['times'];P=archive['reported_positions'];y=archive['bearings'];sides=archive['sides']
                    for rep in range(500):
                        key={'config':qname+'_'+bname,'quality':qname,'geometry':bname,'case_id':tr['case_id'],'replicate':rep,'draw_index':g*500+rep,'side':int(sides[rep])};reason='NONE'
                        try:h,diag=estimate(t,P[rep],y[rep],quality)
                        except (ValueError,FloatingPointError,np.linalg.LinAlgError) as exc:
                            h=np.full(6,np.nan);diag={'accepted':False,'solver_status':-999,'nfev':0,'cost':'NA','optimality':'NA','initializer_fallback':False,'initializer_intersections':0,'initial_h':[float('nan')]*6,'state_at_bound':False,'bias_at_bound':False,'common_at_bound':False,'diff_at_bound':False,'rank_augmented_likelihood_jac':-1,'smallest_largest_jac_ratio':'NA','mean_residual_RMS_standard':'NA'};reason=type(exc).__name__+':'+str(exc)
                        if not diag['accepted'] and reason=='NONE':reason='SOLVER_NONCONVERGENCE_OR_NONFINITE'
                        initial=diag.pop('initial_h');row={**key,'failure_reason':reason,**dict(zip(['r_hat_km','theta_hat_deg','v_hat_mps','psi_hat_deg','common_hat_deg','half_diff_hat_deg'],h)),**dict(zip(['initial_r_km','initial_theta_deg','initial_v_mps','initial_psi_deg','initial_common_deg','initial_diff_deg'],initial)),**diag};writer.writerow(row);completed+=1
                    archive.close();handle.flush();print('FITTED',qname,bname,tr['case_id'],completed,'/48000',flush=True)
                    save(OUT/'MAIN_PROGRESS.json',{'completed':completed,'science_wall_s':time.monotonic()-start,'latest_config':qname+'_'+bname,'latest_case':tr['case_id']})
    save(OUT/'ESTIMATOR_SEAL.json',{'configurations':completed,'files':[meta(runpath)],'one_fit_per_configuration':True,'truth_not_initializer_or_selector':True,'no_global_branch_certificate':True})
    raw=csvread(runpath);errors=evaluate(raw,panel);summary,globalrows,paired=summarize(errors,raw);csvwrite(OUT/'RUN_ERRORS.csv',errors);csvwrite(OUT/'CASE_STATISTICS.csv',summary);csvwrite(OUT/'GLOBAL_STATISTICS.csv',globalrows);csvwrite(OUT/'PAIRED_CHANGES.csv',paired)
    result=decision(summary,globalrows,infos);result.update(status='PAIRED_CONFIRMATION_EXECUTED_ONCE',configurations=completed,base_innovations=6000,information_rows=len(infos),nominal_geometry_checks=len(geom),input_blocks=len(public),scientific_wall_s=time.monotonic()-start,independent_cold_audit='PENDING');save(OUT/'T3_W1_DECISION.json',result);print('MAIN_COMPLETE',completed,flush=True)

if __name__=='__main__':
    try:main()
    except Exception:
        (OUT/'MAIN_EXCEPTION.txt').write_text(traceback.format_exc(),encoding='utf-8')
        save(OUT/'T3_W1_DECISION.json',{'stage':'R4-T3-W1','status':'EXECUTION_INCOMPLETE_NO_RETRY','depth':'NOT_EVALUATED','R4_percent':0,'STOP':True});raise
