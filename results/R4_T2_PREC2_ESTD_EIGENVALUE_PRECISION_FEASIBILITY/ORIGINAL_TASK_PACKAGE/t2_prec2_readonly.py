#!/usr/bin/env python3
"""PREC2 read-only archived PRT/MOD + independent lossless real-k eigenvalue subset.
No acoustic solvers, no waveform/MC, no mutation of source files.
Run: python t2_prec2_readonly.py --repo ROOT --config T2_PREC2_EXPERIMENT_SPEC_20261014_v1.json --out OUTPUT
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, os, struct, time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh_tridiagonal

N2 = 'results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET'
FREQ = (201, 235, 283)

def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()

def read_json(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))

def save_json(path,data):
    Path(path).write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

def save_csv(path,rows,columns):
    with Path(path).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns,lineterminator='\n',extrasaction='raise');w.writeheader();w.writerows(rows)

def parse_mod(path):
    data=Path(path).read_bytes()
    rec_words=struct.unpack_from('<i',data,0)[0]
    rec=rec_words*4
    if not 100<=rec<=20_000_000: raise ValueError('invalid MOD record size')
    n_freq,n_media,n_depths,n_stored=struct.unpack_from('<4i',data,84)
    nmodes=struct.unpack_from('<i',data,5*rec)[0]
    start=(7+nmodes)*rec
    if n_freq!=1 or nmodes<1 or n_depths!=n_stored or start+nmodes*8>len(data):
        raise ValueError('invalid MOD geometry / frequency / wavenumber record')
    k=np.frombuffer(data,dtype='<c8',count=nmodes,offset=start).copy()
    if not np.all(np.isfinite(k)) or np.any(k.real<=0):raise ValueError('invalid MOD wavenumbers')
    return k,{'rec_bytes':rec,'n_depths':n_depths,'nmodes':nmodes}

def parse_prt(path):
    lines=Path(path).read_text(encoding='utf-8',errors='strict').splitlines()
    starts=[i for i,l in enumerate(lines) if 'k (1/m)' in l and 'Phase Speed' in l]
    if len(starts)!=1:raise ValueError('PRT frequency/spectrum header ambiguous')
    rows=[]
    for l in lines[starts[0]+1:]:
        pieces=l.strip().split()
        if len(pieces)<5 or not pieces[0].isdigit():continue
        try:
            v=[float(t.replace('D','E').replace('d','e')) for t in pieces[1:5]]
        except ValueError:continue
        if not all(np.isfinite(vv) for vv in v):raise ValueError('nonfinite PRT row')
        rows.append({'index_1based':int(pieces[0]),'k_re':v[0],'k_im':v[1],'phase_speed':v[2],'group_speed':v[3]})
    if len(rows)<25 or len(set(r['index_1based'] for r in rows))!=len(rows):
        raise ValueError('incomplete or duplicate PRT spectrum')
    return rows

def env_essentials(path):
    lines=Path(path).read_text(encoding='utf-8-sig').splitlines()
    if len(lines)<50:raise ValueError('ENV suspiciously short')
    if not lines[3].strip().strip("'").startswith('CVWT'):
        raise ValueError('expected original CVWT nominal contract, no implicit switch to lossless')
    bottom=next((j for j in range(5,len(lines)) if lines[j].lstrip().startswith("'R'")),None)
    if bottom is None:raise ValueError('rigid bottom R not present')
    arr=np.array([list(map(float,line.split()[:6])) for line in lines[5:bottom]],dtype=float)
    if arr.shape != (251,6):raise ValueError(f'expected 251 SSP nodes, got {arr.shape}')
    z,c=arr[:,0],arr[:,1]
    if not np.all(np.diff(z)>0) or abs(z[0])>1e-9 or abs(z[-1]-5000)>1e-9:
        raise ValueError('SSP depth contract mismatch')
    if not np.allclose(arr[:,2],0) or not np.allclose(arr[:,3],1) or not np.allclose(arr[:,4:],0):
        raise ValueError('not the expected single-fluid, constant density profile')
    return z,c,{'options':lines[3].strip(),'bottom':lines[bottom].strip(),'frequency':float(lines[1]),'mesh':int(lines[4].split()[0]),'clow':float(lines[bottom+1].split()[0]),'chigh':float(lines[bottom+1].split()[1])}

def lossless_lambda_tridiagonal(f,z,c,nodes,mode):
    """2nd-order Dirichlet at z=0, Neumann at z=H; mass-symmetrized last offdiagonal.

    L phi = phi'' + (omega/c(z))**2 phi = k_real**2 phi.
    Lossless independent diagnostic only; does not refine attenuation imaginary part.
    """
    H=float(z[-1]);N=int(nodes);h=H/N
    depths=np.arange(1,N+1,dtype=np.float64)*h
    cz=np.interp(depths,z,c)
    om=2*np.pi*f
    main=(om/cz)**2-2/h**2
    off=np.full(N-1,1/h**2)
    off[-1]*=math.sqrt(2)
    idx=N-mode
    if idx<0:raise ValueError('mode index >= grid dimension')
    lam=eigh_tridiagonal(main,off,eigvals_only=True,select='i',select_range=(idx,idx),lapack_driver='stebz')[0]
    if lam<=0 or not np.isfinite(lam):raise ValueError('nonpropagating/invalid subset mode')
    return float(lam)

def richardson(vals):
    # 4 h levels h, h/2, h/4, h/8; extrapolate in h^2,
    # using lambda = kr**2 (the eigenproblem's actual eigenvalue).
    v=np.asarray(vals,float)
    if v.shape!=(4,):raise ValueError('four grids required')
    r2=np.array([(4*v[i+1]-v[i])/3 for i in range(3)])
    r4=np.array([(16*r2[i+1]-r2[i])/15 for i in range(2)])
    r6=(64*r4[1]-r4[0])/63
    if np.any(r2<=0) or np.any(r4<=0) or r6<=0:raise ValueError('bad richardson eigenvalue')
    return {'lambda_R2':r2.tolist(),'lambda_R4':r4.tolist(),'lambda_R6':float(r6),
            'kR4_coarse':math.sqrt(r4[0]),'kR4_fine':math.sqrt(r4[1]),'kR6':math.sqrt(r6),
            'R4_k_gap':abs(math.sqrt(r4[1])-math.sqrt(r4[0])),
            'R6_vs_R4fine_k_gap':abs(math.sqrt(r6)-math.sqrt(r4[1]))}

def execute(repo:Path,config:dict,out:Path):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=True)
    source=[];printed=[];summary=[];spectrum=[];count=0;diagnostic_k_err=config['phase_target_radians']/config['phase_reference_range_m'];all_zc={}
    for entry in config['archive_sets']:
        basename=entry['stem'];f=int(entry['freq']);mesh=entry['mesh']
        names={ext:repo/entry['paths'][ext] for ext in ('env','mod','prt')}
        for ext,p in names.items():
            if not p.is_file():raise RuntimeError('ARCHIVE_MISSING '+str(p))
            chk=sha(p)
            if chk!=entry['sha256'][ext]:raise RuntimeError('ARCHIVE_HASH_MISMATCH '+str(p))
            source.append({'stem':basename,'type':ext,'path':entry['paths'][ext],'sha256':chk})
        z,c,env=env_essentials(names['env'])
        if f!=env['frequency'] or mesh!=env['mesh'] or env['clow']!=1500 or env['chigh']!=2400:
            raise ValueError('E-STD source contract mismatch')
        if all_zc and (not np.array_equal(z,all_zc['z']) or not np.array_equal(c,all_zc['c'])):
            raise ValueError('SSP differs across frozen jobs; no cross-contract pooling')
        all_zc={'z':z,'c':c}
        km,info=parse_mod(names['mod']);pr=parse_prt(names['prt']);count+=info['nmodes'];
        if len(km)!=entry['modes']:raise ValueError('MOD mode count does not match historical record')
        if 'prt_rows_expected' in entry and len(pr)!=entry['prt_rows_expected']:
            raise ValueError('PRT expected sample count not present; STOP')
        maxdiff=0;sample=[]
        for r in pr:
            i=r['index_1based']-1
            if i<0 or i>=len(km):raise ValueError('PRT mode index not found')
            kfloat=np.float32(r['k_re'])
            if kfloat.view(np.uint32)!=np.float32(km[i].real).view(np.uint32):
                raise RuntimeError(f'PRT_TO_MOD_FLOAT32_IDENTITY_FAIL {basename} {i+1}')
            difference=float(r['k_re'])-float(km[i].real)
            maxdiff=max(maxdiff,abs(difference))
            if not np.isfinite(km[i].imag) or np.sign(r['k_im'])!=np.sign(km[i].imag):
                raise RuntimeError('attenuation sign mismatch')
            printed.append({'stem':basename,'freq_hz':f,'mesh':mesh,'mode':i+1,
                            'k_PRT_real':r['k_re'],'k_MOD_real':float(km[i].real),
                            'delta_k_PRT_minus_MOD':difference,
                            'phase_delta_at_60km_rad':difference*config['phase_reference_range_m'],
                            'k_PRT_imag':r['k_im'],'k_MOD_imag':float(km[i].imag),
                            'K_REAL_MORE_DIGITS_NOT_PROOF_OF_ACCURACY':True})
        summary.append({'stem':basename,'frequency_hz':f,'mesh':mesh,'MOD_modes':len(km),
                        'PRT_rows':len(pr),'PRT_f32_exact':len(pr),
                        'max_abs_delta_k':maxdiff,'max_phase_delta_at60km_rad':maxdiff*config['phase_reference_range_m'],
                        'PROVENANCE':'N2_GIT_ARCHIVES_NOT_NEW_SOLVER'})
        if mesh==160001:
            chosen=[pr[j] for j in config['prt_rank_positions'][:-1]]+[pr[-1]]
            if len(set(t['index_1based'] for t in chosen))!=5:raise RuntimeError('sampled mode collisions')
            for p in chosen:
                j=p['index_1based'];eigs=[lossless_lambda_tridiagonal(f,z,c,N,j) for N in config['FDM_grid_intervals']]
                rc=richardson(eigs)
                row={'freq_hz':f,'mode':j,'K_mod':float(km[j-1].real),'K_prt':p['k_re'],
                     'PRT_MODE_IDENTITY':True,'K_real_lossless_R6':rc['kR6'],
                     'R4_k_gap':rc['R4_k_gap'],'R6_vs_R4fine_k_gap':rc['R6_vs_R4fine_k_gap'],
                     'PRT_minus_R6_k':p['k_re']-rc['kR6'],
                     'MOD_minus_R6_k':float(km[j-1].real)-rc['kR6'],
                     'R4_CONVERGED_AT_PHASE_BUDGET':rc['R4_k_gap']<=diagnostic_k_err,
                     'FDM_LOSSLESS_ONLY_NOT_ACTUAL_THORP_COMPLEX_K':True}
                row.update({f'lambda_N{N}':value for N,value in zip(config['FDM_grid_intervals'],eigs)})
                spectrum.append(row)
        if time.monotonic()-start>config['wall_s']:raise RuntimeError('SCIENCE_WALL_BUDGET_HARD_STOP')
    if 'expected_total_prt_samples' in config and len(printed)!=config['expected_total_prt_samples']:
        raise ValueError('archive print sample total not as frozen')
    cols1=['stem','type','path','sha256'];cols2=list(printed[0]);cols3=list(summary[0]);cols4=list(spectrum[0])
    save_csv(out/'PREC2_SOURCE_IDENTITIES.csv',source,cols1)
    save_csv(out/'PREC2_PRT_MOD_ALL_SAMPLES.csv',printed,cols2)
    save_csv(out/'PREC2_PRT_MOD_SUMMARY.csv',summary,cols3)
    save_csv(out/'PREC2_LOSSLESS_SUBSET_FDM.csv',spectrum,cols4)
    save_json(out/'PREC2_DECISION.json',{
        'stage':'T2-PREC2','parent_sha':config['parent_sha'],
        'status':'SPARSE_PRT_MOD_ROUNDING_REPRODUCED_WITH_LOSSLESS_SUBSET_DIAGNOSTIC',
        'sources_verified':len(source),'samples_verified':len(printed),'MOD_mode_total_duplicated_mesh_views':count,
        'subset_roots':len(spectrum),'FDM_R4_gap_all_within_target':all(r['R4_CONVERGED_AT_PHASE_BUDGET'] for r in spectrum),
        'any_k_full_field_refinement_completed':False,'actual_imaginary_k_refined':False,
        'accuracy_certified_in_actual_lossy_waveguide':False,
        'parent_G0_N2_N3_PREC1_unchanged':True,'A_B_receiver_performance':'NOT_EVALUATED',
        'R4_percent':0,'new_KRAKEN_FIELD_SCOOTER_MC_DFT':0,'state':'STOP_AFTER_COLD_AUDIT',
        'elapsed_s':time.monotonic()-start,
        'note':'FDM is an independent self-adjoint real-k diagnostic, not a high-precision full complex-mode provider.'
    })
    print('PREC2 READ-ONLY COMPLETE: samples',len(printed),'subset',len(spectrum),'elapsed',time.monotonic()-start)

if __name__=='__main__':
    pa=argparse.ArgumentParser();pa.add_argument('--repo',type=Path,required=True);pa.add_argument('--config',type=Path,required=True);pa.add_argument('--out',type=Path,required=True)
    x=pa.parse_args();execute(x.repo.resolve(),read_json(x.config),x.out.resolve())
