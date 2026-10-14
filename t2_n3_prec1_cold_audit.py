#!/usr/bin/env python3
"""Independent T2-N3-PREC1 cold audit; does not import main analysis code.

Independently reads binary MOD with struct, finds analytic eigenroots with brentq,
and reconstructs all four combinations by channelwise scalar loops with cmath.
"""
import argparse,csv,cmath,hashlib,json,math,struct
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

H=5000.;C=1500.;F=201.
BASE='results/R4_T2_N3_BOUNDARY_PHYSICS_ANALYTIC_REFERENCE'
MOD={1800:f'{BASE}/modal_archive/n3_PR_R_f201_h1800_n160001/n3_PR_R_f201_h1800_n160001.mod',2400:f'{BASE}/modal_archive/n3_PR_R_f201_h2400_n160001/n3_PR_R_f201_h2400_n160001.mod'}

def table(file):
 with open(file,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

def binary(path):
 blob=Path(path).read_bytes();rec=4*struct.unpack_from('<i',blob,0)[0]
 n=struct.unpack_from('<i',blob,92)[0];d=[struct.unpack_from('<f',blob,4*rec+4*j)[0] for j in range(n)]
 m=struct.unpack_from('<i',blob,5*rec)[0]
 phi=[];k=[]
 for i in range(m):
  phi.append([complex(*struct.unpack_from('<ff',blob,(7+i)*rec+8*j)) for j in range(n)])
  k.append(complex(*struct.unpack_from('<ff',blob,(7+m)*rec+8*i)))
 return {'n':n,'m':m,'d':d,'phi':phi,'k':k}

def root_spectrum(m):
 # Solve cos(mu H)=0 independently of closed-form index in main implementation.
 mu=[]; k=[]
 for i in range(m):
  a=(i+.25)*math.pi;b=(i+.75)*math.pi
  x=brentq(math.cos,a,b,xtol=1e-13)
  u=x/H
  mu.append(u);k.append(math.sqrt((2*math.pi*F/C)**2-u*u))
 return mu,k

def coordinate_map(root):
 r={}
 for row in table(root/BASE/'package/T2_TILT_G0_ARRAY_GEOMETRY_20261013_v1.csv'):
  if row['resource'] not in ('A0_0','A_P4M','B1_V4','B2_DUP'):continue
  r.setdefault(row['resource'],[]).append((int(row['channel']),*[float(row[x]) for x in ('x_m','y_m','z_m')]))
 return {name:[point[1:] for point in sorted(v)] for name,v in r.items()}

def idx(z, d):
 i=min(range(len(d)),key=lambda i:abs(d[i]-z))
 if abs(z-d[i])>1e-5:raise RuntimeError('COLD_DEPTH_MISMATCH')
 return i

def wave(mod,mu,kk,position,coords,ak,ap):
 _,rr,th,sd=position;x=1000*rr*math.cos(math.radians(th));y=1000*rr*math.sin(math.radians(th))
 source=idx(sd,mod['d']);out=[]
 for rx,ry,rz in coords:
  rd=idx(rz,mod['d']);r=math.hypot(x-rx,y-ry)
  terms=[]
  for i in range(mod['m']):
   k=kk[i] if ak else mod['k'][i]
   if ap:
    p=(math.sqrt(2/H)*math.sin(mu[i]*sd))*(math.sqrt(2/H)*math.sin(mu[i]*rz))
   else:p=mod['phi'][i][source]*mod['phi'][i][rd]
   ker=cmath.sqrt(2*math.pi/(k*r))*cmath.exp(-1j*(k*r+math.pi/4))
   terms.append(p*ker)
  # Accumulate real and imaginary components separately to reduce summation bias.
  out.append(complex(math.fsum(t.real for t in terms),math.fsum(t.imag for t in terms)))
 return out

def qr_error(a,b):
 av=np.array(a,np.complex128);bv=np.array(b,np.complex128)
 if not np.isfinite(av).all() or not np.isfinite(bv).all():raise ValueError('NAN_FIELD')
 av/=np.linalg.norm(av);bv/=np.linalg.norm(bv)
 Q,R=np.linalg.qr(np.column_stack([bv,av]),mode="reduced")
 return float(abs(R[1,1]))

def main(repo,out):
 repo=Path(repo);out=Path(out)
 if (out/'PREC1_COLD_AUDIT_STARTED.json').exists():raise RuntimeError('COLD_AUDIT_ONCE_ONLY_ALREADY_STARTED')
 (out/'PREC1_COLD_AUDIT_STARTED.json').write_text(json.dumps({'stage':'T2_N3_PREC1','independent_cold':True,'solver_calls':0})+'\n',encoding='utf-8')
 items=table(out/'PREC1_ALL_40_CASES.csv');geom=coordinate_map(repo)
 checks=[];count=0;max_delta=0.;max_rounding=0.
 decision=json.loads((out/'PREC1_DECISION.json').read_text(encoding='utf-8'))
 for rel,expected in decision['source_hashes'].items():
  h=hashlib.sha256((repo/rel).read_bytes()).hexdigest();checks.append({'check':'input_seal_identity','key':rel,'deviation':0. if h==expected else 1.,'limit':0.,'PASS':h==expected})
 for high in (1800,2400):
  m=binary(repo/MOD[high]);mus,ks=root_spectrum(m['m'])
  assert all(abs(k.imag)==0 for k in m['k'])
  for i in range(m['m']):
   delta=abs(m['k'][i].real-float(np.float32(ks[i])));max_rounding=max(max_rounding,delta)
   checks.append({'check':'stored_kr_equals_independent_float32_rounding','key':f'{high}_{i}','deviation':delta,'limit':0.0,'PASS':delta==0})
  for x in items:
   if int(x['CHIGH'])!=high:continue
   p=next((n,float(q),float(t),float(z)) for n,q,t,z in [('P01',45.5,0,180),('P02',50,0,200),('P03',55,2.5,200),('P04',58,4.5,210),('P05',59.5,5,220)] if n==x['probe'])
   aa={}
   for ak,ap in [(False,False),(True,False),(False,True),(True,True)]:aa[(ak,ap)]=wave(m,mus,ks,p,geom[x['resource']],ak,ap)
   ref=aa[(True,True)]
   dvals=[qr_error(aa[(False,False)],ref),qr_error(aa[(True,False)],ref),qr_error(aa[(False,True)],ref),qr_error(aa[(False,False)],aa[(True,False)]),qr_error(aa[(False,False)],aa[(False,True)])]
   fields=['C0_stored_k_stored_phi_vs_C3','C1_analytic_k_stored_phi_vs_C3','C2_stored_k_analytic_phi_vs_C3','C0_to_C1_direct','C0_to_C2_direct']
   for field,value in zip(fields,dvals):
    delta=abs(value-float(x[field]));max_delta=max(delta,max_delta)
    checks.append({'check':f'cold_channelwise_{field}','key':f'{high}_{p[0]}_{x["resource"]}','deviation':delta,'limit':5e-10,'PASS':delta<=5e-10})
   count+=1
 assert count==40
 failed=sum(not r['PASS'] for r in checks)
 with (out/'PREC1_COLD_AUDIT_CHECKS.csv').open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(checks[0]),lineterminator='\n');w.writeheader();w.writerows(checks)
 result={'PASS':not failed,'FAIL':failed,'check_count':len(checks),'cold_cases':count,'modal_rounding_checks':sum(1 for c in checks if c['check']=='stored_kr_equals_independent_float32_rounding'),'max_shape_delta':max_delta,'max_k_round_delta':max_rounding,'no_external_solver':True,'original_N3_gate_unchanged':True,'R4_percent':0}
 (out/'PREC1_COLD_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
 print(json.dumps(result,indent=2));
 if failed:raise RuntimeError('COLD_AUDIT_FAIL')

if __name__=='__main__':
 pa=argparse.ArgumentParser();pa.add_argument('--repo',required=True);pa.add_argument('--out',required=True);a=pa.parse_args();main(a.repo,a.out)
