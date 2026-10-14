#!/usr/bin/env python3
"""Generated synthetic MOD records for UNIT TEST ONLY, never actual project evidence."""
import sys,csv,struct,tempfile,shutil
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import t2_n3_prec1 as p
import t2_n3_prec1_cold_audit as cold

def run():
 with tempfile.TemporaryDirectory() as td:
  root=Path(td);b=root/p.REPO_BASE;gdir=b/'package';gdir.mkdir(parents=True)
  shutil.copy(Path(__file__).parent/'T2_TILT_G0_ARRAY_GEOMETRY_20261013_v1.csv',gdir/'T2_TILT_G0_ARRAY_GEOMETRY_20261013_v1.csv')
  geo=p.make_geometry(root);history=[]
  for high,m in [(1800,741),(2400,1046)]:
   stem=f'n3_PR_R_f201_h{high}_n160001';path=b/'modal_archive'/stem/(stem+'.mod');path.parent.mkdir(parents=True)
   z=np.array(sorted(set([180,200,210,220]+[float(v) for xyz in geo.values() for v in xyz[:,2]])),dtype=float);N=len(z)
   rec=max(160,8*N);rec=(rec+3)//4*4
   allbytes=(7+m)*rec+m*8
   blob=bytearray(allbytes)
   struct.pack_into('<i',blob,0,rec//4);struct.pack_into('<4i',blob,84,1,1,N,N);struct.pack_into('<d',blob,3*rec,201.)
   struct.pack_into('<'+str(N)+'f',blob,4*rec,*z);struct.pack_into('<i',blob,5*rec,m)
   mu,k=p.analytic_spectrum(m)
   for i in range(m):
    ff=np.sqrt(2/5000)*np.sin(mu[i]*z)
    for j,val in enumerate(ff):struct.pack_into('<ff',blob,(7+i)*rec+j*8,np.float32(val),0.)
    struct.pack_into('<ff',blob,(7+m)*rec+i*8,np.float32(k[i]),0.)
   path.write_bytes(blob)
   mm=p.mod_read(path)
   for position in p.PROBES:
    for res,xyz in geo.items():
     xx=p.field(mm,position,xyz,False,False)
     zz=p.field(mm,position,xyz,True,True)
     history.append({'stem':stem,'probe':position[0],'resource':res,'spatial_shape_error':p.stable_shape(xx,zz),'PASS':'True'})
  p.save_csv(b/'ANALYTIC_KRAKEN_MODE_AND_FIELD_CROSSCHECK.csv',history)
  out=root/'output'
  d=p.run(root,out,verify_frozen=False)
  cold.main(root,out)
  assert d['all_40_control_rows']==40
  print('TEST_SYNTHETIC_ONLY_PASSED: actual N3 MOD data not present in this test')
if __name__=='__main__':run()
