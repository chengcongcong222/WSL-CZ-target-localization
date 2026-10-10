"""Receiver reads only waveform and public parameters, never private manifest."""
from hla_sp1_core import ROOT,OUT,LOCAL,read,rows,sha,dump,guard,stamp
import numpy as np

def main(kind):
 dest=LOCAL/kind;meta=read(dest/'PUBLIC_METADATA.json');index=rows(dest/'PUBLIC_RECORD_INDEX.csv');path=dest/'RECEIVER_DFT.npy'
 if path.exists():raise RuntimeError('no DFT rerun')
 a=np.lib.format.open_memmap(path,mode='w+',dtype=np.complex128,shape=(len(index),meta['K'],3,8));n=np.arange(meta['N']);weights=np.exp(-2j*np.pi*np.array(meta['analysis_frequencies_hz'])[:,None]*n/meta['fs_hz'])/meta['N']
 for i,r in enumerate(index):
  guard();assert sha(r['path'])==r['sha256'];x=np.load(r['path']);a[i]=np.einsum('knm,fn->kfm',x,weights)
 a.flush();dump(OUT/(kind.upper()+'_DFT_FREEZE.json'),{'time_utc':stamp(),'file':str(path),'sha256':sha(path),'shape':list(a.shape),'waveform_hash_verified':len(index),'private_truth_source_H_gain_reads':0,'public_index_sha256':sha(dest/'PUBLIC_RECORD_INDEX.csv')});print('DFT FROZEN',kind,len(index),flush=True)
if __name__=='__main__':
 import sys
 main(sys.argv[1])
