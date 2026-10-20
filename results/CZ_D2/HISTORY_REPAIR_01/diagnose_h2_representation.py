"""Read-only diagnosis of retained strict H2 float equality check; not a rerun."""
from pathlib import Path
import csv,json,hashlib,zipfile,struct,ast,math
R=Path(__file__).resolve().parents[3];O=Path(__file__).resolve().parent
p=R/'results/R4_SINGLE_HLA_RADIAL_INCREMENT_PRECISION_PILOT/GENERATOR_TRUTH_AND_INNOVATIONS.npz';s=R/'results/R4_SINGLE_HLA_SOURCE_PROFILED_HORIZONTAL_PILOT/OFFGRID_TRUTH_EVALUATION_ONLY.csv'
with zipfile.ZipFile(p) as z:b=z.read('horizontal_truth.npy')
nsize=2 if b[6]==1 else 4;n=struct.unpack('<H' if nsize==2 else '<I',b[8:8+nsize])[0];start=8+nsize+n;v=struct.unpack('<32d',b[start:])
r=list(csv.DictReader(s.open(encoding='utf-8-sig')));keys=['r_km','theta_deg','v_mps','psi_deg'];width=json.loads((O.parent/'DESIGN_FREEZE.json').read_text())['panel']['old_and_new_four_state_neighborhood_half_width'];report=[]
for i,x in enumerate(r):
 for j,k in enumerate(keys):
  a=v[4*i+j];c=float(x[k]);d=abs(a-c);report.append(dict(row=i,coordinate=k,archive_float_hex=a.hex(),csv_float_hex=c.hex(),absolute_difference=d,difference_ULP=d/math.ulp(c),original_neighborhood_width=width[j],within_original_exclusion=d<=width[j]))
manifest=json.loads((p.parent/'OBSERVATION_GENERATION_MANIFEST.json').read_text());hashed=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
result=dict(status='STRICT_FLOAT_EQUALITY_FAIL_RETAINED;EXCLUSION_COVERAGE_UNAFFECTED',original_check_not_overwritten=True,archived_values_source='H2 generator explicitly uses pandas.read_csv then to_numpy; this diagnostic uses Python csv+float and standard-library NPY decoding',max_absolute_difference=max(x['absolute_difference'] for x in report),max_ULP_difference=max(x['difference_ULP'] for x in report),all32_within_original_neighborhood=all(x['within_original_exclusion'] for x in report),source_panel_historical_SHA_match=manifest['source_panel_sha256']==hashed(s),archive_historical_SHA_match=manifest['sha256'][p.name]==hashed(p),no_new_tolerance_or_state=True,formal_execution=False,entries=report)
with (O/'H2_REPRESENTATION_DIAGNOSTIC.json').open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({k:result[k] for k in ['status','max_absolute_difference','max_ULP_difference','all32_within_original_neighborhood','source_panel_historical_SHA_match','archive_historical_SHA_match']}))
