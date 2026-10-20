from pathlib import Path
import csv,json,re,hashlib,collections
R=Path(__file__).resolve().parents[3];O=Path(__file__).resolve().parent;OLD=O.parent/'HISTORY_REPAIR_01'
def rows(p):return list(csv.DictReader(p.open(encoding='utf-8-sig')))
def put(n,x):
 with (O/n).open('x',encoding='utf-8',newline='') as f:
  if n.endswith('.json'):json.dump(x,f,indent=2,ensure_ascii=False)
  else:
   w=csv.DictWriter(f,fieldnames=list(x[0]),lineterminator='\n');w.writeheader();w.writerows(x)
put('REPAIR02_STARTED.json',dict(parent='65f9cec61c05e43b3a628f91b078bfdcc4b5318b',snapshot='a580499475cc4148f1195e6256c4813c6e225726',mode='read-only source semantic review',new_panel=False,formal_execution_authorized=False))
inv=rows(OLD/'SNAPSHOT_INVENTORY.csv');hits=[];stats=collections.Counter()
for x in inv:
 p=R/x['path']
 if p.suffix!='.json':continue
 obj=json.loads(p.read_text(encoding='utf-8-sig'))
 def walk(v,loc):
  if isinstance(v,dict):
   keys=set(v);s='|'.join(keys)
   rkeys=[k for k in keys if re.fullmatch(r'(?:truth_|true_|ref_|alt_)?(?:r|r0|range)(?:_true)?(?:_km|_m)?',k)]
   thkeys=[k for k in keys if re.fullmatch(r'(?:truth_|true_|ref_|alt_)?(?:theta|theta0|bearing)(?:_deg|_rad)?',k)]
   motion=[k for k in keys if re.fullmatch(r'(?:truth_|true_|ref_|alt_)?(?:v|v_mps|speed|psi|psi_deg|heading_deg)',k)]
   if rkeys and (thkeys or motion):
    hits.append(dict(path=x['path'],sha256=x['sha256'],pointer=loc,keys='|'.join(sorted(keys)),preview=json.dumps(v,ensure_ascii=False)[:2500]));stats[x['path'].split('/')[1] if '/' in x['path'] else x['path']]+=1
   for k,a in v.items():walk(a,loc+'/'+k)
  elif isinstance(v,list):
   for i,a in enumerate(v):walk(a,loc+'/'+str(i))
 walk(obj,'')
put('STRUCTURED_STATE_OBJECTS.csv',hits)
put('STRUCTURED_SCAN_SUMMARY.json',dict(objects=len(hits),families=dict(stats)))
print(json.dumps(dict(objects=len(hits),families=dict(stats))))
