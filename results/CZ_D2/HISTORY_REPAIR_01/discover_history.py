"""Frozen-snapshot source discovery. Reads data/headers; never imports research code."""
from pathlib import Path
import subprocess,json,csv,hashlib,collections,re,zipfile,struct,ast,time
R=Path(__file__).resolve().parents[3];O=Path(__file__).resolve().parent
SNAP='a580499475cc4148f1195e6256c4813c6e225726'
def emit(name,rows):
 with (O/name).open('x',encoding='utf-8',newline='') as f:
  if name.endswith('.json'):json.dump(rows,f,ensure_ascii=False,indent=2)
  else:
   w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
with (O/'DISCOVERY_STARTED.json').open('x') as f:json.dump({'historical_snapshot':SNAP,'new_science':False},f)
listing=subprocess.check_output(['git','ls-tree','-r','-z',SNAP],cwd=R).split(b'\0')
inventory=[];arrays=[];headers=[];uses=[];jhits=[];external=[];errors=[]
pattern=re.compile(r'truth|true_state|h_true|htruth|h_truth|calibrat|holdout|offgrid|off_grid|generating_cases|panel',re.I)
for idx,line in enumerate(listing):
 if not line:continue
 metadata,path=line.split(b'\t',1);path=path.decode();p=R/path
 if not p.is_file():errors.append({'path':path,'error':'ABSENT'});continue
 suffix=p.suffix.lower();digest=sha(p)
 inv={'path':path,'parent_git_blob':metadata.decode().split()[-1],'bytes':p.stat().st_size,'sha256':digest,'format':suffix,'discovery':'NON_STATE_ASSET'}
 try:
  if suffix in ['.npz','.npy']:
   inv['discovery']='NUMERIC_ARRAY_HEADERS'
   def header(f):
    if f.read(6)!=b'\x93NUMPY':raise ValueError('not NPY')
    major,minor=f.read(2);size=2 if major==1 else 4
    n=struct.unpack('<H' if size==2 else '<I',f.read(size))[0]
    return ast.literal_eval(f.read(n).decode('latin1'))
   if suffix=='.npz':
    with zipfile.ZipFile(p) as z:
     for member in z.namelist():
      if not member.endswith('.npy'):continue
      with z.open(member) as f:h=header(f)
      arrays.append({'path':path,'member':member,'shape':str(h['shape']),'dtype':str(h['descr']),'candidate_state_key':bool(pattern.search(member))})
   else:
    with p.open('rb') as f:h=header(f)
    arrays.append({'path':path,'member':'','shape':str(h['shape']),'dtype':str(h['descr']),'candidate_state_key':bool(pattern.search(p.name))})
  elif suffix=='.csv':
   with p.open(encoding='utf-8-sig',errors='replace',newline='') as f:fields=next(csv.reader(f),[])
   inv['discovery']='CSV_HEADER'
   headers.append({'path':path,'fields':json.dumps(fields,ensure_ascii=False),'candidate_truth_fields':bool(pattern.search('|'.join(fields))),'named_panel':bool(pattern.search(p.name))})
  elif suffix in ['.py','.m','.yaml','.yml','.toml']:
   inv['discovery']='SOURCE_TEXT_USAGE'
   for n,text in enumerate(p.read_text(encoding='utf-8-sig',errors='replace').splitlines(),1):
    if pattern.search(text):uses.append({'path':path,'line':n,'text':text[:1600]})
  elif suffix=='.json':
   inv['discovery']='JSON_KEYS_AND_EXTERNAL_REFERENCES'
   obj=json.loads(p.read_text(encoding='utf-8-sig'))
   def walk(x,loc):
    if isinstance(x,dict):
     for k,v in x.items():
      here=loc+'/'+str(k)
      if pattern.search(str(k)):
       jhits.append({'path':path,'pointer':here,'type':type(v).__name__,'preview':json.dumps(v,ensure_ascii=False)[:1100]})
      walk(v,here)
    elif isinstance(x,list):
     for i,v in enumerate(x):walk(v,loc+'/'+str(i))
    elif isinstance(x,str) and re.match(r'^[dD]:[\\/]',x) and any(x.lower().endswith(e) for e in ['.npz','.npy','.csv','.json']):
     external.append({'source_manifest':path,'pointer':loc,'external_path':x})
   walk(obj,'')
 except Exception as e:errors.append({'path':path,'error':type(e).__name__+':'+str(e)});inv['discovery']='UNRESOLVED_PARSE'
 inventory.append(inv)
 if idx%1500==0:print('discovery',idx,flush=True)
for name,data in [('SNAPSHOT_INVENTORY.csv',inventory),('ARCHIVE_ARRAY_HEADERS.csv',arrays),('CSV_SCHEMA_DISCOVERY.csv',headers),('GENERATOR_USAGE_DISCOVERY.csv',uses),('JSON_STATE_KEY_DISCOVERY.csv',jhits)]:emit(name,data)
emit('EXTERNAL_ARCHIVE_REFERENCES.json',external);emit('DISCOVERY_ERRORS.json',errors)
emit('DISCOVERY_SUMMARY.json',{'snapshot':SNAP,'files':len(inventory),'array_members':len(arrays),'CSV_files':len(headers),'usage_hits':len(uses),'json_key_hits':len(jhits),'external_refs':len(external),'errors':errors,'new_panel':0,'new_samples':0})
print('DONE',len(inventory),'errors',len(errors),flush=True)
