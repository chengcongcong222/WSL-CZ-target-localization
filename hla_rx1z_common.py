"""RX1Z utilities, no historical imports or import-time writes."""
from pathlib import Path
import json,csv,hashlib,datetime,time
import numpy as np
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_SINGLE_HLA_CARRIER_REFERENCE_CONTRAST_PILOT'
PARENT='41463d9aa71da1ae6ebed9b025de22cb3d679e68'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def rows(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig',newline='')))
def table(p,x):
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(x[0]));w.writeheader();w.writerows(x)
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def guard():
 import psutil
 start=read(OUT/'EXECUTION_STARTED.json');elapsed=(datetime.datetime.now(datetime.timezone.utc)-datetime.datetime.fromisoformat(start['time_utc'])).total_seconds()
 rss=psutil.Process().memory_info().rss
 if elapsed>1800:raise TimeoutError('RX1Z total hard wall budget')
 if rss>4*1024**3:raise MemoryError('RX1Z RSS hard limit')
 return elapsed,rss
