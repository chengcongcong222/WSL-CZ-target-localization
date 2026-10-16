from pathlib import Path
import json,urllib.request,datetime
import netCDF4 as nc
W=Path('/mnt/d/ProjectStorage/WSL-CZ/ENV_CAP'); D=W/'sources/independent_25N140E';D.mkdir(exist_ok=True)
logs=[]
for typ in ['t','s']:
 url=f'https://thredds.rda.ucar.edu/thredds/dodsC/files/d285000/woa23_netcdf/woa23_decav91C0_{typ}00_04.nc.ascii?depth,lat[459:1:460],lon[1279:1:1280],{typ}_an[0:1:0][0:1:101][459:1:460][1279:1:1280]'
 try:
  with urllib.request.urlopen(url,timeout=45) as f:b=f.read(2000000)
  (D/f'{typ}00.ascii').write_bytes(b);logs.append({'url':url,'status':'DOWNLOADED','bytes':len(b)})
 except Exception as e:logs.append({'url':url,'status':'UNAVAILABLE','error':repr(e)})
url='https://coastwatch.pfeg.noaa.gov/erddap/griddap/ETOPO_2022_v1_15s.csv?z[(24.99):(25.01)][(139.99):(140.01)]'
try:
 with urllib.request.urlopen(url,timeout=45) as f:b=f.read(2000000)
 (D/'bathy.csv').write_bytes(b); logs.append({'url':url,'status':'DOWNLOADED','bytes':len(b)})
except Exception as e:logs.append({'url':url,'status':'UNAVAILABLE','error':repr(e)})
(D/'ACCESS.json').write_text(json.dumps(logs,indent=2));print(json.dumps(logs,indent=2))
