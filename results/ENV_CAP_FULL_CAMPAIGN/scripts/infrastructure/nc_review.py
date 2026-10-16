from pathlib import Path
import xarray as xr,gsw,json,numpy as np
W=Path('/mnt/d/ProjectStorage/WSL-CZ/ENV_CAP');P=W/'sources/subsets'
rows=[]
for p in sorted(P.glob('*source17N118E.nc')):
 ds=xr.open_dataset(p,decode_times=False)
 r={'name':p.name,'sizes':dict(ds.sizes),'depth_firstlast':[float(ds.depth[0]),float(ds.depth[-1])],'lat':ds.lat.values.tolist(),'lon':ds.lon.values.tolist(),'title':ds.attrs.get('title'),'history':ds.attrs.get('history'),'time':ds.time.values.tolist(),'time_attributes':dict(ds.time.attrs)}
 typ=p.name.split('_')[2][0]
 ar=ds[typ+'_an'].values
 r['valid_depth_max']=float(ds.depth[np.any(np.isfinite(ar[0]),axis=(1,2))].max())
 rows.append(r);ds.close()
(W/'sources/NETCDF_METADATA_REVIEW.json').write_text(json.dumps({'gsw':gsw.__version__,'rows':rows},indent=2,default=str));print(json.dumps(rows,indent=2,default=str))

