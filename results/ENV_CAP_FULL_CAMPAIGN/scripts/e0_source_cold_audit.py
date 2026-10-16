"""Independent E0 source cold audit, no importing profile preparation/science code."""
from pathlib import Path
import xarray as xr,gsw,numpy as np,csv,json,hashlib,time,resource
R=Path.cwd();O=R/'results/ENV_CAP_FULL_CAMPAIGN';W=Path('/mnt/d/ProjectStorage/WSL-CZ/ENV_CAP');out=O/'E0_SOURCE_INDEPENDENT_COLD_AUDIT.json'
if out.exists():raise RuntimeError('Source cold audit already done')
marker=W/'E0_SOURCE_COLD.SCIENCE_STARTED.json'
if marker.exists():raise RuntimeError('No retry')
marker.write_text(json.dumps({'code_SHA256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'purpose':'independent rawT/S->SSP and ETOPO source provenance cold reconstruction'}));start=time.perf_counter();checks=[];raw={}
def ck(check,ok,error='',tol='',scope='source_data_identity_or_reconstruction'):checks.append({'check':check,'PASS':bool(ok),'error':str(error),'tolerance':str(tol),'scope':scope})
for code in ['01','13','00']:
 pair=[]
 for typ in ['t','s']:
  p=O/'inputs/woa_raw'/f'woa23_decav91C0_{typ}{code}_04_source17N118E.nc';ds=xr.open_dataset(p,decode_times=False);v=ds[typ+'_an'].values[0];z=ds.depth.values
  allvalid=np.all(np.isfinite(v),axis=(1,2));mean=(v[:,0,0]+v[:,0,1]+v[:,1,0]+v[:,1,1])/4
  ck(p.name+'_axes',np.array_equal(ds.lat.values,[16.875,17.125]) and np.array_equal(ds.lon.values,[117.875,118.125]));pair.append((z.copy(),mean.copy(),allvalid.copy()));ds.close()
 ck('paired_depths_'+code,np.array_equal(pair[0][0],pair[1][0]));ok=pair[0][2]&pair[1][2];raw[code]=(pair[0][0][ok],pair[0][1][ok],pair[1][1][ok])
for name,codes in [('SCS17_JW',('01','13')),('SCS17_JA',('01','00')),('SCS17_ANNUAL',('00',))]:
 z,t,s=raw[codes[0]]
 if len(codes)>1:
  zd,td,sd=raw[codes[1]];keep=zd>1500;z=np.concatenate([z,zd[keep]]);t=np.concatenate([t,td[keep]]);s=np.concatenate([s,sd[keep]])
 pr=gsw.p_from_z(-z,17);absolute=gsw.SA_from_SP(s,pr,118,17);con=gsw.CT_from_t(absolute,t,pr);speed=gsw.sound_speed(absolute,con,pr)
 ref=list(csv.DictReader((O/'inputs'/f'{name}.csv').open()));rz=np.array([float(a['depth_m']) for a in ref]);rc=np.array([float(a['sound_speed_mps']) for a in ref]);cold=np.interp(rz,z,speed);err=float(np.max(abs(cold-rc)));ck(name+'_TEOS10_profile',err<=1e-8,err,1e-8);ck(name+'_no_deep_extrapolation',z[-1]>=rz[-1] and rz[-1]==3909.75)
 # Converted originalT/S archive columns were rounded by original exporters; compare new preparation only.
tile=np.load(W/'sources/m45/etopo_tile_C1_LCE_NSCS_17N118E.npy');lat=np.unique(tile[:,0]);lon=np.unique(tile[:,1]);ys=[lat[lat<17][-1],lat[lat>17][0]];xs=[lon[lon<118][-1],lon[lon>118][0]]
a=[]
for y in ys:
 for x in xs:a.append(float(tile[(tile[:,0]==y)&(tile[:,1]==x),2][0]))
wy=(17-ys[0])/(ys[1]-ys[0]);wx=(118-xs[0])/(xs[1]-xs[0]);depth=-((1-wy)*(1-wx)*a[0]+(1-wy)*wx*a[1]+wy*(1-wx)*a[2]+wy*wx*a[3]);ck('ETOPO_four_pixel_source_point',abs(depth-3909.75)<1e-7,depth-3909.75,1e-7,'bathymetry_point_only_not_track_or_sediment')
result={'PASS':all(a['PASS'] for a in checks),'checks':len(checks),'FAIL_count':sum(not a['PASS'] for a in checks),'wall_s':time.perf_counter()-start,'single_process_peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'gsw_version':gsw.__version__,'rows':checks,'scope':'Original raw T/S and actualETOPOtile, independent formula reconstruction; not firstCZ or modal admission'}
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='rows'}));raise SystemExit(0 if result['PASS'] else 2)
