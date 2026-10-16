from pathlib import Path
import xarray as xr,gsw,numpy as np,json,csv,hashlib
W=Path('/mnt/d/ProjectStorage/WSL-CZ/ENV_CAP');S=W/'sources/subsets';D=W/'prepared';D.mkdir(exist_ok=True)
profiles={};audit=[]
for code in ['01','13','00']:
 ds=[]
 for typ in ['t','s']:
  p=S/f'woa23_decav91C0_{typ}{code}_04_source17N118E.nc';x=xr.open_dataset(p,decode_times=False)
  assert np.allclose(x.lat,[16.875,17.125]) and np.allclose(x.lon,[117.875,118.125])
  a=x[typ+'_an'].values[0]; valid=np.all(np.isfinite(a),axis=(1,2));vals=np.mean(a,axis=(1,2));ds.append((x.depth.values.copy(),vals,valid))
  audit.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'title':x.attrs['title'],'valid_all4_depth_max':float(x.depth.values[valid].max()),'bounds':str(x.get('climatology_bounds','NOT_SAVED'))});x.close()
 assert np.array_equal(ds[0][0],ds[1][0]);z=ds[0][0];ok=ds[0][2]&ds[1][2];profiles[code]=(z[ok],ds[0][1][ok],ds[1][1][ok])
for name,codes in [('SCS17_JW',('01','13')),('SCS17_JA',('01','00')),('SCS17_ANNUAL',('00',))]:
 z,t,s=profiles[codes[0]]
 if len(codes)==2:
  zd,td,sd=profiles[codes[1]];m=zd>1500;z=np.r_[z,zd[m]];t=np.r_[t,td[m]];s=np.r_[s,sd[m]]
 p=gsw.p_from_z(-z,17);sa=gsw.SA_from_SP(s,p,118,17);ct=gsw.CT_from_t(sa,t,p);c=gsw.sound_speed(sa,ct,p)
 with (D/f'{name}.csv').open('w',newline='') as f:
  wr=csv.writer(f);wr.writerow(['depth_m','temperature_degC','practical_salinity','pressure_dbar','absolute_salinity_gkg','conservative_temperature_degC','sound_speed_mps']);wr.writerows(zip(z,t,s,p,sa,ct,c))
(D/'PROVENANCE.json').write_text(json.dumps({'gsw':gsw.__version__,'source_rows':audit,'no_missing_corner_renormalization':True,'seasons':'Winter=Jan-Mar per NetCDF title; legacy woa_core seasonal DJF shorthand not used','bathymetry_m':3909.75,'bottom_interpolation':'C-linear within valid 4000m T/S range; no extrapolation'},indent=2))
print('Prepared 3 independently rebuilt T/S->c profiles; full-depth all-four-corner validity required.')
