from common import *
from frontend import extract
import sys
# Hard deny private files in this receiver process; labels opened only later by evaluation.
def guard(event,args):
 if event=='open' and isinstance(args[0],(str,bytes)):
  p=str(args[0]).lower()
  if 'truth_private' in p or 'evaluation_private' in p or 'p0_frozen_panel' in p:raise PermissionError('Receiver attempted private input')
sys.addaudithook(guard)
def main():
 out=RAW/'p2/products';out.mkdir()
 rows=[]
 for ii,p in enumerate(sorted((RAW/'p2/receiver_public').glob('*.npz'))):
  with np.load(p) as d:
   assert set(d.files)=={'samples','background','calibration','calibration_angles_deg','timestamps','reported_positions','nominal_heading','nav_sigma_axis_m','fs_hz','LO_hz'}
   for calibrated in [True,False]:
    products=[extract(y,bg,d['calibration'],d['calibration_angles_deg'],calibrated) for y,bg in zip(d['samples'],d['background'])]
    data=dict(token=p.stem,calibrated=calibrated,source_SHA=sha(p),epochs=products,timestamps=d['timestamps'].tolist(),reported_positions=d['reported_positions'].tolist(),nominal_heading=d['nominal_heading'].tolist(),nav_sigma_axis_m=float(d['nav_sigma_axis_m']),unprovided_orientation='FREE' if not calibrated else 'CALIBRATION_DATA_ESTIMATED',time_offset_bound_s=.25)
    path=out/(p.stem+('_CAL' if calibrated else '_UNCAL')+'.json');path.write_text(json.dumps(data,allow_nan=False)+'\n')
    rows.append(dict(path=str(path),sha256=sha(path),source_SHA=data['source_SHA'],private_reads=0))
  if ii%100==0:print('extracted',ii,flush=True)
 tab('P2_PRODUCT_INDEX.csv',rows)
if __name__=='__main__':main()
