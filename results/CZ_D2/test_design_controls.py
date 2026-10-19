"""CZ-D2 design/arithmetical controls only. No random samples or optimizer calls."""
from pathlib import Path
import csv,hashlib,json,math,unittest,sys
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
D=json.loads((OUT/'DESIGN_FREEZE.json').read_text(encoding='utf-8'))

def rank_quantile(values,q):
    if not 0 < q <= 1: raise ValueError('quantile')
    if not values: return None
    if any(math.isnan(v) or v < 0 for v in values): raise ValueError('invalid error')
    return sorted(values)[math.ceil(q*len(values))-1]

def circular_error(a,b): return abs((a-b+180)%360-180)

def paired(a,b):
    if math.isnan(a) or math.isnan(b): raise ValueError('NaN')
    if not math.isfinite(a) and not math.isfinite(b): return 'both_failed'
    if not math.isfinite(a): return 'baseline_failed_candidate_valid'
    if not math.isfinite(b): return 'candidate_failed_baseline_valid'
    return 'better' if b<a else 'worse' if b>a else 'equal'

def old_neighbor(h,old,widths):
    indices=[i for i,x in enumerate(old) if x is not None]
    return all((circular_error(h[i],old[i]) if i in [1,3] else abs(h[i]-old[i]))<=widths[i] for i in indices)

class DesignControls(unittest.TestCase):
    def test_01_authorization(self):
        self.assertTrue(D['authorization']['Commit_A_design_only'])
        for key in ['new_scientific_samples','formal_execution','Commit_B']:self.assertFalse(D['authorization'][key])
        self.assertFalse(D['panel']['materialized'])
        self.assertTrue(D['STOP']);self.assertEqual(D['R4_percent'],0)
        for name in ['PANEL.csv','OBSERVATIONS.npz','RUN_ERRORS.csv','ESTIMATOR_OUTPUTS.csv']:
            self.assertFalse((OUT/name).exists())
    def test_02_matrix(self):
        self.assertEqual(set(D['geometry']),{'B0','B3'});self.assertEqual(set(D['quality']),{'Q1','Q2'})
        self.assertEqual(D['bias_support_multiplier'],{'P0':1,'P1':2})
        self.assertEqual(12*2*250,D['random']['base_innovations'])
        self.assertEqual(6000*2*2,D['random']['public_records'])
        self.assertEqual(6000*2*2*2,D['random']['method_configurations'])
        self.assertEqual(12*2*2*2,D['information']['nominal_geometries'])
    def test_03_inheritance(self):
        old=json.loads((ROOT/'results/R4_T3_W1_WIDE_BASELINE_DEMO/T3_W1_DESIGN_FREEZE.json').read_text(encoding='utf-8-sig'))
        self.assertEqual(D['quality'],old['qualities'])
        self.assertEqual(D['state_domain'],old['state_domain'])
        for k in ['platform','coordinate','navigation']:self.assertEqual(D[k],old[k])
        for k in old['estimator']:self.assertEqual(D['estimator'][k],old['estimator'][k])
        self.assertNotEqual(D['random']['seed'],old['seed'])
        self.assertEqual(D['information']['scale'],old['information_scale'])
        self.assertEqual(D['budget'],old['budget'])
    def test_04_support(self):
        for q,expected in [('Q1',(0.1,0.05)),('Q2',(0.15,0.075))]:
            x=D['quality'][q];before=dict(x);est=dict(x)
            for k in ['common_bound_deg','half_diff_bound_deg']:est[k]*=2
            self.assertEqual(x,before)
            self.assertEqual((est['common_bound_deg'],est['half_diff_bound_deg']),expected)
            self.assertEqual(est['bearing_sigma_deg'],x['bearing_sigma_deg']);self.assertEqual(est['nav_sigma_m'],x['nav_sigma_m'])
    def test_05_inf_quantiles(self):
        self.assertEqual(rank_quantile([1,2,3,4],.5),2)
        self.assertEqual(rank_quantile(list(range(95))+[math.inf]*5,.95),94)
        self.assertEqual(rank_quantile(list(range(94))+[math.inf]*6,.95),math.inf)
        self.assertEqual(rank_quantile([math.inf]*4,.5),math.inf)
        self.assertIsNone(rank_quantile([],.95))
        with self.assertRaises(ValueError):rank_quantile([math.nan],.5)
    def test_06_angles(self):
        self.assertEqual(circular_error(179,-179),2)
        self.assertEqual(circular_error(359,1),2)
        self.assertEqual(circular_error(-180,180),0)
    def test_07_paired_failure(self):
        self.assertEqual(paired(math.inf,math.inf),'both_failed')
        self.assertEqual(paired(math.inf,1),'baseline_failed_candidate_valid')
        self.assertEqual(paired(1,math.inf),'candidate_failed_baseline_valid')
        self.assertEqual(paired(2,1),'better');self.assertEqual(paired(1,2),'worse');self.assertEqual(paired(1,1),'equal')
        with self.assertRaises(ValueError):paired(math.nan,1)
    def test_08_equal_weight_and_pool_are_different(self):
        a=[0,0,0,100];b=[10,10,10,10]
        self.assertEqual((rank_quantile(a,.5)+rank_quantile(b,.5))/2,5)
        self.assertEqual(rank_quantile(a+b,.5),10)
    def test_09_units_and_exclusion(self):
        # Deliberately synthetic dimensionless fixtures, not physical target samples.
        self.assertTrue(old_neighbor([1,359,3,1],[1,-1,3,1],[0,0,0,0]))
        self.assertFalse(old_neighbor([1,359,4,1],[1,-1,3,1],[0,0,.1,0]))
        self.assertTrue(old_neighbor([1,359,4,1],[1,-1,None,None],[0,0,0,0]))
        self.assertAlmostEqual(math.hypot(3,15),15.297058540778355)
        self.assertNotEqual(1e-6*0.1,1e-6*0.2) # P0/P1 physical bound tolerance differs.
    def test_10_source_hashes(self):
        s=json.loads((OUT/'SOURCE_IDENTITIES.json').read_text(encoding='utf-8'))
        for row in s['sources']:
            with self.subTest(path=row['path']):
                data=(ROOT/row['path']).read_bytes()
                self.assertEqual(len(data),row['bytes']);self.assertEqual(hashlib.sha256(data).hexdigest(),row['sha256'])
        self.assertEqual(hashlib.sha256((OUT/'HISTORICAL_STATE_EXCLUSIONS.csv').read_bytes()).hexdigest(),s['historical_registry_sha256'])
    def test_11_registry_reconstruction(self):
        rr=list(csv.DictReader((OUT/'HISTORICAL_STATE_EXCLUSIONS.csv').open(encoding='utf-8')))
        self.assertEqual(len(rr),216)
        for row in rr:
            p=ROOT/row['source'];i=row['record'];h=[float(row[k]) if row[k] else None for k in ['r_km','theta_deg','v_mps','psi_deg']]
            if p.suffix=='.csv':
                records=list(csv.DictReader(p.open(encoding='utf-8-sig')));old=records[int(i)]
                names=[['r0_km','r_km','truth_r0_km','truth_r_km'],['theta0_deg','theta_deg','truth_theta0_deg','truth_theta_deg'],['v_mps','truth_v'],['psi_deg','truth_psi_deg']]
                wanted=[float(next(old[k] for k in names[j] if k in old)) if h[j] is not None else None for j in range(4)]
            elif 'CAP06' in str(p):
                old=next(x['state'] for x in json.loads(p.read_text())['states'] if x['id']==i)
                wanted=[old[0]/1000,old[1]*180/math.pi,old[3],old[4]*180/math.pi]
            else:wanted=next(x['h'] for x in json.loads(p.read_text(encoding='utf-8-sig'))['generating_cases'] if x['id']==i)
            for a,b in zip(h,wanted):
                if a is None:self.assertIsNone(b)
                else:self.assertAlmostEqual(a,b,places=11)
    def test_12_no_science_entry(self):
        self.assertEqual(D['budget']['main_calls'],1);self.assertEqual(D['budget']['cold_calls'],1);self.assertEqual(D['budget']['retries'],0)
        self.assertEqual(D['budget']['main_wall_s']+D['budget']['cold_wall_s'],7200)
        self.assertEqual(D['information']['rows_with_known_bias_control'],192)
        self.assertEqual(D['depth'],'NOT_EVALUATED')
        self.assertTrue(all(v==0 for v in D['forbidden_calls'].values()))

if __name__=='__main__':
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(DesignControls))
    output={'scope':'COMMIT_A_DESIGN_AND_ARITHMETIC_ONLY_NOT_SCIENTIFIC_COLD_AUDIT','tests_run':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),'PASS':result.wasSuccessful(),'new_panel_samples':0,'MC':0,'optimizer_calls':0,'Fisher_evaluations':0,'propagation_calls':0,'source_identity_records':25,'historical_exclusion_records':216,'formal_execution_authorized':False,'STOP':True}
    (OUT/'DESIGN_CHECKS.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    sys.exit(0 if result.wasSuccessful() else 1)
