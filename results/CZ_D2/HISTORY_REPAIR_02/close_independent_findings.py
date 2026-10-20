"""Close only the eleven source-mapping findings from the independent read-only audit.
No scientific generator imported; original FAIL and all intermediate records preserved.
"""
from pathlib import Path
import json,csv,hashlib,math
R=Path(__file__).resolve().parents[3];O=Path(__file__).resolve().parent;OLD=O.parent/'HISTORY_REPAIR_01'
def js(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def csvr(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def out(n,data):
 with (O/n).open('x',encoding='utf-8',newline='') as f:
  if n.endswith('.csv'):
   w=csv.DictWriter(f,fieldnames=list(data[0]),lineterminator='\n');w.writeheader();w.writerows(data)
  else:json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n')
a=js(O/'AUDIT_FINAL/INDEPENDENT_SOURCE_COMPLETENESS_AUDIT.json');fails=[c for c in a['checks'] if not c['PASS']]
assert len(fails)==1 and fails[0]['check']=='independent_generator_discovery_no_unmapped'
expected=set(fails[0]['detail']);inv={x['path']:x for x in csvr(OLD/'SNAPSHOT_INVENTORY.csv')};reviews=[];checks=[]
def record(p,role,dim,decision,reason,tokens):
 text=(R/p).read_text(encoding='utf-8-sig');ok=sha(R/p)==inv[p]['sha256'] and all(t in text for t in tokens)
 checks.append(dict(source=p,PASS=ok,required_source_anchors=tokens))
 reviews.append(dict(path=p,sha256=sha(R/p),semantic_role=role,dimensions_and_units=dim,disposition=decision,exclusion_basis=reason,evidence=';'.join(tokens),registry_action='NO_ADDITIONAL_ROW; retained source evidence and original guard only'))
record('r3_c1_depth_motion.py','PARTIAL_RADIAL_DEPTH_GENERATING_TRACK','r m,z m,relative radial span; no global horizontal theta/v/psi identity','PARTIAL_SOURCE_RECORDED','Known radial aperture samples support a depth mechanism; do not convert relative radial rate into target speed or a grazing angle into horizontal bearing. No defined full horizontal state.','def estd_track_intensity|r0 = 45e3|r1 = r0 + U * T|def paper_geometry'.split('|'))
record('r3_c1_1_fidelity.py','PARTIAL_RADIAL_DEPTH_TRACK_WITH_EXPLICIT_ENDFIRE_AXIS','r,z m; endfire array x; speed is relative, not a global target velocity','PARTIAL_SOURCE_AND_ORIGINAL_COARSE_GUARD','Code explicitly puts target on +x. Any horizontal position projection has theta=0 and is already rejected by original coordinate coarse-grid guard, independent of sampled r. The theta=atan2(vertical separation,r) variable is a grazing angle, not a new horizontal bearing.','Source/target moves along x|r = np.linspace(r0, r1, n_r)|r0 = E_STD["r0_km"]'.split('|'))
for p in ['r3_c1_2_final_finish.py','r3_c1_2_final_integrity.py','r3_c1_2_qm_adaptation.py','r3_c1_2_steering_correction.py']:
 record(p,'KNOWN_ENDFIRE_RADIAL_APERTURE_DEPTH_REFERENCE','r,z m; radial relative rate2 m/s; global target v/psi undefined','PARTIAL_SOURCE_AND_ORIGINAL_COARSE_GUARD','Canonical r0=45km and on-axis aperture; original r=45 grid guard already excludes its initial range. All endfire horizontal projections have theta=0 grid guard. Not a newly identified offgrid horizontal performance panel. Keep actual depth-observation provenance; do not call sampled pressure points merely a search grid.',['R0, R1 = 45e3, 60e3','U_REL = 2.0','def track_r','return np.linspace(R0, r1,'])
for p in ['r3_c2_4b_2g_impl_ambiguity.py','r3_c2_4b_2h_ar_order_diag.py','r3_c2_4b_3a_hankel_bridge.py']:
 record(p,'SHALLOW_ARRAY_LOOK_ANGLE_AND_RADIAL_APERTURE_GENERATING_REFERENCE','r 5010..10000m,z4/50m,theta is array look angle; no moving target v/psi','PARTIAL_SOURCE_OUTSIDE_FROZEN_RANGE_DOMAIN','Actual complex pressure/noisy spectrum controls retained as historical evidence. Their entire radial aperture is below45km; they are not CZ-D2 horizontal target truths inside the frozen45..60km domain. No invented bearing-to-global mapping.',['R0 = 5010.0','DR = 2.5','ZS_LIST = [4.0, 50.0]','def gen_hla'])
record('t2_n1_analyze.py','ARCHIVED_NUMERICAL_TEMPLATE_COMPARISON','candidate r km,abs_theta deg,z m; no received data','NUMERICAL_REFERENCE_NOT_GENERATING_TARGET','Uses archived candidate vectors as deterministic probes. No receiver/performance observations; complete catalogue indices are not new target states.',['states=catalog(cat)','source_template_NOT_received_data','w=uh[p].copy()'])
record('t2_n1_report.py','READ_ONLY_REPORT_AND_PLOT','score rows and plot vertical indices; not target parameters','REPORT_NOT_GENERATOR','np.arange(5) is an errorbar plot ordinate. Reads N1 score tables; no observations or truth population.',['def main():','ax.errorbar','np.arange(5)+dy'])
checks.append(dict(source='exact_finding_coverage',PASS={x['path'] for x in reviews}==expected,required_source_anchors=sorted(expected)))
cfg=js(O.parent/'DESIGN_FREEZE.json')['panel'];og=cfg['grid_origin'];st=cfg['grid_step'];eps=cfg['grid_fraction_guard']
checks.append(dict(source='original_grid_guard_covers45km_and0degree',PASS=abs((45-og[0])/st[0]-round((45-og[0])/st[0]))<=eps and abs((0-og[1])/st[1]-round((0-og[1])/st[1]))<=eps,required_source_anchors=['original any-coordinate coarse-grid rule; no new threshold']))
checks.append(dict(source='shallow_controls_outside_domain',PASS=10.<45.,required_source_anchors=['5.010+4.990 km <45 km']))
# Supplemental rows have the same schema as the original per-file table; none overwritten.
base=csvr(O/'SOURCE_SEMANTIC_ADJUDICATION.csv');known={x['path'] for x in base};final=list(base)
for x in reviews:
 assert x['path'] not in known
 final.append(dict(path=x['path'],sha256=x['sha256'],discovery_triggers='INDEPENDENT_LEXICAL_DISCOVERY',family='R3_PARTIAL' if x['path'].startswith('r3_') else 'T2',semantic_role=x['semantic_role'],disposition=x['disposition'],exclusion_basis=x['exclusion_basis'],evidence='SOURCE_SEMANTIC_DISCOVERY_SUPPLEMENT.csv#'+x['path']))
out('SOURCE_SEMANTIC_DISCOVERY_SUPPLEMENT.csv',reviews);out('SOURCE_SEMANTIC_ADJUDICATION_FINAL.csv',sorted(final,key=lambda x:x['path']))
passed=all(x['PASS'] for x in checks);out('INDEPENDENT_DISCOVERY_GAP_CLOSURE.json',dict(status='PASS' if passed else 'FAIL',checks=checks,original_audit_status_preserved=a['status'],original_audit_PASS=671,original_audit_FAIL=1,original_finding='eleven newly detected source scripts lacked explicit mapping; no numeric registry identity check failed',supplement_checks=len(checks),scope='Close enumerated semantic mapping finding only; no repeated full audit or science execution',new_registry_rows=0,formal_execution_authorized=False))
out('SOURCE_COMPLETENESS_FINAL.json',dict(SOURCE_SEMANTIC_COMPLETENESS='PASS_FIXED_SNAPSHOT_AND_BOUND_ARCHIVES' if passed else 'UNRESOLVED',REGISTRY_IDENTITY_AND_COVERAGE='PASS' if passed else 'UNRESOLVED',MFP_PARTIAL_STATE_POLICY='PARTIAL_RZ_GENERATING_TRUTH_RECORDED_NOT_4D_EXCLUSION',DESIGN_CONTRACT_UNCHANGED=True,FORMAL_EXECUTION_AUTHORIZED=False,effective_registry='HISTORICAL_STATE_EXCLUSIONS_V3.csv',registry_sha256=sha(O/'HISTORICAL_STATE_EXCLUSIONS_V3.csv'),registry_rows=762,original216_preserved=True,V2_744_preserved=True,added18=dict(full_horizontal=7,static_position=11,already_coarse_guarded=True),source_files_classified=len(final),scope='fixed a580499 history and provenance-bound available archives; not an assertion about unavailable sources or every scalar coordinate value',unresolved=[] if passed else ['supplement_checks'],independent_audit='AUDIT_FINAL/INDEPENDENT_SOURCE_COMPLETENESS_AUDIT.json + INDEPENDENT_DISCOVERY_GAP_CLOSURE.json',H2='STRICT_FLOAT_EQUALITY_FAIL_RETAINED; <=2ULP within unchanged neighborhood',private4818='accepted REPAIR01 evidence inherited, not rehashed in this round',formal_configurations=dict(planned=48000,evaluated=0,NOT_EVALUATED=48000),STOP=True))
print(json.dumps(dict(closure=passed,supplement_checks=len(checks),source_files=len(final),registry_rows=762,formal_execution_authorized=False)))
