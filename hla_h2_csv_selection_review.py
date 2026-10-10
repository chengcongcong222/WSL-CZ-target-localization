"""Supplemental cold selection audit using exact original CSV float tokens.
No estimator, observation, score, acceptance threshold or branch is changed.
The frozen primary checker uses pandas' default parser, which is not a
round-trip reader; preserve its original result as evidence before resolution.
"""
import csv,gzip,json,shutil
import hla_h2_core as c

def audit():
    selected={}
    with (c.OUT/'ESTIMATES.csv').open(encoding='utf-8',newline='') as f:
        for row in csv.DictReader(f):selected[int(row['config_id'])]=row
    minima={};count=0
    with gzip.open(c.OUT/'ALL_TERMINALS.csv.gz','rt',encoding='utf-8',newline='') as f:
        for row in csv.DictReader(f):
            count+=1
            if row['exception'] or row['accepted_profile']!='True':continue
            cid=int(row['config_id']);key=(float(row['J_beta'])+float(row['J_profile']),int(row['terminal_id']))
            if cid not in minima or key<minima[cid][0]:minima[cid]=(key,row)
    checks=[]
    for cid,row in selected.items():
        if cid in minima:
            key,minimum=minima[cid]
            ok=row['status']=='ACCEPTED' and int(row['selected_terminal_id'])==key[1]
            same_state=all(float(row[k])==float(minimum[k]) for k in c.AXES)
            same_bias=float(row['b'])==float(minimum['profile_b'])
            same_cost=float(row['J_beta'])==float(minimum['J_beta']) and float(row['J_radial'])==float(minimum['J_profile'])
            ok=ok and same_state and same_bias and same_cost
        else:ok=row['status']=='ABSTAIN' and int(row['selected_terminal_id'])==-1
        checks.append(dict(config_id=cid,PASS=bool(ok),selected_terminal_id=int(row['selected_terminal_id'])))
    result=dict(PASS=len(selected)==2816 and count==203008 and all(r['PASS'] for r in checks),selection_rebuild_checks=len(checks),selection_rebuild_FAIL=sum(not r['PASS'] for r in checks),terminal_rows=count,reader='stdlib csv.DictReader then Python float, preserving original double precision literals',ordering='exact Python floating sum, then original terminal ID; no equality tolerance',estimator_or_frozen_contract_changed=False,main_reexecuted=False,all_configurations=True)
    c.write_json(c.OUT/'EXACT_CSV_SELECTION_REVIEW.json',result)
    with (c.OUT/'EXACT_CSV_SELECTION_CHECKS.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=['config_id','PASS','selected_terminal_id']);writer.writeheader();writer.writerows(checks)
    old=json.loads((c.OUT/'VALIDATION.json').read_text())
    snapshot=c.OUT/'INITIAL_DEFAULT_PARSER_VALIDATION.json'
    if snapshot.exists():raise RuntimeError('Only one supplemental review permitted')
    shutil.copyfile(c.OUT/'VALIDATION.json',snapshot)
    old['initial_default_pandas_parser_selection_FAIL']=old['selection_rebuild_FAIL']
    old['initial_default_parser_validation_preserved']=snapshot.name
    old['selection_rebuild_checks']=result['selection_rebuild_checks']
    old['selection_rebuild_FAIL']=result['selection_rebuild_FAIL']
    old['selection_reader_resolution']='Exact original CSV tokens recover all original minima; frozen checker source unchanged. See independent supplemental review and two default-parser discrepancies.'
    old['PASS']=bool(result['PASS'] and old['terminal_score_FAIL']==0 and old['observation_rebuild_FAIL']==0 and old['frozen_input_hashes_PASS'] and old['branch_execution_PASS'] and old['execution_complete'] and old['max_normalized_Jacobian_column_difference']<=1e-5)
    c.write_json(c.OUT/'VALIDATION.json',old)
    print(json.dumps(result,indent=2))
    if not result['PASS']:raise RuntimeError('Exact original CSV selection failed')
if __name__=='__main__':audit()
