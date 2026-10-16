#!/usr/bin/env python3
"""Tiny resumable state CLI for multi-session Codex research; no scientific execution.

Usage:
  python stage_state.py init /path/to/results/ENV_CAP_FULL_CAMPAIGN
  python stage_state.py status /path/to/results/ENV_CAP_FULL_CAMPAIGN
  python stage_state.py start ... E1 "source inventory completed"
  python stage_state.py finish ... E1 LIMITED "one CZ2, another uncertified"
  python stage_state.py halt ... "need approval for new area"

State changes keep append-only history; invalid transitions fail without modifying state.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse, json, sys

PHASES=[f'E{i}' for i in range(6)]
END={'PASS','LIMITED','FAIL','NOT_EVALUATED'}

def now():return datetime.now(timezone.utc).replace(microsecond=0).isoformat()

def statefile(folder):return Path(folder)/'RUN_STATE.json'

def write(path,v):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(v,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    tmp.replace(path)

def main():
    pa=argparse.ArgumentParser()
    pa.add_argument('action',choices=['init','status','start','finish','halt'])
    pa.add_argument('directory')
    pa.add_argument('phase_or_note',nargs='?')
    pa.add_argument('status_or_note',nargs='?')
    pa.add_argument('note',nargs='?')
    x=pa.parse_args();src=statefile(x.directory)
    if x.action=='init':
        if src.exists():print('RUN_STATE already exists; not overwriting');return 2
        data={'stage':'ENV_CAP_20261016','parent_SHA':'a9a80191ee4d0f0cecd9db296198606b63a68080','research_branch':'research/env-cap-20261016','phases':{k:'PENDING' for k in PHASES},'STOP':False,'history':[{'time_utc':now(),'action':'init'}]}
        write(src,data);print(json.dumps(data,ensure_ascii=False,indent=2));return 0
    if not src.exists():print('RUN_STATE missing; use init');return 2
    data=json.loads(src.read_text(encoding='utf-8'))
    if x.action=='status':print(json.dumps(data,ensure_ascii=False,indent=2));return 0
    if x.action=='halt':
        note=x.phase_or_note or 'unspecified blocker';data['STOP']=True;data['history'].append({'time_utc':now(),'action':'halt','note':note});write(src,data);print('HALT '+note);return 0
    phase=x.phase_or_note
    if phase not in PHASES:print('Unknown phase');return 2
    current=data['phases'][phase]
    if data.get('STOP'):print('RUN_STATE is STOP; new authorization required');return 2
    if x.action=='start':
        for p in PHASES[:PHASES.index(phase)]:
            if data['phases'][p] not in END:print('Prior phase not terminal: '+p);return 2
        if current!='PENDING':print('Start refused; current '+current);return 2
        data['phases'][phase]='RUNNING';comment=x.status_or_note or ''
    else:
        st=x.status_or_note;comment=x.note or ''
        if current!='RUNNING' or st not in END:
            print('Finish requires RUNNING phase and terminal status PASS/LIMITED/FAIL/NOT_EVALUATED');return 2
        data['phases'][phase]=st
        if phase=='E5':data['STOP']=True
    data['history'].append({'time_utc':now(),'action':x.action,'phase':phase,'new_status':data['phases'][phase],'note':comment})
    write(src,data);print(json.dumps({'phase':phase,'status':data['phases'][phase],'STOP':data['STOP']},ensure_ascii=False));return 0

if __name__=='__main__':sys.exit(main())
