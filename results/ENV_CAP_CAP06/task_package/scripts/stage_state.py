#!/usr/bin/env python3
"""CAP06 lightweight append-only stage register (metadata tool, NOT a scientific Gate)."""
import argparse
import datetime
import json
import os
from pathlib import Path

STAGES=("P0","P1","P2","P3","P4")
FINAL=("PASS","LIMITED","FAIL","NOT_EVALUATED","STOP")


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_atomic(path,data):
    p=Path(path)
    p.parent.mkdir(parents=True,exist_ok=True)
    t=p.with_name(p.name+".tmp")
    t.write_text(json.dumps(data,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    os.replace(t,p)


def event(d,action,**kwargs):
    d.setdefault("history",[]).append({"utc":utc(),"action":action,**kwargs})


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--state",required=True)
    sp=ap.add_subparsers(dest="action",required=True)
    init=sp.add_parser("init");init.add_argument("--template",required=True)
    st=sp.add_parser("start");st.add_argument("phase",choices=STAGES)
    fin=sp.add_parser("finish");fin.add_argument("phase",choices=STAGES);fin.add_argument("status",choices=FINAL);fin.add_argument("--note",default="")
    uni=sp.add_parser("mark-science");uni.add_argument("unit");uni.add_argument("sha256");uni.add_argument("--note",default="")
    ck=sp.add_parser("checkpoint");ck.add_argument("label");ck.add_argument("sha")
    stop=sp.add_parser("stop");stop.add_argument("--reason",required=True)
    sp.add_parser("show")
    a=ap.parse_args(argv)
    p=Path(a.state)
    if a.action=="init":
        if p.exists():raise SystemExit("Refuse overwrite existing RUN_STATE")
        d=read(a.template);event(d,"init");write_atomic(p,d)
        print(json.dumps({"state":str(p),"status":"initialized"}));return
    d=read(p)
    if a.action=="show":print(json.dumps(d,ensure_ascii=False,indent=2));return
    if d.get("STOP"):raise SystemExit("STOP=true: refusing modifications")
    if a.action=="start":
        if d["phases"][a.phase]!="NOT_STARTED":raise SystemExit("phase already started/completed")
        # Branch independence: P1 and P2 can begin in any order once P0 is complete.
        if a.phase!="P0" and d["phases"]["P0"] not in FINAL:raise SystemExit("P0 must be terminal first")
        if a.phase=="P3" and (d["phases"]["P1"] not in FINAL or d["phases"]["P2"] not in FINAL):raise SystemExit("P1 and P2 must be terminal before P3")
        if a.phase=="P4" and d["phases"]["P3"] not in FINAL:raise SystemExit("P3 must be terminal before P4")
        d["phases"][a.phase]="RUNNING";event(d,"phase_start",phase=a.phase)
    elif a.action=="finish":
        if d["phases"][a.phase]!="RUNNING":raise SystemExit("phase must be RUNNING")
        d["phases"][a.phase]=a.status;event(d,"phase_finish",phase=a.phase,status=a.status,note=a.note)
    elif a.action=="mark-science":
        if not d.get("frozen_design_sha"):raise SystemExit("Commit A/frozen_design_sha missing")
        if any(z["unit"]==a.unit for z in d.setdefault("SCIENCE_STARTED_units",[])):raise SystemExit("unit is already frozen/started")
        d["SCIENCE_STARTED_units"].append({"unit":a.unit,"input_config_code_sha256":a.sha256,"utc":utc(),"note":a.note});event(d,"science_started",unit=a.unit)
    elif a.action=="checkpoint":
        if any(x["label"]==a.label for x in d.setdefault("checkpoints_pushed",[])):raise SystemExit("checkpoint already recorded")
        d["checkpoints_pushed"].append({"label":a.label,"sha":a.sha,"utc":utc()});d["latest_git_sha"]=a.sha;event(d,"checkpoint",label=a.label,sha=a.sha)
    elif a.action=="stop":
        if any(v in ("NOT_STARTED","RUNNING") for v in d["phases"].values()):raise SystemExit("all phases must be terminal before STOP")
        d["STOP"]=True;event(d,"stop",reason=a.reason)
    write_atomic(p,d)
    print(json.dumps({"action":a.action,"phase_status":d["phases"],"STOP":d["STOP"]},ensure_ascii=False))

if __name__=="__main__":main()
