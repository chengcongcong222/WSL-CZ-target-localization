#!/usr/bin/env python3
"""Simple dependency, SHA format, resource composition and source index checks. Does NOT authorize science."""
import csv
import json
from pathlib import Path
import re
import sys


def validate(folder):
    folder=Path(folder)
    p=json.loads((folder/"12_STAGE_SPEC.json").read_text(encoding="utf-8"))
    assert re.fullmatch(r"[a-f0-9]{40}",p["parent_sha"])
    assert p["plan_status"]=="DESIGN_SEED_NOT_SCIENTIFIC_EXECUTION"
    assert p["hard_constraints"]["first_CZ_required"] and p["hard_constraints"]["no_partner_staff_or_new_proprietary_data"]
    assert p["hard_constraints"]["original_R4_management_percent"]==0
    ids=[x["id"] for x in p["phases"]];assert ids==["P0","P1","P2","P3","P4"]
    deps={x["id"]:set(x["deps"]) for x in p["phases"]}
    assert deps["P1"]=={"P0"} and deps["P2"]=={"P0"}
    assert deps["P3"]=={"P1","P2"} and deps["P4"]=={"P3"}
    r={x["id"]:x for x in p["resource_configs"]}
    assert [r[k]["total_channels"] for k in ["S0","S1","S2","S3"]]==[8,12,16,20]
    budgets=p["global_budgets"]
    assert budgets["external_acoustic_solver_calls"]==48 and budgets["MC_configurations"]==20000
    for x in p["phases"]:assert x["status"]=="NOT_STARTED"
    with (folder/"11_HISTORICAL_EVIDENCE_INDEX.csv").open(encoding="utf-8",newline="") as f:
        rows=list(csv.DictReader(f))
    assert len(rows)>=10 and len({x["id"] for x in rows})==len(rows)
    with (folder/"13_RESOURCE_CONFIG_MATRIX.csv").open(encoding="utf-8",newline="") as f:
        rr=list(csv.DictReader(f))
    assert len(rr)==4
    print(json.dumps({"PASS":True,"phases":ids,"source_index":len(rows),"resource_configs":list(r),"claim":"static design validation only"},ensure_ascii=False))
    return True

if __name__=="__main__":sys.exit(0 if validate(Path(sys.argv[1] if len(sys.argv)>1 else Path(__file__).resolve().parent.parent)) else 2)
