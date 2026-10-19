#!/usr/bin/env python3
"""Self tests of package orchestration tools only. Not a physics/estimator test."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parent.parent

def import_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m

stage=import_module("cap06_stage_state",ROOT/"scripts/stage_state.py")
plan=import_module("cap06_validate_plan",ROOT/"scripts/validate_plan.py")

class TestCAP06(unittest.TestCase):
    def test_stage_spec(self):
        self.assertTrue(plan.validate(ROOT))
    def test_stage_independence_and_stop(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"state.json"
            template=ROOT/"16_RUN_STATE_TEMPLATE.json"
            stage.main(["--state",str(p),"init","--template",str(template)])
            stage.main(["--state",str(p),"start","P0"])
            stage.main(["--state",str(p),"finish","P0","PASS"])
            stage.main(["--state",str(p),"start","P2"])
            stage.main(["--state",str(p),"finish","P2","LIMITED"])
            stage.main(["--state",str(p),"start","P1"])
            stage.main(["--state",str(p),"finish","P1","NOT_EVALUATED"])
            stage.main(["--state",str(p),"start","P3"])
            stage.main(["--state",str(p),"finish","P3","LIMITED"])
            stage.main(["--state",str(p),"start","P4"])
            stage.main(["--state",str(p),"finish","P4","PASS"])
            stage.main(["--state",str(p),"stop","--reason","delivered"])
            self.assertTrue(json.loads(p.read_text())["STOP"])
            with self.assertRaises(SystemExit):stage.main(["--state",str(p),"start","P0"])
    def test_duplicate_science_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"state.json"
            stage.main(["--state",str(p),"init","--template",str(ROOT/"16_RUN_STATE_TEMPLATE.json")])
            d=stage.read(p);d["frozen_design_sha"]="a"*40;stage.write_atomic(p,d)
            stage.main(["--state",str(p),"mark-science","P2_TEST_A","b"*64])
            with self.assertRaises(SystemExit):stage.main(["--state",str(p),"mark-science","P2_TEST_A","b"*64])

if __name__=="__main__":unittest.main(verbosity=2)
