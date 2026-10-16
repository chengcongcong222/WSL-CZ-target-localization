"""Independent lightweight standard-library tests; no solver or Git mutations."""
import importlib.util, pathlib, tempfile, unittest, json

ROOT=pathlib.Path(__file__).resolve().parents[1]

def load(path, name):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj

class ToolTests(unittest.TestCase):
    def test_manifest_self_check(self):
        mod=load('scripts/verify_manifest.py','pkg_verify')
        self.assertTrue(mod.verify(ROOT)['PASS'])
    def test_stage_state_init_and_reject_duplicate(self):
        import subprocess, sys
        with tempfile.TemporaryDirectory() as d:
            tool=ROOT/'scripts/stage_state.py'
            def call(*args):
                return subprocess.run([sys.executable,str(tool),*args],capture_output=True,text=True)
            self.assertEqual(call('init',d).returncode,0)
            self.assertNotEqual(call('init',d).returncode,0)
            self.assertEqual(call('start',d,'E0','baseline').returncode,0)
            self.assertEqual(call('finish',d,'E0','PASS','reviewed').returncode,0)
            self.assertEqual(call('start',d,'E1','environment candidate').returncode,0)
            self.assertEqual(call('finish',d,'E1','FAIL','no CZ2').returncode,0)
            self.assertEqual(call('start',d,'E2','paper control still possible').returncode,0)
            self.assertEqual(json.loads((pathlib.Path(d)/'RUN_STATE.json').read_text())['phases']['E2'],'RUNNING')
    def test_source_registry_headers(self):
        import csv
        with (ROOT/'10_REFERENCE_SOURCE_INDEX.csv').open(newline='',encoding='utf-8') as f:
            r=list(csv.DictReader(f))
        self.assertGreaterEqual(len(r),20)
        self.assertTrue(all(x['must_exist'] in {'YES','NO'} for x in r))
        self.assertTrue(all(x['evidence_alias'].strip() and '\n' not in x['evidence_alias'] for x in r))

if __name__ == '__main__':unittest.main()
