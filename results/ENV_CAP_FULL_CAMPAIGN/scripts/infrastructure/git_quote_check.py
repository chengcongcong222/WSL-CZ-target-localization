from pathlib import Path
import subprocess,json
p=Path('results/ENV_CAP_FULL_CAMPAIGN/scripts/archive_and_budget.py');s=p.read_text(encoding='utf-8-sig');s=s.replace("'diff','--name-only','a9a80191ee4d0f0cecd9db296198606b63a68080']).decode().splitlines()","'diff','--name-only','-z','a9a80191ee4d0f0cecd9db296198606b63a68080']).decode('utf-8').strip('\\0').split('\\0')");p.write_text(s,encoding='utf-8',newline='\n')
r=subprocess.check_output(['git','-c',f'safe.directory={Path.cwd().as_posix()}','diff','--name-only','-z','a9a80191ee4d0f0cecd9db296198606b63a68080']).decode('utf-8').strip('\0').split('\0');bad=[a for a in r if not a.startswith('results/ENV_CAP_FULL_CAMPAIGN/')];print('Outside stage tracked changes:',bad);assert not bad
print(Path('results/ENV_CAP_FULL_CAMPAIGN/BUDGET_USAGE.json').read_text())
