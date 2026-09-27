"""Run the supported deterministic pipeline and record environment/commands."""
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]


def main():
    for directory in ['results','figures','reports']:(ROOT/directory).mkdir(exist_ok=True)
    report={'baseline_commit':'043a80508e35306f2ddaac6edc9a621c12a25b21',
        'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
        'packages':{p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','statsmodels','matplotlib']},
        'seeds':'No random draws in the supported active pipeline.', 'commands':[],
        'not_performed':['R estimation or covariance recovery','Patient-level clinical fitting',
                         'Historical calibration/utility revalidation','Clinical model validation']}
    commands=[['-m','unittest','discover','-s','tests','-v']]+[['src/'+x] for x in
        ['validate.py','numerical_checks.py','grouped_exposure_response.py','design_scenario.py','figures.py']]
    for index,args in enumerate(commands):
        start=time.monotonic()
        proc=subprocess.run([sys.executable,*args],cwd=ROOT,text=True,capture_output=True)
        log=f'reports/active_command_{index+1}.log'
        (ROOT/log).write_text(proc.stdout+proc.stderr)
        report['commands'].append({'command':[sys.executable,*args],'returncode':proc.returncode,
                                   'elapsed_seconds':time.monotonic()-start,'log':log})
        if proc.returncode:
            (ROOT/'reports/execution_report.json').write_text(json.dumps(report,indent=2)+'\n')
            raise RuntimeError(f'Failed: {args}; inspect {log}')
    inputs=[*sorted((ROOT/'data').glob('*')),*sorted((ROOT/'src').glob('*.py')),*sorted((ROOT/'tests').glob('*.py'))]
    report['input_and_code_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                                    for p in inputs if p.is_file()}
    report['derived_files']=[str(p.relative_to(ROOT)) for folder in ['results','figures'] for p in sorted((ROOT/folder).glob('*')) if p.is_file()]
    (ROOT/'reports/execution_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Supported pipeline completed; reports/execution_report.json records all commands and statuses.')

if __name__=='__main__':main()
