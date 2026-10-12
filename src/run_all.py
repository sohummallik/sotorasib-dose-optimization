"""Run the supported deterministic pipeline and record environment/commands."""
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
PACKAGES = ['numpy', 'pandas', 'scipy', 'statsmodels', 'matplotlib']
SCRIPTS = ['validate.py', 'numerical_checks.py', 'grouped_exposure_response.py',
           'design_scenario.py', 'design_sensitivity.py', 'figures.py']
OUTPUTS = [
    'results/pk_covariate_checks.csv', 'results/pk_check_summary.json',
    'results/grouped_coordinate_sensitivity.csv', 'results/grouped_exposure_response.json',
    'results/design_scenario.json', 'results/design_sensitivity.csv',
    'results/design_precision.csv', 'reports/numerical_checks.json',
] + [f'figures/{name}.{ext}' for name in
     ['pk_covariate_checks', 'grouped_coordinate_sensitivity', 'design_sensitivity']
     for ext in ['png', 'pdf']]


def file_hashes(paths):
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in paths if path.is_file()}


def main():
    for directory in ['results', 'figures', 'reports']:
        (ROOT / directory).mkdir(exist_ok=True)
    report_path = ROOT / 'reports/execution_report.json'
    inputs = [*sorted((ROOT / 'data').glob('*')),
              *sorted((ROOT / 'src').glob('*.py')),
              *sorted((ROOT / 'src').glob('*.R')),
              *sorted((ROOT / 'tests').glob('*.py')),
              *sorted((ROOT / 'tests').glob('*.R')),
              *sorted(ROOT.glob('requirements*.txt')),
              *sorted((ROOT / '.github/workflows').glob('*.yml'))]
    report = {
        'status': 'running',
        'started_at_utc': datetime.now(timezone.utc).isoformat(),
        'python': sys.version, 'executable': sys.executable,
        'platform': platform.platform(), 'packages': {},
        'seeds': 'No random draws in the supported active pipeline.',
        'commands': [], 'input_and_code_sha256': file_hashes(inputs),
        'derived_files': OUTPUTS,
        'not_performed': [
            'R reproduction checks (run separately)',
            'R estimation or covariance recovery', 'Patient-level clinical fitting',
            'Historical calibration/utility revalidation', 'Clinical model validation',
        ],
    }

    def save_report():
        report_path.write_text(json.dumps(report, indent=2) + '\n')

    save_report()
    try:
        for package in PACKAGES:
            report['packages'][package] = importlib.metadata.version(package)
        commands = [['-m', 'unittest', 'discover', '-s', 'tests', '-v']]
        commands += [['src/' + script] for script in SCRIPTS]
        for index, args in enumerate(commands, start=1):
            print(f'[{index}/{len(commands)}] Running {" ".join(args)}', flush=True)
            start = time.monotonic()
            proc = subprocess.run([sys.executable, *args], cwd=ROOT,
                                  text=True, capture_output=True)
            log = f'reports/active_command_{index}.log'
            (ROOT / log).write_text(proc.stdout + proc.stderr)
            report['commands'].append({
                'command': [sys.executable, *args], 'returncode': proc.returncode,
                'elapsed_seconds': time.monotonic() - start, 'log': log,
            })
            save_report()
            if proc.returncode:
                raise RuntimeError(f'Failed: {" ".join(args)}; inspect {log}')
        missing = [name for name in OUTPUTS if not (ROOT / name).is_file()]
        if missing:
            raise RuntimeError(f'Missing expected output files: {", ".join(missing)}')
        report['derived_file_sha256'] = file_hashes(ROOT / name for name in OUTPUTS)
        report['status'] = 'completed'
    except Exception as error:
        report['status'] = 'failed'
        report['error'] = str(error)
        if isinstance(error, importlib.metadata.PackageNotFoundError):
            print('Install the tested dependencies with: '
                  'python -m pip install -r requirements-tested.txt', file=sys.stderr)
        raise
    finally:
        report['finished_at_utc'] = datetime.now(timezone.utc).isoformat()
        save_report()
    print('Supported pipeline completed; reports/execution_report.json records '
          'all commands and statuses.')


if __name__ == '__main__':
    main()
