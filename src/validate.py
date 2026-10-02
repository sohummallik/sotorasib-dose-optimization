"""Reproduce selected published PK covariate point contrasts; not clinical validation."""
import json
from pathlib import Path
import pandas as pd
from model import load_params, individual_params, simulate, interval_metrics, typical_subject

ROOT = Path(__file__).resolve().parents[1]
N_DOSES = 30
GRID_H = 0.03125
TOLERANCE_PCT = 5.0  # retained baseline point-contrast tolerance, not a validity threshold


def contrast_metrics(overrides, dt=GRID_H):
    p = load_params()
    cov = typical_subject(960)
    cov.update(overrides)
    ip = individual_params(p, cov)
    return interval_metrics(simulate(ip, 960, N_DOSES, dt=dt), N_DOSES)


def main():
    targets = json.loads((ROOT / 'data/covariate_targets.json').read_text())['contrasts']
    ref = contrast_metrics({})
    rows = []
    for target in targets:
        result = contrast_metrics(target['overrides'])
        row = {'contrast': target['contrast']}
        for key in ['AUC_tau', 'Cmax']:
            estimate = result[key] / ref[key]
            published = target[key]
            row.update({key+'_ratio_model': estimate, key+'_ratio_published': published,
                        key+'_pct_error': 100*(estimate/published-1)})
        row['pass_5pct'] = all(abs(row[k+'_pct_error']) <= TOLERANCE_PCT for k in ['AUC_tau','Cmax'])
        rows.append(row)
    data = pd.DataFrame(rows)
    (ROOT/'results').mkdir(exist_ok=True)
    data.to_csv(ROOT/'results/pk_covariate_checks.csv', index=False)
    summary = {'contrasts':len(data), 'point_ratios':2*len(data),
        'passing_contrasts':int(data.pass_5pct.sum()),
        'maximum_absolute_percent_error':float(data[['AUC_tau_pct_error','Cmax_pct_error']].abs().to_numpy().max()),
        'grid_h':GRID_H,'dose_mg':960,'day':N_DOSES,'solver':'LSODA',
        'rtol':1e-8,'atol':1e-6,'point_contrast_tolerance_percent':TOLERANCE_PCT,
        'interpretation':'Selected published point-contrast agreement only. Published uncertainty bands and original control stream were not reproduced; clinical validity is not established.'}
    (ROOT/'results/pk_check_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    if not data.pass_5pct.all():
        raise RuntimeError('A published point contrast exceeded the retained tolerance')

if __name__ == '__main__':
    main()
