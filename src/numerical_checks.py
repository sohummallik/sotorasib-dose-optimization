"""Measured numerical errors against analytical identities and refined calculations."""
import json
from pathlib import Path
import numpy as np
from model import load_params,individual_params,typical_subject,simulate,interval_metrics
from validate import contrast_metrics

ROOT=Path(__file__).resolve().parents[1]


def main():
    ip=individual_params(load_params(),typical_subject())
    rows=[]
    for dt,rtol,atol in [(.25,1e-8,1e-6),(.0625,1e-8,1e-6),(.03125,1e-8,1e-6),(.015625,1e-10,1e-8)]:
        m=interval_metrics(simulate(ip,960,8,dt=dt,rtol=rtol,atol=atol),8)
        rows.append({'dt_h':dt,'rtol':rtol,'atol':atol,**m})
    selected=rows[2]; reference=rows[3]
    errors={k:100*(selected[k]/reference[k]-1) for k in ['AUC_tau','Cmax','C_tau']}
    frozen=dict(ip,CLBS=ip['CLSS'],F1BS=ip['F1SS'])
    frozen_result=interval_metrics(simulate(frozen,960,60,dt=.03125),60)
    exact=960*frozen['F1SS']/frozen['CLSS']*1000
    contrasts=json.loads((ROOT/'data/covariate_targets.json').read_text())['contrasts']
    ref1,ref2=contrast_metrics({},.03125),contrast_metrics({},.015625)
    ratios=[]
    for target in contrasts:
        a,b=contrast_metrics(target['overrides'],.03125),contrast_metrics(target['overrides'],.015625)
        for k in ['AUC_tau','Cmax']:
            ratios.append({'contrast':target['contrast'],'metric':k,
                          'percent_change_on_halving_grid':100*((b[k]/ref2[k])/(a[k]/ref1[k])-1)})
    maximum=max(abs(x['percent_change_on_halving_grid']) for x in ratios)
    out={'day8_typical_grid_results':rows,'selected_vs_finer_stricter_percent_difference':errors,
        'frozen_parameter_mass_balance':{'calculated_AUC':frozen_result['AUC_tau'],'exact_AUC':exact,
                                        'relative_error':frozen_result['AUC_tau']/exact-1},
        'covariate_ratio_grid_sensitivity':ratios,'maximum_ratio_percent_change_on_halving_grid':maximum,
        'budgets':{'AUC_relative_error':.0005,'Cmax_relative_grid_change':.001,
                   'solver_relative_change':.00001,'covariate_ratio_percent_grid_change':.1},
        'interpretation':'These are numerical checks of the chosen equations, not tests of clinical model validity.'}
    (ROOT/'reports/numerical_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    if abs(out['frozen_parameter_mass_balance']['relative_error'])>=.0005 or maximum>=.1:
        raise RuntimeError('Numerical error budget exceeded')
    print(json.dumps({k:v for k,v in out.items() if k!='covariate_ratio_grid_sensitivity'},indent=2))

if __name__=='__main__':main()
