"""Descriptive grouped-binomial coordinate sensitivity using FDA Figure 26.

Four grouped counts are fitted at chosen coordinates. Individual exposures and
covariates are unavailable. Conditional model intervals omit coordinate uncertainty,
within-bin variation and confounding; they are not causal exposure-effect intervals.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm

ROOT=Path(__file__).resolve().parents[1]


def fit_grouped(data, coordinates):
    coordinates=np.asarray(coordinates,dtype=float)
    if len(coordinates)!=len(data) or np.any(coordinates<=0):
        raise ValueError('Need a positive coordinate for every group')
    response=np.column_stack([data.responders,data.total-data.responders])
    fit=sm.GLM(response,sm.add_constant(np.log(coordinates)),family=sm.families.Binomial()).fit()
    slope,se=float(fit.params[1]),float(fit.bse[1])
    interval=np.exp((slope+np.array([-1.96,1.96])*se)*np.log(2))
    return {'intercept':float(fit.params[0]),'slope_per_log_AUC':slope,'slope_SE':se,
        'OR_per_doubling':float(np.exp(slope*np.log(2))),
        'conditional_95CI_lower':float(interval[0]),'conditional_95CI_upper':float(interval[1])}


def main():
    data=pd.read_csv(ROOT/'data/fda_figure26_quartiles.csv')
    schemes={'displayed_labels':data.auc_label.to_numpy(),
             'geometric_range_midpoints':np.sqrt(data.auc_lower*data.auc_upper).to_numpy(),
             'arithmetic_range_midpoints':((data.auc_lower+data.auc_upper)/2).to_numpy()}
    rows=[]
    for name,coordinates in schemes.items():
        rows.append({'coordinate_scheme':name,**fit_grouped(data,coordinates),
                     **{f'coordinate_q{i+1}':float(x) for i,x in enumerate(coordinates)}})
    pd.DataFrame(rows).to_csv(ROOT/'results/grouped_coordinate_sensitivity.csv',index=False)
    summary={'figure_population_total':int(data.total.sum()),'number_of_groups':len(data),
        'source':'FDA NDA 214665 multidisciplinary review, Figure 26, PDF page 248',
        'unit':'AUCtau,ss in h*ng/mL', 'models':rows,
        'limitations':['The source text mentions an ER population of 248, while these four bins total 228; not resolved.',
            'The source does not explicitly name the statistic represented by its displayed central labels.',
            'All fits use aggregate counts at chosen coordinates, not patient-level continuous exposure data.',
            'Conditional Wald intervals omit coordinate uncertainty, within-bin exposure variability and confounding.',
            'Alternative midpoints are representation sensitivity analyses, not estimates of the actual within-bin means.']}
    (ROOT/'results/grouped_exposure_response.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(pd.DataFrame(rows)[['coordinate_scheme','OR_per_doubling','conditional_95CI_lower','conditional_95CI_upper']].to_string(index=False))

if __name__=='__main__':
    main()
