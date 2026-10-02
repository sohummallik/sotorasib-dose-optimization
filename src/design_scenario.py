"""One hypothetical equal-allocation ORR superiority design; not sponsor reconstruction."""
import json
from pathlib import Path
import numpy as np
from scipy.stats import norm
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]


def power_normal_approximation(p1,p2,n,alpha=.05):
    """Two-sided z boundary, pooled-null and unpooled-alternative variance.

    Equal n per arm; independent binomial outcomes; normal approximation without
    continuity correction or attrition adjustment.
    """
    pooled=(p1+p2)/2
    se_null=np.sqrt(2*pooled*(1-pooled)/n)
    se_alternative=np.sqrt((p1*(1-p1)+p2*(1-p2))/n)
    z=norm.ppf(1-alpha/2)
    difference=abs(p1-p2)
    return float(norm.cdf((difference-z*se_null)/se_alternative)+
                 norm.cdf((-difference-z*se_null)/se_alternative))


def main():
    p1,p2,alpha,target=.35,.25,.05,.8
    solution=brentq(lambda n:power_normal_approximation(p1,p2,n,alpha)-target,2,100000)
    n=int(np.ceil(solution))
    result={'assumed_ORR_arm_1':p1,'assumed_ORR_arm_2':p2,'difference_percentage_points':10,
        'allocation':'1:1','alpha_two_sided':alpha,'target_power':target,
        'normal_approximation_continuous_n_per_arm':float(solution),
        'n_per_arm':n,'n_total':2*n,
        'power_at_selected_n':power_normal_approximation(p1,p2,n,alpha),
        'power_at_one_fewer_per_arm':power_normal_approximation(p1,p2,n-1,alpha),
        'method':'Two-sided two-proportion z normal approximation, pooled-null rejection boundary and unpooled alternative variance; no continuity correction.',
        'limitations':['Hypothetical rates, not original-trial planning assumptions.',
            'No attrition or unevaluable-patient allowance.',
            'Response-assessment time, meaningful difference and choice of superiority objective require clinical justification.',
            'Not a noninferiority design or a uniquely recommended sample size.'],
        'random_seed':None,'randomness':'None; deterministic calculation.'}
    (ROOT/'results/design_scenario.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
