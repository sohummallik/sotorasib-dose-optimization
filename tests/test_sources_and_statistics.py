"""Source-linked inputs and independently expressed statistical checks."""
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from grouped_exposure_response import fit_grouped
from design_scenario import power_normal_approximation
ROOT=Path(__file__).resolve().parents[1]


class SourceAndStatisticTests(unittest.TestCase):
    def test_fda_figure26_transcription(self):
        # FDA NDA 214665 review, PDF page 248, Figure 26, left panel.
        # https://www.accessdata.fda.gov/drugsatfda_docs/nda/2021/214665Orig1s000MultidisciplineR.pdf
        data=pd.read_csv(ROOT/'data/fda_figure26_quartiles.csv')
        self.assertEqual(data.auc_lower.tolist(),[5000,16000,22000,35000])
        self.assertEqual(data.auc_upper.tolist(),[16000,22000,35000,85000])
        self.assertEqual(data.auc_label.tolist(),[13000,19000,27000,46000])
        self.assertEqual(data.responders.tolist(),[25,26,21,12])
        self.assertEqual(data.total.tolist(),[57]*4)
        self.assertFalse(np.allclose(data.auc_label,np.sqrt(data.auc_lower*data.auc_upper)))

    def test_grouped_fit_equivalent_to_frequency_weighted_bernoulli(self):
        # Expanding counts at fixed coordinates checks the grouped likelihood,
        # but does NOT recover any real patient-level exposures.
        import statsmodels.api as sm
        data=pd.read_csv(ROOT/'data/fda_figure26_quartiles.csv')
        fitted=fit_grouped(data,data.auc_label)
        x=[];y=[]
        for row in data.itertuples():
            x.extend([np.log(row.auc_label)]*row.total)
            y.extend([1]*row.responders+[0]*(row.total-row.responders))
        independent=sm.Logit(y,sm.add_constant(x)).fit(disp=False)
        self.assertAlmostEqual(fitted['slope_per_log_AUC'],independent.params[1],places=8)

    def test_design_integer_boundary(self):
        self.assertLess(power_normal_approximation(.35,.25,328),.8)
        self.assertGreaterEqual(power_normal_approximation(.35,.25,329),.8)
        # Familiar one-tail-in-the-alternative normal planning formula, independent
        # of the numerical root, rounds to the same n for this 10-point scenario.
        pooled=.3
        formula=(norm.ppf(.975)*np.sqrt(2*pooled*(1-pooled))+
                 norm.ppf(.8)*np.sqrt(.35*.65+.25*.75))**2/.1**2
        self.assertEqual(int(np.ceil(formula)),329)

    def test_null_power_matches_alpha(self):
        self.assertAlmostEqual(power_normal_approximation(.3,.3,329),.05,places=12)

if __name__=='__main__':unittest.main()
