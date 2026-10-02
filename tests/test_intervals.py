"""Regression tests for confirmed endpoints/grid defects and oral dose semantics."""
import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from model import load_params, individual_params, typical_subject, simulate, interval_metrics


class IntervalTests(unittest.TestCase):
    def test_exact_endpoints_constant_curve(self):
        data=pd.DataFrame({'time_h':[0,12,24,24,36,48],'conc_ng_mL':[100]*6})
        for day in [1,2]:
            result=interval_metrics(data,day)
            self.assertEqual(result['AUC_tau'],2400)
            self.assertEqual(result['C_tau'],100)

    def test_unaligned_boundaries_linear_interpolation(self):
        # c(t)=2t+3: integral over [24,48] is 1800 exactly.
        times=np.array([0,13,26,39,52])
        data=pd.DataFrame({'time_h':times,'conc_ng_mL':2*times+3})
        result=interval_metrics(data,2)
        self.assertAlmostEqual(result['AUC_tau'],1800,places=12)
        self.assertEqual(result['C_tau'],99)
        self.assertEqual(result['Tmax'],24)

    def test_no_extrapolation(self):
        data=pd.DataFrame({'time_h':[1,12,23],'conc_ng_mL':[1,2,1]})
        with self.assertRaisesRegex(ValueError,'no extrapolation'):
            interval_metrics(data,1)

    def test_discontinuous_duplicate_rejected(self):
        data=pd.DataFrame({'time_h':[0,24,24,48],'conc_ng_mL':[0,10,20,5]})
        with self.assertRaisesRegex(ValueError,'event convention'):
            interval_metrics(data,1)

    def test_invalid_arguments_and_nonfinite_samples(self):
        data=pd.DataFrame({'time_h':[0,24],'conc_ng_mL':[0,1]})
        for tau in [0,-1,np.nan,np.inf]:
            with self.assertRaises(ValueError):interval_metrics(data,1,tau=tau)
        with self.assertRaises(ValueError):interval_metrics(data.assign(conc_ng_mL=[0,np.nan]),1)

    def test_oral_event_continuity_and_integral_additivity(self):
        ip=individual_params(load_params(),typical_subject())
        data=simulate(ip,960,3,dt=.7)
        event=data[data.time_h==24].conc_ng_mL.to_numpy()
        self.assertEqual(len(event),2)
        np.testing.assert_allclose(event[0],event[1],rtol=1e-12,atol=1e-10)
        whole=interval_metrics(data,1,tau=72)['AUC_tau']
        pieces=sum(interval_metrics(data,k)['AUC_tau'] for k in [1,2,3])
        self.assertAlmostEqual(whole/pieces,1,places=12)

    def test_nonbinary_grid_reaches_every_exact_boundary(self):
        ip=individual_params(load_params(),typical_subject())
        data=simulate(ip,960,30,dt=.05)
        self.assertEqual(data.time_h.min(),0)
        self.assertEqual(data.time_h.max(),720)
        for boundary in np.arange(31)*24:
            self.assertIn(boundary,data.time_h.values)

    def test_terminal_concentration_is_not_interval_minimum(self):
        data=pd.DataFrame({'time_h':[0,12,24],'conc_ng_mL':[0,10,2]})
        self.assertEqual(interval_metrics(data,1)['C_tau'],2)

    def test_frozen_parameter_mass_balance(self):
        # At periodic steady state, integral = F*dose/CL for this linear model.
        # An independent mass-balance identity tests more than agreement with old code.
        ip=individual_params(load_params(),typical_subject())
        ip.update(CLBS=ip['CLSS'],F1BS=ip['F1SS'])
        result=interval_metrics(simulate(ip,960,60,dt=.03125),60)
        exact=960*ip['F1SS']/ip['CLSS']*1000
        self.assertLess(abs(result['AUC_tau']/exact-1),5e-4)

    def test_grid_and_solver_tolerance_sensitivity(self):
        ip=individual_params(load_params(),typical_subject())
        coarse=interval_metrics(simulate(ip,960,8,dt=.0625),8)
        fine=interval_metrics(simulate(ip,960,8,dt=.03125),8)
        strict=interval_metrics(simulate(ip,960,8,dt=.03125,rtol=1e-10,atol=1e-8),8)
        # Numerical budgets: <0.05% AUC and <0.1% peak-grid changes; solver <0.001%.
        self.assertLess(abs(coarse['AUC_tau']/fine['AUC_tau']-1),5e-4)
        self.assertLess(abs(coarse['Cmax']/fine['Cmax']-1),1e-3)
        for metric in ['AUC_tau','Cmax','C_tau']:
            self.assertLess(abs(fine[metric]/strict[metric]-1),1e-5)

if __name__=='__main__':unittest.main()
