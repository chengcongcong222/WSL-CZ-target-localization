import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from threadpoolctl import threadpool_limits
import r4_a1_fix2_coverage as fix


class AcousticCoverage(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.limits=threadpool_limits(limits=2);cls.limits.__enter__()
        cls.model=fix.Acoustic(fix.inputs())

    @classmethod
    def tearDownClass(cls):cls.limits.__exit__(None,None,None)

    def test_exact_analytic_features_and_derivatives(self):
        s=np.array([50.123,.234,2.132,4.531]);depth=165.
        value,jac=fix.feature_jacobian(self.model,s,depth,True)
        reference=fix.L.direct_features(self.model.models,fix.L.geometry(s)[1],[depth])[0,0]
        np.testing.assert_allclose(value,reference,atol=1e-10,rtol=0)
        for k,h in enumerate([1e-6,1e-5,1e-6,1e-5]):
            plus=s.copy();minus=s.copy();plus[k]+=h;minus[k]-=h
            numerical=(fix.feature_jacobian(self.model,plus,depth,True)[0]-fix.feature_jacobian(self.model,minus,depth,True)[0])/(2*h)
            relative=np.linalg.norm(numerical-jac[:,:,k])/np.linalg.norm(numerical)
            self.assertLess(relative,1e-4)

    def test_spline_analytic_derivatives(self):
        s=np.array([50.123,.234,2.132,4.531]);_,jac=fix.feature_jacobian(self.model,s,200.,False)
        for k,h in enumerate([1e-6,1e-5,1e-6,1e-5]):
            plus=s.copy();minus=s.copy();plus[k]+=h;minus[k]-=h
            numerical=(fix.feature_jacobian(self.model,plus,200.)[0]-fix.feature_jacobian(self.model,minus,200.)[0])/(2*h)
            self.assertLess(np.linalg.norm(numerical-jac[:,:,k])/np.linalg.norm(numerical),1e-4)

    def test_radial_continuation_matches_independent_fit_without_truth_helpers(self):
        observed=fix.L.geometry([51.123,.432,2.231,-4.321])[0][0]+np.random.default_rng(3777).normal(0,np.radians(.1),121)
        r=50321.;u=2.01
        with patch.object(fix.L,'errors',side_effect=AssertionError('oracle error leaked')),patch.object(fix.L,'neighborhood',side_effect=AssertionError('truth cell leaked')):
            s=fix.radial_profile([r],[u],observed)[0]
        def state(x):return np.array([r/1000,x[0],np.hypot(u,x[1]),x[0]+np.degrees(np.arctan2(x[1],u))])
        fit=least_squares(lambda x:fix.prior.angular_residual(state(x),observed)[0],[.4,-.1],bounds=([-5,-1],[5,1]),gtol=1e-13,xtol=1e-13,ftol=1e-13)
        cost=np.square(fix.prior.angular_residual(s,observed)).sum()
        self.assertLess(abs(cost-np.square(fit.fun).sum()),1e-12)

    def test_cumulative_bank_preserves_earlier_recovery(self):
        frames=pd.DataFrame([dict(budget='B1',r_km=50.,theta_deg=0.,v_mps=2.,psi_deg=5.,J_exact=.0002),
                             dict(budget='B2',r_km=52.,theta_deg=0.,v_mps=1.8,psi_deg=6.,J_exact=1.1),
                             dict(budget='B3',r_km=53.,theta_deg=0.,v_mps=1.7,psi_deg=7.,J_exact=.8)])
        pooled=fix.cumulative_candidates(frames,['B1','B2','B3'])
        self.assertTrue((pooled.groupby('budget').J_exact.min()==.0002).all())
        self.assertEqual(pooled[pooled.budget=='B3'].iloc[0].origin_budget,'B1')


if __name__=='__main__':unittest.main()
