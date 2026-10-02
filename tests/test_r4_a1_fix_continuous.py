import unittest
from unittest.mock import patch
import numpy as np
from scipy.optimize import minimize_scalar
import r4_a1_fix_continuous as fix


class ContinuousRepair(unittest.TestCase):
    def test_noiseless_observation_recovers_without_evaluation_helpers(self):
        s=np.array([49.123,-.417,2.234,-4.345])
        t=np.arange(121)*10.; dt=np.maximum(t-600,0)
        x=1000*s[0]*np.cos(np.radians(s[1]))+s[2]*t*np.cos(np.radians(s[3]))-2*np.minimum(t,600)-2*dt*np.cos(np.pi/12)
        y=1000*s[0]*np.sin(np.radians(s[1]))+s[2]*t*np.sin(np.radians(s[3]))-2*dt*np.sin(np.pi/12)
        bearing=np.arctan2(y,x)
        with patch.object(fix.legacy,'errors',side_effect=AssertionError('evaluation leaked')),patch.object(fix.legacy,'neighborhood',side_effect=AssertionError('truth-cell leaked')):
            cloud,cmin,cutoff,logs,sv=fix.rc2(bearing,0)
        self.assertLess(cmin,1e-20)
        self.assertEqual(len(sv[sv>1e-8]),4)
        np.testing.assert_allclose(cloud[0],s,atol=1e-7,rtol=0)
        self.assertEqual(len(logs),32)

    def test_observation_theta_profile_against_scalar_optimization(self):
        rv=np.array([[48.3,1.7,-7.2],[53.1,2.4,8.4]])
        observed=fix.legacy.geometry([51.12,.32,2.13,5.6])[0][0]+np.random.default_rng(8877).normal(0,np.radians(.1),121)
        got=fix.theta_profile(rv,observed)
        for s in got:
            def cost(theta):
                state=s.copy(); state[1]=theta
                return float(np.square(fix.angular_residual(state,observed)).sum())
            reference=minimize_scalar(cost,bounds=(-5,5),method='bounded',options={'xatol':1e-12})
            self.assertAlmostEqual(s[1],reference.x,places=6)

    def test_geometric_range_and_jacobian_bounds(self):
        # Whole engineering region corners, both trajectory windows.
        corners=np.array(np.meshgrid(*zip(fix.LOW,fix.HIGH))).T.reshape(-1,4)
        bearing,ranges=fix.legacy.geometry(corners)
        self.assertTrue((ranges>=39000).all() and (ranges<=66000).all())
        h=1e-5; shifted=corners.copy(); shifted[:,1]+=h
        derivative=fix.legacy.wrap_rad(fix.legacy.geometry(shifted)[0]-bearing)/np.radians(h)
        self.assertGreater(derivative.min(),45000*39000/66000**2)

    def test_global_population_explores_theta_and_remains_feasible(self):
        observed=fix.legacy.geometry([51.12,.32,2.13,5.6])[0][0]+np.random.default_rng(8877).normal(0,np.radians(.1),121)
        cloud,_,cutoff,_,_=fix.rc2(observed,.1)
        pop=fix.acoustic_population(cloud,observed,.1,cutoff,np.random.default_rng(771001))
        self.assertGreater(np.ptp(pop[:,3]),.1)
        rv=fix.LOW[[0,2,3]]+pop[:,:3]*fix.WIDTH[[0,2,3]]
        s=fix.theta_profile(rv,observed)
        bound=.1*np.sqrt(13.3/121)/(45000*39000/66000**2)
        s[:,1]+=pop[:,3]*bound
        costs=np.square(fix.angular_residual(s,observed)).sum(axis=1)
        self.assertTrue((costs<=cutoff+1e-14).all())


if __name__=='__main__':unittest.main()
