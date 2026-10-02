import itertools
import unittest
import numpy as np
import pandas as pd
import r4_a1_offgrid_bearing as r4


class OffgridRegression(unittest.TestCase):
    def test_circular_angles(self):
        self.assertEqual(float(r4.angle_error_deg(179,-179)),2)
        self.assertEqual(float(r4.angle_error_deg(-179,179)),2)

    def test_midpoint_all_ties(self):
        truth = pd.Series(dict(zip(r4.AXES,[50.5,.25,2.1,5.5])))
        nearest,vertices = r4.neighborhood(truth,r4.grid())
        expected = sorted(np.ravel_multi_index(i,r4.SHAPE) for i in itertools.product((5,6),(10,11),(5,6),(20,21)))
        self.assertEqual(nearest.tolist(),expected)
        self.assertEqual(sorted(vertices.tolist()),expected)

    def test_topology_no_diagonal_or_edge_wrap(self):
        for indices in [((0,0,0,0),(1,1,0,0)),((0,0,0,30),(0,0,1,0))]:
            ids = np.array([np.ravel_multi_index(i,r4.SHAPE) for i in indices])
            self.assertEqual(r4.topology(ids)['n_connected_components'],2)

    def test_costs_wrapping_and_observation_immutability(self):
        est = r4.FrozenEstimator(np.array([[50,0,2,5],[49,.5,1.8,4]]))
        for observation in (est.pred[0]+.001,est.pred[0]+2*np.pi):
            before = observation.copy()
            expected = np.square(np.arctan2(np.sin(est.pred-observation),np.cos(est.pred-observation))).sum(axis=1)
            np.testing.assert_allclose(est.bearing_costs(observation)[:,0],expected,atol=1e-12,rtol=0)
            np.testing.assert_array_equal(observation,before)

    def test_profile_mapping_and_score(self):
        mod = {'depths':np.arange(150,251,2)}
        np.testing.assert_array_equal(r4.effective_depths(mod,[155,205])[1],[154,204])
        f = np.zeros((2,21,3,121)); f[0] = 3; f[1] = 5
        f[0,4] = 1; f[1,10] = 2
        js,z = r4.score_features(f,np.zeros((3,121)))
        np.testing.assert_array_equal(js,[1,2]); np.testing.assert_array_equal(z,[170,200])


if __name__=='__main__':
    unittest.main()
