"""Boundary cases for grid topology, full-state truth and final-only decisions."""
import itertools
import json
import unittest
import numpy as np
import pandas as pd
import r3_final_parameter_contraction_fix as audit


def cloud(states):
    frame = pd.DataFrame(states, columns=audit.AXES)
    idx = np.rint((frame.to_numpy()-audit.ORIGIN)/audit.STEPS).astype(int)
    frame['node_id'] = np.ravel_multi_index(idx.T, audit.SHAPE)
    return audit.normalize(frame)


class AuditBoundaries(unittest.TestCase):
    def test_bbox_includes_missing_bins_and_ranges_split(self):
        frame = cloud([(45,0,2,5),(47,0,2,5)])
        quality = audit.metrics(frame)
        self.assertEqual(quality['n_occupied_r_bins'], 2)
        self.assertEqual(quality['n_bbox_r_bins'], 3)
        self.assertEqual(quality['N_bbox'], 3)
        self.assertAlmostEqual(quality['correlation_sparsity'], 2/3)
        self.assertEqual(audit.component_metrics(frame)['n_range_components'], 2)

    def test_diagonal_nodes_are_not_manhattan_neighbors(self):
        frame = cloud([(50,0,2,4),(51,0,2,5)])
        result = audit.component_metrics(frame)
        self.assertEqual(result['n_connected_components'], 2)
        self.assertEqual(result['n_range_components'], 1)
        self.assertEqual(json.loads(result['component_size_top5']), [1,1])
        bridge = cloud([(50,0,2,4),(50,0,2,5),(51,0,2,5)])
        self.assertEqual(audit.component_metrics(bridge)['largest_component_size'], 3)

    def test_directions_do_not_wrap_at_grid_boundary(self):
        frame = cloud([(50,0,2,-15),(50,0,2,15)])
        self.assertEqual(audit.component_metrics(frame)['n_connected_components'], 2)

    def test_all_four_axes_can_connect(self):
        frame = cloud([(50,0,2,5),(51,0,2,5),(51,.5,2,5),(51,.5,2.2,5),(51,.5,2.2,6)])
        self.assertEqual(audit.component_metrics(frame)['largest_component_size'], 5)

    def test_truth_rank_uses_full_cloud_and_exact_ties(self):
        scored = cloud([(50,0,2,4),(50,0,2,5),(51,0,2,5),(52,0,2,5)])
        scored['J'] = [.1,.2,.2,.3]
        kept = scored.iloc[[1]]
        rank = audit.truth_rank(scored, kept)
        self.assertEqual(rank['n_strictly_better'], 1)
        self.assertEqual(rank['n_tied_with_truth'], 2)
        self.assertEqual(rank['truth_rank_min'], 2)
        self.assertEqual(rank['truth_rank_max'], 3)
        self.assertFalse(rank['top1_truth_exact_in_matched_synthetic_control'])

    def test_truth_is_full_state_not_any_node_at_truth_range(self):
        frame = cloud([(50,0,2,4),(50,0,2,5)])
        self.assertEqual(frame.is_truth.tolist(), [False,True])
        self.assertFalse(audit.metrics(frame.iloc[[0]])['truth_retained'])

    def test_missing_theta_is_joined_and_id_mismatch_is_rejected(self):
        frame = cloud([(50,.5,2,5)])
        recovered = audit.normalize(frame.drop(columns=['theta']), reference=frame)
        self.assertEqual(recovered.theta.iloc[0], .5)
        broken = frame.copy()
        broken['node_id'] += 1
        with self.assertRaises(ValueError):
            audit.normalize(broken)

    def test_off_grid_coordinate_is_rejected(self):
        frame = cloud([(50,0,2,5)])
        frame['v'] = 2.05
        with self.assertRaises(ValueError):
            audit.normalize(frame)

    def test_final_bound_ignores_intermediate_and_requires_strict_lt(self):
        rows = pd.DataFrame([{'stage':s,'z_true':z,'max_rel_err_r':0.,'max_rel_err_v':0.,
                              'max_rel_err_psi':.09,'max_abs_err_theta':0.}
                             for s,z in itertools.product(('S6','S7','S8'),audit.DEPTHS)])
        rows.loc[rows.stage=='S6','max_rel_err_psi'] = .8
        self.assertEqual(audit.final_bound(rows), 'ESTABLISHED')
        rows.loc[rows.stage=='S7','max_rel_err_psi'] = .10
        self.assertEqual(audit.final_bound(rows), 'NOT_ESTABLISHED')
        with self.assertRaises(ValueError):
            audit.final_bound(rows.drop(rows[(rows.stage=='S8') & (rows.z_true==220)].index))


if __name__ == '__main__':
    unittest.main()
