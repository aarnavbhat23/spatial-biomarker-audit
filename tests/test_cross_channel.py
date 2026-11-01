import unittest
import numpy as np
from study_cross_channel import conformal_radius, fit_transform, predict_transform, partition


class CrossChannelTests(unittest.TestCase):
    def test_affine_recovers_heldout_reference_geometry(self):
        x = np.array([[0., 0.], [1., 0.], [0., 1.], [2., 3.], [-1., 4.]])
        matrix = np.array([[1.2, -.2], [.3, .9], [2., -3.]])
        y = np.column_stack([x, np.ones(len(x))]) @ matrix
        estimate = fit_transform(x[:3], y[:3], 'affine')
        np.testing.assert_allclose(predict_transform(x[3:], estimate, 'affine'), y[3:], atol=1e-12)

    def test_conformal_discreteness_and_insufficient_calibration(self):
        self.assertEqual(conformal_radius(np.arange(33.)), 32.)
        self.assertTrue(np.isinf(conformal_radius(np.arange(10.))))
        # Exchangeable rank argument: every test rank <=k is covered.
        n = 33; k = int(np.ceil(.95*(n+1)))
        self.assertGreaterEqual(k/(n+1), .95)

    def test_features_stay_whole_and_extrapolation_is_ordered(self):
        nodes = {f'{i:03}':[(0, float(i), 0.), (1, float(i)+.1, 0.)] for i in range(20)}
        split = partition(nodes, 'spatial_extrapolation')
        self.assertEqual([sum(v == s for v in split.values()) for s in ['train', 'calibration', 'test']], [10, 5, 5])
        self.assertEqual(split['019'], 'test')
        self.assertEqual(split['000'], 'train')


if __name__ == '__main__':
    unittest.main()
