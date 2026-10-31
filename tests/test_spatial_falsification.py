from itertools import product
import unittest
import numpy as np
from study_spatial_falsification import mc_p, nearest_distances, robust_oracle


def direct_score(a, b, radius=3.):
    delta = abs(np.asarray(b)[..., :, None, :]-np.asarray(a))
    distances = np.sqrt(np.sum(np.minimum(delta, 100-delta)**2, axis=-1))
    return np.mean(np.min(distances, axis=-1) <= radius, axis=-1)


class SpatialFalsificationTests(unittest.TestCase):
    def test_periodic_nearest_distances(self):
        a = np.array([[1., 1.], [20., 20.]])
        b = np.array([[99., 1.], [22., 20.], [50., 50.]])
        np.testing.assert_allclose(nearest_distances(a, b), [2, 2, np.sqrt(1800)])

    def test_error_envelope_dominates_all_admissible_toy_p_values(self):
        observed_a = np.array([[10., 10.], [20., 20.]])
        null_b = np.array([[[13.1, 10.], [23.1, 20.]], [[12.6, 10.], [23.5, 20.]], [[12.5, 10.], [22.5, 20.]]])
        for rho, distances in [(.25, [2.8, 2.9]), (.1, [2.6, 2.7])]:
            observed_b = observed_a+np.array(distances)[:, None]*np.array([[1., 0.], [1., 0.]])
            upper, _ = robust_oracle(observed_a, observed_b, null_b, rho)
            choices = np.array([[rho, 0], [-rho, 0], [0, rho], [0, -rho]])
            for offsets in product(range(4), repeat=4):
                a = observed_a+choices[list(offsets[:2])]
                b = observed_b+choices[list(offsets[2:])]
                latent_p = mc_p(float(direct_score(a, b)), direct_score(a, null_b))
                self.assertGreaterEqual(upper+1e-12, latent_p)

    def test_resolution_limit_with_positive_observed_distances(self):
        a = np.array([[10., 10.]])
        b = np.array([[10.1, 10.1]])
        null = np.array([[[80., 80.]], [[15., 15.]]])
        p, lower = robust_oracle(a, b, null, 1.5)
        self.assertEqual(lower, 0.)
        self.assertEqual(p, 1.)


if __name__ == '__main__':
    unittest.main()
