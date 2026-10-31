import unittest
import numpy as np
from scipy.stats import binom, poisson
from study_fiducial_measurement import (design, projection_witness, minimum_round_radius,
                                        paired_decision, detect_centroids)


class FiducialMeasurementTests(unittest.TestCase):
    def test_round_color_cycle_breaks_independent_bit_alias(self):
        rounds, colors = [1, 1, 2, 2], [561, 650, 561, 650]
        matrix = design(rounds, colors, 'rounds_and_colors')
        for difference, feasible in [(np.array([0., 1., 0., 1.]), True),
                                     (np.array([0., 1., 1., 0.]), False)]:
            _, residual = projection_witness(matrix, difference)
            self.assertEqual(abs(residual).max() < 1e-10, feasible)
            radius, coefficients = minimum_round_radius(matrix, difference)
            self.assertEqual(radius is not None, feasible)
            if feasible:
                self.assertAlmostEqual(radius, 0.)
                np.testing.assert_allclose(matrix @ coefficients, difference, atol=1e-12)

    def test_minimum_radius_matches_analytic_same_color_two_rounds(self):
        matrix = design([1, 2], [561, 561], 'rounds_and_colors')
        # A free common color offset centers two drift differences at +/- .5;
        # each of the two worlds then requires radius .25.
        radius, coefficients = minimum_round_radius(matrix, np.array([0., 1.]))
        self.assertAlmostEqual(radius, .25)
        np.testing.assert_allclose(matrix @ coefficients, [0., 1.], atol=1e-12)

    def test_conditional_test_exact_error_at_boundary_probabilities(self):
        for n in [10, 100, 1000]:
            x = np.arange(n+1); y = n-x
            for delta in [0., .25]:
                for eta in [.5, 8.]:
                    decisions = paired_decision(x, y, delta, eta)
                    p0 = (1+delta)/(17+delta)
                    p1 = ((1-delta)*10+eta)/((17-delta)*10+eta)
                    self.assertLessEqual(binom.pmf(x, n, p0)[decisions == 1].sum(), .025+1e-12)
                    self.assertLessEqual(binom.pmf(x, n, p1)[decisions == 0].sum(), .025+1e-12)

    def test_background_mismatch_can_alias_joint_observation(self):
        b, delta, eta = 10., .25, 2.
        positive_target = (1-delta)*b+eta
        null_mismatch = positive_target/b-1
        self.assertLessEqual(abs(null_mismatch), delta)
        np.testing.assert_array_equal(poisson.pmf(np.arange(60), positive_target),
                                      poisson.pmf(np.arange(60), (1+null_mismatch)*b))
        # Both worlds also have exactly the same 16*b background intensity.
        self.assertEqual(16*b, 160)

    def test_centroid_refinement_on_independent_symmetric_spots(self):
        yy, xx = np.indices((120, 120))
        image = 20+1000*np.exp(-((yy-30)**2+(xx-30)**2)/4)
        image += 800*np.exp(-((yy-85)**2+(xx-85)**2)/4)
        centers, count = detect_centroids(image)
        self.assertEqual(count, 2)
        np.testing.assert_allclose(centers, [[30, 30], [85, 85]], atol=1e-10)


if __name__ == '__main__':
    unittest.main()
