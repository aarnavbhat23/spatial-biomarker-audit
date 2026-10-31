import unittest
import numpy as np
from study_fiducial_measurement import design
from study_fiducial_sets import translation_intervals, maximum_separation


class FiducialSetTests(unittest.TestCase):
    def test_intersection_uses_extreme_matches_including_incompatibility(self):
        points = np.array([[0., 1.], [.4, 1.2], [.1, .8]])
        lo, hi = translation_intervals(points, .25)
        np.testing.assert_allclose(lo, [.15, .95])
        np.testing.assert_allclose(hi, [.25, 1.05])
        lo, hi = translation_intervals(points, .1)
        self.assertTrue(np.any(lo > hi))

    def test_fixed_reference_yields_analytic_difference_threshold(self):
        # First round is exact; one common color offset cannot cancel the
        # difference between two rounds. The second interval width is sharp.
        matrix = design([1, 2], [561, 561], 'rounds_and_colors')
        for width in [0., .1, 1., 2.]:
            answer = maximum_separation(matrix, [0., 1.], [0., width]+[0.]*6)
            self.assertAlmostEqual(answer['maximum_alias_separation'], width)

    def test_cycle_incompatibility_cannot_be_repaired_by_larger_allowance(self):
        matrix = design([1, 1, 2, 2], [561, 650, 561, 650], 'rounds_and_colors')
        answer = maximum_separation(matrix, [0., 1., 1., 0.], [0.]+[100.]*7)
        self.assertAlmostEqual(answer['maximum_alias_separation'], 0.)


if __name__ == '__main__':
    unittest.main()
