import unittest
import numpy as np
from scipy.stats import poisson
from study_barcode_collisions import collision_classes, poisson_interval


class BarcodeCollisionTests(unittest.TestCase):
    def test_minimal_four_gene_witness_and_probe(self):
        codes = np.array([[1, 1, 0, 0], [0, 0, 1, 1], [1, 0, 1, 0], [0, 1, 0, 1]])
        worlds = collision_classes(codes)[(1, 1, 1, 1)]
        self.assertEqual(worlds, [(0, 1), (2, 3)])
        self.assertEqual([int(0 in w) for w in worlds], [1, 0])
        self.assertEqual([sum(codes[i, 0] for i in w) for w in worlds], [1, 1])

    def test_additive_sum_not_logical_or(self):
        groups = collision_classes([[1, 0], [0, 1]])
        self.assertEqual(groups[(2, 0)], [(0, 0)])
        self.assertEqual(groups[(1, 1)], [(0, 1)])
        with self.assertRaises(ValueError):
            collision_classes([[1, 0], [1, 0]])

    def test_exact_poisson_coverage_and_tail_inversion(self):
        counts = np.arange(100)
        lower, upper = poisson_interval(counts, 1)
        self.assertEqual(lower[0], 0)
        for rate in [.1, 1., 5., 20.]:
            coverage = poisson.pmf(counts, rate)[(lower <= rate)&(rate <= upper)].sum()
            self.assertGreaterEqual(coverage, .95-1e-12)
        np.testing.assert_allclose(poisson.sf(counts[1:]-1, lower[1:]), .025, atol=1e-12)

    def test_background_artifact_matches_all_repeated_counts(self):
        signal, positive_background, negative_background = 2., 8., 10.
        for exposures in [1, 4, 16, 64]:
            np.testing.assert_array_equal(
                poisson.pmf(np.arange(200), exposures*(positive_background+signal)),
                poisson.pmf(np.arange(200), exposures*negative_background))


if __name__ == '__main__':
    unittest.main()
