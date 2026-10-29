import itertools
import unittest
import numpy as np
from molecular_permutation import evaluate, kernel_bounds, permutations, scores


def direct_score(positions, labels, bandwidth):
    a = positions[np.asarray(labels) == 1]
    b = positions[np.asarray(labels) == 0]
    return np.mean([np.exp(-np.sum((x-y)**2)/(2*bandwidth**2)) for x in a for y in b])


class MolecularPermutationTests(unittest.TestCase):
    def test_scores_equal_direct_cross_pairs(self):
        positions = np.array([[0., 0.], [1., 2.], [3., 1.], [5., 0.]])
        labels = np.array([1, 0, 1, 0])
        kernel = kernel_bounds(positions, 2., .1)[0]
        self.assertAlmostEqual(scores(kernel, labels)[0], direct_score(positions, labels, 2.), places=14)

    def test_stratum_counts_preserved(self):
        labels = np.array([1, 0, 1, 1, 0, 0])
        blocks = np.array([0, 0, 1, 1, 1, 1])
        shuffled = permutations(labels, blocks, 100, np.random.default_rng(3))
        for block in (0, 1):
            np.testing.assert_array_equal(shuffled[:, blocks == block].sum(axis=1), np.full(100, labels[blocks == block].sum()))
        self.assertGreater(len(np.unique(shuffled, axis=0)), 1)

    def test_bound_contains_all_latent_permutation_p_values(self):
        positions = np.c_[np.array([0., .4, 1.5, 2., 3.2, 4.]), np.zeros(6)]
        labels = np.array([1, 0, 1, 0, 1, 0])
        # Complete group of permutations within three fixed pairs.
        actions = []
        for flips in itertools.product((False, True), repeat=3):
            actions.append(np.concatenate([np.array([2*j, 2*j+1])[::-1] if flip else np.array([2*j, 2*j+1]) for j, flip in enumerate(flips)]))
        actions = np.array(actions)
        assignments = labels[actions]
        checked = 0
        for radius in (0., .2):
            for swaps in (0, 1):
                result = evaluate(positions, labels, assignments, 1., radius, swaps)
                for latent_labels in assignments:
                    if np.sum(latent_labels != labels) > 2*swaps:
                        continue
                    for shifts in itertools.product((0.,) if radius == 0 else (-radius, radius), repeat=6):
                        latent_positions = positions.copy()
                        latent_positions[:, 0] += shifts
                        observed = direct_score(latent_positions, latent_labels, 1.)
                        null = np.array([direct_score(latent_positions, latent_labels[action], 1.) for action in actions])
                        latent_p = (1 + np.sum(null >= observed - 1e-14))/(len(actions)+1)
                        self.assertLessEqual(latent_p, result['robust_upper_p'] + 1e-14)
                        checked += 1
        self.assertEqual(checked, 325)

    def test_balanced_global_proximity_is_reverse_mmd(self):
        positions = np.c_[np.arange(6.), np.zeros(6)]
        kernel = kernel_bounds(positions, 1., 0.)[0]
        m = 3
        values = []
        for a in itertools.combinations(range(6), m):
            labels = np.zeros(6, int); labels[list(a)] = 1
            sign = 2*labels-1
            mmd2 = sign @ kernel @ sign / m**2
            score = direct_score(positions, labels, 1.)
            self.assertAlmostEqual(score, kernel.sum()/(4*m*m)-mmd2/4, places=14)
            self.assertGreaterEqual(mmd2, -1e-14)
            values.append(score)
        null_mean = (kernel.sum()-2*m)/(2*m*(2*m-1))
        self.assertAlmostEqual(np.mean(values), null_mean, places=14)
        self.assertLessEqual(max(values)-null_mean, 1/(2*m-1))

    def test_exact_group_null_is_superuniform(self):
        positions = np.c_[np.arange(6.), np.zeros(6)]
        labels = []
        for a in itertools.combinations(range(6), 3):
            y = np.zeros(6, int); y[list(a)] = 1; labels.append(y)
        values = np.array([direct_score(positions, y, 1.) for y in labels])
        pvalues = np.array([np.mean(values >= value-1e-14) for value in values])
        for alpha in (.05, .1, .25, .5):
            self.assertLessEqual(np.mean(pvalues <= alpha), alpha + 1e-14)


if __name__ == '__main__':
    unittest.main()
