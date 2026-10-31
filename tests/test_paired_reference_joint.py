import unittest
import numpy as np
from scipy.stats import chi2, poisson
from study_paired_reference_joint import poisson_bounds, joint_decision


class JointReferenceTests(unittest.TestCase):
    def test_gamma_bounds_agree_with_chi_square_construction(self):
        counts = np.array([0, 1, 4, 100, 10000])
        lo, hi = poisson_bounds(counts, 16)
        expected = np.zeros(5)
        expected[1:] = chi2.ppf(.0125, 2*counts[1:])/32
        np.testing.assert_allclose(lo, expected)
        np.testing.assert_allclose(hi, chi2.ppf(.9875, 2*(counts+1))/32)

    def test_exact_poisson_tail_coverage_including_zero(self):
        for mean in [.01, .1, 1, 5, 30, 100]:
            counts = np.arange(int(poisson.ppf(1-1e-12, mean))+1)
            lo, hi = poisson_bounds(counts, 1)
            probability = poisson.pmf(counts, mean)
            self.assertLessEqual(probability[lo > mean].sum(), .0125+1e-12)
            self.assertLessEqual(probability[hi < mean].sum(), .0125+1e-12)

    def test_joint_error_by_direct_distribution_enumeration(self):
        # Sum the full joint law up to negligible tails, rather than simulate it.
        for truth, mismatch in [(0, .25), (1, -.25)]:
            mx, my = 2*(1+mismatch)+2*truth, 32
            x = np.arange(int(poisson.ppf(1-1e-12, mx))+1)
            y = np.arange(int(poisson.ppf(1-1e-12, my))+1)
            decision = joint_decision(x[:, None], y[None, :], 1, .25, 2)
            prob = poisson.pmf(x, mx)[:, None]*poisson.pmf(y, my)[None, :]
            wrong = (decision >= 0)&(decision != truth)
            self.assertLessEqual(prob[wrong].sum(), .025+2e-12)

    def test_interval_between_absence_and_minimum_signal_is_abstention(self):
        # At very high dose, a contrast of 4 is incompatible with both models
        # (absence <=0 or presence >=8): do not silently choose a hypothesis.
        self.assertEqual(joint_decision(np.array([140000]), np.array([1600000]), 10000, 0, 8)[0], -1)


if __name__ == '__main__':
    unittest.main()
