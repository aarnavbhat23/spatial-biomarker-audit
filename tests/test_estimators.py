import itertools
import unittest

import numpy as np

from estimators import betting_interval, betting_logwealth, difference_estimate, union_area


class TestSampling(unittest.TestCase):
    def test_exact_unbiasedness_and_variance_over_all_samples(self):
        truth = np.array([1.,4.,0.,7.,2.])
        predictions = np.array([2.,3.,3.,5.,1.])
        counts = np.array([3.,6.,4.,9.,3.])
        target = truth.sum()/counts.sum()
        for coefficient in (0.0,0.6,1.0):
            for size in (2,3,5):
                results = [difference_estimate(predictions,truth[list(indices)],indices,counts,coefficient)
                           for indices in itertools.combinations(range(5),size)]
                estimates = np.array([r['estimate'] for r in results])
                self.assertAlmostEqual(estimates.mean(),target,places=13)
                self.assertAlmostEqual(np.mean([r['variance'] for r in results]),estimates.var(),places=13)

    def test_census_recovers_truth_with_zero_width(self):
        result=difference_estimate([0,0,0],[1,2,3],[0,1,2],[2,3,4])
        self.assertAlmostEqual(result['estimate'],6/9)
        self.assertEqual(result['variance'],0)
        self.assertEqual(result['bounded_interval'],[result['estimate']]*2)

    def test_rare_event_exposes_zero_variance_student_failure(self):
        # One positive tile in 100: most size-two samples observe two zeros.
        counts=np.ones(100)
        result=difference_estimate(np.zeros(100),[0,0],[1,2],counts,coefficient=0)
        self.assertEqual(result['student_interval'],[0.,0.])
        self.assertLessEqual(result['bounded_interval'][0],0.01)
        self.assertGreaterEqual(result['bounded_interval'][1],0.01)

    def test_halos_charge_unique_area(self):
        self.assertEqual(union_area([[0,0,2,2],[1,0,3,2]]),6)
        self.assertEqual(union_area([[0,0,2,2],[0,0,2,2]]),4)

    def test_repeated_tile_is_rejected(self):
        with self.assertRaises(ValueError):
            difference_estimate([0,0,0],[0,0],[1,1],[1,1,1])

    def test_betting_wealth_has_unit_expectation_under_random_order(self):
        values=np.array([0.,0.1,0.3,0.6,1.])
        for direction in (True,False):
            wealth=[np.exp(betting_logwealth(values.mean(),values[list(order)],5,direction))
                    for order in itertools.permutations(range(5),3)]
            self.assertAlmostEqual(np.mean(wealth),1.,places=12)

    def test_betting_intervals_have_exact_enumerated_coverage(self):
        truth=np.array([0.,0.,1.,0.,3.])
        predicted=np.array([0.5,1.,1.,0.2,2.])
        counts=np.array([1.,1.,2.,1.,4.])
        target=truth.sum()/counts.sum()
        coverage=[]
        for indices in itertools.permutations(range(5),3):
            interval=betting_interval(predicted,truth[list(indices)],indices,counts,0.7)
            coverage.append(interval is not None and interval[0]-1e-12<=target<=interval[1]+1e-12)
        self.assertGreaterEqual(np.mean(coverage),0.95)


if __name__=='__main__':
    unittest.main()
