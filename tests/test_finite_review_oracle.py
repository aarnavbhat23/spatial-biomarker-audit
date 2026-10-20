import itertools,unittest
import numpy as np
from finite_review_oracle import FiniteReviewOracle


class FiniteReviewTests(unittest.TestCase):
    def test_feasible_set_matches_direct_error_counts_and_bounds(self):
        oracle=FiniteReviewOracle(2,[[(0,1)],[(0,1),(1,0)]],[(0,),(1,)])
        state=((0,1,1,0),(2,1,0,0),(1,1,0,0))
        expected=[]
        for bits in itertools.product((0,1),repeat=4):
            truth=np.array(bits).reshape(2,2);raw=np.array(state[0]).reshape(2,2)
            counts=np.array(state[1]).reshape(2,2);ones=np.array(state[2]).reshape(2,2)
            if np.any((truth!=raw).sum(axis=0)>1):continue
            if np.any(np.where(truth,counts-ones,ones).sum(axis=0)>1):continue
            expected.append(truth)
        self.assertEqual(len(oracle.feasible(state)),len(expected))
        targets=np.array([[t[0,0]*t[1,1],(t[0,0]*t[1,1]+t[1,0]*t[0,1])/2] for t in expected])
        np.testing.assert_allclose(oracle.bounds(state),np.stack([targets.min(axis=0),targets.max(axis=0)],axis=1))

    def test_repeated_noisy_review_can_identify_what_one_review_cannot(self):
        oracle=FiniteReviewOracle(2,[[(0,1)]],[(0,),(1,)],review_cap=3)
        raw=(1,0,0,1)
        once=(raw,(1,1,1,1),raw)
        self.assertEqual(oracle.bounds(once),((0.,1.),))
        twice=(raw,(2,2,2,2),tuple(2*x for x in raw))
        self.assertEqual(oracle.bounds(twice),((1.,1.),))
        with self.assertRaises(ValueError):
            oracle.bounds(((1,0,1,1),(3,0,3,0),(1,0,1,0)))

    def test_systematic_cell_errors_are_not_removed_by_repetition(self):
        oracle=FiniteReviewOracle(2,[[(0,1)]],[(0,),(1,)],review_error_unit='cells')
        raw=(1,0,0,1);state=(raw,(3,3,3,3),tuple(3*x for x in raw))
        self.assertEqual(oracle.bounds(state),((0.,1.),))

    def test_trusted_anchor_breaks_systematic_ambiguity(self):
        oracle=FiniteReviewOracle(2,[[(0,1)]],[(0,),(1,)],review_error_unit='cells',anchor_cost=1.)
        state=((1,1,1,1),(3,3,3,3),(3,3,3,3),(0,-1,-1,-1))
        self.assertEqual(oracle.bounds(state),((0.,0.),))

    def test_pruned_solver_matches_independent_bellman_recursion(self):
        from functools import lru_cache
        oracle=FiniteReviewOracle(2,[[(0,1)],[(0,1),(1,0)]],[(0,),(1,)],review_budgets=(0,0),review_cap=1)
        @lru_cache(None)
        def direct(state):
            if oracle.done(state):return 0.
            actions=oracle.actions(state)
            if not actions:return float('inf')
            return min(oracle.cost(action)+max(direct(child) for _,child in oracle.children(state,action)) for action in actions)
        self.assertAlmostEqual(oracle.optimal_value(oracle.initial),direct(oracle.initial))

    def test_no_error_observations_have_zero_ambiguity_at_census(self):
        oracle=FiniteReviewOracle(2,[[(0,1)]],[(0,),(1,)],budgets=(0,0),review_budgets=(0,0))
        state=((1,0,0,1),(0,)*4,(0,)*4)
        self.assertTrue(oracle.done(state))


if __name__=='__main__':unittest.main()
