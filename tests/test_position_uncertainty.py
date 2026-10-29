from fractions import Fraction
from itertools import product
import unittest
import numpy as np
from position_uncertainty import position_bounds, exact_position_oracle


class PositionUncertaintyTests(unittest.TestCase):
    def test_complete_small_geometry_and_categorical_errors(self):
        raw = np.array([[0, 0], [3, 0], [0, 3], [3, 3]])
        positions = raw[:, None]+np.array([[0,0],[-1,0],[1,0]])
        labels = np.array([1,2,0,2])
        for budget in (0,1):
            bound = position_bounds(positions,labels,3,budget,budget)
            exact = exact_position_oracle(positions,labels,3,budget,budget)
            self.assertLessEqual(bound['outer_bounds'][0],exact['bounds'][0]+1e-10)
            self.assertGreaterEqual(bound['outer_bounds'][1],exact['bounds'][1]-1e-10)

    def test_covariance_is_one_lipschitz_in_total_variation(self):
        distributions = [tuple(Fraction(x,4) for x in counts) for counts in product(range(5),repeat=4) if sum(counts)==4]
        def covariance(p):return p[3]-(p[2]+p[3])*(p[1]+p[3])
        for p in distributions:
            for q in distributions:
                self.assertLessEqual(abs(covariance(p)-covariance(q)),sum(abs(a-b) for a,b in zip(p,q))/2)

    def test_no_edges_report_undefined_possibility(self):
        positions=np.array([[[0,0]],[[10,0]]]); labels=np.array([1,2])
        bound=position_bounds(positions,labels,2,0,0)
        self.assertTrue(bound['may_include_empty_graph'])
        self.assertEqual(bound['outer_bounds'],[-.25,.25])

    def test_exact_positions_and_labels_recover_the_target(self):
        positions=np.array([[[0,0]],[[1,0]],[[3,0]]]); labels=np.array([1,2,0])
        bound=position_bounds(positions,labels,2,0,0)
        exact=exact_position_oracle(positions,labels,2,0,0)
        np.testing.assert_allclose(bound['outer_bounds'],exact['bounds'],atol=1e-10,rtol=0)


if __name__ == '__main__':unittest.main()
