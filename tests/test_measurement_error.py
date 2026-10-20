import itertools
import unittest
import numpy as np
from measurement_error import pair_target,degree_error_margin,sharp_a_bounds,worst_a_flips,enlarge


class MeasurementErrorTests(unittest.TestCase):
    def setUp(self):
        self.edges=np.array([(0,1),(0,2),(0,3),(1,2),(2,0),(2,3),(3,1)])
        self.labels=np.array([[0,1],[1,0],[1,1],[0,1]],dtype=bool)

    def test_sharp_bounds_and_worst_flip_are_attained(self):
        for budget in range(5):
            targets=[]
            for bits in itertools.product((False,True),repeat=4):
                if np.sum(np.asarray(bits)!=self.labels[:,0])>budget:continue
                candidate=self.labels.copy();candidate[:,0]=bits
                targets.append(pair_target(self.edges,candidate))
            bounds=sharp_a_bounds(self.edges,self.labels,budget)
            np.testing.assert_allclose(bounds,[min(targets),max(targets)],atol=1e-15)
            changed,ids,_=worst_a_flips(self.edges,self.labels,budget)
            self.assertLessEqual(len(ids),budget)
            original=pair_target(self.edges,self.labels)
            self.assertAlmostEqual(abs(pair_target(self.edges,changed)-original),max(abs(bounds-original)))

    def test_geometry_margin_covers_all_two_marker_completions(self):
        original=pair_target(self.edges,self.labels)
        for budget_a in range(3):
            for budget_b in range(3):
                margin=degree_error_margin(self.edges,4,budget_a,budget_b)
                interval=enlarge([original,original],margin)
                for bits in itertools.product((False,True),repeat=8):
                    candidate=np.array(bits).reshape(4,2)
                    if np.sum(candidate[:,0]!=self.labels[:,0])>budget_a:continue
                    if np.sum(candidate[:,1]!=self.labels[:,1])>budget_b:continue
                    target=pair_target(self.edges,candidate)
                    self.assertLessEqual(interval[0]-1e-15,target)
                    self.assertLessEqual(target,interval[1]+1e-15)


if __name__=='__main__':unittest.main()
