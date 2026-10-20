import itertools,unittest
import numpy as np
from coupled_marker_bounds import sharp_one_a_bounds


class CoupledMarkerTests(unittest.TestCase):
    def test_all_small_populations_and_allowed_completions_attain_bounds(self):
        edges=np.array([(0,1),(0,2),(1,0),(1,2),(2,0)])
        populations=np.array(list(itertools.product((False,True),repeat=6))).reshape(-1,3,2)
        targets=np.mean(populations[:,edges[:,0],0]&populations[:,edges[:,1],1],axis=1)
        for recorded in populations:
            for budget_b in range(4):
                allowed=(np.sum(populations[:,:,0]!=recorded[:,0],axis=1)<=1)&(np.sum(populations[:,:,1]!=recorded[:,1],axis=1)<=budget_b)
                expected=[np.min(targets[allowed]),np.max(targets[allowed])]
                np.testing.assert_allclose(sharp_one_a_bounds(edges,recorded,budget_b),expected,atol=1e-15)

    def test_individually_zero_marker_sensitivity_misses_joint_change(self):
        edges=np.array([(0,1)]);recorded=np.zeros((2,2),bool)
        np.testing.assert_equal(sharp_one_a_bounds(edges,recorded,0),[0,0])
        np.testing.assert_equal(sharp_one_a_bounds(edges,recorded,1),[0,1])


if __name__=='__main__':unittest.main()
